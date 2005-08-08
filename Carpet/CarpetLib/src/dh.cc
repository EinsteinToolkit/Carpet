#include <cassert>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "vect.hh"

#include "dh.hh"

using namespace std;


// structure to hold a set of grid functions
// which all have the same CCTK vartype
struct gf_set {
  int vartype;                // e.g. CCTK_VARIABLE_REAL, etc.
  vector<ggf*> members;       // members of this set
};



// Constructors
dh::dh (gh& h_,
        const ivect& lghosts_, const ivect& ughosts_,
        const int prolongation_order_space_, const int buffer_width_)
  : h(h_),
    lghosts(lghosts_), ughosts(ughosts_),
    prolongation_order_space(prolongation_order_space_),
    buffer_width(buffer_width_)
{
  assert (all(lghosts>=0 and ughosts>=0));
  assert (prolongation_order_space>=0);
  assert (buffer_width>=0);
  h.add(this);
  CHECKPOINT;
  recompose (false);
}

// Destructors
dh::~dh ()
{
  CHECKPOINT;
  h.remove(this);
}

// Helpers
int dh::prolongation_stencil_size () const
{
  assert (prolongation_order_space>=0);
  return prolongation_order_space/2;
}

// Modifiers
void dh::recompose (const bool do_prolongate)
{
  DECLARE_CCTK_PARAMETERS;
  
  CHECKPOINT;

  boxes.clear();

  allocate_bboxes();
  
  foreach_reflevel_component_mglevel (&dh::setup_sync_and_refine_boxes);
  foreach_reflevel_component_mglevel (&dh::intersect_sync_with_interior);
  foreach_reflevel_component_mglevel (&dh::setup_multigrid_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_refinement_interior_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_refinement_exterior_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_restrict_interior_boxes);
  foreach_reflevel_component_mglevel (&dh::trim_unsynced_boundaries);

  calculate_bases();

  if (output_bboxes) {
    cout << endl << h << endl;
    foreach_reflevel_component_mglevel (&dh::do_output_bboxes);
    output_bases();
  }
  
  foreach_reflevel_component_mglevel (&dh::assert_assert_assert);

  if (! save_memory_during_regridding) {
    save_time(do_prolongate);
  } else {
    save_memory(do_prolongate);
  }
}

void dh::allocate_bboxes ()
{
  boxes.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    boxes.at(ml).resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      boxes.at(ml).at(rl).resize(h.components(rl));
      for (int c=0; c<h.components(rl); ++c) {
        const ibbox intr = h.extents().at(ml).at(rl).at(c);
        dboxes & b = boxes.at(ml).at(rl).at(c);

        // Interior
        // (the interior of the grid has the extent as specified by
        // the user)
        b.interior = intr;

        // Exterior (add ghost zones)
        // (the content of the exterior is completely determined by
        // the interior of this or other components; the content of
        // the exterior is redundant)
        ivect ldist(lghosts), udist(ughosts);
        for (int d=0; d<dim; ++d) {
          if (h.outer_boundaries().at(rl).at(c)[d][0]) ldist[d] = 0;
          if (h.outer_boundaries().at(rl).at(c)[d][1]) udist[d] = 0;
        }
        b.exterior = intr.expand(ldist, udist);

        // Boundaries (ghost zones only)
        // (interior + boundaries = exterior)
        b.boundaries = b.exterior - intr;

      } // for c
    } // for rl
  } // for ml
}

// Loops over each multigrid level, each refinement level, and each
// component, executing the "boxesop" member function argument on the
// corresponding element of the "boxes" member
void dh::foreach_reflevel_component_mglevel (dh::boxesop op)
{
  for (int ml=0; ml<h.mglevels(); ++ml) {
    for (int rl=0; rl<h.reflevels(); ++rl) {
      for (int c=0; c<h.components(rl); ++c) {
        dboxes & box = boxes.at(ml).at(rl).at(c);
        (this->*op)(box, rl, c, ml); // evaluate member function
      }
    }
  }
}

void dh::setup_sync_and_refine_boxes (dh::dboxes & b, int rl, int c, int ml)
{
  // Sync boxes
  const int cs = h.components(rl);
  b.send_sync.resize(cs);
  b.recv_sync.resize(cs);

  // Refinement boxes
  if (rl>0) {
    const int csm1 = h.components(rl-1);
    b.send_ref_coarse.resize(csm1);
    b.recv_ref_coarse.resize(csm1);
    b.recv_ref_bnd_coarse.resize(csm1);
  }
  if (rl<h.reflevels()-1) {
    const int csp1 = h.components(rl+1);
    b.recv_ref_fine.resize(csp1);
    b.send_ref_fine.resize(csp1);
    b.send_ref_bnd_fine.resize(csp1);
  }
}

void dh::intersect_sync_with_interior (dh::dboxes & box, int rl, int c, int ml)
{
  const ibset& bnds = box.boundaries;

  // Sync boxes
  for (int cc=0; cc<h.components(rl); ++cc) {
    dboxes & box1 = boxes.at(ml).at(rl).at(cc);
    // intersect boundaries with interior of that component
    ibset ovlp = bnds & box1.interior;
    ovlp.normalize();
    for (ibset::const_iterator b=ovlp.begin();b!=ovlp.end(); ++b) {
      box .recv_sync.at(cc).push_back(*b);
      box1.send_sync.at(c).push_back(*b);
    }
  }
}

void dh::setup_multigrid_boxes (dh::dboxes & box, int rl, int c, int ml)
{
  const ibbox& intr = box.interior;
  const ibbox& extr = box.exterior;

  // Multigrid boxes
  if (ml>0) {
    dboxes & bbox = boxes.at(ml-1).at(rl).at(c);
    const ibbox intrf = bbox.interior;
    const ibbox extrf = bbox.exterior;
    // Restriction (interior)
    {
      // (the restriction must fill all of the interior of the
      // coarse grid, and may use the exterior of the fine grid)
      const ibbox recv = intr;
      assert (intr.empty() or ! recv.empty());
      const ibbox send = recv.expanded_for(extrf);
      assert (intr.empty() or ! send.empty());
      // TODO: put the check back in, taking outer boundaries
      // into account
#if 0
      assert (send.is_contained_in(extrf));
#endif
      bbox.send_mg_coarse.push_back(send);
      box .recv_mg_fine  .push_back(recv);
    }
    // Prolongation (interior)
    {
      // (the prolongation may use the exterior of the coarse
      // grid, and may fill only the interior of the fine grid,
      // and the bbox must be as large as possible)
      const ibbox recv = extr.contracted_for(intrf) & intrf;
            assert (intr.empty() or ! recv.empty());
      const ibbox send = recv.expanded_for(extr);
            assert (intr.empty() or ! send.empty());
      bbox.recv_mg_coarse.push_back(recv);
      box .send_mg_fine  .push_back(send);
    }
  } // if not finest multigrid level
}

void dh::setup_refinement_interior_boxes (dh::dboxes & box, int rl, int c, int ml)
{
  const ibbox& extr = box.exterior;

  // Refinement boxes
  if (rl<h.reflevels()-1) {
    for (int cc=0; cc<h.components(rl+1); ++cc) {
      dboxes & box1 = boxes.at(ml).at(rl+1).at(cc);
      const ibbox intrf = box1.interior;
      // Prolongation (interior)
      // TODO: prefer boxes from the same processor
      {
        // (the prolongation may use the exterior of the coarse
        // grid, and must fill all of the interior of the fine
        // grid)
        const int pss = prolongation_stencil_size();
        ibset recvs = extr.expand(-pss,-pss).contracted_for(intrf) & intrf;
        const iblistvect& rrc = box1.recv_ref_coarse;
        for (iblistvect::const_iterator lvi=rrc.begin();
             lvi!=rrc.end(); ++lvi)
        {
          for (iblist::const_iterator li=lvi->begin();
               li!=lvi->end(); ++li)
          {
            recvs -= *li;
          }
        }
        recvs.normalize();
        assert (recvs.setsize() <= 1);
        if (recvs.setsize() == 1) {
          const ibbox recv = *recvs.begin();
          const ibbox send = recv.expanded_for(extr);
          assert (! send.empty());
          assert (send.is_contained_in(extr));
          box1.recv_ref_coarse.at(c).push_back(recv);
          box. send_ref_fine  .at(cc).push_back(send);
        }
      }
    } // for cc
  } // if not finest refinement level
}

void dh::setup_refinement_exterior_boxes (dh::dboxes & box, int rl, int c, int ml)
{
  const ibbox& extr = box.exterior;

  // Refinement boxes
  if (rl<h.reflevels()-1) {
    for (int cc=0; cc<h.components(rl+1); ++cc) {
      dboxes & box1 = boxes.at(ml).at(rl+1).at(cc);
      const ibbox intrf = box1.interior;
      const ibbox& extrf = box1.exterior;
      const ibset& bndsf = box1.boundaries;
      // Prolongation (boundaries)
      // TODO: prefer boxes from the same processor
      {
        // (the prolongation may use the exterior of the coarse
        // grid, and must fill all of the boundary of the fine
        // grid)
        const int pss = prolongation_stencil_size();
        // Prolongation boundaries
        ibset pbndsf = bndsf;
        {
          // Do not count what is synced
          const iblistvect& rs = box1.recv_sync;
          for (iblistvect::const_iterator lvi=rs.begin();
               lvi!=rs.end(); ++lvi)
          {
            for (iblist::const_iterator li=lvi->begin();
                 li!=lvi->end(); ++li)
            {
              pbndsf -= *li;
            }
          }
          pbndsf.normalize();
        }
        // Buffer zones
        ibset buffers;
        {
          for (ibset::const_iterator pbi=pbndsf.begin();
               pbi!=pbndsf.end(); ++pbi)
          {
            buffers |= (*pbi).expand(buffer_width, buffer_width) & extrf;
          }
          buffers.normalize();
        }
        // Add boundaries
        const ibbox maxrecvs = extr.expand(-pss,-pss).contracted_for(extrf);
        ibset recvs = buffers & maxrecvs;
        recvs.normalize();
        {
          // Do not prolongate what is already prolongated
          const iblistvect& rrbc = box1.recv_ref_bnd_coarse;
          for (iblistvect::const_iterator lvi=rrbc.begin();
               lvi!=rrbc.end(); ++lvi)
          {
            for (iblist::const_iterator li=lvi->begin();
                 li!=lvi->end(); ++li)
            {
              recvs -= *li;
            }
          }
          recvs.normalize();
        }
        {
          for (ibset::const_iterator ri = recvs.begin();
               ri != recvs.end(); ++ri)
          {
            const ibbox & recv = *ri;
            const ibbox send = recv.expanded_for(extr);
            assert (! send.empty());
            assert (send.is_contained_in(extr));
            assert (send.is_contained_in(extr.expand(-pss,-pss)));
            box1.recv_ref_bnd_coarse.at(c).push_back(recv);
            box .send_ref_bnd_fine  .at(cc).push_back(send);
          }
        }
      }
            
    } // for cc
  } // if not finest refinement level
}

void dh::setup_restrict_interior_boxes (dh::dboxes & box, int rl, int c, int ml)
{
  const ibbox& intr = box.interior;

  // Refinement boxes
  if (rl<h.reflevels()-1) {
    for (int cc=0; cc<h.components(rl+1); ++cc) {
      dboxes & box1 = boxes.at(ml).at(rl+1).at(cc);
      const ibbox intrf = box1.interior;
      // Restriction (interior)
      {
        // (the restriction may fill the interior of the of the
        // coarse grid, and may use the interior of the fine
        // grid, and the bbox must be as large as possible)
        // (the restriction must not use points that are filled
        // by boundary prolongation)
        ibset sends = intrf & intr.expanded_for(intrf);
        // remove what is received during boundary prolongation
        for (int ccc=0; ccc<h.components(rl); ++ccc) {
          const iblist& sendlist = box1.recv_ref_bnd_coarse.at(ccc);
          for (iblist::const_iterator sli = sendlist.begin();
               sli != sendlist.end(); ++sli)
          {
            sends -= *sli;
          }
        }
        sends.normalize();
        for (ibset::const_iterator si = sends.begin();
             si != sends.end(); ++si)
        {
          const ibbox recv = (*si).contracted_for(intr);
          if (! recv.empty()) {
            const ibbox & send = recv.expanded_for(intrf);
            assert (! send.empty());
            box1.send_ref_coarse.at(c).push_back(send);
            box .recv_ref_fine  .at(cc).push_back(recv);
          }
        }
      }
            
    } // for cc
  } // if not finest refinement level
}

void dh::trim_unsynced_boundaries (dh::dboxes & box, int rl, int c, int ml)
{
  // Boundaries that are not synced, or are neither synced nor
  // prolonged to from coarser grids (outer boundaries)
  ibset& sync_not = box.sync_not;
  ibset& recv_not = box.recv_not;
  
  // The whole boundary
  sync_not = box.boundaries;
  recv_not = box.boundaries;
  
  // Subtract boxes received during synchronisation
  const iblistvect& recv_sync = box.recv_sync;
  for (iblistvect::const_iterator lvi=recv_sync.begin();
       lvi!=recv_sync.end(); ++lvi)
  {
    for (iblist::const_iterator li=lvi->begin();
         li!=lvi->end(); ++li)
    {
      sync_not -= *li;
      recv_not -= *li;
    }
  }
  
  // Subtract boxes received during prolongation
  const iblistvect& recv_ref_bnd_coarse = box.recv_ref_bnd_coarse;
  for (iblistvect::const_iterator lvi=recv_ref_bnd_coarse.begin();
       lvi!=recv_ref_bnd_coarse.end(); ++lvi)
  {
    for (iblist::const_iterator li=lvi->begin();
         li!=lvi->end(); ++li)
    {
      recv_not -= *li;
    }
  }
}

void dh::assert_assert_assert (dh::dboxes & box, int rl, int c, int ml)
{
// Assert that all boundaries are synced or received
  {
    const ibset& sync_not = box.sync_not;
#if 0
    const ibset& recv_not = box.recv_not;
#endif
          
    // Check that no boundaries are left over
    if (rl==0) assert (sync_not.empty());
#if 0
    assert (recv_not.empty());
#endif
  }

  // Assert that the interior is received exactly once during
  // prolongation, and that nothing else is received
  {
    if (rl==0) {
      const iblistvect& recv_ref_coarse = box.recv_ref_coarse;
      assert (recv_ref_coarse.empty());
    } else {                    // rl!=0
      const iblistvect& recv_ref_coarse = box.recv_ref_coarse;
      ibset intr = box.interior;
      ibset received;
      for (iblistvect::const_iterator lvi=recv_ref_coarse.begin();
           lvi!=recv_ref_coarse.end(); ++lvi)
      {
        for (iblist::const_iterator li=lvi->begin(); li!=lvi->end(); ++li) {
          const int old_sz = intr.size();
          const int this_sz = li->size();
          intr -= *li;
          const int new_sz = intr.size();
          // TODO
          assert (new_sz + this_sz == old_sz);
          assert ((received & *li).empty());
          received |= *li;
        }
      }
      // TODO
      // This need not be empty at outer boundaries.  Check that
      // those are indeed outer boundaries!  But what size of the
      // boundary region should be used for that?
#if 0
      assert (intr.empty());
#endif
    }
  }
  
  // Assert that the boundaries are received at most once during
  // prolongation and synchronisation, and that nothing else is
  // received
  {
    const iblistvect& recv_sync = box.recv_sync;
    const iblistvect& recv_ref_bnd_coarse = box.recv_ref_bnd_coarse;
    ibset bnds = box.boundaries;
    ibset received;
    for (iblistvect::const_iterator lvi=recv_sync.begin();
         lvi!=recv_sync.end(); ++lvi)
    {
      for (iblist::const_iterator li=lvi->begin();
           li!=lvi->end(); ++li)
      {
        const int old_sz = bnds.size();
        const int this_sz = li->size();
        bnds -= *li;
        const int new_sz = bnds.size();
        assert (new_sz + this_sz == old_sz);
        assert ((received & *li).empty());
        received |= *li;
      }
    }
    for (iblistvect::const_iterator lvi=recv_ref_bnd_coarse.begin();
         lvi!=recv_ref_bnd_coarse.end(); ++lvi)
    {
      for (iblist::const_iterator li=lvi->begin(); li!=lvi->end(); ++li) {
        const int old_sz = bnds.size();
        const int this_sz = li->size();
        bnds -= *li;
        const int new_sz = bnds.size();
        // TODO
        // The new size can be larger if part of the
        // prolongation went into the buffer zone.
//      assert (new_sz + this_sz == old_sz);
        assert (new_sz + this_sz >= old_sz);
#if 0
        assert ((received & *li).empty());
        received |= *li;
#endif
      }
    }
    // TODO
    // This need not be empty at outer boundaries.  Check that
    // those are indeed outer boundaries!  But what size of the
    // boundary region should be used for that?
#if 0
    assert (bnds.empty());
#endif
  }
}

void dh::calculate_bases ()
{
  // Calculate bases
  bases.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    bases.at(ml).resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      dbases & base = bases.at(ml).at(rl);
      base.exterior = ibbox();
      base.interior = ibbox();
      for (int c=0; c<h.components(rl); ++c) {
        dboxes & box = boxes.at(ml).at(rl).at(c);
        base.exterior = (base.exterior.expanded_containing(box.exterior));
        base.interior = (base.interior.expanded_containing(box.interior));
      }
      base.boundaries = base.exterior - base.interior;
    }
  }
}

void dh::save_time (bool do_prolongate)
{
  DECLARE_CCTK_PARAMETERS;

  // sort all grid functions into sets of the same vartype
  vector<gf_set> ggfs;
  for (list<ggf*>::iterator f = gfs.begin(); f != gfs.end(); ++f) {
    gf_set newset;
    newset.vartype = CCTK_VarTypeI ((*f)->varindex);
    assert (newset.vartype >= 0);
    int c;
    for (c = 0; c < ggfs.size(); c++) {
      if (newset.vartype == ggfs[c].vartype) {
        break;
      }
    }
    if (c == ggfs.size()) {
      ggfs.push_back (newset);
    }
    ggfs[c].members.push_back (*f);
  }

  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    (*f)->recompose_crop ();
  }
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c = 0; c < ggfs.size(); c++) {
      for (int g = 0; g < ggfs[c].members.size(); g++) {
        ggfs[c].members[g]->recompose_allocate (rl);
      }
      for (comm_state state(ggfs[c].vartype); ! state.done(); state.step()) {
        for (int g = 0; g < ggfs[c].members.size(); g++) {
          ggfs[c].members[g]->recompose_fill (state, rl, do_prolongate);
        }
      }
      for (int g = 0; g < ggfs[c].members.size(); g++) {
        ggfs[c].members[g]->recompose_free (rl);
      }
      for (comm_state state(ggfs[c].vartype); ! state.done(); state.step()) {
        for (int g = 0; g < ggfs[c].members.size(); g++) {
          ggfs[c].members[g]->recompose_bnd_prolongate (state, rl, do_prolongate);
        }
      }
      for (comm_state state(ggfs[c].vartype); ! state.done(); state.step()) {
        for (int g = 0; g < ggfs[c].members.size(); g++) {
          ggfs[c].members[g]->recompose_sync (state, rl, do_prolongate);
        }
      }
    } // for all grid functions of same vartype
  } // for all refinement levels
}

void dh::save_memory (bool do_prolongate) {
  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    
    (*f)->recompose_crop ();
    for (int rl=0; rl<h.reflevels(); ++rl) {
      (*f)->recompose_allocate (rl);
      for (comm_state state; !state.done(); state.step()) {
        (*f)->recompose_fill (state, rl, do_prolongate);
      }
      (*f)->recompose_free (rl);
      for (comm_state state; !state.done(); state.step()) {
        (*f)->recompose_bnd_prolongate (state, rl, do_prolongate);
      }
      for (comm_state state; !state.done(); state.step()) {
        (*f)->recompose_sync (state, rl, do_prolongate);
      }
    } // for rl
    
  } // for gf
}

// Grid function management
void dh::add (ggf* f)
{
  CHECKPOINT;
  gfs.push_back(f);
}

void dh::remove (ggf* f)
{
  CHECKPOINT;
  gfs.remove(f);
}


// Output
void dh::output (ostream& os) const
{
  os << "dh:"
     << "ghosts=[" << lghosts << "," << ughosts << "],"
     << "gfs={";
  int cnt=0;
  for (list<ggf*>::const_iterator f = gfs.begin();
       f != gfs.end(); ++f) {
    if (cnt++) os << ",";
    (*f)->output(os);
  }
  os << "}";
}

void dh::do_output_bboxes(dh::dboxes & box, int rl, int c, int ml)
{
  cout << endl;
  cout << "dh bboxes:" << endl;
  cout << "ml=" << ml << " rl=" << rl << " c=" << c << endl;
  cout << "exterior=" << box.exterior << endl;
  cout << "interior=" << box.interior << endl;
  cout << "send_mg_fine=" << box.send_mg_fine << endl;
  cout << "send_mg_coarse=" << box.send_mg_coarse << endl;
  cout << "recv_mg_fine=" << box.recv_mg_fine << endl;
  cout << "recv_mg_coarse=" << box.recv_mg_coarse << endl;
  cout << "send_ref_fine=" << box.send_ref_fine << endl;
  cout << "send_ref_coarse=" << box.send_ref_coarse << endl;
  cout << "recv_ref_fine=" << box.recv_ref_fine << endl;
  cout << "recv_ref_coarse=" << box.recv_ref_coarse << endl;
  cout << "send_sync=" << box.send_sync << endl;
  cout << "send_ref_bnd_fine=" << box.send_ref_bnd_fine << endl;
  cout << "boundaries=" << box.boundaries << endl;
  cout << "recv_sync=" << box.recv_sync << endl;
  cout << "recv_ref_bnd_coarse=" << box.recv_ref_bnd_coarse << endl;
  cout << "sync_not=" << box.sync_not << endl;
  cout << "recv_not=" << box.recv_not << endl;
}

void dh::output_bases ()
{
  for (int ml=0; ml<h.mglevels(); ++ml) {
    for (int rl=0; rl<h.reflevels(); ++rl) {
      if (h.components(rl)>0) {
        dbases & base = bases.at(ml).at(rl);
        cout << endl;
        cout << "dh bases:" << endl;
        cout << "ml=" << ml << " rl=" << rl << endl;
        cout << "exterior=" << base.exterior << endl;
        cout << "interior=" << base.interior << endl;
        cout << "boundaries=" << base.boundaries << endl;
      }
    }
  }
}
