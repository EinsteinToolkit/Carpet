#include <cassert>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "vect.hh"

#include "dh.hh"

using namespace std;



// Constructors
dh::dh (gh& h_,
        const ivect& lghosts_, const ivect& ughosts_,
        const int prolongation_order_space_, const int inner_buffer_width_,
        const ivect& lbuffers_, const ivect& ubuffers_)
  : h(h_),
    ghosts(lghosts_, ughosts_),
    prolongation_order_space(prolongation_order_space_),
    inner_buffer_width(inner_buffer_width_),
    buffers(lbuffers_, ubuffers_)
{
  assert (all(ghosts[0]>=0 and ghosts[1]>=0));
  assert (prolongation_order_space>=0);
  assert (inner_buffer_width>=0);
  assert (all(buffers[0]>=0 and buffers[1]>=0));
  h.add(this);
  CHECKPOINT;
  regrid ();
  for (int rl=0; rl<h.reflevels(); ++rl) {
    recompose (rl, false);
  }
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
void dh::regrid ()
{
  DECLARE_CCTK_PARAMETERS;
  
  CHECKPOINT;
  
  boxes.clear();

  allocate_bboxes();
  
  foreach_reflevel_component_mglevel (&dh::setup_allocate);
  foreach_reflevel_component_mglevel (&dh::setup_sync_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_multigrid_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_refinement_prolongation_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_refinement_boundary_prolongation_boxes);
  foreach_reflevel_component_mglevel (&dh::setup_refinement_restriction_boxes);
  foreach_reflevel_component_mglevel (&dh::trim_unsynced_boundaries);

  foreach_reflevel_component_mglevel (&dh::optimise_fields);

  calculate_bases();

  if (output_bboxes) {
    cout << endl << h << endl;
    foreach_reflevel_component_mglevel (&dh::do_output_bboxes);
    output_bases();
  }
  
  foreach_reflevel_component_mglevel (&dh::check_bboxes);
}

void dh::recompose (const int rl, const bool do_prolongate)
{
  assert (rl>=0 and rl<h.reflevels());
  
  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    (*f)->recompose_crop ();
  }
  
  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    (*f)->recompose_allocate (rl);
    for (comm_state state; not state.done(); state.step()) {
      (*f)->recompose_fill (state, rl, do_prolongate);
    }
    (*f)->recompose_free_old (rl);
    for (comm_state state; not state.done(); state.step()) {
      (*f)->recompose_bnd_prolongate (state, rl, do_prolongate);
    }
    for (comm_state state; not state.done(); state.step()) {
      (*f)->recompose_sync (state, rl, do_prolongate);
    }
  } // for all grid functions of same vartype
}

void dh::allocate_bboxes ()
{
  boxes.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    boxes.AT(ml).resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      boxes.AT(ml).AT(rl).resize(h.components(rl));
      for (int c=0; c<h.components(rl); ++c) {
        const ibbox intr = h.extent(ml,rl,c);
        dboxes & b = boxes.AT(ml).AT(rl).AT(c);



        // Interior
        // (the interior of the grid has the extent as specified by
        // the user)
        b.interior = intr;
        
        // Exterior (add ghost zones)
        // (the content of the exterior is completely determined by
        // the interior of this or other components; the content of
        // the exterior is redundant)
        i2vect dist(ghosts);    // for ghosts
        i2vect odist(0);        // for owned regions
        for (int f=0; f<2; ++f) {
          for (int d=0; d<dim; ++d) {
            if (h.outer_boundaries(rl,c)[f][d]) {
              dist[f][d] = 0;
              boxes.AT(ml).AT(rl).AT(c).is_interproc[d][f] = false;
            } else {
              bool const is_empty = intr.lower()[d] > intr.upper()[d];
              if (is_empty) {
                dist[f][d] = 0;
              }
              // Check whether the boundary in this direction is
              // covered by other interiors
              vect<ivect,2> dist1(0,0);
              dist1[f][d] = is_empty ? 0 : dist[f][d];
              ibset bnd = intr.expand(dist1[0], dist1[1]) - intr;
              for (int cc=0; cc<h.components(rl); ++cc) {
                bnd -= h.extent(ml,rl,cc);
              }
              bool const is_interproc = bnd.empty();
              boxes.AT(ml).AT(rl).AT(c).is_interproc[d][f] = is_interproc;
              if (! is_empty and ! is_interproc) {
                dist[f][d] += buffers[f][d];
                odist[f][d] = dist[f][d];
              }
            }
          }
        }
        
        b.exterior = intr.expand(dist[0], dist[1]);
        
        // Owned (interior plus buffer zones)
        // (will be made disjoint below)
        b.owned = intr.expand(odist[0], odist[1]);
        
        // Boundaries (ghost zones only)
        // (interior + boundaries = exterior)
        b.boundaries = b.exterior - intr;

      } // for c
      
      // Make owned regions disjoint
      for (int c=0; c<h.components(rl); ++c) {
        const ibbox intr = h.extent(ml,rl,c);
        dboxes & b = boxes.AT(ml).AT(rl).AT(c);
        
        // 1. Remove all other interiors from this owned region
        for (int cc=0; cc<h.components(rl); ++cc) {
          if (cc != c) {
            b.owned -= boxes.AT(ml).AT(rl).AT(cc).interior;
          }
        }
        b.owned.normalize();
        
        // 2. Make disjoint from all earlier owned regions
        for (int cc=0; cc<c; ++cc) {
          b.owned -= boxes.AT(ml).AT(rl).AT(cc).owned;
        }
        b.owned.normalize();
        
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
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        (this->*op)(box, rl, c, ml); // evaluate member function
      }
    }
  }
}

void dh::setup_allocate (dh::dboxes & b,
                         int const rl, int const c, int const ml)
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

void dh::setup_sync_boxes (dh::dboxes & box,
                           int const rl, int const c, int const ml)
{
  const ibset& bnds = box.boundaries;

  // Sync boxes
  for (int cc=0; cc<h.components(rl); ++cc) {
    dboxes & box1 = boxes.AT(ml).AT(rl).AT(cc);
    ibset ovlp;
    if (cc != c) {
      // Do not sync from the same component
      // Intersect boundaries with owned region of that component
      ovlp = bnds & box1.owned;
      ovlp.normalize();
    }
    for (ibset::const_iterator b=ovlp.begin(); b!=ovlp.end(); ++b) {
      box .recv_sync.AT(cc).push_back(*b);
      box1.send_sync.AT(c).push_back(*b);
    }
  }
}

void dh::setup_multigrid_boxes (dh::dboxes & box,
                                int const rl, int const c, int const ml)
{
  const ibbox& intr = box.interior;
  const ibbox& extr = box.exterior;

  // Multigrid boxes
  if (ml>0) {
    dboxes & bbox = boxes.AT(ml-1).AT(rl).AT(c);
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

void dh::setup_refinement_prolongation_boxes (dh::dboxes & box,
                                              int const rl, int const c, int const ml)
{
  const ibbox& extr = box.exterior;

  // Refinement boxes
  if (rl<h.reflevels()-1) {
    for (int cc=0; cc<h.components(rl+1); ++cc) {
      dboxes & box1 = boxes.AT(ml).AT(rl+1).AT(cc);
      const ibbox intrf = box1.interior;
      // Prolongation (interior)
      // TODO: prefer boxes from the same processor
      {
        // (the prolongation may use the exterior of the coarse
        // grid, and must fill all of the interior of the fine
        // grid)
        const int pss = prolongation_stencil_size();
        ibset recvs = extr.expand(-pss,-pss).contracted_for(intrf) & intrf;
        // Receive only once
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
        //
        for (ibset::const_iterator ri=recvs.begin(); ri!=recvs.end(); ++ri) {
          const ibbox recv = *ri;
          const ibbox send = recv.expanded_for(extr);
          assert (! send.empty());
          assert (send.is_contained_in(extr));
          box1.recv_ref_coarse.AT(c).push_back(recv);
          box. send_ref_fine  .AT(cc).push_back(send);
        }
      }
    } // for cc
  } // if not finest refinement level
}

void dh::setup_refinement_boundary_prolongation_boxes (dh::dboxes & box,
                                                       int const rl, int const c, int const ml)
{
  const ibbox& extr = box.exterior;

  // Refinement boxes
  if (rl<h.reflevels()-1) {
    for (int cc=0; cc<h.components(rl+1); ++cc) {
      dboxes & box1 = boxes.AT(ml).AT(rl+1).AT(cc);
      const ibbox & intrf = box1.interior;
      const ibbox & extrf = box1.exterior;
      const ibset & bndsf = box1.boundaries;
      const ibset & owndf = box1.owned;
      // Prolongation (boundaries)
      // TODO: prefer boxes from the same processor
      {
        // (the boundary prolongation may use the exterior of the
        // coarse grid, and must fill all of the owned boundary of the
        // fine grid)
        const int pss = prolongation_stencil_size();
        ibset recvs = bndsf & owndf;
        recvs.normalize();
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
        // Inner buffer zones
        {
          ibset bufs;
          for (ibset::const_iterator pbi=pbndsf.begin();
               pbi!=pbndsf.end(); ++pbi)
          {
            bufs |= (*pbi).expand(inner_buffer_width, inner_buffer_width);
          }
          bufs &= intrf;
          recvs |= bufs;
          recvs.normalize();
        }
        const ibbox maxrecvs = extr.expand(-pss,-pss).contracted_for(extrf);
        recvs &= maxrecvs;
        recvs.normalize();
        // Receive only once
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
        //
        {
          for (ibset::const_iterator ri = recvs.begin();
               ri != recvs.end(); ++ri)
          {
            const ibbox & recv = *ri;
            const ibbox send = recv.expanded_for(extr);
            assert (! send.empty());
            assert (send.is_contained_in(extr));
            assert (send.is_contained_in(extr.expand(-pss,-pss)));
            box1.recv_ref_bnd_coarse.AT(c).push_back(recv);
            box .send_ref_bnd_fine  .AT(cc).push_back(send);
          }
        }
      }
            
    } // for cc
  } // if not finest refinement level
}

void dh::setup_refinement_restriction_boxes (dh::dboxes & box,
                                             int const rl, int const c, int const ml)
{
  DECLARE_CCTK_PARAMETERS;

  const ibbox& intr = box.interior;

  // Refinement boxes
  if (rl<h.reflevels()-1) {
    for (int cc=0; cc<h.components(rl+1); ++cc) {
      dboxes & box1 = boxes.AT(ml).AT(rl+1).AT(cc);
      const ibbox intrf = box1.interior;
      // Restriction (interior)
      {
        // (the restriction may fill the interior of the of the
        // coarse grid, and may use the interior of the fine
        // grid, and the bbox must be as large as possible)
        // (the restriction must not use points that are filled
        // by boundary prolongation)
        // (the restriction must not fill points that are used for
        // boundary prolongation)
        ibset sends = intrf & intr.expanded_for(intrf);
        // remove what is received during boundary prolongation
        for (iblistvect::const_iterator rlvi = box1.recv_ref_bnd_coarse.begin();
             rlvi != box1.recv_ref_bnd_coarse.end(); ++ rlvi)
        {
          const iblist& recvlist = * rlvi;
          for (iblist::const_iterator rli = recvlist.begin();
               rli != recvlist.end(); ++ rli)
          {
            const ibbox& recv = * rli;
            sends -= recv;
          }
        }
        sends.normalize();
        // coarsify
        ibset recvs;
        for (ibset::const_iterator si = sends.begin();
             si != sends.end(); ++si)
        {
          const ibbox recv = (*si).contracted_for(intr);
          recvs |= recv;
        }
        if (omit_prolongation_points_when_restricting) {
          // remove what is sent during boundary prolongation
          const int pss = prolongation_stencil_size();
          for (int ccc=0; ccc<h.components(rl); ++ccc) {
            const dh::dboxes& box2 = boxes.AT(ml).AT(rl).AT(ccc);
            for (iblistvect::const_iterator slvi =
                   box2.send_ref_bnd_fine.begin();
                 slvi != box2.send_ref_bnd_fine.end(); ++ slvi)
            {
              const iblist& sendlist = * slvi;
              for (iblist::const_iterator sli = sendlist.begin();
                   sli != sendlist.end(); ++sli)
              {
                const ibbox& send = * sli;
                recvs -= send.expand(pss,pss);
              }
            }
          }
        }
        recvs.normalize();
        //
        for (ibset::const_iterator ri = recvs.begin();
             ri != recvs.end(); ++ri)
        {
          const ibbox recv = *ri;
          assert (! recv.empty());
          const ibbox & send = recv.expanded_for(intrf);
          assert (! send.empty());
          box1.send_ref_coarse.AT(c).push_back(send);
          box .recv_ref_fine  .AT(cc).push_back(recv);
        }
      }
            
    } // for cc
  } // if not finest refinement level
}

void dh::trim_unsynced_boundaries (dh::dboxes & box,
                                   int const rl, int const c, int const ml)
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

void
dh::
optimise_field (dboxes & box,
                iblistvect const dboxes::* const field,
                pvect dboxes::* const field_fast,
                int const rl, int const c, int const ml)
{
  size_t num_regions = 0;
  for (size_t cc=0; cc<(box.*field).size(); ++cc) {
    num_regions += (box.*field).AT(cc).size();
  }
  assert ((box.*field_fast).empty());
  (box.*field_fast).reserve (num_regions);
  for (size_t cc=0; cc<(box.*field).size(); ++cc) {
    for (iblist::const_iterator
           r=(box.*field).AT(cc).begin(); r!=(box.*field).AT(cc).end(); ++r)
    {
      pseudoregion pr;
      pr.extent = * r;
      pr.processor = cc;
      (box.*field_fast).push_back (pr);
    }
  }
  assert ((box.*field_fast).size() == num_regions);
}

void
dh::
optimise_fields (dboxes & box,
                 int const rl, int const c, int const ml)
{
  optimise_field (box, &dboxes::recv_ref_fine, &dboxes::recv_ref_fine_fast, rl, c, ml);
  optimise_field (box, &dboxes::recv_ref_coarse, &dboxes::recv_ref_coarse_fast, rl, c, ml);
  optimise_field (box, &dboxes::recv_sync, &dboxes::recv_sync_fast, rl, c, ml);
  optimise_field (box, &dboxes::recv_ref_bnd_coarse, &dboxes::recv_ref_bnd_coarse_fast, rl, c, ml);
}

void dh::check_bboxes (dh::dboxes & box,
                       int const rl, int const c, int const ml)
{
  DECLARE_CCTK_PARAMETERS;

  // Ensure that the bboxes are aligned with the base extent
  {
    if (ml==0) {
      if (rl==0) {
        assert (box.interior.is_aligned_with(h.baseextent));
      } else {
        // TODO: check alignment with next coarser grid
      }
    } else {
      // TODO: check alignment with next finer mglevel
    }
  }

  // Assert that all boundaries are synced or received
  {
    const ibset& sync_not = box.sync_not;
    const ibset& recv_not = box.recv_not;
          
    // Check that no boundaries are left over
    if (rl==0) assert (sync_not.empty());
#if 0
    assert (recv_not.empty());
#endif
    
    // Check that a component does not sync from itself
    assert (box.recv_sync.AT(c).empty());
    assert (box.send_sync.AT(c).empty());
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
    
#if 0
      // TODO
      // This just takes the number of boundary points into account.
      // It should really also take the location of the boundary into
      // account.

      // Check that whatever is left is exactly the outer boundary
      {
        // Determine whether Carpet::domain_from_coordbase is used
        int type;
        void const * ptr =
          CCTK_ParameterGet ("domain_from_coordbase", "Carpet", & type);
        assert (ptr != 0);
        assert (type == PARAMETER_BOOLEAN);
        CCTK_INT const domain_from_coordbase =
          * static_cast <CCTK_INT const *> (ptr);
        
        if (not domain_from_coordbase) {
          CCTK_WARN (2, "Cannot check correctness of grid structure near outer or symmetry boundaries because Carpet::domain_from_coordbase is not used");
        } else {
          
          // Find number of outer boundary points in each direction
          int const size = 2 * dim;
          CCTK_INT nboundaryzones_[size];
          CCTK_INT is_internal_[size];
          CCTK_INT is_staggered_[size];
          CCTK_INT shiftout_[size];
          int const ierr =
            GetBoundarySpecification
            (size, nboundaryzones_, is_internal_, is_staggered_, shiftout_);
          assert (not ierr);
          i2vect nboundaryzones;
          for (int f=0; f<2; ++f) {
            for (int d=0; d<dim; ++d) {
              nboundaryzones[f][d] = nboundaryzones_[2*d+f];
            }
          }
          
          // Calculate boundary bboxes
          b2vect const obs = h.outer_boundaries(rl,c);
          ibset bnd =
            box.interior - 
            box.interior.expand (- (ivect (obs[0]) * nboundaryzones[0]),
                                 - (ivect (obs[1]) * nboundaryzones[1]));
          bnd.normalize();
          
          // TODO: why should not intr == bnd?
          intr -= bnd;
          
          if (! intr.empty()) {
            intr.normalize();
            bnd.normalize();
            cout << "rl " << rl << " c " << c << endl;
            cout << "box.exterior: " << box.exterior << endl;
            cout << "box.interior: " << box.interior << endl;
            cout << "box.boundaries: " << box.boundaries << endl;
            cout << "intr: " << intr << endl;
            cout << "bnd: " << bnd << endl;
          }
          assert (intr.empty());
        }
      }
#endif
    }
  }
  
  // Assert that the owned regions are disjoint, that they together
  // make up the sum of all exterior regions, and that each
  // component's owned region contains its interior region
  if (c==0) {
    ibset combined_ownd;
    ibset combined_extr;
    for (int cc=0; cc<h.components(rl); ++cc) {
      
      dboxes const & box1 = boxes.AT(ml).AT(rl).AT(cc);
      
      ibbox const & intr = box1.interior;
      ibset const & ownd = box1.owned;
      ibbox const & extr = box1.exterior;
      
      ibset const ovlp = combined_ownd & ownd;
      assert (ovlp.empty());
      
      assert ((ownd & intr) == intr);
      
      combined_ownd |= ownd;
      combined_extr |= extr;
    }
    assert (combined_ownd == combined_extr);
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
        // The new size can be larger if part of the prolongation went
        // into the buffer zone.
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
  
  // Assert that points which are used for restricting are not
  // boundary prolongated
  {
    for (iblistvect::const_iterator rlvi = box.recv_ref_bnd_coarse.begin();
         rlvi != box.recv_ref_bnd_coarse.end(); ++ rlvi)
    {
      for (iblist::const_iterator rli = (*rlvi).begin();
           rli != (*rlvi).end(); ++ rli)
      {
        for (iblistvect::const_iterator slvi = box.send_ref_coarse.begin();
             slvi != box.send_ref_coarse.end(); ++ slvi)
        {
          for (iblist::const_iterator sli = (*slvi).begin();
               sli != (*slvi).end(); ++ sli)
          {
            assert ((*rli & *sli).empty());
          }
        }
      }
    }
  }
  
  // Assert that points which are used for boundary prolongation are
  // not restricted
  if (omit_prolongation_points_when_restricting) {
    for (iblistvect::const_iterator slvi = box.send_ref_bnd_fine.begin();
         slvi != box.send_ref_bnd_fine.end(); ++ slvi)
    {
      for (iblist::const_iterator sli = (*slvi).begin();
           sli != (*slvi).end(); ++ sli)
      {
        for (iblistvect::const_iterator rlvi = box.recv_ref_fine.begin();
             rlvi != box.recv_ref_fine.end(); ++ rlvi)
        {
          for (iblist::const_iterator rli = (*rlvi).begin();
               rli != (*rlvi).end(); ++ rli)
          {
            assert ((*sli & *rli).empty());
          }
        }
      }
    }
  }
  
  // Check proper nesting
  // "Proper nesting" means that prolongation from level L to level
  // L+1 does not use any points that are prolongated from level L-1
  // to level L.
  // We extend that notion to require a certain distance D in between.
  if (c == 0) {
    if (rl > 0 and rl < h.reflevels()) {
      // Points that are filled by prolongation from level rl-1
      ibset recvs;
      for (int cc=0; cc<h.components(rl); ++cc) {
        dboxes & box1 = boxes.AT(ml).AT(rl).AT(cc);
        for (iblistvect::const_iterator rlvi = box1.recv_ref_bnd_coarse.begin();
             rlvi != box1.recv_ref_bnd_coarse.end(); ++ rlvi)
        {
          const iblist & recvlist = * rlvi;
          for (iblist::const_iterator rli = recvlist.begin();
               rli != recvlist.end(); ++ rli)
          {
            const ibbox & recv = * rli;
            recvs |= recv;
          }
        }
      }
      recvs.normalize();
      // Extend the received points by the desired distance
      const ivect dist = proper_nesting_distance;
      ibset taboo;
      for (ibset::const_iterator bi = recvs.begin(); bi != recvs.end(); ++ bi) {
        const ibbox & b = * bi;
        const ibbox t = b.expand (dist, dist);
        taboo |= t;
      }
      taboo.normalize();
      // Points that are used for prolongation to level rl+1
      ibset sends;
      for (int cc=0; cc<h.components(rl); ++cc) {
        dboxes & box1 = boxes.AT(ml).AT(rl).AT(cc);
        for (iblistvect::const_iterator slvi = box1.send_ref_bnd_fine.begin();
             slvi != box1.send_ref_bnd_fine.end(); ++ slvi)
        {
          const iblist & sendlist = * slvi;
          for (iblist::const_iterator sli = sendlist.begin();
               sli != sendlist.end(); ++ sli)
          {
            const ibbox & send = * sli;
            sends |= send;
          }
        }
      }
      sends.normalize();
      // Calculate the overlap
      ibset overlap = sends & taboo;
      if (not overlap.empty()) {
        overlap.normalize();
        cout << "Not properly nested: rl=" << rl << ", "
             << "required distance=" << proper_nesting_distance << endl;
        cout << "Received from level " << rl-1 << ": " << recvs << endl;
        cout << "Taboo region on level " << rl << ": " << taboo << endl;
        cout << "Sent to level " << rl+1 << ": " << sends << endl;
        cout << "Overlap on level " << rl << ": " << overlap << endl;
        CCTK_WARN (1, "Not properly nested");
      }
    }
  }
}

void dh::calculate_bases ()
{
  // Calculate bases
  bases.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    bases.AT(ml).resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      dbases & base = bases.AT(ml).AT(rl);
      base.exterior = ibbox();
      base.interior = ibbox();
      for (int c=0; c<h.components(rl); ++c) {
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        base.exterior = (base.exterior.expanded_containing(box.exterior));
        base.interior = (base.interior.expanded_containing(box.interior));
      }
      base.boundaries = base.exterior - base.interior;
    }
  }
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
     << "ghosts=" << ghosts << ","
     << "gfs={";
  int cnt=0;
  for (list<ggf*>::const_iterator f = gfs.begin();
       f != gfs.end(); ++f) {
    if (cnt++) os << ",";
    (*f)->output(os);
  }
  os << "}";
}

void dh::do_output_bboxes(dh::dboxes & box,
                          int const rl, int const c, int const ml)
{
  cout << endl;
  cout << "dh bboxes:" << endl;
  cout << "ml=" << ml << " rl=" << rl << " c=" << c << endl;
  cout << "exterior=" << box.exterior << endl;
  cout << "is_interproc=" << box.is_interproc << endl;
  cout << "interior=" << box.interior << endl;
  cout << "owned=" << box.owned << endl;
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
        dbases & base = bases.AT(ml).AT(rl);
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
