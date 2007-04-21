#include <cassert>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "timestat.hh"
#include "vect.hh"

#include "dh.hh"

using namespace std;

using namespace CarpetLib;



// Constructors
dh::
dh (gh & h_,
    i2vect const & ghost_width_, i2vect const & buffer_width_,
    int const prolongation_order_space_)
  : h(h_),
    ghost_width(ghost_width_), buffer_width(buffer_width_),
    prolongation_order_space(prolongation_order_space_)
{
  assert (all (all (ghost_width >= 0)));
  assert (all (all (buffer_width >= 0)));
  assert (prolongation_order_space >= 0);
  h.add (this);
  CHECKPOINT;
  regrid ();
  for (int rl = 0; rl < h.reflevels(); ++ rl) {
    recompose (rl, false);
  }
}



// Destructors
dh::
~dh ()
{
  CHECKPOINT;
  h.remove (this);
}



// Helpers
int
dh::
prolongation_stencil_size ()
  const
{
  assert (prolongation_order_space >= 0);
  return prolongation_order_space / 2;
}



// Modifiers
void
dh::
regrid ()
{
  DECLARE_CCTK_PARAMETERS;
  
  CHECKPOINT;
  
  static Timer total ("dh::regrid");
  total.start ();
  
  oldboxes.clear();
  swap (boxes, oldboxes);
  
  
  
  boxes.resize (h.mglevels());
  for (int ml = 0; ml < h.mglevels(); ++ ml) {
    boxes.AT(ml).resize (h.reflevels());
    for (int rl = 0; rl < h.reflevels(); ++ rl) {
      boxes.AT(ml).AT(rl).resize (h.components(rl));
      
      cboxes & level = boxes.AT(ml).AT(rl);
      
      
      
      // Domain:
      
      ibbox const & domain_exterior = h.baseextent(ml,rl);
      assert (not domain_exterior.empty());
      
      i2vect const & boundary_width = h.boundary_width;
      assert (all (all (boundary_width >= 0)));
      
      ibbox const domain_active = domain_exterior.expand (- boundary_width);
      assert (not domain_active.empty());
      
      ibset domain_boundary = domain_exterior - domain_active;
      domain_boundary.normalize();
      
      
      
      for (int c = 0; c < h.components(rl); ++ c) {
        
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        
        
        
        // Interior:
        
        ibbox & intr = box.interior;
        
        // The interior of the grid has the extent as specified by the
        // regridding thorn
        intr = h.extent (ml,rl,c);
        
        // The interior must not be empty
        assert (not intr.empty());
        
        // The interior must be contained in the domain
        assert (intr <= h.baseextent(ml,rl));
        
        // All interiors must be disjunct
        for (int cc = 0; cc < c; ++ cc) {
          assert (not intr.intersects (level.AT(cc).interior));
        }
        
        
        
        // Outer boundary faces:
        
        b2vect & is_outer_boundary = box.is_outer_boundary;
        
        // The outer boundary faces are where the interior extends up
        // to the outer boundary of the domain.  It is not possible to
        // check whether it extends past the active part of the
        // domain, since this would be wrong when the outer boundary
        // width is zero.
        is_outer_boundary[0] = intr.lower() == domain_exterior.lower(); 
        is_outer_boundary[1] = intr.upper() == domain_exterior.upper(); 
        
        
        
        // Exterior:
        
        ibbox & extr = box.exterior;
        
        assert (all (all (ghost_width >= 0)));
        extr = intr.expand (i2vect (not is_outer_boundary) * ghost_width);
        
        // The exterior must not be empty
        assert (not extr.empty());
        
        // The exterior must be contained in the domain
        assert (extr <= domain_exterior);
        
        
        
        // Cactus ghost zones (which include outer boundaries):
        
        ibset & ghosts = box.ghosts;
        
        ghosts = extr - intr;
        ghosts.normalize();
        
        // The ghosts must be contained in the domain.  Different from
        // the boundaries, the ghost can include part of the outer
        // boundary of the domain.
        assert (ghosts <= domain_exterior);
        
        
        
        // Communicated region:
        
        ibbox & comm = box.communicated;
        
        comm = extr.expand (i2vect (is_outer_boundary) * (- boundary_width));
        
        // The communicated region must not be empty
        assert (not comm.empty());
        
        // The communicated region must be contained in the active
        // part of the domain
        assert (comm <= domain_active);
        
        
        
        // Outer boundary:
        
        ibset & outer_boundaries = box.outer_boundaries;
        
        outer_boundaries = extr - comm;
        outer_boundaries.normalize();
        
        // The outer boundary must be contained in the outer boundary
        // of the domain
        assert (outer_boundaries <= domain_boundary);
        
        
        
        // Owned region:
        
        ibbox & owned = box.owned;
        
        owned = intr.expand (i2vect (is_outer_boundary) * (- boundary_width));
        
        // The owned region must not be empty
        assert (not owned.empty());
        
        // The owned region must be contained in the active part of
        // the domain
        assert (owned <= domain_active);
        
        // All owned regions must be disjunct
        for (int cc = 0; cc < c; ++ cc) {
          assert (not owned.intersects (level.AT(cc).owned));
        }
        
        
        
        // Boundary (Carpet ghost zones, which do not include outer
        // boundaries):
        
        ibset & boundaries = box.boundaries;
        
        boundaries = comm - owned;
        boundaries.normalize();
        
        // The boundary must be contained in the active part of the
        // domain.  This prevents that a region is too close to the
        // outer boundary, so that it has ghost zones overlapping with
        // the outer boundary.
        assert (boundaries <= domain_active);
        
      } // for c
      
      
      
      // Conjunction of all buffer zones:
      
      // Enlarge active part of domain
      i2vect const safedist = i2vect (0);
      ibbox const domain_enlarged = domain_active.expand (safedist);
      
      // All owned regions
      ibset allowned;
      for (int c = 0; c < h.components(rl); ++ c) {
        dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
        allowned += box.owned;
      }
      allowned.normalize();
      assert (allowned <= domain_active);
      
      // All not-owned regions
      ibset notowned = domain_enlarged - allowned;
      notowned.normalize();
      
      // All not-active points
      ibset notactive;
      for (ibset::const_iterator
             ri = notowned.begin(); ri != notowned.end(); ++ ri)
      {
        ibbox const & r = * ri;
        ibbox const r_enlarged = r.expand (buffer_width);
        notactive |= r_enlarged;
      }
      notactive.normalize();
      
      // All buffer zones
      ibset allbuffers = allowned & notactive;
      allbuffers.normalize();
      
      // Buffer zones must be in the active part of the domain
      assert (allbuffers <= domain_active);
      
      
      
      for (int c = 0; c < h.components(rl); ++ c) {
        
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        
        
        
        // Buffer zones:
        
        box.buffers = box.owned & allbuffers;
        box.buffers.normalize();
        
        
        
        // Active region:
        
        box.active = box.owned - box.buffers;
        box.active.normalize();
        
      } // for c
      
      
      
      // The conjunction of all buffer zones must equal allbuffers
      
      ibset allbuffers1;
      for (int c = 0; c < h.components(rl); ++ c) {
        dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
        allbuffers1 += box.buffers;
      }
      allbuffers1.normalize();
      assert (allbuffers1 == allbuffers);
      
      
      
      // Test constituency relations:
      
      for (int c = 0; c < h.components(rl); ++ c) {
        dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
        
        assert ((box.active & box.buffers).empty());
        assert ((box.active | box.buffers) == box.owned);
        
        assert ((box.owned & box.boundaries).empty());
        assert ((box.owned | box.boundaries) == box.communicated);
        
        assert ((box.communicated & box.outer_boundaries).empty());
        assert ((box.communicated | box.outer_boundaries) == box.exterior);
        
        assert (box.boundaries <= box.ghosts);
        
        assert ((box.interior & box.ghosts).empty());
        assert ((box.interior | box.ghosts) == box.exterior);
        
      } // for c
      
      
      
      // Communication schedule:
      
      for (int c = 0; c < h.components(rl); ++ c) {
        
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        
        
        
        // Multigrid restriction:
        
        if (ml > 0) {
          int const oml = ml - 1;
          
          // Multigrid restriction must fill all active points
          
          dboxes const & obox = boxes.AT(oml).AT(rl).AT(c);
          
          ibset needrecv = box.active;
          
          ibset contracted_oactive;
          for (ibset::const_iterator
                 ai = obox.active.begin(); ai != obox.active.end(); ++ ai)
          {
            ibbox const & oactive = * ai;
            contracted_oactive += oactive.contracted_for (box.interior);
          }
          contracted_oactive.normalize();
          
          ibset ovlp = needrecv & contracted_oactive;
          ovlp.normalize();
          
          for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            box.mg_rest_recv.push_back (recv);
            ibbox const send = recv.expanded_for (obox.interior);
            assert (send <= obox.exterior);
            box.mg_rest_send.push_back (send);
          }
          
          needrecv -= ovlp;
          needrecv.normalize();
            
          // All points must have been received
          assert (needrecv.empty());
          
        } // if ml > 0
        
        
        
        // Multigrid prolongation:
        
        if (ml > 0) {
          int const oml = ml - 1;
          
          // Multigrid prolongation must fill all active points (this
          // could probably be relaxed)
          
          dboxes const & obox = boxes.AT(oml).AT(rl).AT(c);
          
          ibset oneedrecv = obox.active;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size());
          
          ibset expanded_active;
          for (ibset::const_iterator
                 ai = box.active.begin(); ai != box.active.end(); ++ ai)
          {
            ibbox const & active = * ai;
            expanded_active += active.expanded_for (obox.interior);
          }
          expanded_active.normalize();
          
          ibset ovlp = oneedrecv & expanded_active;
          ovlp.normalize();
            
          for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            box.mg_prol_recv.push_back (recv);
            ibbox const send =
              recv.expanded_for (box.interior).expand (stencil_size);
            assert (send <= box.exterior);
            box.mg_prol_send.push_back (send);
          }
          
          oneedrecv -= ovlp;
          oneedrecv.normalize();
          
          // All points must have been received
          assert (oneedrecv.empty());
          
        } // if ml > 0
        
        
        
        // Refinement prolongation:
        
        if (rl > 0) {
          int const orl = rl - 1;
          
          // Refinement prolongation must fill all active points
          
          ibset needrecv = box.active;
          
          box.ref_prol_recv.resize (h.components(orl));
          box.ref_prol_send.resize (h.components(orl));
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size());
          
          assert (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0));
          i2vect const reffact =
            i2vect (h.reffacts.at(rl) / h.reffacts.at(orl));
          
          for (int cc = 0; cc < h.components(orl); ++ cc) {
            dboxes const & obox = boxes.AT(ml).AT(orl).AT(cc);
            
            ibset contracted_oactive;
            for (ibset::const_iterator
                   ai = obox.active.begin(); ai != obox.active.end(); ++ ai)
            {
              ibbox const & oactive = * ai;
              // untested for cell centering
              contracted_oactive +=
                oactive.contracted_for (box.interior).expand (reffact);
            }
            contracted_oactive.normalize();
            
            ibset ovlp = needrecv & contracted_oactive;
            ovlp.normalize();
            
            for (ibset::const_iterator ri =
                   ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              box.ref_prol_recv.AT(cc).push_back (recv);
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              assert (send <= obox.exterior);
              box.ref_prol_send.AT(cc).push_back (send);
            }
            
            needrecv -= ovlp;
            
          } // for cc
          
          needrecv.normalize();
          
          // All points must have been received
          assert (needrecv.empty());
          
        } // if rl > 0
        
        
        
        // Synchronisation:
        
        // Synchronisation should fill as many boundary points as
        // possible
        
#if 0
        // Outer boundaries are not synchronised, since they cannot be
        // filled by boundary prolongation either, and therefore the
        // user code must set them anyway.
        ibset needrecv = box.boundaries;
#else
        // Outer boundaries are synchronised for backward
        // compatibility.
        ibset needrecv = box.ghosts;
#endif
        
        box.sync_recv.resize (h.components(rl));
        box.sync_send.resize (h.components(rl));
        
        ibset & sync = box.sync;
        
        for (int cc = 0; cc < h.components(rl); ++ cc) {
          dboxes const & obox = boxes.AT(ml).AT(rl).AT(cc);
          
#if 0
          ibset ovlp = needrecv & obox.owned;
#else
          ibset ovlp = needrecv & obox.interior;
#endif
          ovlp.normalize();
          
          if (cc == c) {
            assert (ovlp.empty());
          }
          
          for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            box.sync_recv.AT(cc).push_back (recv);
            ibbox const & send = recv;
            box.sync_send.AT(cc).push_back (send);
          }
          
          needrecv -= ovlp;
          sync += ovlp;
          
        } // for cc
        
        needrecv.normalize();
        sync.normalize();
        
        
        
        // Boundary prolongation:
        
        if (rl > 0) {
          int const orl = rl - 1;
          
          // Outer boundary points cannot be boundary prolongated
          needrecv &= box.communicated;
          
          // Prolongation must fill what cannot be synchronised, and
          // also all buffer zones
          needrecv += box.buffers;
          needrecv.normalize();
          
          box.ref_bnd_prol_recv.resize (h.components(orl));
          box.ref_bnd_prol_send.resize (h.components(orl));
          
          ibset & bndref = box.bndref;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size());
          
          assert (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0));
          i2vect const reffact =
            i2vect (h.reffacts.at(rl) / h.reffacts.at(orl));
          
          for (int cc = 0; cc < h.components(orl); ++ cc) {
            dboxes const & obox = boxes.AT(ml).AT(orl).AT(cc);
            
            ibset contracted_oactive;
            for (ibset::const_iterator
                   ai = obox.active.begin(); ai != obox.active.end(); ++ ai)
            {
              ibbox const & oactive = * ai;
              // untested for cell centering
              contracted_oactive +=
                oactive.contracted_for (box.interior).expand (reffact);
            }
            contracted_oactive.normalize();
            
            ibset ovlp = needrecv & contracted_oactive;
            ovlp.normalize();
            
            for (ibset::const_iterator ri =
                   ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              box.ref_bnd_prol_recv.AT(cc).push_back (recv);
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              assert (send <= obox.exterior);
              box.ref_bnd_prol_send.AT(cc).push_back (send);
            }
            
            needrecv -= ovlp;
            bndref += ovlp;
            
          } // for cc
          
          needrecv.normalize();
          bndref.normalize();
          
        } // if rl > 0
        
        // All points must now have been received, either through
        // synchronisation or through boundary prolongation
        assert (needrecv.empty());
        
        
        
        // Optimise fields:
        
        optimise_field
          (box, c, & dboxes::mg_rest_recv,     & dboxes::fast_mg_rest_recv);
        optimise_field
          (box, c, & dboxes::mg_rest_send,     & dboxes::fast_mg_rest_send);
        optimise_field
          (box, c, & dboxes::mg_prol_recv,     & dboxes::fast_mg_prol_recv);
        optimise_field
          (box, c, & dboxes::mg_prol_send,     & dboxes::fast_mg_prol_send);
        
        optimise_field
          (box, & dboxes::ref_prol_recv,     & dboxes::fast_ref_prol_recv);
        optimise_field
          (box, & dboxes::ref_prol_send,     & dboxes::fast_ref_prol_send);
        optimise_field
          (box, & dboxes::sync_recv,         & dboxes::fast_sync_recv);
        optimise_field
          (box, & dboxes::sync_send,         & dboxes::fast_sync_send);
        optimise_field
          (box, & dboxes::ref_bnd_prol_recv, & dboxes::fast_ref_bnd_prol_recv);
        optimise_field
          (box, & dboxes::ref_bnd_prol_send, & dboxes::fast_ref_bnd_prol_send);
        
        
        
      } // for c
      
      
      
      // Refinement restriction:
      
      if (rl > 0) {
        int const orl = rl - 1;
        
        ibset needrecv;
        for (int c = 0; c < h.components(rl); ++ c) {
          dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
          dboxes const & obox0 = boxes.AT(ml).AT(orl).AT(0);
          
          // Refinement restriction may fill all active points, and
          // must use all active points
          
          for (ibset::const_iterator
                 ai = box.active.begin(); ai != box.active.end(); ++ ai)
          {
            ibbox const & active = * ai;
            needrecv += active.contracted_for (obox0.interior);
          }
          needrecv.normalize();
        } // for c
        
        for (int cc = 0; cc < h.components(orl); ++ cc) {
          dboxes & obox = boxes.AT(ml).AT(orl).AT(cc);
          
          obox.ref_rest_recv.resize (h.components(rl));
          obox.ref_rest_send.resize (h.components(rl));
          
          for (int c = 0; c < h.components(rl); ++ c) {
            dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
            
            ibset contracted_active;
            for (ibset::const_iterator
                   ai = box.active.begin(); ai != box.active.end(); ++ ai)
            {
              ibbox const & active = * ai;
              contracted_active += active.contracted_for (obox.interior);
            }
            contracted_active.normalize();
            
            ibset ovlp = obox.active & contracted_active;
            ovlp.normalize();
            
            for (ibset::const_iterator ri =
                   ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              obox.ref_rest_recv.AT(c).push_back (recv);
              ibbox const send = recv.expanded_for (box.interior);
              assert (send <= box.active);
              obox.ref_rest_send.AT(c).push_back (send);
            }
            
            needrecv -= ovlp;
            
          } // for c
          
          optimise_field
            (obox, & dboxes::ref_rest_recv, & dboxes::fast_ref_rest_recv);
          optimise_field
            (obox, & dboxes::ref_rest_send, & dboxes::fast_ref_rest_send);
          
        } // for cc
        
        needrecv.normalize();
        
        // All points must have been received
        assert (needrecv.empty());
        
      } // if rl > 0
      
      
      
      // Regridding schedule:
      
      if (int (oldboxes.size()) > ml and int (oldboxes.AT(ml).size()) > rl) {
        
        int const oldcomponents = oldboxes.AT(ml).AT(rl).size();
        
        for (int c = 0; c < h.components(rl); ++ c) {
          
          dboxes & box = boxes.AT(ml).AT(rl).AT(c);
          
          
          
          // Synchronisation:
          
          // Synchronisation should fill as many active points as
          // possible
          
          ibset needrecv = box.active;
          
          box.old2new_sync_recv.resize (oldcomponents);
          box.old2new_sync_send.resize (oldcomponents);
          
          for (int cc = 0; cc < oldcomponents; ++ cc) {
            dboxes const & obox = oldboxes.AT(ml).AT(rl).AT(cc);
            
            ibset ovlp = needrecv & obox.owned;
            ovlp.normalize();
            
            for (ibset::const_iterator ri =
                   ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              box.old2new_sync_recv.AT(cc).push_back (recv);
              ibbox const & send = recv;
              box.old2new_sync_send.AT(cc).push_back (send);
            }
            
            needrecv -= ovlp;
            
          } // for cc
          
          needrecv.normalize();
          
          
          
          // Prolongation:
          
          if (rl > 0) {
            int const orl = rl - 1;
            
            // Prolongation must fill what cannot be synchronised
            
            box.old2new_ref_prol_recv.resize (h.components(orl));
            box.old2new_ref_prol_send.resize (h.components(orl));
            
            i2vect const stencil_size = i2vect (prolongation_stencil_size());
            
            assert (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0));
            i2vect const reffact =
              i2vect (h.reffacts.at(rl) / h.reffacts.at(orl));
            
            for (int cc = 0; cc < h.components(orl); ++ cc) {
              dboxes const & obox = boxes.AT(ml).AT(orl).AT(cc);
              
              ibset contracted_oactive;
              for (ibset::const_iterator
                     ai = obox.active.begin(); ai != obox.active.end(); ++ ai)
              {
                ibbox const & oactive = * ai;
                // untested for cell centering
                contracted_oactive +=
                  oactive.contracted_for (box.interior).expand (reffact);
              }
              contracted_oactive.normalize();
              
              ibset ovlp = needrecv & contracted_oactive;
              ovlp.normalize();
              
              for (ibset::const_iterator ri =
                     ovlp.begin(); ri != ovlp.end(); ++ ri)
              {
                ibbox const & recv = * ri;
                box.old2new_ref_prol_recv.AT(cc).push_back (recv);
                ibbox const send =
                  recv.expanded_for (obox.interior).expand (stencil_size);
                assert (send <= obox.exterior);
                box.old2new_ref_prol_send.AT(cc).push_back (send);
              }
              
              needrecv -= ovlp;
              
            } // for cc
            
            needrecv.normalize();
            
          } // if rl > 0
          
          // All points must now have been received, either through
          // synchronisation or through prolongation
          assert (needrecv.empty());
          
    
          
          // Optimise fields:
          
          optimise_field
            (box,
             & dboxes::old2new_sync_recv,
             & dboxes::fast_old2new_sync_recv);
          optimise_field
            (box,
             & dboxes::old2new_sync_send,
             & dboxes::fast_old2new_sync_send);
          optimise_field
            (box,
             & dboxes::old2new_ref_prol_recv,
             & dboxes::fast_old2new_ref_prol_recv);
          optimise_field
            (box,
             & dboxes::old2new_ref_prol_send,
             & dboxes::fast_old2new_ref_prol_send);
          
        } // for c
        
      } // if not oldboxes.empty
      
      
      
      // Output:
      
      if (output_bboxes) {
        
        for (int c = 0; c < h.components(rl); ++ c) {
          dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
          
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << box;
          cout << endl;
          
        } // for c
        
      } // if output_bboxes
      
      
      
    } // for rl
  }   // for m
  
  
  
  total.stop (0);
}



void
dh::
optimise_field (dboxes & box,
                iblistvect const dboxes::* const field,
                pvect dboxes::* const fast_field)
{
  size_t num_regions = 0;
  for (iblistvect::const_iterator
         li = (box.*field).begin(); li != (box.*field).end(); ++ li)
  {
    iblist const & l = * li;
    num_regions += l.size();
  }
  
  assert ((box.*fast_field).empty());
  (box.*fast_field).reserve (num_regions);
  
  {
    int p = 0;
    for (iblistvect::const_iterator
           li = (box.*field).begin(); li != (box.*field).end(); ++ li, ++ p)
    {
      iblist const & l = * li;
      for (iblist::const_iterator bi = l.begin(); bi != l.end(); ++ bi) {
        ibbox const & b = * bi;
        
        pseudoregion pr;
        pr.extent = b;
        pr.processor = p;
        
        (box.*fast_field).push_back (pr);
        
      }
    }
  }
  
  assert ((box.*fast_field).size() == num_regions);
}



void
dh::
optimise_field (dboxes & box,
                int const proc,
                iblist const dboxes::* const field,
                pvect dboxes::* const fast_field)
{
  size_t const num_regions = (box.*field).size();
  
  assert ((box.*fast_field).empty());
  (box.*fast_field).reserve (num_regions);
  
  iblist const & l = box.*field;
  for (iblist::const_iterator bi = l.begin(); bi != l.end(); ++ bi) {
    ibbox const & b = * bi;
    
    pseudoregion pr;
    pr.extent = b;
    pr.processor = proc;
    
    (box.*fast_field).push_back (pr);
  }
  
  assert ((box.*fast_field).size() == num_regions);
}



void
dh::
recompose (int const rl, bool const do_prolongate)
{
  assert (rl>=0 and rl<h.reflevels());
  
  static Timer timer ("dh::recompose");
  timer.start ();
  
  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    (*f)->recompose_crop ();
  }
  
  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    (*f)->recompose_allocate (rl);
    for (comm_state state; not state.done(); state.step()) {
      (*f)->recompose_fill (state, rl, do_prolongate);
    }
    (*f)->recompose_free_old (rl);
  } // for all grid functions of same vartype
  
  timer.stop (0);
}



// Grid function management
void
dh::
add (ggf * const f)
{
  CHECKPOINT;
  gfs.push_back (f);
}

void
dh::
remove (ggf * const f)
{
  CHECKPOINT;
  gfs.remove (f);
}



// Output
ostream &
dh::
output (ostream & os)
  const
{
  os << "dh:"
     << "ghost_width=" << ghost_width << ","
     << "buffer_width=" << buffer_width << ","
     << "prolongation_order_space=" << prolongation_order_space << ","
     << "boxes=" << boxes << ","
     << "gfs={";
  {
    bool isfirst = true;
    for (list<ggf*>::const_iterator
           f = gfs.begin(); f != gfs.end(); ++ f, isfirst = false)
    {
      if (not isfirst) os << ",";
      os << *f;
    }
  }
  os << "}";
  return os;
}

ostream &
dh::dboxes::
output (ostream & os)
  const
{
  os << "dh::dboxes:" << eol;
  
  // Regions:
  
  os << "regions:" << eol;
  os << "exterior:" << exterior << eol;
  os << "is_outer_boundary:" << is_outer_boundary << eol;
  os << "outer_boundaries:" << outer_boundaries << eol;
  os << "communicated:" << communicated << eol;
  os << "boundaries:" << boundaries << eol;
  os << "owned:" << owned << eol;
  os << "buffers:" << buffers << eol;
  os << "active:" << active << eol;
  os << "sync:" << sync << eol;
  os << "bndref:" << bndref << eol;
  os << "ghosts:" << ghosts << eol;
  os << "interior:" << interior << eol;
  
  // Communication schedule:
  
  os << "communication:" << eol;
  os << "mg_rest_recv: " << mg_rest_recv << eol;
  os << "mg_rest_send: " << mg_rest_send << eol;
  os << "mg_prol_recv: " << mg_prol_recv << eol;
  os << "mg_prol_send: " << mg_prol_send << eol;
  os << "ref_prol_recv:" << ref_prol_recv << eol;
  os << "ref_prol_send:" << ref_prol_send << eol;
  os << "ref_rest_recv:" << ref_rest_recv << eol;
  os << "ref_rest_send:" << ref_rest_send << eol;
  os << "sync_recv:" << sync_recv << eol;
  os << "sync_send:" << sync_send << eol;
  os << "ref_bnd_prol_recv:" << ref_bnd_prol_recv << eol;
  os << "ref_bnd_prol_send:" << ref_bnd_prol_send << eol;
  
  // Regridding schedule:
  
  os << "old2new_sync_recv:" << old2new_sync_recv << eol;
  os << "old2new_sync_send:" << old2new_sync_send << eol;
  os << "old2new_prol_recv:" << old2new_ref_prol_recv << eol;
  os << "old2new_prol_send:" << old2new_ref_prol_send << eol;
  
  return os;
}
