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

// Calculate this quantity on this processor?  It does not need to be
// calculated if it won't be used later on.
inline
int
dh::this_proc (int const rl, int const c) const
{
  return h.processor (rl, c);
}

inline
bool
dh::on_this_proc (int const rl, int const c) const
{
  return this_proc (rl, c) == dist::rank();
}

inline
int
dh::this_oldproc (int const rl, int const c) const
{
  return h.old_processor (rl, c);
}

inline
bool
dh::on_this_oldproc (int const rl, int const c) const
{
  return this_oldproc (rl, c) == dist::rank();
}



bool there_was_an_error = false;

static
void
assert_error (char const * restrict const checkstring,
              int const ml, int const rl,
              char const * restrict const message)
{
  CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "[ml=%d rl=%d] The following grid structure consistency check failed:\n   %s\n   %s",
              ml, rl, message, checkstring);
  there_was_an_error = true;
}

static
void
assert_error (char const * restrict const checkstring,
              int const ml, int const rl, int const c,
              char const * restrict const message)
{
  CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "[ml=%d rl=%d c=%d] The following grid structure consistency check failed:\n   %s\n   %s",
              ml, rl, c, message, checkstring);
  there_was_an_error = true;
}

static
void
assert_error (char const * restrict const checkstring,
              int const ml, int const rl, int const c, int const cc,
              char const * restrict const message)
{
  CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "[ml=%d rl=%d c=%d cc=%d] The following grid structure consistency check failed:\n   %s\n   %s",
              ml, rl, c, cc, message, checkstring);
  there_was_an_error = true;
}

#define ASSERT_rl(check, message)                       \
  do {                                                  \
    if (not (check)) {                                  \
      assert_error (#check, ml, rl, message);           \
    }                                                   \
  } while (false)

#define ASSERT_c(check, message)                        \
  do {                                                  \
    if (not (check)) {                                  \
      assert_error (#check, ml, rl, c, message);        \
    }                                                   \
  } while (false)

#define ASSERT_cc(check, message)                       \
  do {                                                  \
    if (not (check)) {                                  \
      assert_error (#check, ml, rl, c, cc, message);    \
    }                                                   \
  } while (false)



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
  fast_oldboxes.clear();
  swap (fast_boxes, fast_oldboxes);
  
  
  
  boxes.resize (h.mglevels());
  fast_boxes.resize (h.mglevels());
  for (int ml = 0; ml < h.mglevels(); ++ ml) {
    boxes.AT(ml).resize (h.reflevels());
    fast_boxes.AT(ml).resize (h.reflevels());
    for (int rl = 0; rl < h.reflevels(); ++ rl) {
      boxes.AT(ml).AT(rl).resize (h.components(rl));
      fast_boxes.AT(ml).AT(rl).resize (dist::size());
      
      cboxes & level = boxes.AT(ml).AT(rl);
      fast_cboxes & fast_level = fast_boxes.AT(ml).AT(rl);
      
      
      
      // Domain:
      
      ibbox const & domain_exterior = h.baseextent(ml,rl);
      // Variables may have size zero
      // ASSERT_rl (not domain_exterior.empty(),
      //            "The exterior of the domain must not be empty");
      
      i2vect const & boundary_width = h.boundary_width;
      ASSERT_rl (all (all (boundary_width >= 0)),
                 "The gh boundary widths must not be negative");
      
      ibbox const domain_active = domain_exterior.expand (- boundary_width);
      // Variables may have size zero
      // ASSERT_rl (not domain_active.empty(),
      //            "The active part of the domain must not be empty");
      ASSERT_rl (domain_active <= domain_exterior,
                 "The active part of the domain must be contained in the exterior part of the domain");
      
      ibset domain_boundary = domain_exterior - domain_active;
      domain_boundary.normalize();
      
      
      
      for (int c = 0; c < h.components(rl); ++ c) {
        
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        
        
        
        // Interior:
        
        ibbox & intr = box.interior;
        
        // The interior of the grid has the extent as specified by the
        // regridding thorn
        intr = h.extent (ml,rl,c);
        
        // (The interior must not be empty)
        // Variables may have size zero
        // ASSERT_c (not intr.empty(),
        //           "The interior must not be empty");
        
        // The interior must be contained in the domain
        ASSERT_c (intr <= h.baseextent(ml,rl),
                  "The interior must be contained in the domain");
        
        // All interiors must be disjunct
        for (int cc = 0; cc < c; ++ cc) {
          ASSERT_cc (not intr.intersects (level.AT(cc).interior),
                     "All interiors must be disjunct");
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
        
        ASSERT_c (all (all (ghost_width >= 0)),
                  "The gh ghost widths must not be negative");
        extr = intr.expand (i2vect (not is_outer_boundary) * ghost_width);
        
        // (The exterior must not be empty)
        // Variables may have size zero
        // ASSERT_c (not extr.empty(),
        //           "The experior must not be empty");
        
        // The exterior must be contained in the domain
        ASSERT_c (extr <= domain_exterior,
                  "The exterior must be contained in the domain");
        
        
        
        // Cactus ghost zones (which include outer boundaries):
        
        ibset & ghosts = box.ghosts;
        
        ghosts = extr - intr;
        ghosts.normalize();
        
        // The ghosts must be contained in the domain.  Different from
        // the boundaries, the ghost can include part of the outer
        // boundary of the domain.
        ASSERT_c (ghosts <= domain_exterior,
                  "The ghost zones must be contained in the domain");
        
        
        
        // Communicated region:
        
        ibbox & comm = box.communicated;
        
        comm = extr.expand (i2vect (is_outer_boundary) * (- boundary_width));
        
        // (The communicated region must not be empty)
        // Variables may have size zero
        // ASSERT_c (not comm.empty(),
        //           "The communicated region must not be empty");
        
        // The communicated region must be contained in the active
        // part of the domain
        ASSERT_c (comm <= domain_active,
                  "The communicated region must be contained in the active part of the domain");
        
        
        
        // Outer boundary:
        
        ibset & outer_boundaries = box.outer_boundaries;
        
        outer_boundaries = extr - comm;
        outer_boundaries.normalize();
        
        // The outer boundary must be contained in the outer boundary
        // of the domain
        ASSERT_c (outer_boundaries <= domain_boundary,
                  "The outer boundary must be contained in the outer boundary of the domain");
        
        
        
        // Owned region:
        
        ibbox & owned = box.owned;
        
        owned = intr.expand (i2vect (is_outer_boundary) * (- boundary_width));
        
        // (The owned region must not be empty)
        // Variables may have size zero
        // ASSERT_c (not owned.empty(),
        //           "The owned region must not be empty");
        
        // The owned region must be contained in the active part of
        // the domain
        ASSERT_c (owned <= domain_active,
                  "The owned region must be contained in the active part of the domain");
        
        // All owned regions must be disjunct
        for (int cc = 0; cc < c; ++ cc) {
          ASSERT_cc (not owned.intersects (level.AT(cc).owned),
                     "All owned regions must be disjunct");
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
        ASSERT_c (boundaries <= domain_active,
                  "The boundary must be contained in the active part of the domain");
        
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
      ASSERT_rl (allowned <= domain_active,
                 "The owned regions must be contained in the active part of the domain");
      
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
      ASSERT_rl (allbuffers <= domain_active,
                 "The buffer zones must be in the active part of the domain");
      
      
      
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
      ASSERT_rl (allbuffers1 == allbuffers,
                 "Buffer zone consistency check");
      
      
      
      // Test constituency relations:
      
      for (int c = 0; c < h.components(rl); ++ c) {
        dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
        
        ASSERT_c ((box.active & box.buffers).empty(),
                  "Consistency check");
        ASSERT_c ((box.active | box.buffers) == box.owned,
                  "Consistency check");
        
        ASSERT_c ((box.owned & box.boundaries).empty(),
                  "Consistency check");
        ASSERT_c ((box.owned | box.boundaries) == box.communicated,
                  "Consistency check");
        
        ASSERT_c ((box.communicated & box.outer_boundaries).empty(),
                  "Consistency check");
        ASSERT_c ((box.communicated | box.outer_boundaries) == box.exterior,
                  "Consistency check");
        
        ASSERT_c (box.boundaries <= box.ghosts,
                  "Consistency check");
        
        ASSERT_c ((box.interior & box.ghosts).empty(),
                  "Consistency check");
        ASSERT_c ((box.interior | box.ghosts) == box.exterior,
                  "Consistency check");
        
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
          
          for (ibset::const_iterator
                 ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            ibbox const send = recv.expanded_for (obox.interior);
            ASSERT_c (send <= obox.exterior,
                      "Multigrid restriction: Send region must be contained in exterior");
            if (on_this_proc (rl, c)) {
              int const p = dist::rank();
              fast_level.AT(p).fast_mg_rest_sendrecv.push_back
                (sendrecv_pseudoregion_t (send, c, recv, c));
            }
          }
          
          needrecv -= ovlp;
          needrecv.normalize();
          
          // All points must have been received
          ASSERT_c (needrecv.empty(),
                    "Multigrid restriction: All points must have been received");
          
        } // if ml > 0
        
        
        
        // Multigrid prolongation:
        
        if (ml > 0) {
          int const oml = ml - 1;
          
          // Multigrid prolongation must fill all active points
          // (this could probably be relaxed)
          
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
          
          for (ibset::const_iterator
                 ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            ibbox const send =
              recv.expanded_for (box.interior).expand (stencil_size);
            ASSERT_c (send <= box.exterior,
                      "Multigrid prolongation: Send region must be contained in exterior");
            if (on_this_proc (rl, c)) {
              int const p = dist::rank();
              fast_level.AT(p).fast_mg_prol_sendrecv.push_back
                (sendrecv_pseudoregion_t (send, c, recv, c));
            }
          }
          
          oneedrecv -= ovlp;
          oneedrecv.normalize();
          
          // All points must have been received
          ASSERT_c (oneedrecv.empty(),
                    "Multigrid prolongation: All points must have been received");
          
        } // if ml > 0
        
        
        
        // Refinement prolongation:
        
        if (rl > 0) {
          int const orl = rl - 1;
          
          // Refinement prolongation must fill all active points
          
          ibset needrecv = box.active;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size());
          
          ASSERT_c (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                    "Refinement factors must be integer multiples of each other");
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
            
            for (ibset::const_iterator
                   ri =ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              ASSERT_c (send <= obox.exterior,
                        "Refinement prolongation: Send region must be contained in exterior");
              if (on_this_proc (rl, c) or on_this_proc (orl, cc)) {
                int const p = dist::rank();
                fast_level.AT(p).fast_ref_prol_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            
          } // for cc
          
          needrecv.normalize();
          
          // All points must have been received
          ASSERT_c (needrecv.empty(),
                    "Refinement prolongation: All points must have been received");
          
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
            ASSERT_cc (ovlp.empty(),
                       "A region may not synchronise from itself");
          }
          
          for (ibset::const_iterator
                 ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            ibbox const & send = recv;
            if (on_this_proc (rl, c) or on_this_proc (rl, cc)) {
              int const p = dist::rank();
              fast_level.AT(p).fast_sync_sendrecv.push_back
                (sendrecv_pseudoregion_t (send, cc, recv, c));
            }
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
          
          ibset & bndref = box.bndref;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size());
          
          ASSERT_c (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                    "Refinement factors must be integer multiples of each other");
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
            
            for (ibset::const_iterator
                   ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              ASSERT_c (send <= obox.exterior,
                        "Boundary prolongation: Send region must be contained in exterior");
              if (on_this_proc (rl, c) or on_this_proc (orl, cc)) {
                int const p = dist::rank();
                fast_level.AT(p).fast_ref_bnd_prol_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            bndref += ovlp;
              
          } // for cc
          
          needrecv.normalize();
          bndref.normalize();
          
        } // if rl > 0
        
        // All points must now have been received, either through
        // synchronisation or through boundary prolongation
        ASSERT_c (needrecv.empty(),
                  "Synchronisation and boundary prolongation: All points must have been received");
        
      } // for c
      
      
      
      // Refinement restriction:
      
      if (rl > 0) {
        int const orl = rl - 1;
        fast_cboxes & fast_olevel = fast_boxes.AT(ml).AT(orl);
        
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
            
            for (ibset::const_iterator
                   ri =ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const send = recv.expanded_for (box.interior);
              ASSERT_c (send <= box.active,
                        "Refinement restriction: Send region must be contained in active part");
              if (on_this_proc (rl, c) or on_this_proc (orl, cc)) {
                int const p = dist::rank();
                fast_olevel.AT(p).fast_ref_rest_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, c, recv, cc));
              }
            }
            
            needrecv -= ovlp;
              
          } // for c
          
        } // for cc
        
        needrecv.normalize();
        
        // All points must have been received
        ASSERT_rl (needrecv.empty(),
                   "Refinement restriction: All points must have been received");
        
      } // if rl > 0
      
      
      
      // Regridding schedule:
      
      for (int c = 0; c < h.components(rl); ++ c) {
        
        dboxes & box = boxes.AT(ml).AT(rl).AT(c);
        
        ibset needrecv = box.active;
        
        
          
        // Synchronisation:
        
        if (int (oldboxes.size()) > ml and int (oldboxes.AT(ml).size()) > rl) {
          
          int const oldcomponents = oldboxes.AT(ml).AT(rl).size();
          
          // Synchronisation copies from the same level of the old
          // grid structure.  It should fill as many active points as
          // possible
          
          for (int cc = 0; cc < oldcomponents; ++ cc) {
            dboxes const & obox = oldboxes.AT(ml).AT(rl).AT(cc);
            
            ibset ovlp = needrecv & obox.owned;
            ovlp.normalize();
            
            for (ibset::const_iterator
                   ri =ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const & send = recv;
              if (on_this_proc (rl, c) or on_this_oldproc (rl, cc)) {
                int const p = dist::rank();
                fast_level.AT(p).fast_old2new_sync_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            
          } // for cc
          
          needrecv.normalize();
        
        } // if not oldboxes.empty
        
        
        
        // Prolongation:
        
        if (rl > 0) {
          int const orl = rl - 1;
          
          // Prolongation interpolates from the next coarser level of
          // the new grid structure.  It must fill what cannot be
          // synchronised
            
          i2vect const stencil_size = i2vect (prolongation_stencil_size());
          
          ASSERT_c (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                    "Refinement factors must be integer multiples of each other");
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
            
            for (ibset::const_iterator
                   ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              ASSERT_c (send <= obox.exterior,
                        "Regridding prolongation: Send region must be contained in exterior");
              if (on_this_proc (rl, c) or on_this_proc (orl, cc)) {
                int const p = dist::rank();
                fast_level.AT(p).fast_old2new_ref_prol_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            
          } // for cc
          
          needrecv.normalize();
          
        } // if rl > 0
        
        if (int (oldboxes.size()) > ml and int (oldboxes.AT(ml).size()) > 0) {
          // All points must now have been received, either through
          // synchronisation or through prolongation
          ASSERT_c (needrecv.empty(),
                        "Regridding prolongation: All points must have been received");
        }
        
      } // for c
      
    } // for rl
  }   // for m
  
  
  
  // Output:
  if (output_bboxes or there_was_an_error) {
    
    for (int ml = 0; ml < h.mglevels(); ++ ml) {
      for (int rl = 0; rl < h.reflevels(); ++ rl) {
        for (int c = 0; c < h.components(rl); ++ c) {
          dboxes const & box = boxes.AT(ml).AT(rl).AT(c);
          fast_dboxes const & fast_box = fast_boxes.AT(ml).AT(rl).AT(c);
          
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << box;
          cout << fast_box;
          cout << endl;
          
        } // for c
      }   // for rl
    }     // for m
    
  } // if output_bboxes
  
  if (there_was_an_error) {
    CCTK_WARN (CCTK_WARN_ABORT,
               "The grid structure is inconsistent.  "
               "It is impossible to continue.");
  }
  
  
  
  total.stop (0);
}



void
dh::
recompose (int const rl, bool const do_prolongate)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (rl>=0 and rl<h.reflevels());
  
  static Timer timer ("dh::recompose");
  timer.start ();
  
  for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
    (*f)->recompose_crop ();
  }
  
  if (combine_recompose) {
    // Recompose all grid functions of this refinement levels at once.
    // This may be faster, but requires more memory.
    for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
      (*f)->recompose_allocate (rl);
    }
    for (comm_state state; not state.done(); state.step()) {
      for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
        (*f)->recompose_fill (state, rl, do_prolongate);
      }
    }
    for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
      (*f)->recompose_free_old (rl);
    }
  } else {
    // Recompose the grid functions sequentially.  This may be slower,
    // but requires less memory.  This is the default.
    for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
      (*f)->recompose_allocate (rl);
      for (comm_state state; not state.done(); state.step()) {
        (*f)->recompose_fill (state, rl, do_prolongate);
      }
      (*f)->recompose_free_old (rl);
    }
  }
  
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



// Memory usage

size_t
dh::
memory ()
  const
{
  return
    memoryof (ghost_width) +
    memoryof (buffer_width) +
    memoryof (prolongation_order_space) +
    memoryof (boxes) +
    memoryof (fast_boxes) +
    memoryof (fast_oldboxes) +
    memoryof (gfs);
}

size_t
dh::dboxes::
memory ()
  const
{
  return
    memoryof (exterior) +
    memoryof (is_outer_boundary) +
    memoryof (outer_boundaries) +
    memoryof (communicated) +
    memoryof (boundaries) +
    memoryof (owned) +
    memoryof (buffers) +
    memoryof (active) +
    memoryof (sync) +
    memoryof (bndref) +
    memoryof (ghosts) +
    memoryof (interior);
}

size_t
dh::fast_dboxes::
memory ()
  const
{
  return
    memoryof (fast_mg_rest_sendrecv) +
    memoryof (fast_mg_prol_sendrecv) +
    memoryof (fast_ref_prol_sendrecv) +
    memoryof (fast_ref_rest_sendrecv) +
    memoryof (fast_sync_sendrecv) +
    memoryof (fast_ref_bnd_prol_sendrecv) +
    memoryof (fast_old2new_sync_sendrecv) +
    memoryof (fast_old2new_ref_prol_sendrecv);
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
     << "fast_boxes=" << fast_boxes << ","
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
  // Regions:
  os << "dh::dboxes:" << eol;
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
  return os;
}

ostream &
dh::fast_dboxes::
output (ostream & os)
  const
{
  // Communication schedule:
  os << "dh::fast_dboxes:" << eol;
  os << "fast_mg_rest_sendrecv: " << fast_mg_rest_sendrecv << eol;
  os << "fast_mg_prol_sendrecv: " << fast_mg_prol_sendrecv << eol;
  os << "fast_ref_prol_sendrecv: " << fast_ref_prol_sendrecv << eol;
  os << "fast_ref_rest_sendrecv: " << fast_ref_rest_sendrecv << eol;
  os << "fast_sync_sendrecv: " << fast_sync_sendrecv << eol;
  os << "fast_ref_bnd_prol_sendrecv: " << fast_ref_bnd_prol_sendrecv << eol;
  os << "fast_old2new_sync_sendrecv:" << fast_old2new_sync_sendrecv << eol;
  os << "fast_old2new_ref_prol_sendrecv:" << fast_old2new_ref_prol_sendrecv << eol;
  return os;
}
