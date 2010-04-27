#include <cassert>
#include <cstddef>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "CarpetTimers.hh"

#include "mpi_string.hh"
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



list<dh*> dh::alldh;



// Constructors
dh::
dh (gh & h_,
    vector<i2vect> const & ghost_widths_, vector<i2vect> const & buffer_widths_,
    vector<int> const & prolongation_orders_space_)
  : h(h_),
    ghost_widths(ghost_widths_), buffer_widths(buffer_widths_),
    prolongation_orders_space(prolongation_orders_space_)
{
  size_t const maxreflevels = h.reffacts.size();
  assert (ghost_widths.size() >= maxreflevels);
  assert (buffer_widths.size() >= maxreflevels);
  assert (prolongation_orders_space.size() >= maxreflevels);
  for (size_t rl=0; rl<maxreflevels; ++rl) {
    assert (all (all (ghost_widths.AT(rl) >= 0)));
    assert (all (all (buffer_widths.AT(rl) >= 0)));
    assert (prolongation_orders_space.AT(rl) >= 0);
  }
  
  alldhi = alldh.insert(alldh.end(), this);
  gh_handle = h.add (this);
  CHECKPOINT;
  regrid (false);
  for (int rl = 0; rl < h.reflevels(); ++ rl) {
    recompose (rl, false);
  }
  regrid_free (false);
}



// Destructors
dh::
~dh ()
{
  CHECKPOINT;
  h.erase (gh_handle);
  alldh.erase(alldhi);
}



// Helpers
int
dh::
prolongation_stencil_size (int const rl)
  const
{
  assert (prolongation_orders_space.AT(rl) >= 0);
  return prolongation_orders_space.AT(rl) / 2;
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
              char const * restrict const file, int const line,
              int const ml, int const rl,
              char const * restrict const message)
{
  CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "\n%s:%d:\n   [ml=%d rl=%d] The following grid structure consistency check failed:\n   %s\n   %s",
              file, line, ml, rl, message, checkstring);
  there_was_an_error = true;
}

static
void
assert_error (char const * restrict const checkstring,
              char const * restrict const file, int const line,
              int const ml, int const rl, int const c,
              char const * restrict const message)
{
  CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "\n%s:%d:\n   [ml=%d rl=%d c=%d] The following grid structure consistency check failed:\n   %s\n   %s",
              file, line, ml, rl, c, message, checkstring);
  there_was_an_error = true;
}

static
void
assert_error (char const * restrict const checkstring,
              char const * restrict const file, int const line,
              int const ml, int const rl, int const c, int const cc,
              char const * restrict const message)
{
  CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "\n%s:%d:\n   [ml=%d rl=%d c=%d cc=%d] The following grid structure consistency check failed:\n   %s\n   %s",
              file, line, ml, rl, c, cc, message, checkstring);
  there_was_an_error = true;
}

#ifdef CARPET_OPTIMISE

// For highest efficiency, omit all self-checks
#define ASSERT_rl(check, message)
#define ASSERT_c(check, message)
#define ASSERT_cc(check, message)

#else

#define ASSERT_rl(check, message)                                       \
  do {                                                                  \
    if (not (check)) {                                                  \
      assert_error (#check, __FILE__, __LINE__, ml, rl, message);       \
    }                                                                   \
  } while (false)

#define ASSERT_c(check, message)                                        \
  do {                                                                  \
    if (not (check)) {                                                  \
      assert_error (#check, __FILE__, __LINE__, ml, rl, c, message);    \
    }                                                                   \
  } while (false)

#define ASSERT_cc(check, message)                                       \
  do {                                                                  \
    if (not (check)) {                                                  \
      assert_error (#check, __FILE__, __LINE__, ml, rl, c, cc, message); \
    }                                                                   \
  } while (false)

#endif



void
dh::
regrid (bool const do_init)
{
  DECLARE_CCTK_PARAMETERS;
    
  static Carpet::Timer timer ("CarpetLib::dh::regrid");
  timer.start();
  
  CHECKPOINT;
  
  static Timer total ("CarpetLib::dh::regrid");
  total.start ();
  
  light_mboxes old_light_boxes;
  swap (light_boxes, old_light_boxes);
  
  full_mboxes full_boxes;
  
  fast_boxes.clear();
  
  
  
  light_boxes.resize (h.mglevels());
  local_boxes.resize (h.mglevels());
  full_boxes.resize (h.mglevels());
  fast_boxes.resize (h.mglevels());
  for (int ml = 0; ml < h.mglevels(); ++ ml) {
    light_boxes.AT(ml).resize (h.reflevels());
    local_boxes.AT(ml).resize (h.reflevels());
    full_boxes.AT(ml).resize (h.reflevels());
    fast_boxes.AT(ml).resize (h.reflevels());
    for (int rl = 0; rl < h.reflevels(); ++ rl) {
      light_boxes.AT(ml).AT(rl).resize (h.components(rl));
      local_boxes.AT(ml).AT(rl).resize (h.local_components(rl));
      full_boxes.AT(ml).AT(rl).resize (h.components(rl));
      
      light_cboxes & light_level = light_boxes.AT(ml).AT(rl);
      local_cboxes & local_level = local_boxes.AT(ml).AT(rl);
      full_cboxes & full_level = full_boxes.AT(ml).AT(rl);
      fast_dboxes & fast_level = fast_boxes.AT(ml).AT(rl);
      
      vector<fast_dboxes> fast_level_otherprocs (dist::size());
      
      
      
      i2vect const& ghost_width = ghost_widths.AT(rl);
      i2vect const& buffer_width = buffer_widths.AT(rl);
      
      
      
      // Domain:
      
      static Carpet::Timer timer_domain ("CarpetLib::dh::regrid::domain");
      timer_domain.start();
      
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
      
      ibset const domain_boundary = domain_exterior - domain_active;
      
      timer_domain.stop();
      
      
      
      static Carpet::Timer timer_region ("CarpetLib::dh::regrid::region");
      timer_region.start();
      
      for (int c = 0; c < h.components(rl); ++ c) {
        
        full_dboxes & box = full_level.AT(c);
        
        
        
        // Interior:
        
        ibbox & intr = box.interior;
        intr = ibbox::poison();
        
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
#ifdef CARPET_DEBUG
        for (int cc = 0; cc < c; ++ cc) {
          ASSERT_cc (not intr.intersects (full_level.AT(cc).interior),
                     "All interiors must be disjunct");
        }
#endif
        
        
        
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
        extr = ibbox::poison();
        
        ASSERT_c (all (all (ghost_width >= 0)),
                  "The gh ghost widths must not be negative");
        extr = intr.expand (i2vect (not is_outer_boundary) * ghost_width);
        
        // (The exterior must not be empty)
        // Variables may have size zero
        // ASSERT_c (not extr.empty(),
        //           "The exterior must not be empty");
        
        // The exterior must be contained in the domain
        ASSERT_c (extr <= domain_exterior,
                  "The exterior must be contained in the domain");
        
        
        
        // Cactus ghost zones (which include outer boundaries):
        
        ibset & ghosts = box.ghosts;
        ghosts = ibset::poison();
        
        ghosts = extr - intr;
        
        // The ghosts must be contained in the domain.  Different from
        // the boundaries, the ghost can include part of the outer
        // boundary of the domain.
        ASSERT_c (ghosts <= domain_exterior,
                  "The ghost zones must be contained in the domain");
        
        
        
        // Communicated region:
        
        ibbox & comm = box.communicated;
        comm = ibbox::poison();
        
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
        outer_boundaries = ibset::poison();
        
        outer_boundaries = extr - comm;
        
        // The outer boundary must be contained in the outer boundary
        // of the domain
        ASSERT_c (outer_boundaries <= domain_boundary,
                  "The outer boundary must be contained in the outer boundary of the domain");
        
        
        
        // Owned region:
        
        ibbox & owned = box.owned;
        owned = ibbox::poison();
        
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
#ifdef CARPET_DEBUG
        for (int cc = 0; cc < c; ++ cc) {
          ASSERT_cc (not owned.intersects (full_level.AT(cc).owned),
                     "All owned regions must be disjunct");
        }
#endif
        
        
        
        // Boundary (Carpet ghost zones, which do not include outer
        // boundaries):
        
        ibset & boundaries = box.boundaries;
        boundaries = ibset::poison();
        
        boundaries = comm - owned;
        
        // The boundary must be contained in the active part of the
        // domain.  This prevents that a region is too close to the
        // outer boundary, so that it has ghost zones overlapping with
        // the outer boundary.
        ASSERT_c (boundaries <= domain_active,
                  "The boundary must be contained in the active part of the domain");
        
      } // for c
      
      timer_region.stop();
      
      
      
      // Conjunction of all buffer zones:
      
      static Carpet::Timer timer_buffers ("CarpetLib::dh::regrid::buffers");
      timer_buffers.start();
      
      // Enlarge active part of domain
      i2vect const safedist = i2vect (0);
      ibbox const domain_enlarged = domain_active.expand (safedist);
      
      // All owned regions
      ibset const allowned (full_level, & full_dboxes::owned);
      ASSERT_rl (allowned <= domain_active,
                 "The owned regions must be contained in the active part of the domain");
      
      // All not-owned regions
      ibset const notowned = domain_enlarged - allowned;
      
      // All not-active points
      ibset const notactive = notowned.expand (buffer_width);
      
      // All not-active points, in stages
      int const num_substeps =
        any (any (ghost_width == 0)) ?
        0 :
        minval (minval (buffer_width / ghost_width));
      if (not all (all (buffer_width == num_substeps * ghost_width))) {
        ostringstream buf;
        buf << "The buffer width " << buffer_width << " is not a multiple of the ghost width " << ghost_width << " on level " << rl;
        CCTK_WARN (CCTK_WARN_COMPLAIN, buf.str().c_str());
      }
      
      vector<ibset> notactive_stepped (num_substeps+1);
      notactive_stepped.AT(0) = notowned;
      for (int substep = 1; substep <= num_substeps; ++ substep) {
        notactive_stepped.AT(substep) =
          notactive_stepped.AT(substep-1).expand (ghost_width);
      }
      if (all (all (buffer_width == num_substeps * ghost_width))) {
        ASSERT_rl (notactive_stepped.AT(num_substeps) == notactive,
                   "The stepped not-active region must be equal to the not-active region");
      }
      
      // All buffer zones
      ibset const allbuffers = allowned & notactive;
      
      // All active points
      ibset const allactive = allowned - notactive;
      
      // All stepped buffer zones
      vector<ibset> allbuffers_stepped (num_substeps);
      ibset allbuffers_stepped_combined;
      for (int substep = 0; substep < num_substeps; ++ substep) {
        allbuffers_stepped.AT(substep) = 
          allowned &
          (notactive_stepped.AT(substep+1) - notactive_stepped.AT(substep));
        allbuffers_stepped_combined += allbuffers_stepped.AT(substep);
      }
      if (all (all (buffer_width == num_substeps * ghost_width))) {
        ASSERT_rl (allbuffers_stepped_combined == allbuffers,
                   "The stepped buffer zones must be equal to the buffer zones");
      }
      
      // Buffer zones must be in the active part of the domain
      ASSERT_rl (allactive <= domain_active,
                 "The active region must be in the active part of the domain");
      ASSERT_rl (allbuffers <= domain_active,
                 "The buffer zones must be in the active part of the domain");
      ASSERT_rl ((allactive & allbuffers).empty(), 
                 "The active points and the buffer zones cannot overlap");
      ASSERT_rl (allactive + allbuffers == allowned,
                 "The active points and the buffer points together must be exactly the owned region");
      
      
      
      for (int c = 0; c < h.components(rl); ++ c) {
        full_dboxes & box = full_level.AT(c);
        
        // Buffer zones:
        box.buffers = box.owned & allbuffers;
        
        // Active region:
        box.active = box.owned  & allactive;
        ASSERT_c (box.active == box.owned - box.buffers,
                  "The active region must equal the owned region minus the buffer zones");
      } // for c
      
      for (int lc = 0; lc < h.local_components(rl); ++ lc) {
        int const c = h.get_component (rl, lc);
        local_dboxes & local_box = local_level.AT(lc);
        full_dboxes const& box = full_level.AT(c);
        
        // Stepped buffer zones:
        local_box.buffers = box.buffers;
        
        local_box.buffers_stepped.resize (num_substeps);
        for (int substep = 0; substep < num_substeps; ++ substep) {
          local_box.buffers_stepped.AT(substep) =
            box.owned & allbuffers_stepped.AT(substep);
        }
        
        local_box.active = box.active;
      } // for lc
      
      
      
      // The conjunction of all buffer zones must equal allbuffers
      ibset const allbuffers1 (full_level, & full_dboxes::buffers);
      ASSERT_rl (allbuffers1 == allbuffers,
                 "Buffer zone consistency check");
      
      timer_buffers.stop();
      
      
      
      // Test constituency relations:
      
      static Carpet::Timer timer_test ("CarpetLib::dh::regrid::test");
      timer_test.start();
      
      for (int c = 0; c < h.components(rl); ++ c) {
        full_dboxes const & box = full_level.AT(c);
        
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
      
      timer_test.stop();
       
      
      
      // Communication schedule:
      
      static Carpet::Timer timer_comm ("CarpetLib::dh::regrid::comm");
      timer_comm.start();
      
      for (int lc = 0; lc < h.local_components(rl); ++ lc) {
        int const c = h.get_component (rl, lc);
        
        full_dboxes & box = full_level.AT(c);
        
        
        
        // Multigrid restriction:
        
        static Carpet::Timer timer_comm_mgrest
          ("CarpetLib::dh::regrid::comm::mgrest");
        timer_comm_mgrest.start();
        
        if (ml > 0) {
          int const oml = ml - 1;
          
          // Multigrid restriction must fill all active points
          
          full_dboxes const & obox = full_boxes.AT(oml).AT(rl).AT(c);
          
          ibset needrecv = box.active;
          
          ibset const contracted_oactive
            (obox.active.contracted_for (box.interior));
          ibset const ovlp = needrecv & contracted_oactive;
          
          for (ibset::const_iterator
                 ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            ibbox const send = recv.expanded_for (obox.interior);
            ASSERT_c (send <= obox.exterior,
                      "Multigrid restriction: Send region must be contained in exterior");
            fast_level.fast_mg_rest_sendrecv.push_back
              (sendrecv_pseudoregion_t (send, c, recv, c));
          }
          
          needrecv -= ovlp;
          
          // All points must have been received
          ASSERT_c (needrecv.empty(),
                    "Multigrid restriction: All points must have been received");
          
        } // if ml > 0
        
        timer_comm_mgrest.stop();
        
        
        
        // Multigrid prolongation:
        
        static Carpet::Timer timer_comm_mgprol
          ("CarpetLib::dh::regrid::comm::mprol");
        timer_comm_mgprol.start();
        
        if (ml > 0) {
          int const oml = ml - 1;
          
          // Multigrid prolongation must fill all active points
          // (this could probably be relaxed)
          
          full_dboxes const & obox = full_boxes.AT(oml).AT(rl).AT(c);
          
          ibset oneedrecv = obox.active;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size(rl));
          
          ibset const expanded_active (box.active.expanded_for (obox.interior));
          ibset const ovlp = oneedrecv & expanded_active;
          
          for (ibset::const_iterator
                 ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
          {
            ibbox const & recv = * ri;
            ibbox const send =
              recv.expanded_for (box.interior).expand (stencil_size);
            ASSERT_c (send <= box.exterior,
                      "Multigrid prolongation: Send region must be contained in exterior");
            fast_level.fast_mg_prol_sendrecv.push_back
              (sendrecv_pseudoregion_t (send, c, recv, c));
          }
          
          oneedrecv -= ovlp;
          
          // All points must have been received
          ASSERT_c (oneedrecv.empty(),
                    "Multigrid prolongation: All points must have been received");
          
        } // if ml > 0
        
        timer_comm_mgprol.stop();
        
        
        
        // Refinement prolongation:
        
        static Carpet::Timer timer_comm_refprol
          ("CarpetLib::dh::regrid::comm::refprol");
        timer_comm_refprol.start();
        
        if (rl > 0) {
          int const orl = rl - 1;
          
          // Refinement prolongation must fill all active points
          
          ibset needrecv = box.active;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size(rl));
          
          ASSERT_c (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                    "Refinement factors must be integer multiples of each other");
          i2vect const reffact =
            i2vect (h.reffacts.at(rl) / h.reffacts.at(orl));
          
          for (int cc = 0; cc < h.components(orl); ++ cc) {
            full_dboxes const & obox = full_boxes.AT(ml).AT(orl).AT(cc);
            
#if 0
            // untested for cell centering
            ibset const expanded_oactive
              (obox.active.contracted_for (box.interior).expand (reffact));
#else
            ibset const expanded_oactive
              (obox.active.expanded_for (box.interior).expand (reffact));
#endif
            ibset const ovlp = needrecv & expanded_oactive;
            
            for (ibset::const_iterator
                   ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              ASSERT_c (send <= obox.exterior,
                        "Refinement prolongation: Send region must be contained in exterior");
              fast_level.fast_ref_prol_sendrecv.push_back
                (sendrecv_pseudoregion_t (send, cc, recv, c));
              if (not on_this_proc (orl, cc)) {
                fast_dboxes & fast_level_otherproc =
                  fast_level_otherprocs.AT(this_proc(orl, cc));
                fast_level_otherproc.fast_ref_prol_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            
          } // for cc
          
          // All points must have been received
          ASSERT_c (needrecv.empty(),
                    "Refinement prolongation: All points must have been received");
          
        } // if rl > 0
        
        timer_comm_refprol.stop();
        
        
        
        // Synchronisation:
        
        static Carpet::Timer timer_comm_sync
          ("CarpetLib::dh::regrid::comm::sync");
        timer_comm_sync.start();
        
        {
          
          // Synchronisation should fill as many boundary points as
          // possible
          
#if 0
          // Outer boundaries are not synchronised, since they cannot
          // be filled by boundary prolongation either, and therefore
          // the user code must set them anyway.
          ibset needrecv = box.boundaries;
#else
          // Outer boundaries are synchronised for backward
          // compatibility.
          ibset needrecv = box.ghosts;
#endif
          ibset const needrecv_orig = needrecv;
          
          ibset & sync = box.sync;
          
          for (int cc = 0; cc < h.components(rl); ++ cc) {
            full_dboxes const & obox = full_level.AT(cc);
            
#if 0
            ibset ovlp = needrecv & obox.owned;
#else
            ibset ovlp = needrecv & obox.interior;
#endif
            
            if (cc == c) {
              ASSERT_cc (ovlp.empty(),
                         "A region may not synchronise from itself");
            }
            
            for (ibset::const_iterator
                   ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const & send = recv;
              fast_level.fast_sync_sendrecv.push_back
                (sendrecv_pseudoregion_t (send, cc, recv, c));
              if (not on_this_proc (rl, cc)) {
                fast_dboxes & fast_level_otherproc =
                  fast_level_otherprocs.AT(this_proc(rl, cc));
                fast_level_otherproc.fast_sync_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            sync += ovlp;
            
          } // for cc
          
          
        }
        
        timer_comm_sync.stop();
        
        
        
        // Boundary prolongation:
        
        static Carpet::Timer timer_comm_refbndprol
          ("CarpetLib::dh::regrid::comm::refbndprol");
        timer_comm_refbndprol.start();
        
        if (rl > 0) {
          int const orl = rl - 1;
          
#if 0
          // Outer boundaries are not synchronised, since they cannot
          // be filled by boundary prolongation either, and therefore
          // the user code must set them anyway.
          ibset needrecv = box.boundaries;
#else
          // Outer boundaries are synchronised for backward
          // compatibility.
          ibset needrecv = box.ghosts;
#endif
          
          // Points which are synchronised need not be boundary
          // prolongated
          needrecv -= box.sync;
          
          // Outer boundary points cannot be boundary prolongated
          needrecv &= box.communicated;
          
          // Prolongation must fill what cannot be synchronised, and
          // also all buffer zones
          needrecv += box.buffers;
          
          ibset & bndref = box.bndref;
          
          i2vect const stencil_size = i2vect (prolongation_stencil_size(rl));
          
          ASSERT_c (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                    "Refinement factors must be integer multiples of each other");
          i2vect const reffact =
            i2vect (h.reffacts.at(rl) / h.reffacts.at(orl));
          ivect const reffact1 = h.reffacts.at(rl) / h.reffacts.at(orl);
          
          for (int cc = 0; cc < h.components(orl); ++ cc) {
            full_dboxes const & obox = full_boxes.AT(ml).AT(orl).AT(cc);
            
#if 0
            // untested for cell centering
            ibset const expanded_oactive
              (obox.active.contracted_for (box.interior).expand (reffact));
#else
            ibset const expanded_oactive
              (obox.active.expanded_for (box.interior).expand (reffact));
#endif
            ibset const ovlp = needrecv & expanded_oactive;
            
            for (ibset::const_iterator
                   ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
            {
              ibbox const & recv = * ri;
              ibbox const send =
                recv.expanded_for (obox.interior).expand (stencil_size);
              ASSERT_c (send <= obox.exterior,
                        "Boundary prolongation: Send region must be contained in exterior");
              fast_level.fast_ref_bnd_prol_sendrecv.push_back
                (sendrecv_pseudoregion_t (send, cc, recv, c));
              if (not on_this_proc (orl, cc)) {
                fast_dboxes & fast_level_otherproc =
                  fast_level_otherprocs.AT(this_proc(orl, cc));
                fast_level_otherproc.fast_ref_bnd_prol_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
              }
            }
            
            needrecv -= ovlp;
            bndref += ovlp;
            
          } // for cc
          
          // All points must now have been received, either through
          // synchronisation or through boundary prolongation
          ASSERT_c (needrecv.empty(),
                    "Synchronisation and boundary prolongation: All points must have been received");
          
        } // if rl > 0
        
        timer_comm_refbndprol.stop();
        
      } // for lc
      
      
      
      // Refinement restriction:
      
      static Carpet::Timer timer_comm_refrest
        ("CarpetLib::dh::regrid::comm::refrest");
      timer_comm_refrest.start();
      
      if (rl > 0) {
        int const orl = rl - 1;
        fast_dboxes & fast_olevel = fast_boxes.AT(ml).AT(orl);
        
        if (h.components(orl) > 0) {
          for (int lc = 0; lc < h.local_components(rl); ++ lc) {
            int const c = h.get_component (rl, lc);
            
            full_dboxes const & box = full_level.AT(c);
            full_dboxes const & obox0 = full_boxes.AT(ml).AT(orl).AT(0);
            
            // Refinement restriction may fill all active points, and
            // must use all active points
            
            ibset needrecv (box.active.contracted_for (obox0.interior));
            
            for (int cc = 0; cc < h.components(orl); ++ cc) {
              full_dboxes & obox = full_boxes.AT(ml).AT(orl).AT(cc);
              
              ibset const contracted_active
                (box.active.contracted_for (obox0.interior));
              ibset const ovlp = obox.active & contracted_active;
              
              for (ibset::const_iterator
                     ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
              {
                ibbox const & recv = * ri;
                ibbox const send = recv.expanded_for (box.interior);
                ASSERT_c (send <= box.active,
                          "Refinement restriction: Send region must be contained in active part");
                fast_olevel.fast_ref_rest_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, c, recv, cc));
                if (not on_this_proc (orl, cc)) {
                  fast_dboxes & fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(orl, cc));
                  fast_level_otherproc.fast_ref_rest_sendrecv.push_back
                    (sendrecv_pseudoregion_t (send, c, recv, cc));
                }
              }
              
              needrecv -= ovlp;
              
            } // for cc
            
            // All points must have been received
            ASSERT_rl (needrecv.empty(),
                       "Refinement restriction: All points must have been received");
            
          } // for lc
        }   // if orl not empty
        
      } // if rl > 0
      
      timer_comm_refrest.stop();
        
      timer_comm.stop();
      
      
      
      // Regridding schedule:
      
      fast_level.do_init = do_init;
      if (do_init) {
        
        static Carpet::Timer timer_regrid ("CarpetLib::dh::regrid::regrid");
        timer_regrid.start();
        
        for (int lc = 0; lc < h.local_components(rl); ++ lc) {
          int const c = h.get_component (rl, lc);
          
          full_dboxes & box = full_level.AT(c);
          
          ibset needrecv = box.active;
          
          
          
          // Regridding synchronisation:
          
          static Carpet::Timer timer_regrid_sync
            ("CarpetLib::dh::regrid::regrid::sync");
          timer_regrid_sync.start();
          
          if (int (old_light_boxes.size()) > ml and
              int (old_light_boxes.AT(ml).size()) > rl)
          {
            
            int const oldcomponents = old_light_boxes.AT(ml).AT(rl).size();
            
            // Synchronisation copies from the same level of the old
            // grid structure.  It should fill as many active points
            // as possible.
            
            for (int cc = 0; cc < oldcomponents; ++ cc) {
              light_dboxes const & obox = old_light_boxes.AT(ml).AT(rl).AT(cc);
              
              ibset const ovlp = needrecv & obox.owned;
              
              for (ibset::const_iterator
                     ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
              {
                ibbox const & recv = * ri;
                ibbox const & send = recv;
                fast_level.fast_old2new_sync_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
                if (not on_this_oldproc (rl, cc)) {
                  fast_dboxes & fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(rl, cc));
                  fast_level_otherproc.fast_old2new_sync_sendrecv.push_back
                    (sendrecv_pseudoregion_t (send, cc, recv, c));
                }
              }
              
              needrecv -= ovlp;
              
            } // for cc
            
          } // if not old_light_boxes.empty
          
          timer_regrid_sync.stop();
          
          
          
          // Regridding prolongation:
          
          static Carpet::Timer timer_regrid_prolongate
            ("CarpetLib::dh::regrid::regrid::prolongate");
          timer_regrid_prolongate.start();
          
          if (rl > 0) {
            int const orl = rl - 1;
            
            // Prolongation interpolates from the next coarser level
            // of the new grid structure.  It must fill what cannot be
            // synchronised.
            
            i2vect const stencil_size = i2vect (prolongation_stencil_size(rl));
            
            ASSERT_c (all (h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                      "Refinement factors must be integer multiples of each other");
            i2vect const reffact =
              i2vect (h.reffacts.at(rl) / h.reffacts.at(orl));
            
            for (int cc = 0; cc < h.components(orl); ++ cc) {
              full_dboxes const & obox = full_boxes.AT(ml).AT(orl).AT(cc);
              
#if 0
              // untested for cell centering
              ibset const expanded_oactive
                (obox.active.contracted_for (box.interior).expand (reffact));
#else
              ibset const expanded_oactive
                (obox.active.expanded_for (box.interior).expand (reffact));
#endif
              ibset const ovlp = needrecv & expanded_oactive;
              
              for (ibset::const_iterator
                     ri = ovlp.begin(); ri != ovlp.end(); ++ ri)
              {
                ibbox const & recv = * ri;
                ibbox const send =
                  recv.expanded_for (obox.interior).expand (stencil_size);
                ASSERT_c (send <= obox.exterior,
                          "Regridding prolongation: Send region must be contained in exterior");
                fast_level.fast_old2new_ref_prol_sendrecv.push_back
                  (sendrecv_pseudoregion_t (send, cc, recv, c));
                if (not on_this_proc (orl, cc)) {
                  fast_dboxes & fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(orl, cc));
                  fast_level_otherproc.fast_old2new_ref_prol_sendrecv.
                    push_back (sendrecv_pseudoregion_t (send, cc, recv, c));
                }
              }
              
              needrecv -= ovlp;
              
            } // for cc
            
          } // if rl > 0
          
          if (int (old_light_boxes.size()) > ml and
              int (old_light_boxes.AT(ml).size()) > 0)
          {
            // All points must now have been received, either through
            // synchronisation or through prolongation
            ASSERT_c (needrecv.empty(),
                      "Regridding prolongation: All points must have been received");
          }
          
          timer_regrid_prolongate.stop();
          
        } // for lc
        
        timer_regrid.stop();
        
      } // if do_init
      
      
      
      for (int lc = 0; lc < h.local_components(rl); ++ lc) {
        int const c = h.get_component (rl, lc);
        
        light_level.AT(c).exterior = full_level.AT(c).exterior;
        light_level.AT(c).owned    = full_level.AT(c).owned;
        light_level.AT(c).interior = full_level.AT(c).interior;
#if 0
        dh::dboxes::ibset2ibboxs (full_level.AT(c).active,
                                  light_level.AT(c).active,
                                  light_level.AT(c).numactive);
#endif
        
        light_level.AT(c).exterior_size = full_level.AT(c).exterior.size();
        light_level.AT(c).owned_size    = full_level.AT(c).owned.size();
        light_level.AT(c).active_size   = full_level.AT(c).active.size();
        
      } // for lc
      
      
      
      // Broadcast grid structure and communication schedule
      
      {
        
        static Carpet::Timer timer_bcast_boxes
          ("CarpetLib::dh::regrid::bcast_boxes");
        timer_bcast_boxes.start();
        
        int const count_send = h.local_components(rl);
        vector<light_dboxes> level_send (count_send);
        for (int lc = 0; lc < h.local_components(rl); ++ lc) {
          int const c = h.get_component (rl, lc);
          level_send.AT(lc) = light_level.AT(c);
        }
        vector<vector<light_dboxes> > const level_recv =
          allgatherv (dist::comm(), level_send);
        vector<int> count_recv (dist::size(), 0);
        for (int c = 0; c < h.components(rl); ++ c) {
          int const p = this_proc (rl, c);
          if (p != dist::rank()) {
            light_level.AT(c) = level_recv.AT(p).AT(count_recv.AT(p));
            ++ count_recv.AT(p);
          }
        }
        for (int p = 0; p < dist::size(); ++ p) {
          if (p != dist::rank()) {
            assert (count_recv.AT(p) == int(level_recv.AT(p).size()));
          }
        }
        
        timer_bcast_boxes.stop();
        
      }
      
      {
        
        static Carpet::Timer timer_bcast_comm
          ("CarpetLib::dh::regrid::bcast_comm");
        timer_bcast_comm.start();
        
        static Carpet::Timer timer_bcast_comm_ref_prol
          ("CarpetLib::dh::regrid::bcast_comm::ref_prol");
        timer_bcast_comm_ref_prol.start();
        broadcast_schedule (fast_level_otherprocs, fast_level,
                            & fast_dboxes::fast_ref_prol_sendrecv);
        timer_bcast_comm_ref_prol.stop();
        
        static Carpet::Timer timer_bcast_comm_sync
          ("CarpetLib::dh::regrid::bcast_comm::sync");
        timer_bcast_comm_sync.start();
        broadcast_schedule (fast_level_otherprocs, fast_level,
                            & fast_dboxes::fast_sync_sendrecv);
        timer_bcast_comm_sync.stop();
        
        static Carpet::Timer timer_bcast_comm_ref_bnd_prol
          ("CarpetLib::dh::regrid::bcast_comm::ref_bnd_prol");
        timer_bcast_comm_ref_bnd_prol.start();
        broadcast_schedule (fast_level_otherprocs, fast_level,
                            & fast_dboxes::fast_ref_bnd_prol_sendrecv);
        timer_bcast_comm_ref_bnd_prol.stop();
        
        if (rl > 0) {
          int const orl = rl - 1;
          fast_dboxes & fast_olevel = fast_boxes.AT(ml).AT(orl);
          static Carpet::Timer timer_bcast_comm_ref_rest
            ("CarpetLib::dh::regrid::bcast_comm::ref_rest");
          timer_bcast_comm_ref_rest.start();
          broadcast_schedule (fast_level_otherprocs, fast_olevel,
                              & fast_dboxes::fast_ref_rest_sendrecv);
        timer_bcast_comm_ref_rest.stop();
        }
        
        // TODO: Maybe broadcast old2new schedule only if do_init is
        // set
        static Carpet::Timer timer_bcast_comm_old2new_sync
          ("CarpetLib::dh::regrid::bcast_comm::old2new_sync");
        timer_bcast_comm_old2new_sync.start();
        broadcast_schedule (fast_level_otherprocs, fast_level,
                            & fast_dboxes::fast_old2new_sync_sendrecv);
        timer_bcast_comm_old2new_sync.stop();
        
        static Carpet::Timer timer_bcast_comm_old2new_ref_prol
          ("CarpetLib::dh::regrid::bcast_comm::old2new_ref_prol");
        timer_bcast_comm_old2new_ref_prol.start();
        broadcast_schedule (fast_level_otherprocs, fast_level,
                            & fast_dboxes::fast_old2new_ref_prol_sendrecv);
        timer_bcast_comm_old2new_ref_prol.stop();
        
        timer_bcast_comm.stop();
        
      }
      
      
      
      // Output:
      if (output_bboxes or there_was_an_error) {
        
        cout << eol;
        cout << "ml=" << ml << " rl=" << rl << eol;
        cout << "baseextent=" << h.baseextent(ml,rl) << eol;
        
        for (int c = 0; c < h.components(rl); ++ c) {
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << "extent=" << h.extent(ml,rl,c) << eol;
          cout << "outer_boundaries=" << h.outer_boundaries(rl,c) << eol;
          cout << "processor=" << h.outer_boundaries(rl,c) << eol;
        } // for c
        
        for (int c = 0; c < h.components(rl); ++ c) {
          full_dboxes const & box = full_boxes.AT(ml).AT(rl).AT(c);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << box;
        } // for c
        
        for (int c = 0; c < h.components(rl); ++ c) {
          light_dboxes const & box = light_boxes.AT(ml).AT(rl).AT(c);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << box;
        } // for c
        
        for (int lc = 0; lc < h.local_components(rl); ++ lc) {
          int const c = h.get_component (rl, lc);
          local_dboxes const & box = local_boxes.AT(ml).AT(rl).AT(lc);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " lc=" << lc << " c=" << c << eol;
          cout << box;
        } // for lc
        
        fast_dboxes const & fast_box = fast_boxes.AT(ml).AT(rl);
        cout << eol;
        cout << "ml=" << ml << " rl=" << rl << eol;
        cout << fast_box;
        
      } // if output_bboxes
      
      
      
      // Free memory early to save space
      if (int (old_light_boxes.size()) > ml and
          int (old_light_boxes.AT(ml).size()) > rl)
      {
        old_light_boxes.AT(ml).AT(rl).clear();
      }
      
      if (ml > 0) {
        if (rl > 0) {
          full_boxes.AT(ml-1).AT(rl-1).clear();
        }
        if (rl == h.reflevels()-1) {
          full_boxes.AT(ml-1).AT(rl).clear();
        }
      }
      if (ml == h.mglevels()-1) {
        if (rl > 0) {
          full_boxes.AT(ml).AT(rl-1).clear();
        }
        if (rl == h.reflevels()-1) {
          full_boxes.AT(ml).AT(rl).clear();
        }
      }
      
    } // for rl
    
    if (ml > 0) {
      full_boxes.AT(ml-1).clear();
    }
    if (ml == h.mglevels()-1) {
      full_boxes.AT(ml).clear();
    }
    
  } // for ml
  
  
  
  // Output:
  if (output_bboxes or there_was_an_error) {
    
    cout << eol;
    cout << "memoryof(gh)=" << memoryof(h) << eol;
    cout << "memoryof(dh)=" << memoryof(*this) << eol;
    cout << "memoryof(dh.light_boxes)=" << memoryof(light_boxes) << eol;
    cout << "memoryof(dh.local_boxes)=" << memoryof(local_boxes) << eol;
    cout << "memoryof(dh.fast_boxes)=" << memoryof(fast_boxes) << eol;
    int gfcount = 0;
    size_t gfmemory = 0;
    for (list<ggf*>::const_iterator
           gfi = gfs.begin(); gfi != gfs.end(); ++ gfi)
    {
      ++ gfcount;
      gfmemory += memoryof(**gfi);
    }
    cout << "#gfs=" << gfcount << eol;
    cout << "memoryof(gfs)=" << gfmemory << eol;
    
  } // if output_bboxes
  
  if (there_was_an_error) {
    CCTK_WARN (CCTK_WARN_ABORT,
               "The grid structure is inconsistent.  It is impossible to continue.");
  }
  
  
  
  total.stop (0);
  timer.stop();
}



void
dh::
broadcast_schedule (vector<fast_dboxes> & fast_level_otherprocs,
                    fast_dboxes & fast_level,
                    srpvect fast_dboxes::* const schedule_item)
{
  static Carpet::Timer timer_bs1 ("CarpetLib::dh::bs1");
  timer_bs1.start();
  vector <srpvect> send (dist::size());
  for (int p=0; p<dist::size(); ++p) {
    swap (send.AT(p), fast_level_otherprocs.AT(p).*schedule_item);
  }
  timer_bs1.stop();
  
  static Carpet::Timer timer_bs2 ("CarpetLib::dh::bs2");
  timer_bs2.start();
  srpvect const recv = alltoallv1 (dist::comm(), send);
  timer_bs2.stop();
  
  static Carpet::Timer timer_bs3 ("CarpetLib::dh::bs3");
  timer_bs3.start();
  (fast_level.*schedule_item).insert
    ((fast_level.*schedule_item).end(), recv.begin(), recv.end());
  timer_bs3.stop();
}



void
dh::
regrid_free (bool const do_init)
{
  if (do_init) {
    for (int ml = 0; ml < h.mglevels(); ++ ml) {
      for (int rl = 0; rl < h.reflevels(); ++ rl) {
        fast_boxes.AT(ml).AT(rl).fast_old2new_sync_sendrecv.clear();
        fast_boxes.AT(ml).AT(rl).fast_old2new_ref_prol_sendrecv.clear();
      }
    }
  } else {
    for (int ml = 0; ml < h.mglevels(); ++ ml) {
      for (int rl = 0; rl < h.reflevels(); ++ rl) {
        assert (fast_boxes.AT(ml).AT(rl).fast_old2new_sync_sendrecv.empty());
        assert (fast_boxes.AT(ml).AT(rl).fast_old2new_ref_prol_sendrecv.empty());
      }
    }
  }
}



void
dh::
recompose (int const rl, bool const do_prolongate)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (rl>=0 and rl<h.reflevels());
  
  static Carpet::Timer timer ("CarpetLib::dh::recompose");
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
#warning "TODO: If this works, rename do_prolongate to do_init here, and remove the do_prolongate parameter from ggf::recompose_fill"
#if 0
    for (comm_state state; not state.done(); state.step()) {
      for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
        (*f)->recompose_fill (state, rl, do_prolongate);
      }
    }
#endif
    if (do_prolongate) {
      for (comm_state state; not state.done(); state.step()) {
        for (list<ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
          (*f)->recompose_fill (state, rl, true);
        }
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
#if 0
      for (comm_state state; not state.done(); state.step()) {
        (*f)->recompose_fill (state, rl, do_prolongate);
      }
#endif
      if (do_prolongate) {
        for (comm_state state; not state.done(); state.step()) {
          (*f)->recompose_fill (state, rl, true);
        }
      }
      (*f)->recompose_free_old (rl);
    }
  }
  
  timer.stop ();
}



// Grid function management
dh::ggf_handle
dh::
add (ggf * const f)
{
  CHECKPOINT;
  return gfs.insert (gfs.end(), f);
}

void
dh::
erase (ggf_handle const fi)
{
  CHECKPOINT;
  gfs.erase (fi);
}



// Equality

bool
dh::full_dboxes::
operator== (full_dboxes const & b) const
{
  return
    exterior         == b.exterior         and
    all(all(is_outer_boundary == b.is_outer_boundary)) and
    outer_boundaries == b.outer_boundaries and
    communicated     == b.communicated     and
    boundaries       == b.boundaries       and
    owned            == b.owned            and
    buffers          == b.buffers          and
    active           == b.active           and
    sync             == b.sync             and
    bndref           == b.bndref           and
    ghosts           == b.ghosts           and
    interior         == b.interior;
}



// MPI datatypes

MPI_Datatype
mpi_datatype (dh::light_dboxes const &)
{
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static dh::light_dboxes s;
#define ENTRY(type, name)                                               \
    {                                                                   \
      sizeof s.name / sizeof(type), /* count elements */                \
        (char*)&s.name - (char*)&s, /* offsetof doesn't work (why?) */  \
        dist::mpi_datatype<type>(), /* find MPI datatype */             \
        STRINGIFY(name),    /* field name */                            \
        STRINGIFY(type),    /* type name */                             \
        }
    dist::mpi_struct_descr_t const descr[] = {
      ENTRY(int, exterior),
      ENTRY(int, owned),
      ENTRY(int, interior),
#if 0
      ENTRY(int, numactive),
      ENTRY(int, active),
#endif
      ENTRY(dh::light_dboxes::size_type, exterior_size),
      ENTRY(dh::light_dboxes::size_type, owned_size),
      ENTRY(dh::light_dboxes::size_type, active_size),
      {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}
    };
#undef ENTRY
    newtype =
      dist::create_mpi_datatype (sizeof descr / sizeof descr[0], descr,
                                 "dh::light::dboxes", sizeof s);
#if 0
    int type_size;
    MPI_Type_size (newtype, & type_size);
    assert (type_size <= sizeof s);
    MPI_Aint type_lb, type_ub;
    MPI_Type_lb (newtype, & type_lb);
    MPI_Type_ub (newtype, & type_ub);
    assert (type_ub - type_lb == sizeof s);
#endif
    initialised = true;
  }
  return newtype;
}

MPI_Datatype
mpi_datatype (dh::fast_dboxes const &)
{
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static dh::fast_dboxes s;
#define ENTRY(type, name)                                               \
    {                                                                   \
      sizeof s.name / sizeof(type), /* count elements */                \
        (char*)&s.name - (char*)&s, /* offsetof doesn't work (why?) */  \
        dist::mpi_datatype<type>(), /* find MPI datatype */             \
        STRINGIFY(name),    /* field name */                            \
        STRINGIFY(type),    /* type name */                             \
        }
    dist::mpi_struct_descr_t const descr[] = {
      ENTRY (dh::srpvect, fast_mg_rest_sendrecv),
      ENTRY (dh::srpvect, fast_mg_prol_sendrecv),
      ENTRY (dh::srpvect, fast_ref_prol_sendrecv),
      ENTRY (dh::srpvect, fast_ref_rest_sendrecv),
      ENTRY (dh::srpvect, fast_sync_sendrecv),
      ENTRY (dh::srpvect, fast_ref_bnd_prol_sendrecv),
      ENTRY (dh::srpvect, fast_old2new_sync_sendrecv),
      ENTRY (dh::srpvect, fast_old2new_ref_prol_sendrecv),
      {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}
    };
#undef ENTRY
    newtype =
      dist::create_mpi_datatype (sizeof descr / sizeof descr[0], descr,
                                 "dh::fast_dboxes", sizeof s);
    initialised = true;
  }
  return newtype;
}



// Memory usage

size_t
dh::
memory ()
  const
{
  return
    sizeof alldhi +             // memoryof (alldhi) +
    sizeof & h +                // memoryof (& h) +
    sizeof gh_handle +          // memoryof (gh_handle) +
    memoryof (ghost_widths) +
    memoryof (buffer_widths) +
    memoryof (prolongation_orders_space) +
    memoryof (light_boxes) +
    memoryof (fast_boxes) +
    memoryof (gfs);
}

size_t
dh::
allmemory ()
{
  size_t mem = memoryof(alldh);
  for (list<dh*>::const_iterator
         dhi = alldh.begin(); dhi != alldh.end(); ++ dhi)
  {
    mem += memoryof(**dhi);
  }
  return mem;
}

#if 0
int const dh::light::dboxes::maxactive;

void
dh::light_dboxes::
ibset2ibboxs (ibset const& s, ibbox* const bs, int& nbs)
{
  assert (s.setsize() <= maxactive);
  nbs = 0;
  for (ibset::const_iterator si = s.begin(); si != s.end(); ++si) {
    bs[nbs++] = *si;
  }
}

void
dh::light_dboxes::
ibboxs2ibset (ibbox const* const bs, int const& nbs, ibset& s)
{
  s = ibset();
  assert (nbs>=0 and nbs<=maxactive);
  for (int n=0; n<nbs; ++n) {
    s += bs[n];
  }
}
#endif

size_t
dh::light_dboxes::
memory ()
  const
{
  return
    memoryof (exterior) +
    memoryof (owned) +
    memoryof (interior) +
#if 0
    memoryof (numactive) +
    memoryof (active) +
#endif
    memoryof (exterior_size) +
    memoryof (owned_size) +
    memoryof (active_size);
}

size_t
dh::local_dboxes::
memory ()
  const
{
  return
    memoryof (buffers) +
    memoryof (buffers_stepped) +
    memoryof (active) +
    memoryof (restricted_region) +
    memoryof (restriction_boundaries) +
    memoryof (prolongation_boundaries) +
    memoryof (coarse_boundary) +
    memoryof (fine_boundary);
}

size_t
dh::full_dboxes::
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



// Input

istream &
dh::light_dboxes::
input (istream & is)
{
  // Regions:
  try {
    skipws (is);
    consume (is, "dh::light_dboxes:{");
    skipws (is);
    consume (is, "exterior:");
    is >> exterior;
    exterior_size = exterior.size();
    skipws (is);
    consume (is, "owned:");
    is >> owned;
    owned_size = owned.size();
    skipws (is);
    consume (is, "interior:");
    is >> interior;
    skipws (is);
    consume (is, "active_size:");
    is >> active_size;
    skipws (is);
    consume (is, "}");
  } catch (input_error & err) {
    cout << "Input error while reading a dh::light_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &
dh::local_dboxes::
input (istream & is)
{
  // Regions:
  try {
    skipws (is);
    consume (is, "dh::local_dboxes:{");
    skipws (is);
    consume (is, "buffers:");
    is >> buffers;
    skipws (is);
    consume (is, "buffers_stepped:");
    is >> buffers_stepped;
    skipws (is);
    consume (is, "active:");
    is >> active;
    skipws (is);
    consume (is, "restricted_region:");
    is >> restricted_region;
    skipws (is);
    consume (is, "restriction_boundaries:");
    is >> restriction_boundaries;
    skipws (is);
    consume (is, "prolongation_boundaries:");
    is >> prolongation_boundaries;
    skipws (is);
    consume (is, "coarse_boundary:");
    is >> coarse_boundary;
    skipws (is);
    consume (is, "fine_boundary:");
    is >> fine_boundary;
    skipws (is);
    consume (is, "}");
  } catch (input_error & err) {
    cout << "Input error while reading a dh::local_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &
dh::full_dboxes::
input (istream & is)
{
  // Regions:
  try {
    skipws (is);
    consume (is, "dh::full_dboxes:{");
    skipws (is);
    consume (is, "exterior:");
    is >> exterior;
    skipws (is);
    consume (is, "is_outer_boundary:");
    is >> is_outer_boundary;
    skipws (is);
    consume (is, "outer_boundaries:");
    is >> outer_boundaries;
    skipws (is);
    consume (is, "communicated:");
    is >> communicated;
    skipws (is);
    consume (is, "boundaries:");
    is >> boundaries;
    skipws (is);
    consume (is, "owned:");
    is >> owned;
    skipws (is);
    consume (is, "buffers:");
    is >> buffers;
    skipws (is);
    consume (is, "active:");
    is >> active;
    skipws (is);
    consume (is, "sync:");
    is >> sync;
    skipws (is);
    consume (is, "bndref:");
    is >> bndref;
    skipws (is);
    consume (is, "ghosts:");
    is >> ghosts;
    skipws (is);
    consume (is, "interior:");
    is >> interior;
    skipws (is);
    consume (is, "}");
  } catch (input_error & err) {
    cout << "Input error while reading a dh::full_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &
dh::fast_dboxes::
input (istream & is)
{
  // Communication schedule:
  try {
    skipws (is);
    consume (is, "dh::fast_dboxes:{");
    skipws (is);
    consume (is, "fast_mg_rest_sendrecv:");
    is >> fast_mg_rest_sendrecv;
    skipws (is);
    consume (is, "fast_mg_prol_sendrecv:");
    is >> fast_mg_prol_sendrecv;
    skipws (is);
    consume (is, "fast_ref_prol_sendrecv:");
    is >> fast_ref_prol_sendrecv;
    skipws (is);
    consume (is, "fast_ref_rest_sendrecv:");
    is >> fast_ref_rest_sendrecv;
    skipws (is);
    consume (is, "fast_sync_sendrecv:");
    is >> fast_sync_sendrecv;
    skipws (is);
    consume (is, "fast_ref_bnd_prol_sendrecv:");
    is >> fast_ref_bnd_prol_sendrecv;
    skipws (is);
    consume (is, "fast_old2new_sync_sendrecv:");
    is >> fast_old2new_sync_sendrecv;
    skipws (is);
    consume (is, "fast_old2new_ref_prol_sendrecv:");
    is >> fast_old2new_ref_prol_sendrecv;
    skipws (is);
    consume (is, "}");
  } catch (input_error & err) {
    cout << "Input error while reading a dh::fast_dboxes" << endl;
    throw err;
  }
  return is;
}



// Output

ostream &
dh::
output (ostream & os)
  const
{
  os << "dh:"
     << "ghost_widths=" << ghost_widths << ","
     << "buffer_widths=" << buffer_widths << ","
     << "prolongation_orders_space=" << prolongation_orders_space << ","
     << "light_boxes=" << light_boxes << ","
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
dh::light_dboxes::
output (ostream & os)
  const
{
  // Regions:
  os << "dh::light_dboxes:{" << eol
     << "   exterior: " << exterior << eol
     << "   owned: " << owned << eol
     << "   interior: " << interior << eol
     << "   active_size: " << active_size << eol
     << "}" << eol;
  return os;
}

ostream &
dh::local_dboxes::
output (ostream & os)
  const
{
  // Regions:
  os << "dh::local_dboxes:{" << eol
     << "   buffers: " << buffers << eol
     << "   buffers_stepped: " << buffers_stepped << eol
     << "   active: " << active << eol
     << "   restricted_region: " << restricted_region << eol
     << "   restriction_boundaries: " << restriction_boundaries << eol
     << "   prolongation_boundaries: " << prolongation_boundaries << eol
     << "   coarse_boundary: " << coarse_boundary << eol
     << "   fine_boundary: " << fine_boundary << eol
     << "}" << eol;
  return os;
}

ostream &
dh::full_dboxes::
output (ostream & os)
  const
{
  // Regions:
  os << "dh::full_dboxes:{" << eol
     << "   exterior: " << exterior << eol
     << "   is_outer_boundary: " << is_outer_boundary << eol
     << "   outer_boundaries: " << outer_boundaries << eol
     << "   communicated: " << communicated << eol
     << "   boundaries: " << boundaries << eol
     << "   owned: " << owned << eol
     << "   buffers: " << buffers << eol
     << "   active: " << active << eol
     << "   sync: " << sync << eol
     << "   bndref: " << bndref << eol
     << "   ghosts: " << ghosts << eol
     << "   interior: " << interior << eol
     << "}" << eol;
  return os;
}

ostream &
dh::fast_dboxes::
output (ostream & os)
  const
{
  // Communication schedule:
  os << "dh::fast_dboxes:{" << eol
     << "   fast_mg_rest_sendrecv: " << fast_mg_rest_sendrecv << eol
     << "   fast_mg_prol_sendrecv: " << fast_mg_prol_sendrecv << eol
     << "   fast_ref_prol_sendrecv: " << fast_ref_prol_sendrecv << eol
     << "   fast_ref_rest_sendrecv: " << fast_ref_rest_sendrecv << eol
     << "   fast_sync_sendrecv: " << fast_sync_sendrecv << eol
     << "   fast_ref_bnd_prol_sendrecv: " << fast_ref_bnd_prol_sendrecv << eol
     << "   fast_old2new_sync_sendrecv: " << fast_old2new_sync_sendrecv << eol
     << "   fast_old2new_ref_prol_sendrecv: " << fast_old2new_ref_prol_sendrecv << eol
     << "}" << eol;
  return os;
}
