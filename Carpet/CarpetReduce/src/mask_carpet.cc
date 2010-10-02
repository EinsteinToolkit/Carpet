#include <cassert>
#include <sstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <dh.hh>

#include <carpet.hh>

#include <loopcontrol.h>

#include "bits.h"
#include "mask_carpet.hh"
  
  
  
#define SWITCH_TO_LEVEL(cctkGH, rl)                     \
  do {                                                  \
    bool switch_to_level_ = true;                       \
    assert (is_singlemap_mode());                       \
    int const rl_ = (rl);                               \
    int const m_ = Carpet::map;                         \
    BEGIN_GLOBAL_MODE (cctkGH) {                        \
      ENTER_LEVEL_MODE (cctkGH, rl_) {                  \
        ENTER_SINGLEMAP_MODE (cctkGH, m_, CCTK_GF) {
#define END_SWITCH_TO_LEVEL                     \
        } LEAVE_SINGLEMAP_MODE;                 \
      } LEAVE_LEVEL_MODE;                       \
    } END_GLOBAL_MODE;                          \
    assert (switch_to_level_);                  \
    switch_to_level_ = false;                   \
  } while (false)



namespace CarpetMask {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  /**
   * Reduce the weight on the current and the next coarser level to
   * make things consistent.  Set the weight to 0 inside the
   * restriction region of the next coarser level, maybe to 1/2 near
   * the boundary of that region, and also to 1/2 near the
   * prolongation boundary of this level.
   */
  
  void
  CarpetMaskSetup (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (not is_singlemap_mode()) {
      CCTK_WARN (0, "This routine may only be called in singlemap mode");
    }
    
    if (reflevel > 0) {
      
      ivect const ione  = ivect(1);
      
      dh const & dd = *vdd.AT(Carpet::map);
      
      ivect const reffact =
        spacereffacts.AT(reflevel) / spacereffacts.AT(reflevel-1);
      assert (all (reffact == 2));
      
      
      
      // Set prolongation boundaries of this level
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        
        ibbox const & ext =
          dd.light_boxes.AT(mglevel).AT(reflevel).AT(component).exterior;
        ibset const & active = 
          dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).active;
        
        vect<vect<ibset,2>,dim> const & boundaries =
          dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).prolongation_boundaries;
        
        ibset const notactive = ext - active;
        
        for (int d=0; d<dim; ++d) {
          for (int f=0; f<2; ++f) {
            assert ((notactive & boundaries[d][f]).empty());
          }
        }
        
        LOOP_OVER_BSET (cctkGH, notactive, box, imin, imax) {
          
          if (verbose) {
            ostringstream buf;
            buf << "Setting buffer region on level " << reflevel << " to weight 0: " << imin << ":" << imax-ione;
            CCTK_INFO (buf.str().c_str());
          }
          
          // Set weight on the boundary to 0
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(CarpetMaskSetup_buffers,
                   i,j,k,
                   imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            iweight[ind] = 0;
          } LC_ENDLOOP3(CarpetMaskSetup_buffers);
          
        } END_LOOP_OVER_BSET;
        
        for (int d=0; d<dim; ++d) {
          for (int f=0; f<2; ++f) {
            LOOP_OVER_BSET (cctkGH, boundaries[d][f], box, imin, imax) {
              
              if (verbose) {
                ostringstream buf;
                buf << "Setting prolongation boundary on level " << reflevel << " direction " << d << " face " << f << " to weight 1/2: " << imin << ":" << imax-ione;
                CCTK_INFO (buf.str().c_str());
              }
              
              // Set weight on the boundary to 1/2
              int bmask = 0;
              for (unsigned b=0; b<BMSK(dim); ++b) {
                if (BGET(b,d)==f) {
                  bmask = BSET(bmask, b);
                }
              }
              assert (dim == 3);
#pragma omp parallel
              LC_LOOP3(CarpetMaskSetup_prolongation,
                       i,j,k,
                       imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                iweight[ind] &= ~bmask;
              } LC_ENDLOOP3(CarpetMaskSetup_prolongation);
              
            } END_LOOP_OVER_BSET;
          } // for f
        } // for d
        
      } END_LOCAL_COMPONENT_LOOP;
      
      
      
      // Set restriction region on next coarser level
      SWITCH_TO_LEVEL (cctkGH, reflevel-1) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS;
          
          ibset const & refined =
            dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).restricted_region;
          vect<vect<ibset,2>,dim> const & boundaries =
            dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).restriction_boundaries;
          
          LOOP_OVER_BSET (cctkGH, refined, box, imin, imax) {
            
            if (verbose) {
              ostringstream buf;
              buf << "Setting restricted region on level " << reflevel << " to weight 0: " << imin << ":" << imax-ione;
              CCTK_INFO (buf.str().c_str());
            }
            
            // Set weight in the restricted region to 0
            assert (dim == 3);
#pragma omp parallel
            LC_LOOP3(CarpetMaskSetup_restriction,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              iweight[ind] = 0;
            } LC_ENDLOOP3(CarpetMaskSetup_restriction);
            
          } END_LOOP_OVER_BSET;
          
          vector<int> imask (prod(ivect::ref(cctk_lsh)));
          vector<int> mask (prod(ivect::ref(cctk_lsh)));
          
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(CarpetMaskSetup_restriction_boundary_init,
                   i,j,k,
                   0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            imask[ind] = 0;
            mask[ind] = 0;
          } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_init);
          
          for (int d=0; d<dim; ++d) {
            for (int f=0; f<2; ++f) {
              LOOP_OVER_BSET (cctkGH, boundaries[d][f], box, imin, imax) {
                
                if (verbose) {
                  ostringstream buf;
                  buf << "Setting restriction boundary on level " << reflevel << " direction " << d << " face " << f << " to weight 1/2: " << imin << ":" << imax-ione;
                  CCTK_INFO (buf.str().c_str());
                }
                
                // Set weight on the boundary to 1/2
                int bmask = 0;
                for (unsigned b=0; b<BMSK(dim); ++b) {
                  if (BGET(b,d)==f) {
                    bmask = BSET(bmask, b);
                  }
                }
                assert (dim == 3);
#pragma omp parallel
                LC_LOOP3(CarpetMaskSetup_restriction_boundary_partial,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
                {
                  int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                  imask[ind] |= bmask;
                  if (mask[ind] == 0) {
                    mask[ind] = 1;
                  }
                  mask[ind] *= 2;
                } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_partial);
                
              } END_LOOP_OVER_BSET;
            } // for f
          } // for d
          
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(CarpetMaskSetup_restriction_boundary_apply,
                   i,j,k,
                   0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (mask[ind] > 0) {
              iweight[ind] &= imask[ind];
            }
          } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_apply);
          
        } END_LOCAL_COMPONENT_LOOP;
      } END_SWITCH_TO_LEVEL;
      
    } // if reflevel>0
  }
  
} // namespace CarpetMask
