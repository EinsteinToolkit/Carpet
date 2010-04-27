#include <cassert>
#include <sstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <dh.hh>

#include <carpet.hh>

#include <loopcontrol.h>

#include "mask_carpet.hh"



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
      
      gh const & hh = *vhh.AT(Carpet::map);
      dh const & dd = *vdd.AT(Carpet::map);
      
      ivect const reffact =
        spacereffacts.AT(reflevel) / spacereffacts.AT(reflevel-1);
      assert (all (reffact == 2));
      
      
      
      // Set prolongation boundaries of this level
      {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          
          DECLARE_CCTK_ARGUMENTS;
          
          ibbox const & ext =
            dd.light_boxes.AT(mglevel).AT(reflevel).AT(component).exterior;
          ibset const & active = 
            dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).active;
          
          vect<ibset,dim> const & boundaries =
            dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).prolongation_boundaries;
          
          ibset const notactive = ext - active;
          
          for (int d=0; d<dim; ++d) {
            assert ((notactive & boundaries[d]).empty());
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
              weight[ind] = 0.0;
            } LC_ENDLOOP3(CarpetMaskSetup_buffers);
            
          } END_LOOP_OVER_BSET;
          
          for (int d=0; d<dim; ++d) {
            LOOP_OVER_BSET (cctkGH, boundaries[d], box, imin, imax) {
              
              if (verbose) {
                ostringstream buf;
                buf << "Setting prolongation boundary on level " << reflevel << " direction " << d << " to weight 1/2: " << imin << ":" << imax-ione;
                CCTK_INFO (buf.str().c_str());
              }
              
              // Set weight on the boundary to 1/2
              assert (dim == 3);
#pragma omp parallel
              LC_LOOP3(CarpetMaskSetup_prolongation,
                       i,j,k,
                       imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                weight[ind] *= 0.5;
              } LC_ENDLOOP3(CarpetMaskSetup_prolongation);
              
              
            } END_LOOP_OVER_BSET;
          } // for d
          
        } END_LOCAL_COMPONENT_LOOP;
      }
      
      
      
      // Set restriction region on next coarser level
      {
        int const oldreflevel = reflevel;
        int const oldgrouptype = mc_grouptype;
        int const oldmap = Carpet::map;
        leave_singlemap_mode (cctkGH);
        leave_level_mode (cctkGH);
        enter_level_mode (cctkGH, oldreflevel-1);
        enter_singlemap_mode (cctkGH, oldmap, oldgrouptype);
        
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          
          DECLARE_CCTK_ARGUMENTS;
          
          ibset const & refined =
            dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).restricted_region;
          vect<ibset,dim> const & boundaries =
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
              weight[ind] = 0;
            } LC_ENDLOOP3(CarpetMaskSetup_restriction);
            
          } END_LOOP_OVER_BSET;
          
          assert (dim == 3);
          vector<int> mask (cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]);
          
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(CarpetMaskSetup_restriction_boundary_init,
                   i,j,k,
                   0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            mask[ind] = 0;
          } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_init);
          
          for (int d=0; d<dim; ++d) {
            LOOP_OVER_BSET (cctkGH, boundaries[d], box, imin, imax) {
              
              if (verbose) {
                ostringstream buf;
                buf << "Setting restriction boundary on level " << reflevel << " direction " << d << " to weight 1/2: " << imin << ":" << imax-ione;
                CCTK_INFO (buf.str().c_str());
              }
              
              // Set weight on the boundary to 1/2
              assert (dim == 3);
#pragma omp parallel
              LC_LOOP3(CarpetMaskSetup_restriction_boundary_partial,
                       i,j,k,
                       imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                if (mask[ind] == 0) {
                  mask[ind] = 1;
                }
                mask[ind] *= 2;
              } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_partial);
              
            } END_LOOP_OVER_BSET;
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
              weight[ind] *= 1.0 - 1.0 / mask[ind];
            }
          } LC_ENDLOOP3(CarpetMaskSetup_restriction_boundary_apply);
          
        } END_LOCAL_COMPONENT_LOOP;
        
        leave_singlemap_mode (cctkGH);
        leave_level_mode (cctkGH);
        enter_level_mode (cctkGH, oldreflevel);
        enter_singlemap_mode (cctkGH, oldmap, oldgrouptype);
      }
      
    } // if reflevel>0
  }
  
} // namespace CarpetMask
