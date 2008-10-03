#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <util_Table.h>

#include <loopcontrol.h>



void
MaskBase_InitMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
#if 0
  CCTK_INT symtable;
  CCTK_INT symmetry_handle[6];
  CCTK_INT symmetry_zone_width[6];
  
  int imin[3], imax[3];         /* boundary extent */
  
  int istat;
#endif
  
  
  
  /* Initialise the weight to 1 everywhere */
  if (verbose) {
    CCTK_INFO ("Initialising to weight 1");
  }
  
#pragma omp parallel
  LC_LOOP3(MaskBase_InitMask_interior,
           i,j,k,
           0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
           cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    
    int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
    weight[ind] = 1.0;
    
  } LC_ENDLOOP3(MaskBase_InitMask_interior);
  
  
  
#if 0
  /* Set the weight to 0 on inter-processor boundaries */
  if (verbose) {
    CCTK_INFO ("Setting inter-processor boundaries to weight 0");
  }
  
  /* Loop over all dimensions and faces */
  for (int d=0; d<3; ++d) {
    for (int f=0; f<2; ++f) {
      /* If this is an inter-processor boundary */
      if (! cctk_bbox[2*d+f]) {
        
        /* Calculate the extent of the boundary hyperslab */
        for (int dd=0; dd<3; ++dd) {
          imin[dd] = 0;
          imax[dd] = cctk_lsh[dd];
        }
        if (f==0) {
          /* lower face */
          imax[d] = imin[d] + cctk_nghostzones[d];
        } else {
          /* upper face */
          imin[d] = imax[d] - cctk_nghostzones[d];
        }
        
        /* Loop over the boundary slab */
#pragma omp parallel
        LC_LOOP3(MaskBase_InitMask_ghosts,
                 i,j,k,
                 imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                 cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          weight[ind] = 0.0;
          
        } LC_ENDLOOP3(MaskBase_InitMask_ghosts);
        
      }
    }
  }
#endif
  
  
  
#if 0
  
  /* Take the symmetry boundaries into account */
  if (verbose) {
    CCTK_INFO ("Setting symmetry boundaries to weight 0");
  }
  
  /* Get symmetry information */
  symtable = SymmetryTableHandleForGrid (cctkGH);
  if (symtable < 0) {
    CCTK_WARN (0, "Could not get symmetry table");
  }
  istat = Util_TableGetIntArray
    (symtable, 6, symmetry_handle, "symmetry_handle");
  if (istat != 6) {
    CCTK_WARN (0, "Could not get \"symmetry_handle\" entry from symmetry table");
  }
  istat = Util_TableGetIntArray
    (symtable, 6, symmetry_zone_width, "symmetry_zone_width");
  if (istat != 6) {
    CCTK_WARN (0, "Could not get \"symmetry_zone_width\" entry from symmetry table");
  }
  
  /* Loop over all dimensions and faces */
  for (int d=0; d<3; ++d) {
    for (int f=0; f<2; ++f) {
      /* If this is a symmetry face */
      if (symmetry_handle[2*d+f] >= 0) {
        /* If this processor has the outer boundary */
        if (cctk_bbox[2*d+f]) {
          
          if (symmetry_zone_width[2*d+f] < 0
              || symmetry_zone_width[2*d+f] > cctk_lsh[d]) {
            CCTK_WARN (0, "The symmetry table entry \"symmetry_zone_width\" contains illegal values");
          }
          
          /* Calculate the extent of the boundary hyperslab */
          for (int dd=0; dd<3; ++dd) {
            imin[dd] = 0;
            imax[dd] = cctk_lsh[dd];
          }
          if (f==0) {
            /* lower face */
            imax[d] = imin[d] + symmetry_zone_width[2*d+f];
          } else {
            /* upper face */
            imin[d] = imax[d] - symmetry_zone_width[2*d+f];
          }
          
          /* Loop over the boundary slab */
#pragma omp parallel
          LC_LOOP3(MaskBase_InitMask_boundary,
                   i,j,k,
                   imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            weight[ind] = 0.0;
            
          } LC_ENDLOOP3(MaskBase_InitMask_boundary);
          
        }
      }
    }
  }
#endif
}
