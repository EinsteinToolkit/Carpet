/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/mask_coords.c,v 1.2 2004/08/02 11:43:35 schnetter Exp $ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



void
CoordBase_SetupMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int imin[3], imax[3];         /* boundary extent */
  
  int i, j, k;
  int d, f;
  int dd;
  
  int ierr;
  
  
  
  ierr = GetBoundarySpecification
    (6, nboundaryzones, is_internal, is_staggered, shiftout);
  if (ierr != 0) {
    CCTK_WARN (0, "Could not get boundary specification");
  }
  
  /* Loop over all dimensions and faces */
  for (d=0; d<3; ++d) {
    for (f=0; f<2; ++f) {
      /* If this processor has the outer boundary */
      if (cctk_bbox[2*d+f]) {
        
        int npoints;
        
        if (is_internal[2*d+f]) {
          /* The boundary extends inwards */
          npoints = - nboundaryzones[2*d+f] + shiftout[2*d+f];
        } else {
          /* The boundary extends outwards */
          npoints = nboundaryzones[2*d+f] + shiftout[2*d+f] - 1;
        }
        
        /* If there are boundary that should be ignored */
        if (npoints >= 0) {
          
          if (npoints < 0 || npoints > cctk_lsh[d]) {
            CCTK_WARN (0, "Illegal number of boundary points");
          }
          
          /* Calculate the extent of the boundary */
          for (dd=0; dd<3; ++dd) {
            imin[dd] = 0;
            imax[dd] = cctk_lsh[dd];
          }
          
          if (f==0) {
            /* lower face */
            imax[d] = imin[d] + npoints;
          } else {
            /* upper face */
            imin[d] = imax[d] - npoints;
          }
          
          /* Loop over the boundary */
          if (verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Setting symmetry boundary in direction %d face %d to weight 0", d, f);
          }
          for (k=imin[2]; k<imax[2]; ++k) {
            for (j=imin[1]; j<imax[1]; ++j) {
              for (i=imin[0]; i<imax[0]; ++i) {
                
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                weight[ind] = 0.0;
                
              }
            }
          }
          
          /* When the boundary is not staggered, then give the points
             on the boundary the weight 1/2 */
          if (! is_staggered[2*d+f]) {
            
            /* Adapt the extent of the boundary */
            if (f==0) {
              /* lower face */
              imin[d] = imax[d];
              imax[d] = imin[d] + 1;
            } else {
              /* upper face */
              imax[d] = imin[d];
              imin[d] = imax[d] - 1;
            }
          
            /* Loop over the points next to boundary */
            if (verbose) {
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Setting symmetry axis in direction %d face %d to weight 1/2", d, f);
            }
            for (k=imin[2]; k<imax[2]; ++k) {
              for (j=imin[1]; j<imax[1]; ++j) {
                for (i=imin[0]; i<imax[0]; ++i) {
                  
                  int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                  weight[ind] *= 0.5;
                  
                }
              }
            }
            
          } /* if the boundary is not staggered */
          
        } /* if there are boundary points */
        
      } /* if is outer boundary */
    } /* loop over faces */
  } /* loop over directions */
  
}
