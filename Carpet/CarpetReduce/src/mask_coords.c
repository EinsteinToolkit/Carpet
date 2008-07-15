#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



void
CoordBase_SetupMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int imin[3], imax[3];         /* domain extent */
  int bmin[3], bmax[3];         /* boundary extent */
  
  int ierr;
  
  
  
  if (CCTK_IsFunctionAliased ("MultiPatch_GetBoundarySpecification")) {
    int const m = MultiPatch_GetMap (cctkGH);
    assert (m >= 0);
    ierr = MultiPatch_GetBoundarySpecification
      (m, 6, nboundaryzones, is_internal, is_staggered, shiftout);
  } else {
    ierr = GetBoundarySpecification
      (6, nboundaryzones, is_internal, is_staggered, shiftout);
  }
  if (ierr != 0) {
    CCTK_WARN (0, "Could not get boundary specification");
  }
  
  /* Loop over all dimensions and faces */
  for (int d=0; d<3; ++d) {
    for (int f=0; f<2; ++f) {
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
          
          /* Calculate the extent of the domain */
          for (int dd=0; dd<3; ++dd) {
            imin[dd] = 0;
            imax[dd] = cctk_lsh[dd];
          }
          
          /* Calculate the extent of the boundary */
          for (int dd=0; dd<3; ++dd) {
            bmin[dd] = imin[dd];
            bmax[dd] = imax[dd];
          }
          if (f==0) {
            /* lower face */
            bmax[d] = imin[d] + npoints;
          } else {
            /* upper face */
            bmin[d] = imax[d] - npoints;
          }
          
          /* Loop over the boundary */
          if (verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Setting boundary points in direction %d face %d to weight 0", d, f);
          }
#pragma omp parallel for
          for (int k=bmin[2]; k<bmax[2]; ++k) {
            for (int j=bmin[1]; j<bmax[1]; ++j) {
              for (int i=bmin[0]; i<bmax[0]; ++i) {
                
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                weight[ind] = 0.0;
                
              }
            }
          }
          
          /* When the boundary is not staggered, then give the points
             on the boundary the weight 1/2 */
          if (! is_staggered[2*d+f]) {
            
            /* Check whether the domain is empty */
            /* TODO: This check is flawed, and therefore disabled */
            if (0 && imin[d] == imax[d] - 1) {
              
              /* The domain is empty.  The correct thing to do would
                 be to set the weights to 0.  But this is boring,
                 because then there are no interior points left, and
                 all reductions become trivial.  Instead, leave the
                 non-staggered boundary points alone, and assume that
                 the user wants to perform a reduction in one
                 dimension less.  */
              /* do nothing */
              
            } else {
            
              /* Adapt the extent of the boundary */
              if (f==0) {
                /* lower face */
                bmin[d] = bmax[d];
                bmax[d] = bmin[d] + 1;
              } else {
                /* upper face */
                bmax[d] = bmin[d];
                bmin[d] = bmax[d] - 1;
              }
          
              /* Loop over the points next to boundary */
              if (verbose) {
                CCTK_VInfo (CCTK_THORNSTRING,
                            "Setting non-staggered boundary points in direction %d face %d to weight 1/2", d, f);
              }
#pragma omp parallel for
              for (int k=bmin[2]; k<bmax[2]; ++k) {
                for (int j=bmin[1]; j<bmax[1]; ++j) {
                  for (int i=bmin[0]; i<bmax[0]; ++i) {
                  
                    int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                    weight[ind] *= 0.5;
                  
                  }
                }
              }

            } /* if the domain is not empty */
            
          } /* if the boundary is not staggered */
          
        } /* if there are boundary points */
        
      } /* if is outer boundary */
    } /* loop over faces */
  } /* loop over directions */
  
}
