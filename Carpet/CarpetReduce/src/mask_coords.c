#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>



void
CoordBase_SetupMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int bnd_points[6];            /* points outside the domain */
  int int_points[3];            /* global interior points */
  
  int gmin[3], gmax[3];         /* global domain extent */
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
  
  
  
  /* Calculate the number of boundary points */
  for (int d=0; d<3; ++d) {
    for (int f=0; f<2; ++f) {
      
      if (is_internal[2*d+f]) {
        /* The boundary extends inwards */
        bnd_points[2*d+f] = shiftout[2*d+f] - is_staggered[2*d+f];
      } else {
        /* The boundary extends outwards */
        bnd_points[2*d+f] =
          nboundaryzones[2*d+f] + shiftout[2*d+f] + is_staggered[2*d+f] - 1;
      }
      
    }
  }
  
  /* Calculate the global extent of the domain */
  for (int d=0; d<3; ++d) {
    gmin[d] = 0;
    gmax[d] = cctk_gsh[d];
  }
  
  /* Ensure that the domain specification is consistent */
  for (int d=0; d<3; ++d) {
    int_points[d] =
      (gmax[d] - bnd_points[2*d]) - (gmin[d] + bnd_points[2*d+1]);
    if (int_points[d] < 0) {
      CCTK_WARN (0, "Number of internal grid points is negative");
    }
  }
  
  
  
  /* Loop over all dimensions and faces */
  for (int d=0; d<3; ++d) {
    for (int f=0; f<2; ++f) {
      /* If this processor has the outer boundary */
      if (cctk_bbox[2*d+f]) {
        
        /* If there are boundary that should be ignored */
        if (bnd_points[2*d+f] >= 0) {
          
          if (bnd_points[2*d+f] < 0 || bnd_points[2*d+f] > cctk_lsh[d]) {
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
            bmax[d] = imin[d] + bnd_points[2*d+f];
          } else {
            /* upper face */
            bmin[d] = imax[d] - bnd_points[2*d+f];
          }
          
          /* Loop over the boundary */
          if (verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Setting boundary points in direction %d face %d to weight 0", d, f);
          }
#pragma omp parallel
          LC_LOOP3(CoordBase_SetupMask_boundary,
                   i,j,k,
                   bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            weight[ind] = 0.0;
            
          } LC_ENDLOOP3(CoordBase_SetupMask);
          
          /* When the boundary is not staggered, then give the points
             on the boundary the weight 1/2 */
          if (! is_staggered[2*d+f]) {
            
            /* Check whether the domain is empty */
            if (int_points[d] == 0) {
              
              /* The domain is empty.  The correct thing to do would
                 be to set the weights to 0.  But this is boring,
                 because then there are no interior points left, and
                 all reductions become trivial.  Instead, leave the
                 non-staggered boundary points alone, and assume that
                 the user wants to perform a reduction in one fewer
                 dimensions.  */
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
#pragma omp parallel
              LC_LOOP3(CoordBase_SetupMask_boundary2,
                       i,j,k,
                       bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                weight[ind] *= 0.5;
                
              } LC_ENDLOOP3(CoordBase_SetupMask_boundary2);
              
            } /* if the domain is not empty */
            
          } /* if the boundary is not staggered */
          
        } /* if there are boundary points */
        
      } /* if is outer boundary */
    } /* loop over faces */
  } /* loop over directions */
  
}
