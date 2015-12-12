#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "interpolate.h"


void
PeriodicCarpet_RegisterBC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT faces[6];
  faces[0] = faces[1] = periodic || periodic_x;
  faces[2] = faces[3] = periodic || periodic_y;
  faces[4] = faces[5] = periodic || periodic_z;
  
  CCTK_INT width[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  int ierr = GetBoundarySpecification
    (6, width, is_internal, is_staggered, shiftout);
  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "Could not get the boundary specification");
  }
  
  CCTK_INT const handle = SymmetryRegister("periodic");
  if (handle < 0) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Could not register periodicity boundary condition");
  }
  
  ierr = SymmetryRegisterGrid(cctkGH, handle, faces, width);
  if (ierr < 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "Could not register the periodic boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }

  ierr = SymmetryRegisterGridInterpolator(cctkGH, handle, PeriodicCarpet_Interpolate);
  if (ierr < 0)
  {
    CCTK_WARN (0, "Could not register the symmetry interpolator");
  }
}
