#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include "interpolate.h"

CCTK_INT
PeriodicCarpet_Interpolate(
    CCTK_POINTER_TO_CONST restrict const cctkGH, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[],
    CCTK_INT const faces) {
  DECLARE_CCTK_PARAMETERS;

  int dim = 3;

  int do_periodic[dim];

  CCTK_POINTER_TO_CONST restrict new_interp_coords[dim];
  CCTK_INT newfaces;

  int dir;

  int iret;

  int m;
  int n;
  int d;

  int ierr;

  CCTK_REAL physical_min[dim];
  CCTK_REAL physical_max[dim];
  CCTK_REAL tmp[dim];
  ierr = GetDomainSpecification(dim, physical_min, physical_max, tmp, tmp, tmp,
                                tmp, tmp);
  assert(!ierr);

  CCTK_REAL periodic_bbox[6] = {physical_min[0], physical_max[0],
                                physical_min[1], physical_max[1],
                                physical_min[2], physical_max[2]};

  /* Get symmetry information */
  do_periodic[0] = periodic || periodic_x;
  do_periodic[1] = periodic || periodic_y;
  do_periodic[2] = periodic || periodic_z;

  newfaces = faces;
  for (d = 0; d < dim; ++d) {
    if (do_periodic[d]) {
      /* Lower face */
      CCTK_INT dl = 2 * d;
      assert(newfaces & (1 << dl));
      newfaces &= ~(1 << dl);

      /* Upper face */
      CCTK_INT du = 2 * d + 1;
      assert(newfaces & (1 << du));
      newfaces &= ~(1 << du);
    }
  }

  /* Fold coordinates */
  assert(interp_coords_type == CCTK_VARIABLE_REAL);
  for (dir = 0; dir < dim; ++dir) {
    if (do_periodic[dir]) {
      CCTK_REAL *restrict coords_in_dir =
          malloc(N_interp_points * sizeof(CCTK_REAL));
      assert(N_interp_points == 0 || coords_in_dir);

      for (n = 0; n < N_interp_points; ++n) {
        CCTK_REAL const pos = ((CCTK_REAL const *)interp_coords[dir])[n];
        CCTK_REAL newpos = pos;
        while (newpos < periodic_bbox[2 * dir]) {
          newpos = newpos + periodic_bbox[2 * dir + 1] - periodic_bbox[2 * dir];
        }
        while (newpos > periodic_bbox[2 * dir + 1]) {
          newpos = newpos - periodic_bbox[2 * dir + 1] + periodic_bbox[2 * dir];
        }
        coords_in_dir[n] = newpos;
      }
      new_interp_coords[dir] = coords_in_dir;
    } else {
      new_interp_coords[dir] = interp_coords[dir];
    }
  }

  /* Recursive call */
  iret = SymmetryInterpolateFaces(
      cctkGH, N_dims, local_interp_handle, param_table_handle,
      coord_system_handle, N_interp_points, interp_coords_type,
      (CCTK_POINTER_TO_CONST)new_interp_coords, N_input_arrays,
      input_array_indices, N_output_arrays, output_array_types,
      (CCTK_POINTER_TO_CONST)output_arrays, newfaces);

  /* Free coordinates */
  for (dir = 0; dir < dim; ++dir) {
    if (do_periodic[dir]) {
      free((void *)new_interp_coords[dir]);
    }
  }

  /* Return */
  return iret;
}
