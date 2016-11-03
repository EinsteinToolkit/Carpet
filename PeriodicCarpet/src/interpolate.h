#ifndef PERIODICCARPET_H
#define PERIODICCARPET_H

#include <cctk.h>
#include <cctk_Arguments.h>

#ifdef __cplusplus
extern "C" {
#endif

CCTK_INT
PeriodicCarpet_Interpolate(
    CCTK_POINTER_TO_CONST restrict const cctkGH, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[],
    CCTK_INT const faces);

#ifdef __cplusplus
}
#endif

#endif /* ! defined PERIODICCARPET_H */
