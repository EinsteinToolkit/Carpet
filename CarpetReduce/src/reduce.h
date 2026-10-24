#ifndef CARPETREDUCE_H
#define CARPETREDUCE_H

#include <cctk.h>

#ifdef __cplusplus
namespace CarpetReduce {
extern "C" {
#endif

/* Scheduled functions */
int CarpetReduceStartup(void);

/* One global array reduction. EllSOR sums a 0-d residual through this. */
int Carpet_ReduceArraysGlobally(
    const cGH *GH, int dest_proc, int local_reduce_handle,
    int param_table_handle, int N_input_arrays,
    const void *const input_arrays[], int input_dims,
    const CCTK_INT input_array_dims[], const CCTK_INT input_array_type_codes[],
    int M_output_values, const CCTK_INT output_value_type_codes[],
    void *const output_values[]);

#ifdef __cplusplus
} /* extern "C" */
} /* namespace CarpetReduce */
#endif

#endif /* !defined(CARPETREDUCE_H) */
