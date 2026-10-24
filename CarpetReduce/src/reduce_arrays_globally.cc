#include <cctk.h>
#include <util_Table.h>

#include <cstring>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <carpet.hh>

#include "reduce.hh"
#include "reduce_arrays_globally.hh"

// Other thorns include the CST header configs/<config>/bindings/include/carpet.hh,
// which is not Carpet/src/carpet.hh. Declare the helper defined there.
namespace Carpet {
MPI_Comm csamr_owner_reduction_comm();
}

namespace CarpetReduce {

using namespace std;
using namespace Carpet;

namespace {

int read_optional_int(int table, const char *key, int fallback) {
  int value = fallback;
  int const status = Util_TableGetInt(table, &value, key);
  if (status == 1 || status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
    return value;
  CCTK_VERROR("CarpetReduce global array reduction could not read table "
              "key %s (status %d).",
              key, status);
  return value;
}

void set_int(int table, const char *key, int value) {
  int const status = Util_TableSetInt(table, value, key);
  if (status < 0) {
    CCTK_VERROR("CarpetReduce global array reduction could not set table "
                "key %s (status %d).",
                key, status);
  }
}

bool table_has(int table, const char *key) {
  int const status = Util_TableQueryValueInfo(table, NULL, NULL, key);
  if (status < 0) {
    CCTK_VERROR("CarpetReduce global array reduction could not query table "
                "key %s (status %d).",
                key, status);
  }
  return status > 0;
}

bool is_complex_type(int type_code) {
  switch (type_code) {
  case CCTK_VARIABLE_COMPLEX:
    return true;
#ifdef HAVE_CCTK_COMPLEX8
  case CCTK_VARIABLE_COMPLEX8:
    return true;
#endif
#ifdef HAVE_CCTK_COMPLEX16
  case CCTK_VARIABLE_COMPLEX16:
    return true;
#endif
#ifdef HAVE_CCTK_COMPLEX32
  case CCTK_VARIABLE_COMPLEX32:
    return true;
#endif
  default:
    return false;
  }
}

MPI_Op mpi_op_for(int operation) {
  switch (mpi_combine_for(operation)) {
  case kMpiCombineMaximum:
    return MPI_MAX;
  case kMpiCombineMinimum:
    return MPI_MIN;
  case kMpiCombineSum:
    return MPI_SUM;
  default:
    CCTK_VERROR("CarpetReduce global array reduction was given "
                "global_operation %d. LocalReduce uses 1 (maximum), "
                "2 (minimum), or 3 (sum).",
                operation);
  }
  return MPI_SUM;
}

} // namespace

extern "C" int Carpet_ReduceArraysGlobally(
    const cGH *GH, int dest_proc, int local_reduce_handle,
    int param_table_handle, int N_input_arrays,
    const void *const input_arrays[], int input_dims,
    const CCTK_INT input_array_dims[], const CCTK_INT input_array_type_codes[],
    int M_output_values, const CCTK_INT output_value_type_codes[],
    void *const output_values[]) {
  if (GH == NULL)
    CCTK_ERROR("CarpetReduce global array reduction was given no grid "
               "hierarchy.");
  if (local_reduce_handle < 0) {
    CCTK_VERROR("CarpetReduce global array reduction was given local "
                "reduction handle %d.",
                local_reduce_handle);
  }
  if (N_input_arrays < 0 || M_output_values < 0 || input_dims < 0) {
    CCTK_VERROR("CarpetReduce global array reduction was given a negative "
                "count (arrays %d, outputs %d, dimensions %d).",
                N_input_arrays, M_output_values, input_dims);
  }
  if (N_input_arrays > 0 &&
      (input_arrays == NULL || input_array_type_codes == NULL)) {
    CCTK_ERROR("CarpetReduce global array reduction is missing its inputs.");
  }
  if (input_dims > 0 && input_array_dims == NULL) {
    CCTK_ERROR("CarpetReduce global array reduction is missing its "
               "dimensions.");
  }
  if (M_output_values > 0 &&
      (output_values == NULL || output_value_type_codes == NULL)) {
    CCTK_ERROR("CarpetReduce global array reduction is missing its outputs.");
  }

  bool created_table = false;
  if (Util_TableQueryNKeys(param_table_handle) < 0) {
    param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
    if (param_table_handle < 0) {
      CCTK_VERROR("CarpetReduce global array reduction could not create a "
                  "parameter table (status %d).",
                  param_table_handle);
    }
    created_table = true;
  }

  // LocalReduce writes global_operation only when this flag is set.
  set_int(param_table_handle, "global_calling", 1);
  // A 0-d sum never stores weight_sum. LocalReduce still reads it, so an
  // absent key would be an uninitialised value. Zero is the unused default.
  if (!table_has(param_table_handle, "weight_sum")) {
    int const status = Util_TableSetReal(param_table_handle, 0.0, "weight_sum");
    if (status < 0) {
      CCTK_VERROR("CarpetReduce global array reduction could not set "
                  "weight_sum (status %d).",
                  status);
    }
  }

  if (input_dims > 0) {
    vector<CCTK_INT> minimums(input_dims, 0);
    vector<CCTK_INT> maximums(input_dims);
    for (int dim = 0; dim < input_dims; ++dim)
      maximums[dim] = input_array_dims[dim];
    if (!table_has(param_table_handle, "input_array_min_subscripts")) {
      int const status = Util_TableSetGenericArray(
          param_table_handle, CCTK_VARIABLE_INT, input_dims, &minimums[0],
          "input_array_min_subscripts");
      if (status < 0) {
        CCTK_VERROR("CarpetReduce global array reduction could not set the "
                    "minimum subscripts (status %d).",
                    status);
      }
    }
    if (!table_has(param_table_handle, "input_array_max_subscripts")) {
      int const status = Util_TableSetGenericArray(
          param_table_handle, CCTK_VARIABLE_INT, input_dims, &maximums[0],
          "input_array_max_subscripts");
      if (status < 0) {
        CCTK_VERROR("CarpetReduce global array reduction could not set the "
                    "maximum subscripts (status %d).",
                    status);
      }
    }
  }

  int const local_status = CCTK_ReduceLocalArrays(
      input_dims, local_reduce_handle, param_table_handle, N_input_arrays,
      input_array_dims, input_array_type_codes, input_arrays, M_output_values,
      output_value_type_codes, output_values);
  if (local_status < 0) {
    CCTK_VERROR("CarpetReduce global array reduction: the local reduction "
                "failed (status %d).",
                local_status);
  }

  int const operation =
      read_optional_int(param_table_handle, "global_operation", -1);
  if (!local_reduce_op_is_known(operation)) {
    CCTK_VERROR("CarpetReduce global array reduction did not receive a "
                "global_operation from the local reduction (got %d).",
                operation);
  }
  int const perform_division = read_optional_int(
      param_table_handle, "perform_division", kLocalReduceSkipNormalisation);
  int const perform_2_root = read_optional_int(
      param_table_handle, "perform_2_root", kLocalReduceSkipNormalisation);
  int const perform_3_root = read_optional_int(
      param_table_handle, "perform_3_root", kLocalReduceSkipNormalisation);
  int const perform_4_root = read_optional_int(
      param_table_handle, "perform_4_root", kLocalReduceSkipNormalisation);
  if (perform_division == kLocalReduceApplyNormalisation ||
      perform_2_root == kLocalReduceApplyNormalisation ||
      perform_3_root == kLocalReduceApplyNormalisation ||
      perform_4_root == kLocalReduceApplyNormalisation) {
    CCTK_VERROR("CarpetReduce global array reduction was asked to divide "
                "or take a root (division %d, roots %d %d %d). That "
                "normalisation is not implemented.",
                perform_division, perform_2_root, perform_3_root,
                perform_4_root);
  }

  MPI_Comm comm = CarpetMPIComm();
  MPI_Comm const owner_comm = csamr_owner_reduction_comm();
  if (owner_comm != MPI_COMM_NULL)
    comm = owner_comm;
  if (comm == MPI_COMM_NULL) {
    CCTK_ERROR("CarpetReduce global array reduction ran before the MPI "
               "communicator existed.");
  }
  int nprocs = 0;
  int const size_status = MPI_Comm_size(comm, &nprocs);
  if (size_status != MPI_SUCCESS) {
    CCTK_VERROR("CarpetReduce global array reduction could not read the "
                "communicator size (status %d).",
                size_status);
  }
  bool allreduce = false;
  if (dest_proc == -1) {
    allreduce = true;
  } else if (dest_proc >= 0 && dest_proc < nprocs) {
    allreduce = false;
  } else {
    CCTK_VERROR("CarpetReduce global array reduction destination processor "
                "%d is outside 0..%d and is not -1.",
                dest_proc, nprocs - 1);
  }

  MPI_Op const op = mpi_op_for(operation);
  for (int out = 0; out < M_output_values; ++out) {
    int const type_code = output_value_type_codes[out];
    if (is_complex_type(type_code)) {
      CCTK_VERROR("CarpetReduce global array reduction does not implement "
                  "complex type %d.",
                  type_code);
    }
    if (output_values[out] == NULL) {
      CCTK_VERROR("CarpetReduce global array reduction output %d is null.",
                  out);
    }
    int const type_size = CCTK_VarTypeSize(type_code);
    if (type_size <= 0) {
      CCTK_VERROR("CarpetReduce global array reduction does not know the "
                  "size of type %d.",
                  type_code);
    }
    vector<char> local(type_size);
    memcpy(&local[0], output_values[out], type_size);
    MPI_Datatype const mpitype = CarpetMPIDatatype(type_code);
    int mpi_status = MPI_SUCCESS;
    if (allreduce) {
      mpi_status = MPI_Allreduce(&local[0], output_values[out], 1, mpitype, op,
                                 comm);
    } else {
      mpi_status = MPI_Reduce(&local[0], output_values[out], 1, mpitype, op,
                              dest_proc, comm);
    }
    if (mpi_status != MPI_SUCCESS) {
      CCTK_VERROR("CarpetReduce global array reduction MPI call failed "
                  "(status %d, output %d).",
                  mpi_status, out);
    }
  }

  if (created_table) {
    int const status = Util_TableDestroy(param_table_handle);
    if (status < 0) {
      CCTK_VERROR("CarpetReduce global array reduction could not destroy "
                  "its parameter table (status %d).",
                  status);
    }
  }
  return 0;
}

} // namespace CarpetReduce
