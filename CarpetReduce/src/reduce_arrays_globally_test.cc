#include <cassert>
#include <cstdio>

#include "reduce_arrays_globally.hh"

int main() {
  assert(kLocalReduceOpMaximum == 1);
  assert(kLocalReduceOpMinimum == 2);
  assert(kLocalReduceOpSum == 3);
  assert(kLocalReduceApplyNormalisation == 0);
  assert(kLocalReduceSkipNormalisation == 1);
  assert(local_reduce_op_is_known(kLocalReduceOpMaximum));
  assert(local_reduce_op_is_known(kLocalReduceOpMinimum));
  assert(local_reduce_op_is_known(kLocalReduceOpSum));
  assert(!local_reduce_op_is_known(0));
  assert(!local_reduce_op_is_known(4));
  assert(mpi_combine_for(kLocalReduceOpMaximum) == kMpiCombineMaximum);
  assert(mpi_combine_for(kLocalReduceOpMinimum) == kMpiCombineMinimum);
  assert(mpi_combine_for(kLocalReduceOpSum) == kMpiCombineSum);
  assert(mpi_combine_for(0) == kMpiCombineUnknown);
  assert(mpi_combine_for(4) == kMpiCombineUnknown);
  std::printf("reduce arrays globally checks passed\n");
  return 0;
}
