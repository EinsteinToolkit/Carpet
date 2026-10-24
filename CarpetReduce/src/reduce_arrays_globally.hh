#ifndef CARPETREDUCE_ARRAYS_GLOBALLY_HH
#define CARPETREDUCE_ARRAYS_GLOBALLY_HH

/* LocalReduce writes these integers into the table key "global_operation".
   ReductionMax.c stores 1, ReductionMin.c stores 2, and ReductionSum.c
   (and the norms and the average) stores 3. */
enum {
  kLocalReduceOpMaximum = 1,
  kLocalReduceOpMinimum = 2,
  kLocalReduceOpSum = 3
};

/* LocalReduce stores 0 when the global operator must still divide or take
   a root, and 1 when the reduced value is already final. */
enum {
  kLocalReduceApplyNormalisation = 0,
  kLocalReduceSkipNormalisation = 1
};

inline int local_reduce_op_is_known(int operation) {
  return operation == kLocalReduceOpMaximum ||
         operation == kLocalReduceOpMinimum || operation == kLocalReduceOpSum;
}

/* Which MPI collective the global operator must use. Kept out of mpi.h so
   the standalone test can lock the mapping. */
enum {
  kMpiCombineMaximum = 1,
  kMpiCombineMinimum = 2,
  kMpiCombineSum = 3,
  kMpiCombineUnknown = -1
};

inline int mpi_combine_for(int operation) {
  switch (operation) {
  case kLocalReduceOpMaximum:
    return kMpiCombineMaximum;
  case kLocalReduceOpMinimum:
    return kMpiCombineMinimum;
  case kLocalReduceOpSum:
    return kMpiCombineSum;
  default:
    return kMpiCombineUnknown;
  }
}

#endif
