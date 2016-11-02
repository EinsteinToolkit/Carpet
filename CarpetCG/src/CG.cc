#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "TATelliptic.h"
#include "carpet.hh"

namespace Carpet {
// TODO: fix this
void Restrict(const cGH *cgh);
};

namespace CarpetCG {

using namespace std;
using namespace Carpet;

extern "C" {
int CarpetCG_register();
}

int CarpetCG_solve(const cGH *const cctkGH, const int *const var,
                   const int *const res, const int nvars,
                   const int options_table, const calcfunc calcres,
                   const calcfunc applybnds, void *const userdata);

namespace common {

const int *var;
const int *res;
int nvars;
int options_table;
calcfunc calcres;
calcfunc applybnds;
void *userdata;

int ierr;

vector<CCTK_INT> nboundaryzones(2 * dim);

CCTK_REAL factor;
vector<CCTK_REAL> factors;

vector<int> fromindex; // Can't pass things through CallLocalFunction.
vector<int> toindex;   // Instead assume everything is going from a GF to a GF.

CCTK_REAL realconstant;     // With only one required constant
CCTK_REAL realoutput;       // And one output
CCTK_REAL realoutput_count; // And one _more_ output

void call_calcres(cGH *const cctkGH) {
  if (ierr)
    return;
  ierr = calcres(cctkGH, options_table, userdata);
}

void call_applybnds(cGH *const cctkGH) {
  if (ierr)
    return;
  ierr = applybnds(cctkGH, options_table, userdata);
}

void global_sum(cGH const *restrict const cctkGH,
                CCTK_REAL *restrict const var) {
  static int initialised = 0;
  static int sum_handle;

  CCTK_REAL local, global;

  int ierr;

  if (!initialised) {
    initialised = 1;
    sum_handle = CCTK_ReductionArrayHandle("sum");
    assert(sum_handle >= 0);
  }

  local = *var;
  ierr = CCTK_ReduceLocScalar(cctkGH, -1, sum_handle, &local, &global,
                              CCTK_VARIABLE_REAL);
  assert(!ierr);
  *var = global;
}

// Copy
void call_copy(cGH *const cctkGH) {

  cGroup groupdata;
  ierr = CCTK_GroupData(CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert(!ierr);
  cGroupDynamicData groupdyndata;
  ierr = CCTK_GroupDynamicData(cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                               &groupdyndata);
  assert(!ierr);
  const int thedim = groupdata.dim;
  assert(thedim >= 0 && thedim <= dim);
  assert(thedim == groupdyndata.dim);
  int lsh[dim], nghostzones[dim];
  for (int d = 0; d < thedim; ++d) {
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
  }
  for (int d = 0; d < dim; ++d) {
    assert(lsh[d] >= 0);
    assert(nghostzones[d] >= 0 && 2 * nghostzones[d] <= lsh[d]);
  }

  int const lsize = lsh[0] * lsh[1] * lsh[2];
  int n;
  double chksum = 0;
  double chksum2 = 0;

  for (n = 0; n < nvars; ++n) {
    CCTK_REAL *dst = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, toindex[n]);
    CCTK_REAL *src = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, fromindex[n]);
    //        cout << "Copy indices " << toindex[n] << " (" <<
    //        CCTK_VarName(toindex[n]) << ") " << fromindex[n] << " (" <<
    //        CCTK_VarName(fromindex[n]) << ") " << endl;

    memcpy(dst, src, lsize * sizeof(*dst));
    for (int k = 0; k < lsh[2]; ++k) {
      for (int j = 0; j < lsh[2]; ++j) {
        for (int i = 0; i < lsh[2]; ++i) {
          int ind = i + lsh[0] * (j + lsh[1] * k);
          chksum += dst[ind];
          chksum2 += src[ind];
        }
      }
    }
    //        cout << "copied " << chksum << " " << chksum2 << endl;
  }
}

// Correct residual sign
void call_correct_residual_sign(cGH *const cctkGH) {

  cGroup groupdata;
  ierr = CCTK_GroupData(CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert(!ierr);
  cGroupDynamicData groupdyndata;
  ierr = CCTK_GroupDynamicData(cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                               &groupdyndata);
  assert(!ierr);
  const int thedim = groupdata.dim;
  assert(thedim >= 0 && thedim <= dim);
  assert(thedim == groupdyndata.dim);
  int lsh[dim], nghostzones[dim], bbox[2 * dim];
  int lbnd[dim], ubnd[dim];
  for (int d = 0; d < thedim; ++d) {
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = 0; d < dim; ++d) {
    assert(lsh[d] >= 0);
    assert(nghostzones[d] >= 0 && 2 * nghostzones[d] <= lsh[d]);
  }
  for (int d = 0; d < dim; ++d) {
    lbnd[d] =
        (bbox[2 * d] ? ((nboundaryzones[2 * d] >= 0) ? nboundaryzones[2 * d]
                                                     : nghostzones[d])
                     : nghostzones[d]);
    ubnd[d] = lsh[d] - (bbox[2 * d + 1] ? ((nboundaryzones[2 * d + 1] >= 0)
                                               ? nboundaryzones[2 * d + 1]
                                               : nghostzones[d])
                                        : nghostzones[d]);
  }

  int const lsize = lsh[0] * lsh[1] * lsh[2];
  int n;
  int i, j, k;
  int ind;

  for (n = 0; n < nvars; ++n) {
    CCTK_REAL *src = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, fromindex[n]);
    for (k = lbnd[2]; k < ubnd[2]; ++k) {
      for (j = lbnd[1]; j < ubnd[1]; ++j) {
        for (i = lbnd[0]; i < ubnd[0]; ++i) {
          ind = i + lsh[0] * (j + lsh[1] * k);
          src[ind] *= copysign(1.0, factors[n]);
        }
      }
    }
  }
}

// ``Precondition''. Ha ha ha.
void call_apply_preconditioner(cGH *const cctkGH) {

  cGroup groupdata;
  ierr = CCTK_GroupData(CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert(!ierr);
  cGroupDynamicData groupdyndata;
  ierr = CCTK_GroupDynamicData(cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                               &groupdyndata);
  assert(!ierr);
  const int thedim = groupdata.dim;
  assert(thedim >= 0 && thedim <= dim);
  assert(thedim == groupdyndata.dim);
  int lsh[dim], nghostzones[dim], bbox[2 * dim];
  int lbnd[dim], ubnd[dim];
  for (int d = 0; d < thedim; ++d) {
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = 0; d < dim; ++d) {
    assert(lsh[d] >= 0);
    assert(nghostzones[d] >= 0 && 2 * nghostzones[d] <= lsh[d]);
  }
  for (int d = 0; d < dim; ++d) {
    lbnd[d] =
        (bbox[2 * d] ? ((nboundaryzones[2 * d] >= 0) ? nboundaryzones[2 * d]
                                                     : nghostzones[d])
                     : nghostzones[d]);
    ubnd[d] = lsh[d] - (bbox[2 * d + 1] ? ((nboundaryzones[2 * d + 1] >= 0)
                                               ? nboundaryzones[2 * d + 1]
                                               : nghostzones[d])
                                        : nghostzones[d]);
  }

  int const lsize = lsh[0] * lsh[1] * lsh[2];
  int n;
  int i, j, k;
  int ind;

  for (n = 0; n < nvars; ++n) {
    CCTK_REAL *dst = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, toindex[n]);
    CCTK_REAL *src = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, fromindex[n]);
    for (k = lbnd[2]; k < ubnd[2]; ++k) {
      for (j = lbnd[1]; j < ubnd[1]; ++j) {
        for (i = lbnd[0]; i < ubnd[0]; ++i) {
          ind = i + lsh[0] * (j + lsh[1] * k);
          dst[ind] = factor * fabs(factors[n]) * src[ind];
        }
      }
    }
  }
}

// Add scaled
void call_add_scaled(cGH *const cctkGH) {

  cGroup groupdata;
  ierr = CCTK_GroupData(CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert(!ierr);
  cGroupDynamicData groupdyndata;
  ierr = CCTK_GroupDynamicData(cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                               &groupdyndata);
  assert(!ierr);
  const int thedim = groupdata.dim;
  assert(thedim >= 0 && thedim <= dim);
  assert(thedim == groupdyndata.dim);
  int lsh[dim], nghostzones[dim], bbox[2 * dim];
  int lbnd[dim], ubnd[dim];
  for (int d = 0; d < thedim; ++d) {
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = 0; d < dim; ++d) {
    assert(lsh[d] >= 0);
    assert(nghostzones[d] >= 0 && 2 * nghostzones[d] <= lsh[d]);
  }
  for (int d = 0; d < dim; ++d) {
    lbnd[d] =
        (bbox[2 * d] ? ((nboundaryzones[2 * d] >= 0) ? nboundaryzones[2 * d]
                                                     : nghostzones[d])
                     : nghostzones[d]);
    ubnd[d] = lsh[d] - (bbox[2 * d + 1] ? ((nboundaryzones[2 * d + 1] >= 0)
                                               ? nboundaryzones[2 * d + 1]
                                               : nghostzones[d])
                                        : nghostzones[d]);
  }

  int n;
  int i, j, k;
  int ind;

  CCTK_REAL const alpha = realconstant;
  double chksum = 0;
  double chksum2 = 0;

  for (n = 0; n < nvars; ++n) {
    CCTK_REAL *src = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, fromindex[n]);
    CCTK_REAL *dst = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, toindex[n]);
    //        cout << "Indices " << toindex[n] << " (" <<
    //        CCTK_VarName(toindex[n]) << ") " << fromindex[n] << " (" <<
    //        CCTK_VarName(fromindex[n]) << ") " << endl;
    for (k = lbnd[2]; k < ubnd[2]; ++k) {
      for (j = lbnd[1]; j < ubnd[1]; ++j) {
        for (i = lbnd[0]; i < ubnd[0]; ++i) {
          ind = i + lsh[0] * (j + lsh[1] * k);
          dst[ind] += alpha * src[ind];
          chksum += dst[ind];
          chksum2 += src[ind];
        }
      }
    }
  }
  //      cout << "add scaled " << alpha << " " << chksum << " " << chksum2 <<
  //      endl;
}

// Scale and add
void call_scale_and_add(cGH *const cctkGH) {

  cGroup groupdata;
  ierr = CCTK_GroupData(CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert(!ierr);
  cGroupDynamicData groupdyndata;
  ierr = CCTK_GroupDynamicData(cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                               &groupdyndata);
  assert(!ierr);
  const int thedim = groupdata.dim;
  assert(thedim >= 0 && thedim <= dim);
  assert(thedim == groupdyndata.dim);
  int lsh[dim], nghostzones[dim], bbox[2 * dim];
  int lbnd[dim], ubnd[dim];
  for (int d = 0; d < thedim; ++d) {
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = 0; d < dim; ++d) {
    assert(lsh[d] >= 0);
    assert(nghostzones[d] >= 0 && 2 * nghostzones[d] <= lsh[d]);
  }
  for (int d = 0; d < dim; ++d) {
    lbnd[d] =
        (bbox[2 * d] ? ((nboundaryzones[2 * d] >= 0) ? nboundaryzones[2 * d]
                                                     : nghostzones[d])
                     : nghostzones[d]);
    ubnd[d] = lsh[d] - (bbox[2 * d + 1] ? ((nboundaryzones[2 * d + 1] >= 0)
                                               ? nboundaryzones[2 * d + 1]
                                               : nghostzones[d])
                                        : nghostzones[d]);
  }

  int n;
  int i, j, k;
  int ind;

  CCTK_REAL const alpha = realconstant;

  for (n = 0; n < nvars; ++n) {
    CCTK_REAL *src = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, fromindex[n]);
    CCTK_REAL *dst = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, toindex[n]);
    for (k = lbnd[2]; k < ubnd[2]; ++k) {
      for (j = lbnd[1]; j < ubnd[1]; ++j) {
        for (i = lbnd[0]; i < ubnd[0]; ++i) {
          ind = i + lsh[0] * (j + lsh[1] * k);
          dst[ind] = alpha * dst[ind] + src[ind];
        }
      }
    }
  }
}

// Dot product
void call_dot_product(cGH *const cctkGH) {

  cGroup groupdata;
  ierr = CCTK_GroupData(CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert(!ierr);
  cGroupDynamicData groupdyndata;
  ierr = CCTK_GroupDynamicData(cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                               &groupdyndata);
  assert(!ierr);
  const int thedim = groupdata.dim;
  assert(thedim >= 0 && thedim <= dim);
  assert(thedim == groupdyndata.dim);
  int lsh[dim], nghostzones[dim], bbox[2 * dim];
  int lbnd[dim], ubnd[dim];
  for (int d = 0; d < thedim; ++d) {
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = thedim; d < dim; ++d) {
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2 * d] = groupdyndata.bbox[2 * d];
    bbox[2 * d + 1] = groupdyndata.bbox[2 * d + 1];
  }
  for (int d = 0; d < dim; ++d) {
    assert(lsh[d] >= 0);
    assert(nghostzones[d] >= 0 && 2 * nghostzones[d] <= lsh[d]);
  }
  for (int d = 0; d < dim; ++d) {
    lbnd[d] =
        (bbox[2 * d] ? ((nboundaryzones[2 * d] >= 0) ? nboundaryzones[2 * d]
                                                     : nghostzones[d])
                     : nghostzones[d]);
    ubnd[d] = lsh[d] - (bbox[2 * d + 1] ? ((nboundaryzones[2 * d + 1] >= 0)
                                               ? nboundaryzones[2 * d + 1]
                                               : nghostzones[d])
                                        : nghostzones[d]);
  }

  int n;
  int i, j, k;
  int ind;
  CCTK_REAL res;

  double chksum = 0;
  double chksum2 = 0;

  res = 0;
  realoutput_count = 0;

  // Shift the count to here instead - what a waste of space
  for (k = 0; k < lsh[2]; ++k) {
    for (j = 0; j < lsh[1]; ++j) {
      for (i = 0; i < lsh[0]; ++i) {
        realoutput_count++;
      }
    }
  }

  for (n = 0; n < nvars; ++n) {
    CCTK_REAL *a = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, toindex[n]);
    CCTK_REAL *b = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, fromindex[n]);
    //        cout << "Indices " << toindex[n] << " (" <<
    //        CCTK_VarName(toindex[n]) << ") " << fromindex[n] << " (" <<
    //        CCTK_VarName(fromindex[n]) << ") " << endl;
    //        cout << "Pointers " << a << " " << b << endl;

    for (k = lbnd[2]; k < ubnd[2]; ++k) {
      for (j = lbnd[1]; j < ubnd[1]; ++j) {
        for (i = lbnd[0]; i < ubnd[0]; ++i) {
          ind = i + lsh[0] * (j + lsh[1] * k);
          res += (a[ind]) * (b[ind]);
          chksum += a[ind];
          chksum2 += b[ind];
        }
      }
    }
  }
  //       cout << "before ("
  //            << lbnd[0] << ","
  //            << lbnd[1] << ","
  //            << lbnd[2] << "), ("
  //            << ubnd[0] << ","
  //            << ubnd[1] << ","
  //            << ubnd[2] << "): "
  //            << res << "   " << chksum << "," << chksum2 << endl;
  global_sum(cctkGH, &res);
  //      cout << "after: " << res << endl;
  realoutput = res;
}

} // namespace common

// Register this solver
int CarpetCG_register() {
  int const ierr = TATelliptic_RegisterSolver(CarpetCG_solve, "CarpetCG");
  assert(!ierr);
  return 0;
}

// Solve
int CarpetCG_solve(cGH const *restrict const cctkGH,
                   int const *restrict const var, int const *restrict const res,
                   int const nvars, int const options_table,
                   calcfunc const calcres, calcfunc const applybounds,
                   void *const userdata) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Domain descriptors
  cGroup groupdata;
  cGroupDynamicData groupdyndata;
  int gsh[dim], lsh[dim], nghostzones[dim], bbox[2 * dim];
  int lbnd[dim], ubnd[dim];

  // Options

  CCTK_INT maxiters;
  CCTK_INT maxiters2;
  CCTK_INT smaxiters;
  CCTK_REAL stepsize;
  CCTK_REAL maxstepsize;
  CCTK_REAL minerror;
  CCTK_REAL sminerror;

  int iter;
  int iter2;
  int siter;

  int do_abort;

  CCTK_REAL alpha; // secant step size
  CCTK_REAL alpha_old, sum_alpha;
  CCTK_REAL beta;    // CG step size
  CCTK_REAL delta_0; // initial A-norm of residual
  CCTK_REAL delta_new, delta_d, delta_old, delta_mid;
  CCTK_REAL epsilon; // norm of residual
  CCTK_REAL eta, eta_prev;

  time_t nexttime, currenttime;

  CCTK_REAL gsize; // global grid size (no. points as real)

  int n, nn;
  int d, f;

  int nelems;
  int ierr;

  // Security check
  if (!CCTK_IsThornActive(CCTK_THORNSTRING)) {
    CCTK_WARN(0, "Thorn " CCTK_THORNSTRING " has not been activated.  It is "
                 "therefore not possible to call "
                 "CarpetCG_solve.");
  }

  // Check arguments
  assert(cctkGH);

  assert(nvars > 0);
  assert(var);
  assert(res);
  for (n = 0; n < nvars; ++n) {
    assert(var[n] >= 0 && var[n] < CCTK_NumVars());
    assert(CCTK_GroupTypeFromVarI(var[n]) == CCTK_GF ||
           CCTK_GroupTypeFromVarI(var[n]) == CCTK_ARRAY);
    assert(CCTK_GroupDimFromVarI(var[n]) <= dim);
    assert(CCTK_VarTypeI(var[n]) == CCTK_VARIABLE_REAL);
    assert(CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(var[n])));
  }
  for (n = 0; n < nvars; ++n) {
    assert(res[n] >= 0 && res[n] < CCTK_NumVars());
    assert(CCTK_GroupTypeFromVarI(res[n]) == CCTK_GF ||
           CCTK_GroupTypeFromVarI(res[n]) == CCTK_ARRAY);
    assert(CCTK_GroupDimFromVarI(res[n]) <= dim);
    assert(CCTK_VarTypeI(res[n]) == CCTK_VARIABLE_REAL);
    assert(CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(res[n])));
  }
  for (n = 0; n < nvars; ++n) {
    assert(var[n] != res[n]);
    for (nn = 0; nn < n; ++nn) {
      assert(var[nn] != var[n]);
      assert(var[nn] != res[n]);
      assert(res[nn] != var[n]);
      assert(res[nn] != res[n]);
    }
  }

  assert(options_table >= 0);

  assert(calcres);
  assert(applybounds);

  nelems = Util_TableGetIntArray(
      options_table, 2 * dim, &(common::nboundaryzones[0]), "nboundaryzones");
  if (nelems == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (d = 0; d < dim; ++d) {
      common::nboundaryzones[2 * d] = -1;
      common::nboundaryzones[2 * d + 1] = -1;
    }
  } else if (nelems != 2 * dim) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Options table key \"nboundaryzones\" is not an integer array "
               "of length %d",
               2 * dim);
    return -1;
  }

  (common::fromindex).reserve(cg_maxsolvevars);
  (common::toindex).reserve(cg_maxsolvevars);
  (common::factors).reserve(nvars);
  common::factor = 0.5;
  for (int i = 0; i < nvars; i++) {
    common::factors[i] = 1.0;
  }

  maxiters = 1000;
  maxiters2 = 1000;
  minerror = 1.0e-8;
  smaxiters = 10;
  sminerror = 1.0e-4;
  stepsize = 1.0e-6;
  maxstepsize = 1.0;

  ierr = Util_TableGetInt(options_table, &maxiters, "maxiters");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetInt(options_table, &maxiters2, "maxiters2");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetInt(options_table, &smaxiters, "smaxiters");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &minerror, "minerror");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &sminerror, "sminerror");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &stepsize, "stepsize");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &maxstepsize, "maxstepsize");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);

  assert(maxiters >= 0);
  assert(maxiters2 >= 0);
  assert(smaxiters >= 0);
  assert(minerror >= 0);
  assert(sminerror >= 0);
  assert(stepsize > 0);

  // Switch to global mode
  BEGIN_GLOBAL_MODE(cctkGH) {

    // Fill common block
    common::var = var;
    common::res = res;
    common::nvars = nvars;
    common::options_table = options_table;
    common::calcres = calcres;
    common::applybnds = applybounds;
    common::userdata = userdata;

    if (verbose || veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Solving through nonlinear conjugate gradient iterations");
    }

    /*
     * Literature:
     *
     * Jonathan Richard Shewchuk, "An Introduction to the Conjugate
     * Gradient Method Without the Agonizing Pain", Technical Report
     * CMU-CS-94-125, Aug. 1994
     *
     * Available from http://www-2.cs.cmu.edu/~jrs/jrspapers.html
     */

    /*
     * Notation:
     *
     *    here        there   meaning
     *
     *    iter        i       current cg iteration
     *    siter       j       current secant iteration
     *    iter2       k       current cg iteration since last restart
     *    maxiters2   n       PARA cg iterations before restart
     *                alpha   secant step size
     *                beta    cg step size
     *                delta   cg distance
     *                epsilon PARA secant error tolerance (straight epsilon)
     *    minerror    vareps  PARA cg error tolerance (curved epsilon)
     *                eta     secant distance
     *    stepsize    sigma   PARA secand method step parameter
     *
     *    dir         d       cg direction
     *   -calcres     f'      PARA nonlinear operator
     *    res         r       cg residual
     *    prs         s       preconditioned cg residual
     *    var         x       INIT unknown variable
     *
     *   -wgt         M       preconditioning matrix (not in this version)
     */

    /*
     * Algorithm:
     * (Preconditioned nonlinear conjugate gradients with secant and
     * Polak-Ribière)
     *
     *    01. i <= 0
     *    02. k <= 0
     *    03. r <= - f'(x)
     *    04. calculate M \approx f''(x)
     *    05. s <= M^-1 r
     *    06. d <= s
     *    07. delta_new <= r^T d
     *    08. delta_0 <= delta_new
     *    09. WHILE i < i_max AND delta_new > vareps^2 delta_0 DO
     *    10.    j <= 0
     *    11.    delta_d <= d^T d
     *    12.    alpha <= - sigma_0
     *    13.    eta_prev <= [f'(x + sigma_0 d)]^T d
     *    14.    DO
     *    15.       eta <= [f'(x)]^T d
     *    16.       alpha <= alpha   (eta) / (eta_prev - eta)
     *    17.       x <= x + alpha d
     *    18.       eta_prev <= eta
     *    19.       j <= j + 1
     *    20.    WHILE j < j_max AND alpha^2 delta_d > epsilon^2
     *    21.    r <= - f'(x)
     *    22.    delta_old <= delta_new
     *    23.    delta_mid <= r^T s
     *    24.    calculate M \approx f''(x)
     *    25.    s <= M^-1 r
     *    26.    delta_new <= r^T s
     *    27.    beta <= (delta_new - delta_mid) / (delta_old)
     *    28.    k <= k + 1
     *    29.    IF k = n OR beta <= 0 THEN
     *    30.       d <= s
     *    31.       k <= 0
     *    32.    ELSE
     *    33.       d <= s + beta d
     *    34.    i <= i + 1
     */

    // Pointers to grid variables and temporary storage
    vector<int> varptrs(nvars);
    vector<int> resptrs(nvars);
    //    CCTK_REAL * restrict * restrict prsptrs;
    //    CCTK_REAL * restrict * restrict dirptrs;

    // Get storage pointers
    for (n = 0; n < nvars; ++n) {
      varptrs[n] = var[n];
      resptrs[n] = res[n];
    }

    // Allocate temporary memory
    //    prsptrs = malloc (nvars * sizeof *prsptrs);
    //    dirptrs = malloc (nvars * sizeof *prsptrs);
    assert(cg_maxsolvevars >= nvars);
    vector<int> prsptrs(cg_maxsolvevars);
    vector<int> dirptrs(cg_maxsolvevars);
    for (n = 0; n < nvars; ++n) {
      char varname[1000];
      sprintf(varname, "CarpetCG::prsvars[%d]", n);
      prsptrs[n] = CCTK_VarIndex(varname);
      assert((prsptrs[n] >= 0) && (prsptrs[n] < CCTK_NumVars()));
      sprintf(varname, "CarpetCG::dirvars[%d]", n);
      dirptrs[n] = CCTK_VarIndex(varname);
      assert((dirptrs[n] >= 0) && (dirptrs[n] < CCTK_NumVars()));
    }

    nexttime = time(0);

    /* 01. i <= 0 */
    iter = 0;

    /* 02. k <= 0 */
    iter2 = 0;

    /* 03. r <= - f'(x) */
    /* 04. calculate M \approx f''(x) */
    common::ierr = 0;
    CallLocalFunction((cGH *)cctkGH, common::call_calcres);
    ierr = common::ierr;
    assert(!ierr);

    for (int n = 0; n < nvars; n++) {
      common::fromindex[n] = resptrs[n];
    }
    for (int n = nvars; n < cg_maxsolvevars; n++) {
      common::fromindex[n] = -1;
    }
    CallLocalFunction((cGH *)cctkGH, common::call_correct_residual_sign);

    /* 05. s <= M^-1 r */
    // No preconditioning
    for (int n = 0; n < nvars; n++) {
      common::fromindex[n] = resptrs[n];
      common::toindex[n] = prsptrs[n];
    }
    for (int n = nvars; n < cg_maxsolvevars; n++) {
      common::fromindex[n] = -1;
      common::toindex[n] = -1;
    }

    //    cout << "dirptrs " << dirptrs[0] << endl;

    CallLocalFunction((cGH *)cctkGH, common::call_apply_preconditioner);

    //    cout << "dirptrs " << dirptrs[0] << endl;

    /* 06. d <= s */
    for (int n = 0; n < nvars; n++) {
      common::fromindex[n] = prsptrs[n];
      common::toindex[n] = dirptrs[n];
    }
    for (int n = nvars; n < cg_maxsolvevars; n++) {
      common::fromindex[n] = -1;
      common::toindex[n] = -1;
    }

    //    cout << "dirptrs " << dirptrs[0] << endl;

    CallLocalFunction((cGH *)cctkGH, common::call_copy);

    //    cout << "dirptrs " << dirptrs[0] << endl;

    /* 07. delta_new <= r^T d */
    //    output (cctkGH, var, res, wgt, nvars, iter);
    common::realoutput_count = 0;
    for (int n = 0; n < nvars; n++) {
      common::fromindex[n] = prsptrs[n];
      common::toindex[n] = resptrs[n];
    }
    for (int n = nvars; n < cg_maxsolvevars; n++) {
      common::fromindex[n] = -1;
      common::toindex[n] = -1;
    }

    //    cout << "dirptrs " << dirptrs[0] << endl;

    CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
    delta_new = common::realoutput;
    gsize = common::realoutput_count;
    //    cout << "delta_new " << delta_new << " gsize " << gsize << endl;

    for (int n = 0; n < nvars; n++) {
      common::fromindex[n] = resptrs[n];
      common::toindex[n] = resptrs[n];
    }
    for (int n = nvars; n < cg_maxsolvevars; n++) {
      common::fromindex[n] = -1;
      common::toindex[n] = -1;
    }

    CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
    epsilon = common::realoutput;
    //     cout << "epsilon " << epsilon << endl;

    /* 08. delta_0 <= delta_new */
    delta_0 = delta_new;

    /* 09. WHILE i < i_max AND delta_new > vareps^2 delta_0 DO */
    while (iter < maxiters && epsilon / (nvars * gsize) > ipow(minerror, 2)) {

      if (verbose || veryverbose) {
        currenttime = time(0);
        if (veryverbose || (iter % outevery == 0 && currenttime >= nexttime)) {
          CCTK_VInfo(
              CCTK_THORNSTRING,
              "Iteration %d (%d since restart): residual is %g (%g,%d,%g)",
              iter, iter2, (double)sqrt(epsilon / (nvars * gsize)),
              (double)epsilon, nvars, (double)gsize);
          if (outeveryseconds > 0) {
            while (nexttime <= currenttime)
              nexttime += outeveryseconds;
          }
        }
      }

      /* 10. j <= 0 */
      siter = 0;

      /* 11. delta_d <= d^T d */
      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = dirptrs[n];
        common::toindex[n] = dirptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }

      CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
      delta_d = common::realoutput;

      //       cout << "delta_d " << delta_d << endl;

      /* 12. alpha <= - sigma_0 */
      alpha = -stepsize;
      sum_alpha = alpha;
      do_abort = 0;

      /* 13. eta_prev <= [f'(x + sigma_0 d)]^T d */
      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = dirptrs[n];
        common::toindex[n] = varptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }
      common::realconstant = -alpha;
      //       cout << "Add scaled, const " << -alpha << endl;

      CallLocalFunction((cGH *)cctkGH, common::call_add_scaled);
      common::ierr = 0;
      CallLocalFunction((cGH *)cctkGH, common::call_applybnds);
      assert(!ierr);

      ierr = 0;
      CallLocalFunction((cGH *)cctkGH, common::call_calcres);
      assert(!ierr);

      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = resptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
      }
      CallLocalFunction((cGH *)cctkGH, common::call_correct_residual_sign);

      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = dirptrs[n];
        common::toindex[n] = resptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }
      CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
      eta = -common::realoutput;

      //       cout << "eta " << eta << endl;

      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = dirptrs[n];
        common::toindex[n] = varptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }
      common::realconstant = alpha;
      //       cout << "Add scaled, const " << alpha << endl;

      CallLocalFunction((cGH *)cctkGH, common::call_add_scaled);
      ierr = 0;
      CallLocalFunction((cGH *)cctkGH, common::call_applybnds);
      assert(!ierr);

      eta_prev = eta;

      /* 14. DO */
      do {

        if (veryverbose) {
          CCTK_VInfo(
              CCTK_THORNSTRING,
              "   Secant iteration %d: step size is %g, orthogonality is %g",
              siter, (double)alpha,
              (double)sqrt(fabs(eta_prev) / (nvars * gsize)));
        }

        /* 15. eta <= [f'(x)]^T d */
        ierr = 0;
        CallLocalFunction((cGH *)cctkGH, common::call_calcres);
        assert(!ierr);

        for (int n = 0; n < nvars; n++) {
          common::fromindex[n] = resptrs[n];
        }
        for (int n = nvars; n < cg_maxsolvevars; n++) {
          common::fromindex[n] = -1;
        }
        CallLocalFunction((cGH *)cctkGH, common::call_correct_residual_sign);

        for (int n = 0; n < nvars; n++) {
          common::fromindex[n] = dirptrs[n];
          common::toindex[n] = resptrs[n];
        }
        for (int n = nvars; n < cg_maxsolvevars; n++) {
          common::fromindex[n] = -1;
          common::toindex[n] = -1;
        }

        CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
        eta = -common::realoutput;

        //         cout << "eta " << eta << endl;

        /* 16. alpha <= alpha   (eta) / (eta_prev - eta) */
        alpha_old = alpha;
        alpha *= eta / (eta_prev - eta);
        sum_alpha += alpha;
        if (veryverbose) {
          CCTK_VInfo(
              CCTK_THORNSTRING,
              "   Changing step size, iteration %d: was %g, now %g (%g, %g)",
              siter, (double)alpha_old, (double)alpha, (double)eta,
              (double)eta_prev);
        }
        assert(sum_alpha > -maxstepsize);
        if (sum_alpha > maxstepsize) {
          alpha = maxstepsize - alpha_old;
          do_abort = 1;
          if (verbose) {
            CCTK_VInfo(CCTK_THORNSTRING,
                       "   Secant iteration %d: limiting total step size",
                       siter);
          }
        }

        /* 17. x <= x + alpha d */

        for (int n = 0; n < nvars; n++) {
          common::fromindex[n] = dirptrs[n];
          common::toindex[n] = varptrs[n];
        }
        for (int n = nvars; n < cg_maxsolvevars; n++) {
          common::fromindex[n] = -1;
          common::toindex[n] = -1;
        }
        common::realconstant = alpha;
        //         cout << "Add scaled, const " << alpha << endl;

        CallLocalFunction((cGH *)cctkGH, common::call_add_scaled);
        ierr = 0;
        CallLocalFunction((cGH *)cctkGH, common::call_applybnds);
        assert(!ierr);

        /* 18. eta_prev <= eta */
        eta_prev = eta;

        /* 19. j <= j + 1 */
        ++siter;

        /* 20. WHILE j < j_max AND alpha^2 delta_d > epsilon^2 */
      } while (siter < smaxiters &&
               ipow(alpha, 2) * delta_d > ipow(sminerror, 2) && !do_abort);

      if (veryverbose) {
        CCTK_VInfo(
            CCTK_THORNSTRING,
            "   Secant iteration %d: step size is %g, orthogonality is %g",
            siter, (double)alpha,
            (double)sqrt(fabs(eta_prev) / (nvars * gsize)));
      }

      /* 21. r <= - f'(x) */
      /* 24. calculate M \approx f''(x) */
      ierr = 0;
      CallLocalFunction((cGH *)cctkGH, common::call_calcres);
      assert(!ierr);

      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = resptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
      }
      CallLocalFunction((cGH *)cctkGH, common::call_correct_residual_sign);

      /* 23. delta_mid <= r^T s */
      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = prsptrs[n];
        common::toindex[n] = resptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }

      CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
      delta_mid = common::realoutput;
      //       cout << "delta_mid " << delta_mid << endl;

      /* 25. s <= M^-1 r */
      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = resptrs[n];
        common::toindex[n] = prsptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }

      CallLocalFunction((cGH *)cctkGH, common::call_apply_preconditioner);

      /* 22. delta_old <= delta_new */
      delta_old = delta_new;

      /* 26. delta_new <= r^T s */
      //       output (cctkGH, var, res, wgt, nvars, iter+1);
      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = prsptrs[n];
        common::toindex[n] = resptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }

      CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
      delta_new = common::realoutput;
      //       cout << "delta_new " << delta_new << endl;

      for (int n = 0; n < nvars; n++) {
        common::fromindex[n] = resptrs[n];
        common::toindex[n] = resptrs[n];
      }
      for (int n = nvars; n < cg_maxsolvevars; n++) {
        common::fromindex[n] = -1;
        common::toindex[n] = -1;
      }

      CallLocalFunction((cGH *)cctkGH, common::call_dot_product);
      epsilon = common::realoutput;
      //       cout << "epsilon " << epsilon << endl;

      /* 27. beta <= (delta_new - delta_mid) / (delta_old) */
      beta = (delta_new - delta_mid) / delta_old;

      /* 28. k <= k + 1 */
      ++iter2;

      /* 29. IF k = n OR beta <= 0 THEN */
      if (iter2 >= maxiters2 || beta <= 0 || do_abort) {
        /* Restart */

        /* 30. d <= s */
        for (int n = 0; n < nvars; n++) {
          common::fromindex[n] = prsptrs[n];
          common::toindex[n] = dirptrs[n];
        }
        for (int n = nvars; n < cg_maxsolvevars; n++) {
          common::fromindex[n] = -1;
          common::toindex[n] = -1;
        }
        CallLocalFunction((cGH *)cctkGH, common::call_copy);

        /* 31. k <= 0 */
        iter2 = 0;

        /* 32. ELSE */
      } else {
        /* Continue with CG iterations */

        /* 33. d <= s + beta d */
        for (int n = 0; n < nvars; n++) {
          common::fromindex[n] = prsptrs[n];
          common::toindex[n] = dirptrs[n];
        }
        for (int n = nvars; n < cg_maxsolvevars; n++) {
          common::fromindex[n] = -1;
          common::toindex[n] = -1;
        }
        common::realconstant = beta;
        //        cout << "beta " << beta << endl;

        CallLocalFunction((cGH *)cctkGH, common::call_scale_and_add);

      } /* if iter2 */

      /* 34. i <= i + 1 */
      ++iter;

    } /* for iter */

    if (verbose || veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Iteration %d (%d since restart): residual is %g", iter, iter2,
                 (double)sqrt(epsilon / (nvars * gsize)));
    }

    /* Free temporary memory */
    //     for (n=0; n<nvars; ++n) {
    //       free (prsptrs[n]);
    //       free (dirptrs[n]);
    //     }

    //     free (varptrs);
    //     free (resptrs);
    //     free (prsptrs);
    //     free (dirptrs);
  }
  END_GLOBAL_MODE;

  return 0;
}

} // namespace CarpetCG
