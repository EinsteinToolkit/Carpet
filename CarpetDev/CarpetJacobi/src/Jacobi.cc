// $Header: /home/eschnett/C/carpet/Carpet/CarpetDev/CarpetJacobi/src/Jacobi.cc,v 1.2 2003/11/19 14:05:14 schnetter Exp $

#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "carpet.hh"
#include "TATelliptic.h"



namespace Carpet {
  // TODO: fix this
  void Restrict (const cGH* cgh);
};



namespace CarpetJacobi {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  extern "C" {
    void CarpetJacobi_register ();
  }
  
  int CarpetJacobi_solve (const cGH * const cctkGH,
                          const int * const var,
                          const int * const res,
                          const int nvars,
                          const int options_table,
                          const calcfunc calcres,
                          const calcfunc applybnds,
                          void * const userdata);
  
  
  
  namespace common {
    
    const int * var;
    const int * res;
    int nvars;
    int options_table;
    calcfunc calcres;
    calcfunc applybnds;
    void * userdata;
    
    CCTK_REAL factor;
    CCTK_REAL * factors;
    
    CCTK_REAL norm_count;
    CCTK_REAL norm_l2;
    
    int ierr;
    
    void call_calcres (cGH * const cctkGH)
    {
      if (ierr) return;
      ierr = calcres (cctkGH, options_table, userdata);
    }
    
    void norm (cGH * const cctkGH)
    {
      int ierr;
      
      cGroup groupdata;
      ierr = CCTK_GroupData (CCTK_GroupIndexFromVarI(var[0]), &groupdata);
      assert (!ierr);
      cGroupDynamicData groupdyndata;
      ierr = CCTK_GroupDynamicData (cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                                    &groupdyndata);
      assert (!ierr);
      const int thedim = groupdata.dim;
      assert (thedim>=0 && thedim<=dim);
      assert (thedim == groupdyndata.dim);
      int lsh[dim], nghostzones[dim];
      for (int d=0; d<thedim; ++d) {
        lsh[d] = groupdyndata.lsh[d];
        nghostzones[d] = groupdyndata.nghostzones[d];
      }
      for (int d=thedim; d<dim; ++d) {
        lsh[d] = 1;
        nghostzones[d] = 0;
      }
      for (int d=0; d<dim; ++d) {
        assert (lsh[d]>=0);
        assert (nghostzones[d] >= 0 && 2*nghostzones[d] <= lsh[d]);
      }
      
      for (int n=0; n<nvars; ++n) {
        const CCTK_REAL * restrict resptr
          = (const CCTK_REAL *)CCTK_VarDataPtrI (cctkGH, 0, res[n]);
        assert (resptr);
        const CCTK_REAL fac = factors[n];
        for (int k=nghostzones[2]; k<lsh[2]-nghostzones[2]; ++k) {
          for (int j=nghostzones[1]; j<lsh[1]-nghostzones[1]; ++j) {
            for (int i=nghostzones[0]; i<lsh[0]-nghostzones[0]; ++i) {
              const int ind = i + lsh[0] * (j + lsh[1] * k);
              ++norm_count;
              // TODO: scale the norm by the resolution?
              norm_l2 += pow(fac * resptr[ind], 2);
            }
          }
        }
      }
    }
    
    void smooth (cGH * const cctkGH)
    {
      int ierr;
      
      cGroup groupdata;
      ierr = CCTK_GroupData (CCTK_GroupIndexFromVarI(var[0]), &groupdata);
      assert (!ierr);
      cGroupDynamicData groupdyndata;
      ierr = CCTK_GroupDynamicData (cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                                    &groupdyndata);
      assert (!ierr);
      const int thedim = groupdata.dim;
      assert (thedim>=0 && thedim<=dim);
      assert (thedim == groupdyndata.dim);
      int lsh[dim], nghostzones[dim];
      for (int d=0; d<thedim; ++d) {
        lsh[d] = groupdyndata.lsh[d];
        nghostzones[d] = groupdyndata.nghostzones[d];
      }
      for (int d=thedim; d<dim; ++d) {
        lsh[d] = 1;
        nghostzones[d] = 0;
      }
      for (int d=0; d<dim; ++d) {
        assert (lsh[d]>=0);
        assert (nghostzones[d] >= 0 && 2*nghostzones[d] <= lsh[d]);
      }
      
      for (int n=0; n<nvars; ++n) {
        CCTK_REAL * restrict varptr
          = (CCTK_REAL *)CCTK_VarDataPtrI (cctkGH, 0, var[n]);
        assert (varptr);
        const CCTK_REAL * restrict resptr
          = (const CCTK_REAL *)CCTK_VarDataPtrI (cctkGH, 0, res[n]);
        assert (resptr);
        assert (resptr != varptr);
        const CCTK_REAL fac = factor * factors[n];
        for (int k=nghostzones[2]; k<lsh[2]-nghostzones[2]; ++k) {
          for (int j=nghostzones[1]; j<lsh[1]-nghostzones[1]; ++j) {
            for (int i=nghostzones[0]; i<lsh[0]-nghostzones[0]; ++i) {
              const int ind = i + lsh[0] * (j + lsh[1] * k);
              varptr[ind] += fac * resptr[ind];
            }
          }
        }
      }
    }
    
    void call_applybnds (cGH * const cctkGH)
    {
      if (ierr) return;
      ierr = applybnds (cctkGH, options_table, userdata);
    }
    
  } // namespace common
  
  
  
  void CarpetJacobi_register ()
  {
    int const ierr = TATelliptic_RegisterSolver
      (CarpetJacobi_solve, "CarpetJacobi");
    assert (!ierr);
  }
  
  
  
  int CarpetJacobi_solve (const cGH * const cctkGH,
                          const int * const var,
                          const int * const res,
                          const int nvars,
                          const int options_table,
                          const calcfunc calcres,
                          const calcfunc applybnds,
                          void * const userdata)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Error control
    int ierr;
    
    if (! CCTK_IsThornActive(CCTK_THORNSTRING)) {
      CCTK_WARN (0, "Thorn " CCTK_THORNSTRING " has not been activated.  It is therefore not possible to call CarpetJacobi_solve.");
    }
    
    // Check arguments
    assert (cctkGH);
    
    assert (nvars > 0);
    assert (var);
    assert (res);
    for (int n=0; n<nvars; ++n) {
      assert (var[n] >= 0 && var[n] < CCTK_NumVars());
      assert (CCTK_GroupTypeFromVarI(var[n]) == CCTK_GF
              || CCTK_GroupTypeFromVarI(var[n]) == CCTK_ARRAY);
      assert (CCTK_GroupDimFromVarI(var[n]) <= dim);
      assert (CCTK_VarTypeI(var[n]) == CCTK_VARIABLE_REAL);
      assert (CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(var[n])));
    }
    for (int n=0; n<nvars; ++n) {
      assert (res[n] >= 0 && res[n] < CCTK_NumVars());
      assert (CCTK_GroupTypeFromVarI(res[n]) == CCTK_GF
              || CCTK_GroupTypeFromVarI(res[n]) == CCTK_ARRAY);
      assert (CCTK_GroupDimFromVarI(res[n]) <= dim);
      assert (CCTK_VarTypeI(res[n]) == CCTK_VARIABLE_REAL);
      assert (CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(res[n])));
    }
    for (int n=0; n<nvars; ++n) {
      assert (var[n] != res[n]);
      for (int nn=0; nn<n; ++nn) {
        assert (var[nn] != var[n]);
        assert (var[nn] != res[n]);
        assert (res[nn] != var[n]);
        assert (res[nn] != res[n]);
      }
    }
    
#if 0
    // Get domain description
    cGroup groupdata;
    ierr = CCTK_GroupData (CCTK_GroupIndexFromVarI(var[0]), &groupdata);
    assert (!ierr);
    cGroupDynamicData groupdyndata;
    ierr = CCTK_GroupDynamicData (cctkGH, CCTK_GroupIndexFromVarI(var[0]),
                                  &groupdyndata);
    assert (!ierr);
    const int thedim = groupdata.dim;
    assert (thedim>=0 && thedim<=dim);
    assert (thedim == groupdyndata.dim);
    int lsh[dim], nghostzones[dim];
    for (int d=0; d<thedim; ++d) {
      lsh[d] = groupdyndata.lsh[d];
      nghostzones[d] = groupdyndata.nghostzones[d];
    }
    for (int d=thedim; d<dim; ++d) {
      lsh[d] = 1;
      nghostzones[d] = 0;
    }
    for (int d=0; d<dim; ++d) {
      assert (lsh[d]>=0);
      assert (nghostzones[d] >= 0 && 2*nghostzones[d] <= lsh[d]);
    }
    
    // Check all variables
    for (int n=0; n<nvars; ++n) {
      ierr = CCTK_GroupData (CCTK_GroupIndexFromVarI(var[n]), &groupdata);
      assert (!ierr);
      ierr = CCTK_GroupDynamicData
        (cctkGH, CCTK_GroupIndexFromVarI(var[n]), &groupdyndata);
      assert (!ierr);
      assert (groupdata.dim == thedim);
      assert (groupdyndata.dim == thedim);
      for (int d=0; d<thedim; ++d) {
        assert (groupdyndata.lsh[d] == lsh[d]);
        assert (groupdyndata.nghostzones[d] == nghostzones[d]);
      }
    }
#endif
    
    assert (options_table >= 0);
    
    assert (calcres);
    assert (applybnds);
    
    // Get options from options table
    CCTK_INT maxiters = 10000;
    CCTK_REAL minerror = 1e-8;
    CCTK_REAL factor = 1.0e-2;    // slow start
    CCTK_REAL minfactor = 1.0e-8;
    CCTK_REAL maxfactor = 1.0;
    CCTK_REAL acceleration = 1.2; // slow acceleration
    CCTK_REAL deceleration = 0.1; // emergency brake
    vector<CCTK_REAL> factors (nvars);
    for (int n=0; n<nvars; ++n) {
      factors[n] = 1.0;
    }
    
    ierr = Util_TableGetInt (options_table, &maxiters, "maxiters");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetReal (options_table, &minerror, "minerror");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetReal (options_table, &factor, "factor");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetReal (options_table, &minfactor, "minfactor");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetReal (options_table, &maxfactor, "maxfactor");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetReal (options_table, &acceleration, "acceleration");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetReal (options_table, &deceleration, "deceleration");
    assert (ierr==1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    ierr = Util_TableGetRealArray (options_table, nvars, &factors[0], "factors");
    assert (ierr==nvars || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    
    assert (maxiters >= 0);
    assert (minerror >= 0);
    assert (factor > 0);
    assert (minfactor > 0 && minfactor <= factor);
    assert (maxfactor >= factor);
    assert (acceleration >= 1);
    assert (deceleration > 0 && deceleration <= 1);
    for (int n=0; n<nvars; ++n) {
      assert (factors[n] != 0);
    }
    
    // Check state
    // (want global mode, or at least level mode)
    assert (reflevel==-1 || component==-1);
    const int saved_reflevel = reflevel;
    const int saved_mglevel = mglevel;
    if (reflevel!=-1) {
      set_mglevel ((cGH *)cctkGH, -1);
      set_reflevel ((cGH *)cctkGH, -1);
    }
    
    // Fill common block
    common::var = var;
    common::res = res;
    common::nvars = nvars;
    common::options_table = options_table;
    common::calcres = calcres;
    common::applybnds = applybnds;
    common::userdata = userdata;
    common::factor = factor;
    common::factors = &factors[0];
    
    // Init statistics
    CCTK_REAL norm2 = HUGE_VAL;
    CCTK_REAL speed = 0.0;
    CCTK_REAL throttle = 1.0;
    bool did_hit_min = false;
    bool did_hit_max = false;
    const time_t starttime = time(0);
    time_t nexttime = starttime;
    
    if (verbose || veryverbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "Solving through adaptive Jacobi iterations");
    }
    
    int iter;
    for (iter=0; iter<=maxiters; ++iter) {
      
      // Calculate residual
      common::ierr = 0;
      CallLocalFunction ((cGH *)cctkGH, common::call_calcres);
      ierr = common::ierr;
      if (ierr != 0) {
        if (verbose || veryverbose) {
          CCTK_VInfo (CCTK_THORNSTRING, "Residual calculation reported error; aborting solver.");
        }
        Util_TableSetInt (options_table, iter, "iters");
        Util_TableDeleteKey (options_table, "error");
        goto done;
      }
      
      // Save old norm
      const CCTK_REAL old_norm2 = norm2;
      
      // Calculate norm
      common::norm_count = 0;
      common::norm_l2 = 0;
      CallLocalFunction ((cGH *)cctkGH, common::norm);
      const int sum_handle = CCTK_ReductionArrayHandle ("sum");
      assert (sum_handle >= 0);
      CCTK_REAL reduce_in[2], reduce_out[2];
      reduce_in[0] = common::norm_count;
      reduce_in[1] = common::norm_l2;
      ierr = CCTK_ReduceLocArrayToArray1D (cctkGH, -1, sum_handle,
                                           reduce_in, reduce_out, 2,
                                           CCTK_VARIABLE_REAL);
      norm2 = sqrt(reduce_out[1] / reduce_out[0]);
      
      // Calculate convergence speed
      speed = old_norm2 / norm2;
      
      // Log output
      if (verbose || veryverbose) {
        const time_t currenttime = time(0);
        if (iter == maxiters
#if HAVE_ISNAN
            || isnan(norm2)
#endif
            || norm2 <= minerror
            || (iter % outevery == 0 && currenttime >= nexttime)) {
          if (verbose || veryverbose) {
            if (did_hit_min) {
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Convergence factor became too small: artificially kept at minfactor (may lead to instability)");
            }
            if (did_hit_max) {
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Convergence factor became too large: artificially kept at maxfactor (may lead to slower convergence)");
            }
            did_hit_min = false;
            did_hit_max = false;
          }
          if (veryverbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Iteration %d, time %g, factor %g, throttle %g, residual %g, speed %g",
                        iter, (double)(currenttime-starttime), (double)factor, (double)throttle, (double)norm2, (double)speed);
          } else {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Iteration %d, time %g, residual %g",
                        iter, (double)(currenttime-starttime), (double)norm2);
          }
          if (outeveryseconds > 0) {
            while (nexttime <= currenttime) nexttime += outeveryseconds;
          }
        }
      }
      
      if (iter == maxiters) break;
      
#if HAVE_ISNAN
      // Exit if things go bad
      if (isnan(norm2)) {
        if (verbose || veryverbose) {
          CCTK_VInfo (CCTK_THORNSTRING, "Encountered NaN.  Aborting solver.");
        }
        break;
      }
#endif
      
      // Return when desired accuracy is reached
      if (norm2 <= minerror) {
        if (verbose || veryverbose) {
          CCTK_VInfo (CCTK_THORNSTRING, "Done solving.");
        }
        Util_TableSetInt (options_table, iter, "iters");
        Util_TableSetReal (options_table, norm2, "error");
        ierr = 0;
        goto done;
      }
      
      // Adjust speed
      if (speed > 1.0) {
        throttle = acceleration;
      } else {
        throttle = deceleration;
      }
      factor *= throttle;
      if (factor < minfactor) {
        did_hit_min = true;
        factor = minfactor;
      }
      if (factor > maxfactor) {
        did_hit_max = true;
        factor = maxfactor;
      }
      
      // Smooth
      CallLocalFunction ((cGH *)cctkGH, common::smooth);
      
      // Apply boundary conditions
      common::ierr = 0;
      CallLocalFunction ((cGH *)cctkGH, common::call_applybnds);
      ierr = common::ierr;
      if (ierr != 0) {
        if (verbose || veryverbose) {
          CCTK_VInfo (CCTK_THORNSTRING, "Boundary conditions reported error; aborting solver.");
        }
        Util_TableSetInt (options_table, iter, "iters");
        Util_TableDeleteKey (options_table, "error");
        goto done;
      }
      
#if 0
      // Restrict
      BEGIN_REVERSE_REFLEVEL_LOOP(cctkGH) {
        BEGIN_MGLEVEL_LOOP(cctkGH) {
          Restrict (cctkGH);
          // TODO: do something here
//           CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", (cGH *)cctkGH, CallFunction);
          CallLocalFunction ((cGH *)cctkGH, common::call_applybnds);
          assert (! common::ierr);
        } END_MGLEVEL_LOOP;
      } END_REVERSE_REFLEVEL_LOOP;
#endif
      
    } // for iter
    
    // Did not solve
    if (verbose || veryverbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "Did not converge.");
    }
    Util_TableSetInt (options_table, iter, "iters");
    Util_TableSetReal (options_table, norm2, "error");
    ierr = -1;
    
  done:
    
    // Restore state
    if (reflevel!=saved_reflevel) {
      set_reflevel ((cGH *)cctkGH, saved_reflevel);
      set_mglevel ((cGH *)cctkGH, saved_mglevel);
    }
    
    return ierr;
  }

} // namespace CarpetJacobi
