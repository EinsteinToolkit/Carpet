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

#include "th.hh"

namespace Carpet {
// TODO: fix this
void CycleTimeLevels(cGH *cctkGH);
void Restrict(const cGH *cctkGH);
};

#define DEBUG 0 // either 0 or 1

namespace CarpetJacobi {

using namespace std;
using namespace Carpet;

extern "C" {
int CarpetJacobi_register();
}

int CarpetJacobi_solve(const cGH *const cctkGH, const int *const var,
                       const int *const res, const int nvars,
                       const int options_table, const calcfunc calcres,
                       const calcfunc applybnds, void *const userdata);

int CarpetJacobi_register() {
  int const ierr =
      TATelliptic_RegisterSolver(CarpetJacobi_solve, "CarpetJacobi");
  assert(!ierr);
  return 0;
}

int CarpetJacobi_solve(const cGH *const cctkGH, const int *const var,
                       const int *const res, const int nvars,
                       const int options_table, const calcfunc calcres,
                       const calcfunc applybnds, void *const userdata) {
  DECLARE_CCTK_PARAMETERS;

  // Error control
  int ierr;

  if (!CCTK_IsThornActive(CCTK_THORNSTRING)) {
    CCTK_WARN(0, "Thorn " CCTK_THORNSTRING " has not been activated.  It is "
                 "therefore not possible to call "
                 "CarpetJacobi_solve.");
  }

  // Check arguments
  assert(cctkGH);

  assert(nvars > 0);
  assert(var);
  assert(res);
  for (int n = 0; n < nvars; ++n) {
    assert(var[n] >= 0 && var[n] < CCTK_NumVars());
    assert(CCTK_GroupTypeFromVarI(var[n]) == CCTK_GF);
    assert(CCTK_GroupDimFromVarI(var[n]) == dim);
    assert(CCTK_VarTypeI(var[n]) == CCTK_VARIABLE_REAL);
    assert(CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(var[n])));
  }
  for (int n = 0; n < nvars; ++n) {
    assert(res[n] >= 0 && res[n] < CCTK_NumVars());
    assert(CCTK_GroupTypeFromVarI(res[n]) == CCTK_GF);
    assert(CCTK_GroupDimFromVarI(res[n]) == dim);
    assert(CCTK_VarTypeI(res[n]) == CCTK_VARIABLE_REAL);
    assert(CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(res[n])));
  }
  for (int n = 0; n < nvars; ++n) {
    assert(var[n] != res[n]);
    for (int nn = 0; nn < n; ++nn) {
      assert(var[nn] != var[n]);
      assert(var[nn] != res[n]);
      assert(res[nn] != var[n]);
      assert(res[nn] != res[n]);
    }
  }

  assert(options_table >= 0);

  assert(calcres);
  assert(applybnds);

  // Get options from options table
  CCTK_INT maxiters = 10000;
  CCTK_REAL minerror = 1e-8;
  CCTK_REAL factor = 1.0e-2; // slow start
  CCTK_REAL minfactor = 1.0e-8;
  CCTK_REAL maxfactor = 1.0;
  CCTK_REAL acceleration = 1.2; // slow acceleration
  CCTK_REAL deceleration = 0.1; // emergency brake
  vector<CCTK_REAL> factors(nvars);
  for (int n = 0; n < nvars; ++n) {
    factors[n] = 1.0;
  }

  ierr = Util_TableGetInt(options_table, &maxiters, "maxiters");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &minerror, "minerror");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &factor, "factor");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &minfactor, "minfactor");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &maxfactor, "maxfactor");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &acceleration, "acceleration");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetReal(options_table, &deceleration, "deceleration");
  assert(ierr == 1 || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  ierr = Util_TableGetRealArray(options_table, nvars, &factors[0], "factors");
  assert(ierr == nvars || ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY);

  assert(maxiters >= 0);
  assert(minerror >= 0);
  assert(factor > 0);
  assert(minfactor > 0 && minfactor <= factor);
  assert(maxfactor >= factor);
  assert(acceleration >= 1);
  assert(deceleration > 0 && deceleration <= 1);
  for (int n = 0; n < nvars; ++n) {
    assert(factors[n] != 0);
  }

  assert(is_global_mode() || is_level_mode());

  // The level to solve on
  const int solve_level = is_level_mode() ? reflevel : reflevels - 1;

  if (solve_level < solve_minlevel) {
    if (verbose || veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Did not solve because the level is too coarse.");
    }
    Util_TableSetInt(options_table, 0, "iters");
    Util_TableSetReal(options_table, -1.0, "error");

    // Return as if the solving had not converged
    return -1;
  }

  // TODO: assert that all levels are at the same time

  // Switch to global mode
  BEGIN_GLOBAL_MODE(cctkGH) {

    // Reset the initial data time
    global_time = 0;
    delta_time = 1;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = global_time;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = delta_time;
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = global_time;
      for (int tl = 0; tl < timelevels; ++tl) {
        tt->set_time(mglevel, reflevel, tl, global_time - tl * delta_time);
      }
      tt->set_delta(mglevel, reflevel,
                    1.0 / ipow(maxval(spacereflevelfact), reflevelpower));
    }
    END_REFLEVEL_LOOP;

    // Init statistics
    vector<CCTK_REAL> norm_counts(solve_level + 1);
    vector<CCTK_REAL> norm_l2s(solve_level + 1);
    CCTK_REAL norm2 = HUGE_VAL;
    CCTK_REAL speed = 0.0;
    CCTK_REAL throttle = 1.0;
    bool did_hit_min = false;
    bool did_hit_max = false;
    const time_t starttime = time(0);
    time_t nexttime = starttime;

    if (verbose || veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Solving through adaptive Jacobi iterations");
    }

    const int icoarsestep =
        ipow(mglevelfact * maxval(maxspacereflevelfact), reflevelpower);
    const int istep = ipow(mglevelfact * maxval(maxspacereflevelfact /
                                                spacereffacts.at(solve_level)),
                           reflevelpower);
    int iter = istep;
    for (;; iter += istep) {

      global_time =
          1.0 * iter / ipow(maxval(maxspacereflevelfact), reflevelpower);
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = global_time;
      if (DEBUG)
        cout << "CJ global iter " << iter << " time " << global_time << flush
             << endl;

      // Calculate residual, smooth, and apply boundary conditions
      ierr = 0;
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        const int do_every =
            ipow(mglevelfact * maxval(maxspacereflevelfact / spacereflevelfact),
                 reflevelpower);
        if (reflevel <= solve_level && (iter - istep) % do_every == 0) {

          // Advance time
          tt->advance_time(mglevel, reflevel);
          *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) =
              (1.0 * (iter - istep + do_every) /
               ipow(maxval(maxspacereflevelfact), reflevelpower));
          CycleTimeLevels(const_cast<cGH *>(cctkGH));
          if (DEBUG)
            cout << "CJ residual iter " << iter << " reflevel " << reflevel
                 << " time " << cctkGH->cctk_time << flush << endl;

          // Advance time levels
          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              for (int n = 0; n < nvars; ++n) {
                if (CCTK_ActiveTimeLevelsVI(cctkGH, var[n]) > 1) {
                  size_t npoints = 1;
                  for (int d = 0; d < cctkGH->cctk_dim; ++d) {
                    npoints *= cctkGH->cctk_lsh[d];
                  }
                  memcpy(CCTK_VarDataPtrI(cctkGH, 0, var[n]),
                         CCTK_VarDataPtrI(cctkGH, 1, var[n]),
                         npoints * sizeof(CCTK_REAL));
                }
              } // for n
            }
            END_LOCAL_COMPONENT_LOOP;
          }
          END_MAP_LOOP;

          // Calculate residual
          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              if (!ierr) {
                ierr = calcres(cctkGH, options_table, userdata);
              }
            }
            END_LOCAL_COMPONENT_LOOP;
          }
          END_MAP_LOOP;

          // Smooth and apply boundary conditions
          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {

            CCTK_REAL levfac = 0;
            for (int d = 0; d < dim; ++d) {
              CCTK_REAL const delta =
                  cctkGH->cctk_delta_space[d] / cctkGH->cctk_levfac[d];
              levfac += 1 / ipow(delta, 2);
            }
            levfac = 1 / levfac;

            BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              if (!ierr) {

                for (int n = 0; n < nvars; ++n) {
                  CCTK_REAL *restrict varptr =
                      (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, var[n]);
                  assert(varptr);
                  const CCTK_REAL *restrict resptr =
                      (const CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, res[n]);
                  assert(resptr);
                  assert(resptr != varptr);
                  const CCTK_REAL fac = factor * factors[n] * levfac;
                  if (DEBUG)
                    cout << "CJ    smoothing fac=" << fac << endl;
                  assert(cctkGH->cctk_dim == 3);
                  for (int k = cctkGH->cctk_nghostzones[2];
                       k < cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];
                       ++k) {
                    for (int j = cctkGH->cctk_nghostzones[1];
                         j < cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
                         ++j) {
                      for (int i = cctkGH->cctk_nghostzones[0];
                           i <
                           cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
                           ++i) {
                        const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        varptr[ind] += fac * resptr[ind];
                      }
                    }
                  }
                } // for n

#if 1
                // Apply boundary conditions
                ierr = applybnds(cctkGH, options_table, userdata);
#endif

              } // if ! ierr
            }
            END_LOCAL_COMPONENT_LOOP;
          }
          END_MAP_LOOP;

#if 0
            // Apply boundary conditions
            ierr = applybnds (cctkGH, options_table, userdata);
            int const ierr2 = CallScheduleGroup (cctkGH, "ApplyBCs");

#endif

        } // if do_every
      }
      END_REFLEVEL_LOOP;

      if (ierr) {
        if (verbose || veryverbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "Residual calculation or boundary "
                                       "conditions reported error; aborting "
                                       "solver.");
        }
        Util_TableSetInt(options_table, iter, "iters");
        Util_TableDeleteKey(options_table, "error");
        goto done;
      }

      // Restrict and calculate norm
      ierr = 0;
      BEGIN_REVERSE_REFLEVEL_LOOP(cctkGH) {
        const int do_every =
            ipow(mglevelfact * maxval(maxspacereflevelfact / spacereflevelfact),
                 reflevelpower);
        if (reflevel <= solve_level && iter % do_every == 0) {

          if (DEBUG)
            cout << "CJ restrict iter " << iter << " reflevel " << reflevel
                 << " time " << cctkGH->cctk_time << flush << endl;

          if (!ierr) {
            if (reflevel < solve_level) {
              Restrict(cctkGH);
            }
          }

          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              if (!ierr) {
                ierr = applybnds(cctkGH, options_table, userdata);
              }
            }
            END_LOCAL_COMPONENT_LOOP;
          }
          END_MAP_LOOP;

          // Initialise norm
          CCTK_REAL norm_count = 0;
          CCTK_REAL norm_l2 = 0;

          BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              if (!ierr) {

                // Calculate norm
                for (int n = 0; n < nvars; ++n) {
                  const CCTK_REAL *restrict resptr =
                      (const CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, res[n]);
                  assert(resptr);
                  const CCTK_REAL fac = factors[n];
                  if (DEBUG)
                    cout << "CJ    norm fac=" << fac << endl;
                  assert(cctkGH->cctk_dim == 3);
                  for (int k = cctkGH->cctk_nghostzones[2];
                       k < cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];
                       ++k) {
                    for (int j = cctkGH->cctk_nghostzones[1];
                         j < cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
                         ++j) {
                      for (int i = cctkGH->cctk_nghostzones[0];
                           i <
                           cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
                           ++i) {
                        const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        ++norm_count;
                        // TODO: scale the norm by the resolution?
                        norm_l2 += ipow(fac * resptr[ind], 2);
                      }
                    }
                  }
                  if (DEBUG)
                    cout << "CJ    norm=" << norm_l2 << endl;
                } // for n
              }
            }
            END_LOCAL_COMPONENT_LOOP;
          }
          END_MAP_LOOP;

          const int sum_handle = CCTK_ReductionArrayHandle("sum");
          assert(sum_handle >= 0);
          CCTK_REAL reduce_in[2], reduce_out[2];
          reduce_in[0] = norm_count;
          reduce_in[1] = norm_l2;
          ierr =
              CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum_handle, reduce_in,
                                           reduce_out, 2, CCTK_VARIABLE_REAL);
          norm_counts.at(reflevel) = reduce_out[0];
          norm_l2s.at(reflevel) = reduce_out[1];

// TODO
#if 0
            CCTK_OutputVarAs (cctkGH, "WaveToyFO::phi", "phi-ell");
            CCTK_OutputVarAs (cctkGH, "IDSWTEsimpleFO::residual", "residual-ell");
#endif
        }
      }
      END_REVERSE_REFLEVEL_LOOP;

      if (ierr) {
        if (verbose || veryverbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "Residual calculation or boundary "
                                       "conditions reported error; aborting "
                                       "solver.");
        }
        Util_TableSetInt(options_table, iter, "iters");
        Util_TableDeleteKey(options_table, "error");
        goto done;
      }

      // Save old norm
      const CCTK_REAL old_norm2 = norm2;

      // Calculate new norm
      CCTK_REAL norm_count = 0;
      CCTK_REAL norm_l2 = 0;
      for (int rl = 0; rl <= solve_level; ++rl) {
        norm_count += norm_counts[rl];
        norm_l2 += norm_l2s[rl];
      }
      norm2 = sqrt(norm_l2 / norm_count);

      // Calculate convergence speed
      speed = old_norm2 / norm2;

      // Log output
      if (verbose || veryverbose) {
        const time_t currenttime = time(0);
        if ((iter % icoarsestep == 0 && iter >= maxiters) || isnan(norm2) ||
            norm2 <= minerror ||
            (iter % outevery == 0 && currenttime >= nexttime)) {
          if (verbose || veryverbose) {
            if (did_hit_min) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Convergence factor became too small: artificially "
                         "kept at minfactor (may lead to instability)");
            }
            if (did_hit_max) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Convergence factor became too large: artificially "
                         "kept at maxfactor (may lead to slower convergence)");
            }
            did_hit_min = false;
            did_hit_max = false;
          }
          if (veryverbose) {
            CCTK_VInfo(CCTK_THORNSTRING, "Iteration %d, time %g, factor %g, "
                                         "throttle %g, residual %g, speed %g",
                       iter, double(currenttime - starttime), double(factor),
                       double(throttle), double(norm2), double(speed));
          } else {
            CCTK_VInfo(CCTK_THORNSTRING, "Iteration %d, time %g, residual %g",
                       iter, double(currenttime - starttime), double(norm2));
          }
          if (outeveryseconds > 0) {
            while (nexttime <= currenttime)
              nexttime += outeveryseconds;
          }
        }
      }

      // Exit if things go bad
      if (isnan(norm2)) {
        if (verbose || veryverbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "Encountered NaN.  Aborting solver.");
        }
        break;
      }

      // Return when desired accuracy is reached
      if (norm2 <= minerror) {
        if (verbose || veryverbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "Done solving.");
        }
        Util_TableSetInt(options_table, iter, "iters");
        Util_TableSetReal(options_table, norm2, "error");
        ierr = 0;
        goto done;
      }

      // Exit if the maximum number of iterations has been reached
      if (iter % icoarsestep == 0 && iter >= maxiters) {
        break;
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

    } // for iter

    // Did not solve
    if (verbose || veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Did not converge.");
    }
    Util_TableSetInt(options_table, iter, "iters");
    Util_TableSetReal(options_table, norm2, "error");
    ierr = -1;

  done:;

    // Reset the initial time
    // TODO: reset the initial time a bit more intelligently
    global_time = 0;
    delta_time = 1;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = global_time;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = delta_time;
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = global_time;
      for (int tl = 0; tl < timelevels; ++tl) {
        tt->set_time(mglevel, reflevel, tl, global_time - tl * delta_time);
      }
      tt->set_delta(mglevel, reflevel, 1.0 / timereflevelfact);
    }
    END_REFLEVEL_LOOP;
  }
  END_GLOBAL_MODE;

  return ierr;
}

} // namespace CarpetJacobi
