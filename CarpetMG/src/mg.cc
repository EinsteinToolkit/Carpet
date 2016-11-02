#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "TATelliptic.h"
#include "carpet.hh"

namespace CarpetMG {

using namespace std;
using namespace Carpet;

char const *const solver_name = "CarpetMG";

struct options_t {

  // Width of outer boundary
  CCTK_INT bndwidth;

  // Symmetry information
  CCTK_INT symmetry_handle[2 * dim];
  CCTK_INT symmetry_zone_width[2 * dim];

  // Coarsest level on which the full multigrid algorithm should start
  CCTK_INT firstlevel;

  // Direct solver
  CCTK_INT maxiters;    // Maximum number of iterations
  CCTK_REAL sor_factor; // SOR factor

  // Weight factors for the equations
  CCTK_REAL factor;
  vector<CCTK_REAL> factors;

  // Multigrid iterations
  CCTK_INT mgiters;   // maximum number of iterations (per level)
  CCTK_INT presteps;  // number of presmoothing steps
  CCTK_INT poststeps; // number of postsmoothing steps
};

extern "C" {

int CarpetMG_register();

void CarpetMG_paramcheck(CCTK_ARGUMENTS);
}

int solve(cGH const *restrict const cctkGH, int const *restrict const var,
          int const *restrict const res, int const nvar,
          int const options_table, calcfunc const calcres,
          calcfunc const applybnds, void *const userdata);
int multigrid(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
              vector<CCTK_INT> const &res, vector<CCTK_INT> const &rhs,
              vector<CCTK_INT> const &sav, vector<CCTK_INT> const &wgt,
              vector<CCTK_INT> const &aux, int const options_table,
              calcfunc const calcres, calcfunc const applybnds,
              void *const userdata, options_t const &options,
              CCTK_REAL const minerror);

int direct_solve(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
                 vector<CCTK_INT> const &res, vector<CCTK_INT> const &rhs,
                 vector<CCTK_INT> const &wgt, vector<CCTK_INT> const &aux,
                 int const options_table, calcfunc const calcres,
                 calcfunc const applybnds, void *const userdata,
                 options_t const &options, CCTK_REAL const minerror);

void smooth(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
            vector<CCTK_INT> const &res, vector<CCTK_INT> const &rhs,
            vector<CCTK_INT> const &wgt, options_t const &options,
            CCTK_REAL const old_error, CCTK_REAL &error,
            bool const use_sor = false);

void norm(cGH const *restrict const cctkGH, vector<CCTK_INT> const &res,
          vector<CCTK_INT> const &rhs, options_t const &options,
          CCTK_REAL &error);

void residual(cGH const *restrict const cctkGH, int const options_table,
              calcfunc const calcres, void *const userdata) throw(char const *);

void boundary(cGH const *restrict const cctkGH, int const options_table,
              calcfunc const applybnds,
              void *const userdata) throw(char const *);

void restrict_var(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
                  options_t const &options);

void prolongate_var(cGH const *restrict const cctkGH,
                    vector<CCTK_INT> const &var, options_t const &options);

void zero(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
          options_t const &options);

void copy(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
          vector<CCTK_INT> const &src, options_t const &options);

void add(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
         vector<CCTK_INT> const &src, options_t const &options);

void subtract(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
              vector<CCTK_INT> const &src, options_t const &options);

void getvars(int const options_table, char const *restrict const name,
             vector<CCTK_INT> &vars);

void checkvar(int const vi, vector<int> &usedvars);

void interior_shape(cGH const *restrict const cctkGH, options_t const &options,
                    int *restrict const imin, int *restrict const imax);

int indwidth();

// Register this solver
int CarpetMG_register() {
  int const ierr = TATelliptic_RegisterSolver(solve, solver_name);
  assert(!ierr);

  return 0;
}

// Check parameters
void CarpetMG_paramcheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(direct_solver, solver_name)) {
    CCTK_VParamWarn(CCTK_THORNSTRING, "It is not possible to use \"%s\" as "
                                      "direct solver for \"%s\" -- this would "
                                      "lead to an infinite recursion",
                    solver_name, solver_name);
  }
}

// Solve the system
int solve(cGH const *restrict const cctkGH, int const *restrict const var_,
          int const *restrict const res_, int const nvar,
          int const options_table, calcfunc const calcres,
          calcfunc const applybnds, void *const userdata) {
  DECLARE_CCTK_PARAMETERS;

  // Check arguments
  assert(cctkGH);
  assert(var_);
  assert(res_);
  assert(nvar >= 0);
  assert(calcres);
  assert(applybnds);

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "%*s[%d] beginning to solve", indwidth(), "",
               reflevel);
  }

  // Solve the equation on this and some coarser levels, leaving all
  // finer levels untouched.
  assert(is_level_mode());

  // Copy solution variables from arguments
  vector<CCTK_INT> var(nvar);
  for (int n = 0; n < nvar; ++n) {
    var.at(n) = var_[n];
  }

  // Copy residual variables from arguments
  vector<CCTK_INT> res(nvar);
  for (int n = 0; n < nvar; ++n) {
    res.at(n) = res_[n];
  }

  // Get RHS variables from options
  vector<CCTK_INT> rhs(nvar);
  {
    int const icnt =
        Util_TableGetIntArray(options_table, nvar, &rhs.front(), "rhs");
    assert(icnt == nvar);
  }

  // Get save variables from options
  vector<CCTK_INT> sav(nvar);
  {
    int const icnt =
        Util_TableGetIntArray(options_table, nvar, &sav.front(), "sav");
    assert(icnt == nvar);
  }

  // Get weight variables from options
  vector<CCTK_INT> wgt;
  getvars(options_table, "wgt", wgt);
  int const nwgt = wgt.size();
  assert(nwgt == 0 or nwgt == nvar);

  // Get auxiliary variables from options
  vector<CCTK_INT> aux;
  getvars(options_table, "aux", aux);
  int const naux = aux.size();
  assert(naux >= 0);

  // Check variables
  {
    vector<int> usedvars;

    for (int n = 0; n < nvar; ++n) {
      checkvar(var.at(n), usedvars);
      checkvar(res.at(n), usedvars);
      checkvar(rhs.at(n), usedvars);
      checkvar(sav.at(n), usedvars);
    }

    for (int n = 0; n < nwgt; ++n) {
      checkvar(wgt.at(n), usedvars);
    }

    for (int n = 0; n < naux; ++n) {
      checkvar(aux.at(n), usedvars);
    }
  }

  options_t options;

  // Get outer boundary information
  {
    options.bndwidth = 1;
    int const icnt =
        Util_TableGetInt(options_table, &options.bndwidth, "bndwidth");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }

  // Get symmetry information
  {
    int const symtable = SymmetryTableHandleForGrid(cctkGH);
    assert(symtable >= 0);
    {
      int const icnt = (Util_TableGetIntArray(
          symtable, 6, options.symmetry_handle, "symmetry_handle"));
      assert(icnt == 6);
    }
    {
      int const icnt = (Util_TableGetIntArray(
          symtable, 6, options.symmetry_zone_width, "symmetry_zone_width"));
      assert(icnt == 6);
    }
  }

  // Get full multigrid information
  {
    options.firstlevel = 0;
    int const icnt =
        Util_TableGetInt(options_table, &options.firstlevel, "firstlevel");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }

  // Get direct solver information
  {
    options.maxiters = 100;
    int const icnt =
        Util_TableGetInt(options_table, &options.maxiters, "maxiters");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }
  {
    options.sor_factor = 1.0;
    int const icnt =
        Util_TableGetReal(options_table, &options.sor_factor, "sor_factor");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }

  // Get equation information
  {
    options.factor = 0.5;
    int const icnt =
        Util_TableGetReal(options_table, &options.factor, "factor");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }
  {
    options.factors.resize(var.size());
    for (int n = 0; n < options.factors.size(); ++n) {
      options.factors.at(n) = 1.0;
    }
    int const icnt =
        (Util_TableGetRealArray(options_table, options.factors.size(),
                                &options.factors.front(), "factors"));
    assert(icnt == options.factors.size() or
           icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }

  // Get multigrid information
  {
    options.mgiters = 10;
    int const icnt =
        Util_TableGetInt(options_table, &options.mgiters, "mgiters");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }
  {
    options.presteps = 2;
    int const icnt =
        Util_TableGetInt(options_table, &options.presteps, "presteps");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }
  {
    options.poststeps = 2;
    int const icnt =
        Util_TableGetInt(options_table, &options.poststeps, "poststeps");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }

  CCTK_REAL minerror = 1.0e-8;
  {
    int const icnt = Util_TableGetReal(options_table, &minerror, "minerror");
    assert(icnt == 1 or icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  }

  if (reflevel < options.firstlevel)
    return +1;

  // Initialise RHS to zero
  BEGIN_GLOBAL_MODE(cctkGH) {
    BEGIN_REFLEVEL_LOOP(cctkGH) { zero(cctkGH, rhs, options); }
    END_REFLEVEL_LOOP;
  }
  END_GLOBAL_MODE;

  CCTK_REAL error;
  int const ierr =
      multigrid(cctkGH, var, res, rhs, sav, wgt, aux, options_table, calcres,
                applybnds, userdata, options, minerror);

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "%*s[%d] finished solving", indwidth(), "",
               reflevel);
  }

  return ierr;
}

// Perform recursive multigrid cycles until converged
// Return values:
// positive: did not converge, continue
// zero:     did converge, exit
// negative: error, abort
int multigrid(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
              vector<CCTK_INT> const &res, vector<CCTK_INT> const &rhs,
              vector<CCTK_INT> const &sav, vector<CCTK_INT> const &wgt,
              vector<CCTK_INT> const &aux, int const options_table,
              calcfunc const calcres, calcfunc const applybnds,
              void *const userdata, options_t const &options,
              CCTK_REAL const minerror) {
  DECLARE_CCTK_PARAMETERS;

  // Solve on this and some coarser levels
  assert(is_level_mode());

  // Solve directly when on the coarsest level
  if (reflevel == 0) {
    return direct_solve(cctkGH, var, res, rhs, wgt, aux, options_table, calcres,
                        applybnds, userdata, options, minerror);
  }

  try {

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "%*s[%d] beginning multigrid solve, desired residual %g",
                 indwidth(), "", reflevel, double(minerror));
    }

    // Loop until converged
    CCTK_REAL error = HUGE_VAL;
    int iter = 0;
    while (error > minerror) {

      if (iter >= options.mgiters) {
        // Did not converge
        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "%*s[%d] finished multigrid solve after "
                                       "%d iterations, residual %g",
                     indwidth(), "", reflevel, iter, double(error));
        }
        return 1;
      }
      ++iter;

      // Presmooth
      {
        int step = 0;
        while (error > minerror) {

          ++step;
          if (step > options.presteps)
            break;

          residual(cctkGH, options_table, calcres, userdata);
          CCTK_REAL const old_error = error;
          smooth(cctkGH, var, res, rhs, wgt, options, old_error, error);
          boundary(cctkGH, options_table, applybnds, userdata);

          if (veryverbose) {
            CCTK_VInfo(
                CCTK_THORNSTRING,
                "%*s[%d] iteration %d, presmoothing step %d, residual %g",
                indwidth(), "", reflevel, iter, step, double(error));
          }

        } // while error > minerror
      }

      // Restrict

      residual(cctkGH, options_table, calcres, userdata);
      norm(cctkGH, res, rhs, options, error);
      subtract(cctkGH, res, rhs, options);
      // TODO: restrict and fixup aux
      assert(aux.empty());

      int const coarse_reflevel = reflevel - 1;
      BEGIN_GLOBAL_MODE(cctkGH) {
        enter_level_mode(const_cast<cGH *>(cctkGH), coarse_reflevel);
        try {

          // Restrict variable
          restrict_var(cctkGH, var, options);
          boundary(cctkGH, options_table, applybnds, userdata);
          copy(cctkGH, sav, var, options);

          // Restrict residual
          zero(cctkGH, res, options);
          restrict_var(cctkGH, res, options);
          copy(cctkGH, rhs, res, options);
          residual(cctkGH, options_table, calcres, userdata);
          subtract(cctkGH, rhs, res, options);

          CCTK_REAL coarse_error;
          norm(cctkGH, res, rhs, options, coarse_error);
          if (veryverbose) {
            CCTK_VInfo(CCTK_THORNSTRING,
                       "%*s[%d] iteration %d, initial coarse residual %g",
                       indwidth(), "", reflevel, iter, double(coarse_error));
          }

          // Recurse
          ivect const reffact =
              (spacereffacts.at(reflevel) / spacereffacts.at(coarse_reflevel));
          CCTK_REAL const coarse_minerror = error / prod(reffact);
          int const ierr =
              multigrid(cctkGH, var, res, rhs, sav, wgt, aux, options_table,
                        calcres, applybnds, userdata, options, coarse_minerror);
          if (ierr < 0)
            throw "multigrid";

          // Prolongate
          copy(cctkGH, res, var, options);
          subtract(cctkGH, res, sav, options);

        } catch (char const *) {
          // TODO
          assert(0);
        }

        leave_level_mode(const_cast<cGH *>(cctkGH));
      }
      END_GLOBAL_MODE;

      // TODO
      // save old solution
      copy(cctkGH, sav, var, options);

      zero(cctkGH, res, options);
      prolongate_var(cctkGH, res, options);
      add(cctkGH, var, res, options);
      boundary(cctkGH, options_table, applybnds, userdata);

      CCTK_REAL const old_error = error;

      try {
        residual(cctkGH, options_table, calcres, userdata);
      } catch (char const *) {
        assert(0);
      }
      norm(cctkGH, res, rhs, options, error);
      if (error > old_error) {
        CCTK_VWarn(
            1, __LINE__, __FILE__, CCTK_THORNSTRING,
            "Residual increased during recursion at level %d from %g to %g",
            reflevel, double(old_error), double(error));
      }

      if (veryverbose) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "%*s[%d] iteration %d, after recursion, residual %g",
                   indwidth(), "", reflevel, iter, double(error));
      }

      // Postsmooth
      {
        int step = 0;
        while (error > minerror) {

          ++step;
          if (step > options.poststeps)
            break;

          residual(cctkGH, options_table, calcres, userdata);
          CCTK_REAL const old_error = error;
          smooth(cctkGH, var, res, rhs, wgt, options, old_error, error);
          boundary(cctkGH, options_table, applybnds, userdata);

          if (veryverbose) {
            CCTK_VInfo(
                CCTK_THORNSTRING,
                "%*s[%d] iteration %d, postsmoothing step %d, residual %g",
                indwidth(), "", reflevel, iter, step, double(error));
          }

        } // while error > minerror
      }

    } // while error > minerror

    if (verbose) {
      CCTK_VInfo(
          CCTK_THORNSTRING,
          "%*s[%d] finished multigrid solve after %d iterations, residual %g",
          indwidth(), "", reflevel, iter, double(error));
    }

    // Everything went fine
    return 0;

  } catch (char const *) {

    // There was an error
    return -1;
  }
}

// Solve directly
// This assumes that the grid covers the whole domain
// Return values:
// positive: did not converge, continue
// zero:     did converge, exit
// negative: error, abort
int direct_solve(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
                 vector<CCTK_INT> const &res, vector<CCTK_INT> const &rhs,
                 vector<CCTK_INT> const &wgt, vector<CCTK_INT> const &aux,
                 int const options_table, calcfunc const calcres,
                 calcfunc const applybnds, void *const userdata,
                 options_t const &options, CCTK_REAL const minerror) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "%*s[%d] beginning direct solve, desired residual %g",
               indwidth(), "", reflevel, double(minerror));
  }

  assert(is_level_mode());

  try {

    // Loop until converged
    CCTK_REAL error = HUGE_VAL;
    int iter = 0;
    while (error > minerror) {

      if (iter >= options.maxiters) {
        // Did not converge
        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "%*s[%d] finished direct solve after %d "
                                       "iterations, residual is %g",
                     indwidth(), "", reflevel, iter, double(error));
        }
        return 1;
      }
      ++iter;

      residual(cctkGH, options_table, calcres, userdata);
      CCTK_REAL const old_error = error;
      smooth(cctkGH, var, res, rhs, wgt, options, old_error, error, true);
      boundary(cctkGH, options_table, applybnds, userdata);

      if (veryverbose) {
        CCTK_VInfo(CCTK_THORNSTRING, "%*s[%d] iteration %d, residual %g",
                   indwidth(), "", reflevel, iter, double(error));
      }

    } // while error > min_error

    if (verbose) {
      CCTK_VInfo(
          CCTK_THORNSTRING,
          "%*s[%d] finished direct solve after %d iterations, residual is %g",
          indwidth(), "", reflevel, iter, double(error));
    }

    // Everything went fine
    return 0;

  } catch (void *) {

    // There was an error
    return -1;
  }
}

// Smooth
void smooth(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
            vector<CCTK_INT> const &res, vector<CCTK_INT> const &rhs,
            vector<CCTK_INT> const &wgt, options_t const &options,
            CCTK_REAL const old_error, CCTK_REAL &error, bool const use_sor) {
  DECLARE_CCTK_ARGUMENTS;

  // Initialise errors
  CCTK_REAL count = 0.0;
  CCTK_REAL error2 = 0.0;

  // Calculate grid spacings
  CCTK_REAL dxinv2 = 0.0;
  CCTK_REAL dx[dim];
  for (int d = 0; d < dim; ++d) {
    // TODO: correct this for solving on grid arrays instead of grid
    // functions
    dx[d] = CCTK_DELTA_SPACE(d);
    dxinv2 += 1.0 / ipow(dx[d], 2);
  }
  CCTK_REAL const mdxinv2inv = 1.0 / (-2.0 * dxinv2);

  // Smooth and calculate errors
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      for (int n = 0; n < var.size(); ++n) {

        CCTK_REAL *restrict const varptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, var.at(n))));
        assert(varptr);
        CCTK_REAL const *restrict const resptr =
            (static_cast<CCTK_REAL const *>(
                CCTK_VarDataPtrI(cctkGH, 0, res.at(n))));
        assert(resptr);
        CCTK_REAL const *restrict const rhsptr =
            (static_cast<CCTK_REAL const *>(
                CCTK_VarDataPtrI(cctkGH, 0, rhs.at(n))));
        assert(rhsptr);
        CCTK_REAL const *restrict const wgtptr =
            (wgt.size() > 0 ? (static_cast<CCTK_REAL const *>(
                                  CCTK_VarDataPtrI(cctkGH, 0, wgt.at(n))))
                            : NULL);
        assert(wgt.empty() or wgtptr);

        CCTK_REAL const fac =
            ((use_sor ? options.sor_factor : options.factor) *
             (wgtptr ? 1.0 : options.factors.at(n) * mdxinv2inv));

        int imin[3], imax[3];
        interior_shape(cctkGH, options, imin, imax);

        for (int k = imin[2]; k < imax[2]; ++k) {
          for (int j = imin[1]; j < imax[1]; ++j) {
            for (int i = imin[0]; i < imax[0]; ++i) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

              CCTK_REAL const diff = resptr[ind] - rhsptr[ind];
              CCTK_REAL const w = wgtptr ? fac / wgtptr[ind] : fac;

              varptr[ind] -= w * diff;

              ++count;
              error2 += ipow(diff, 2);
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;

  // Reduce errors
  int const sum_handle = CCTK_ReductionArrayHandle("sum");
  assert(sum_handle >= 0);
  CCTK_REAL reduce_in[2], reduce_out[2];
  reduce_in[0] = count;
  reduce_in[1] = error2;
  int const ierr = CCTK_ReduceLocArrayToArray1D(
      cctkGH, -1, sum_handle, reduce_in, reduce_out, 2, CCTK_VARIABLE_REAL);
  count = reduce_out[0];
  error2 = reduce_out[1];
  if (count > 0) {
    error = sqrt(error2 / count);
  } else {
    error = 0.0;
  }

  // Sanity check
  if (error > old_error) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Residual increased during smoothing at level %d from %g to %g",
               reflevel, double(old_error), double(error));
  }
}

// Calculate the norm of the residual without smoothing
void norm(cGH const *restrict const cctkGH, vector<CCTK_INT> const &res,
          vector<CCTK_INT> const &rhs, options_t const &options,
          CCTK_REAL &error) {
  DECLARE_CCTK_ARGUMENTS;

  // Initialise errors
  CCTK_REAL count = 0.0;
  CCTK_REAL error2 = 0.0;

  // Calculate errors
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      for (int n = 0; n < res.size(); ++n) {

        CCTK_REAL const *restrict const resptr =
            (static_cast<CCTK_REAL const *>(
                CCTK_VarDataPtrI(cctkGH, 0, res.at(n))));
        assert(resptr);
        CCTK_REAL const *restrict const rhsptr =
            (static_cast<CCTK_REAL const *>(
                CCTK_VarDataPtrI(cctkGH, 0, rhs.at(n))));
        assert(rhsptr);

        int imin[3], imax[3];
        interior_shape(cctkGH, options, imin, imax);

        for (int k = imin[2]; k < imax[2]; ++k) {
          for (int j = imin[1]; j < imax[1]; ++j) {
            for (int i = imin[0]; i < imax[0]; ++i) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

              CCTK_REAL const diff = resptr[ind] - rhsptr[ind];

              ++count;
              error2 += ipow(diff, 2);
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;

  // Reduce errors
  int const sum_handle = CCTK_ReductionArrayHandle("sum");
  assert(sum_handle >= 0);
  CCTK_REAL reduce_in[2], reduce_out[2];
  reduce_in[0] = count;
  reduce_in[1] = error2;
  int const ierr = CCTK_ReduceLocArrayToArray1D(
      cctkGH, -1, sum_handle, reduce_in, reduce_out, 2, CCTK_VARIABLE_REAL);
  count = reduce_out[0];
  error2 = reduce_out[1];
  if (count > 0) {
    error = sqrt(error2 / count);
  } else {
    error = 0.0;
  }
}

// Calculate the residual
void residual(cGH const *restrict const cctkGH, int const options_table,
              calcfunc const calcres,
              void *const userdata) throw(char const *) {
  int ierr = 0;
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      if (!ierr) {
        ierr = calcres(cctkGH, options_table, userdata);
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
  if (ierr)
    throw "calcres";
}

// Apply the boundary condition
void boundary(cGH const *restrict const cctkGH, int const options_table,
              calcfunc const applybnds,
              void *const userdata) throw(char const *) {
  int ierr = 0;
  ierr = applybnds(cctkGH, options_table, userdata);
  if (ierr)
    throw "applybnds";
  ierr = CallScheduleGroup(const_cast<cGH *>(cctkGH), "ApplyBCs");
  assert(!ierr);
}

// Restrict the solution variables from the next finer level
void restrict_var(cGH const *restrict const cctkGH, vector<CCTK_INT> const &var,
                  options_t const &options) {
  assert(reflevel < reflevels - 1);
  int const tl = 0;
  for (comm_state state; !state.done(); state.step()) {
    for (int m = 0; m < maps; ++m) {
      for (int n = 0; n < var.size(); ++n) {
        int const vi = var.at(n);
        assert(vi >= 0);
        int const gi = CCTK_GroupIndexFromVarI(vi);
        assert(gi >= 0);
        int const v0 = CCTK_FirstVarIndexI(gi);
        assert(v0 >= 0);
        arrdata.at(gi).at(m).data.at(vi - v0)->ref_restrict_all(
            state, tl, reflevel, mglevel);
      } // for n
    }   // for m
  }     // for state
}

// Prolongate the solution variables from the next coarser level
void prolongate_var(cGH const *restrict const cctkGH,
                    vector<CCTK_INT> const &var, options_t const &options) {
  assert(reflevel > 0);
  int const tl = 0;
  CCTK_REAL const time = tt->get_time(mglevel, reflevel, tl);
  for (comm_state state; !state.done(); state.step()) {
    for (int m = 0; m < maps; ++m) {
      for (int n = 0; n < var.size(); ++n) {
        int const vi = var.at(n);
        assert(vi >= 0);
        int const gi = CCTK_GroupIndexFromVarI(vi);
        assert(gi >= 0);
        int const v0 = CCTK_FirstVarIndexI(gi);
        assert(v0 >= 0);
        arrdata.at(gi).at(m).data.at(vi - v0)->ref_prolongate_all(
            state, tl, reflevel, mglevel, time);
      } // for n
    }   // for m
  }     // for state
}

void zero(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
          options_t const &options) {
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      for (int n = 0; n < dst.size(); ++n) {
        CCTK_REAL *restrict const dstptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, dst.at(n))));
        assert(dstptr);
        for (int k = 0; k < cctk_lsh[2]; ++k) {
          for (int j = 0; j < cctk_lsh[1]; ++j) {
            for (int i = 0; i < cctk_lsh[0]; ++i) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
              dstptr[ind] = 0.0;
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
}

void copy(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
          vector<CCTK_INT> const &src, options_t const &options) {
  assert(dst.size() == src.size());
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      for (int n = 0; n < dst.size(); ++n) {
        CCTK_REAL *restrict const dstptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, dst.at(n))));
        assert(dstptr);
        CCTK_REAL const *restrict const srcptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, src.at(n))));
        assert(srcptr);
        for (int k = 0; k < cctk_lsh[2]; ++k) {
          for (int j = 0; j < cctk_lsh[1]; ++j) {
            for (int i = 0; i < cctk_lsh[0]; ++i) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
              dstptr[ind] = srcptr[ind];
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
}

void add(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
         vector<CCTK_INT> const &src, options_t const &options) {
  assert(dst.size() == src.size());
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      for (int n = 0; n < dst.size(); ++n) {
        CCTK_REAL *restrict const dstptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, dst.at(n))));
        assert(dstptr);
        CCTK_REAL const *restrict const srcptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, src.at(n))));
        assert(srcptr);
        for (int k = 0; k < cctk_lsh[2]; ++k) {
          for (int j = 0; j < cctk_lsh[1]; ++j) {
            for (int i = 0; i < cctk_lsh[0]; ++i) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
              dstptr[ind] += srcptr[ind];
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
}

void subtract(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
              vector<CCTK_INT> const &src, options_t const &options) {
  assert(dst.size() == src.size());
  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      for (int n = 0; n < dst.size(); ++n) {
        CCTK_REAL *restrict const dstptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, dst.at(n))));
        assert(dstptr);
        CCTK_REAL const *restrict const srcptr =
            (static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, src.at(n))));
        assert(srcptr);
        for (int k = 0; k < cctk_lsh[2]; ++k) {
          for (int j = 0; j < cctk_lsh[1]; ++j) {
            for (int i = 0; i < cctk_lsh[0]; ++i) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
              dstptr[ind] -= srcptr[ind];
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
}

// Get a variable list from an options table
void getvars(int const options_table, char const *restrict const name,
             vector<CCTK_INT> &vars) {
  int const nvars = Util_TableGetIntArray(options_table, 0, NULL, name);
  if (nvars == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    vars.resize(0);
  } else if (nvars >= 0) {
    assert(nvars >= 0);
    vars.resize(nvars);
    int const icnt =
        Util_TableGetIntArray(options_table, nvars, &vars.front(), name);
    assert(icnt == nvars);
  } else {
    assert(0);
  }
}

// Check a grid variable index
void checkvar(int const vi, vector<int> &usedvars) {
  int ierr;

  assert(vi >= 0 and vi < CCTK_NumVars());

  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);

  cGroup group;
  ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);

  assert(group.grouptype == CCTK_GF);
  assert(group.vartype == CCTK_VARIABLE_REAL);
  assert(group.dim == 3);

  for (int n = 0; n < usedvars.size(); ++n) {
    assert(vi != usedvars.at(n));
  }
  usedvars.push_back(vi);
}

// Calculate the shape of the interior
void interior_shape(cGH const *restrict const cctkGH, options_t const &options,
                    int *restrict const imin, int *restrict const imax) {
  DECLARE_CCTK_ARGUMENTS;

  for (int d = 0; d < dim; ++d) {

    imin[d] = 0;
    if (cctk_bbox[2 * d]) {
      if (options.symmetry_handle[2 * d] >= 0) {
        imin[d] += options.symmetry_zone_width[2 * d];
      } else {
        imin[d] += options.bndwidth;
      }
    } else {
      imin[d] += cctk_nghostzones[d];
    }

    imax[d] = cctk_lsh[d];
    if (cctk_bbox[2 * d + 1]) {
      if (options.symmetry_handle[2 * d + 1] >= 0) {
        imax[d] -= options.symmetry_zone_width[2 * d + 1];
      } else {
        imax[d] -= options.bndwidth;
      }
    } else {
      imax[d] -= cctk_nghostzones[d];
    }

    assert(imin[d] <= imax[d]);

  } // for d
}

// Determine indentation
int indwidth() {
  if (reflevel == -1)
    return 0;
  return 3 * (reflevels - reflevel - 1);
}

} // namespace CarpetMG
