#include "Requirements.hh"
#include "all_clauses.hh"
#include "all_state.hh"
#include "gridpoint.hh"
#include "util.hh"

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>
#include <util_String.h>

#include <cstdlib>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

namespace Requirements {

// Rules:
//
// 1. Everything that is required by a routine must be provided by
//    another routine which is scheduled earlier.
//
// 2. Things can be provided only once, not multiple times.
//    Except when they are also required.

inline ostream &operator<<(ostream &os, const all_clauses_t &a) {
  a.output(os);
  return os;
}

all_clauses_t all_clauses;
all_state_t all_state;

void Setup(int const maps) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Setup maps=%d", maps);
    }
    all_state.setup(maps);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void ChangeStorage(vector<int> const &groups, vector<int> const &timelevels,
                   int const reflevel) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      std::ostringstream stream;
      stream << "groups: " << groups << " timelevels: " << timelevels;
      CCTK_VInfo(CCTK_THORNSTRING, "ChangeStorage reflevel=%d %s", reflevel,
                 stream.str().c_str());
    }
    all_state.change_storage(groups, timelevels, reflevel);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void Regrid(int const reflevels) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Regrid reflevels=%d", reflevels);
    }
    all_state.regrid(reflevels);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void Recompose(int const iteration, int const reflevel,
               valid::valid_t const where) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Recompose reflevel=%d where=%s", reflevel,
                 where == valid::nowhere
                     ? "nowhere"
                     : where == valid::interior
                           ? "interior"
                           : where == valid::everywhere ? "everywhere" : NULL);
    }
    all_state.recompose(iteration, reflevel, where);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void RegridFree() {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "RegridFree");
    }
    all_state.regrid_free();
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void Cycle(int const reflevel) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Cycle reflevel=%d", reflevel);
    }
    all_state.cycle(reflevel);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void BeforeRoutine(cFunctionData const *const function_data,
                   int const iteration, int const reflevel, int const map,
                   int const timelevel, int const timelevel_offset) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    all_state.before_routine(function_data, all_clauses, iteration, reflevel,
                             map, timelevel, timelevel_offset);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void AfterRoutine(cFunctionData const *const function_data, int const iteration,
                  int const reflevel, int const map, int const timelevel,
                  int const timelevel_offset) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    all_state.after_routine(function_data, all_clauses, iteration, reflevel,
                            map, timelevel, timelevel_offset);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void Sync(cFunctionData const *const function_data, vector<int> const &groups,
          int const iteration, int const reflevel, int const timelevel) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Sync reflevel=%d timelevel=%d", reflevel,
                 timelevel);
    }
    all_state.sync(function_data, groups, iteration, reflevel, timelevel);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

void Restrict(vector<int> const &groups, int const iteration,
              int const reflevel) {
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Restrict reflevel=%d", reflevel);
    }
    all_state.restrict1(groups, iteration, reflevel);
  }
  if (inconsistencies_are_fatal and gridpoint_t::there_was_an_error) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Aborting because schedule clauses were not satisfied");
  }
}

////////////////////////////////////////////////////////////////////////////

// Check that the grid is in the correct state, i.e. all necessary
// parts are valid, for the "current" function.
extern "C" void Carpet_Requirements_CheckReads(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const nvars,
    CCTK_INT const *const varinds, char const *const clause) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    // TODO: come up with a scheme to avoid constructing and destroying clauses
    cFunctionData const *const function_data =
        CCTK_ScheduleQueryCurrentFunction(cctkGH);
    int const reflevel = GetRefinementLevel(cctkGH);
    int const map = GetMap(cctkGH);
    int const timelevel = GetTimeLevel(cctkGH);
    int const timelevel_offset = GetTimeLevelOffset(cctkGH);
    // TODO: design an interface to all_state.before_routine that operates
    //       on indices and clauses directly
    for (int v = 0; v < nvars; ++v) {
      cFunctionData temp_function_data = *function_data;
      char *const fullname = CCTK_FullName(varinds[v]);
      char *reads;
      int const len_written = Util_asprintf(&reads, "%s(%s)", fullname, clause);
      assert(len_written > 0);
      temp_function_data.n_WritesClauses = 0;
      temp_function_data.WritesClauses = NULL;
      temp_function_data.n_ReadsClauses = 1;
      temp_function_data.ReadsClauses = (char const **)&reads;
      all_clauses.get_clauses(&temp_function_data);
      BeforeRoutine(&temp_function_data, cctkGH->cctk_iteration, reflevel, map,
                    timelevel, timelevel_offset);
      all_clauses.remove_clauses(&temp_function_data);
      free(fullname);
      free(reads);
    }
  }
}

// Register the fact that certain parts of the grid have been
// written in certain variables due to executing the "current"
// function.
extern "C" void Carpet_Requirements_NotifyWrites(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const nvars,
    CCTK_INT const *const varinds, char const *const clause) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    // TODO: come up with a scheme to avoid constructing and
    // destroying clauses
    cFunctionData const *const function_data =
        CCTK_ScheduleQueryCurrentFunction(cctkGH);
    int const reflevel = GetRefinementLevel(cctkGH);
    int const map = GetMap(cctkGH);
    int const timelevel = GetTimeLevel(cctkGH);
    int const timelevel_offset = GetTimeLevelOffset(cctkGH);
    // TODO: design an interface to all_state.before_routine that
    //       operates on indices and claues directly
    for (int v = 0; v < nvars; ++v) {
      cFunctionData temp_function_data = *function_data;
      char *const fullname = CCTK_FullName(varinds[v]);
      char *writes;
      int const len_written =
          Util_asprintf(&writes, "%s(%s)", fullname, clause);
      assert(len_written > 0);
      temp_function_data.n_WritesClauses = 1;
      temp_function_data.WritesClauses = (char const **)&writes;
      temp_function_data.n_ReadsClauses = 0;
      temp_function_data.ReadsClauses = NULL;
      all_clauses.get_clauses(&temp_function_data);
      AfterRoutine(&temp_function_data, cctkGH->cctk_iteration, reflevel, map,
                   timelevel, timelevel_offset);
      all_clauses.remove_clauses(&temp_function_data);
      free(fullname);
      free(writes);
    }
  }
}

extern "C" void
Carpet_Requirements_Invalidate(CCTK_POINTER_TO_CONST const cctkGH_,
                               CCTK_INT const nvars,
                               CCTK_INT const *const varinds) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  DECLARE_CCTK_PARAMETERS;
  if (check_requirements) {
    vector<int> vars(nvars);
    for (int v = 0; v < nvars; ++v) {
      vars.AT(v) = varinds[v];
    }
    int const reflevel = GetRefinementLevel(cctkGH);
    int const map = GetMap(cctkGH);
    int const timelevel = GetTimeLevel(cctkGH);
    all_state.invalidate(vars, reflevel, map, timelevel);
  }
}

////////////////////////////////////////////////////////////////////////////

// scheduled routines to handle boundary and symmetry conditions
extern "C" void CarpetCheckReadsBeforeBoundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  int num_vars, err;
  vector<CCTK_INT> vars, faces, widths, tables;

  num_vars = Boundary_SelectedGVs(cctkGH, 0, NULL, NULL, NULL, NULL, NULL);
  if (num_vars < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error retrieving number of selected GVs: %d", num_vars);
  }
  vars.resize(num_vars);
  faces.resize(num_vars);
  widths.resize(num_vars);
  tables.resize(num_vars);

  /* get selected vars for all bc */
  err = Boundary_SelectedGVs(cctkGH, num_vars, &vars[0], &faces[0], &widths[0],
                             &tables[0], NULL);
  if (err < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error in Boundary_SelectedGVs for all boundary conditions");
  } else if (err != num_vars) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Boundary_SelectedGVs returned %d selected variables for "
               "all boundary conditions, but %d expected\n",
               err, num_vars);
  }

  Requirements_CheckReads(cctkGH, num_vars, &vars[0], "interior");
}

extern "C" void CarpetNotifyWritesAfterBoundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  int num_vars, err;
  vector<CCTK_INT> vars, faces, widths, tables;

  num_vars = Boundary_SelectedGVs(cctkGH, 0, NULL, NULL, NULL, NULL, NULL);
  if (num_vars < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error retrieving number of selected GVs: %d", num_vars);
  }
  vars.resize(num_vars);
  faces.resize(num_vars);
  widths.resize(num_vars);
  tables.resize(num_vars);

  /* get selected vars for all bc */
  err = Boundary_SelectedGVs(cctkGH, num_vars, &vars[0], &faces[0], &widths[0],
                             &tables[0], NULL);
  if (err < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error in Boundary_SelectedGVs for all boundary conditions");
  } else if (err != num_vars) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Boundary_SelectedGVs returned %d selected variables for "
               "all boundary conditions, but %d expected\n",
               err, num_vars);
  }

  Requirements_NotifyWrites(cctkGH, num_vars, &vars[0],
                            "boundary;boundary_ghostzones");
}

} // namespace Requirements
