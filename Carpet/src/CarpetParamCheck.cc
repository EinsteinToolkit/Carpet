#include <cassert>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

namespace Carpet {

using namespace std;

/** Ensure that the parameters have legal values.
 *
 * Note that this checking happens only after most of Carpet has
 * already been set up.
 */
void CarpetParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CarpetParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if (poison_taper and
      not CCTK_IsThornActive("MoL")) {
    CCTK_PARAMWARN("Carpet::poison_taper requires MoL");
  }

  if ((CCTK_ParameterQueryTimesSet("periodic", "Carpet") and periodic) or
      (CCTK_ParameterQueryTimesSet("periodic_x", "Carpet") and periodic_x) or
      (CCTK_ParameterQueryTimesSet("periodic_y", "Carpet") and periodic_y) or
      (CCTK_ParameterQueryTimesSet("periodic_z", "Carpet") and periodic_z)) {
    CCTK_PARAMWARN("Some of the parameters \"Carpet::periodic*\" have been set "
                   "to \"yes\".  These parameters are there for compatibility "
                   "reasons only and must not be used.");
  }

  if (adaptive_stepsize and max_refinement_levels > 1) {
    CCTK_PARAMWARN("Adaptive time step sizes do not work with mesh refinement "
                   "yet.  Please use only a single level, and set "
                   "max_refinement_levels=1.");
  }

  // These fire for a single level as well. Buffer zones are named in
  // the tapered-grid message because a tapered grid without them has
  // only one ghost-width, which does not cover a multi-stage step.
  if (use_csamr and not use_tapered_grids) {
    CCTK_PARAMWARN(
        "Carpet::use_csamr requires Carpet::use_tapered_grids. Tapered "
        "grids and buffer zones stop an RK stage from prolongating out "
        "of the parent during the step.");
  }
  if (use_csamr and not use_buffer_zones) {
    CCTK_PARAMWARN(
        "Carpet::use_csamr requires Carpet::use_buffer_zones. Without "
        "buffer zones a refined grid has only one ghost-width, which "
        "does not cover a multi-stage Runge-Kutta step.");
  }
  if (use_csamr and adaptive_stepsize) {
    CCTK_PARAMWARN("Carpet::use_csamr does not support adaptive timesteps.");
  }

  // The serial schedule has no tapered step for Hermite to replace.
  // The CSAMR schedule keeps the buffer on the first 2:1 fine step.
  if (hermite_dense_output and not use_csamr and use_buffer_zones) {
    CCTK_PARAMWARN("Carpet::hermite_dense_output on the serial schedule "
                   "replaces buffer zones. Set Carpet::use_buffer_zones = no.");
  }
  if (hermite_dense_output and not use_csamr and use_tapered_grids) {
    CCTK_PARAMWARN("Carpet::hermite_dense_output on the serial schedule "
                   "replaces tapered buffers. Set Carpet::use_tapered_grids "
                   "= no.");
  }
  if (hermite_dense_output and adaptive_stepsize) {
    CCTK_PARAMWARN("Carpet::hermite_dense_output does not support adaptive "
                   "timesteps.");
  }
  if (hermite_dense_output) {
    int parameter_type = 0;
    void const *const parameter =
        CCTK_ParameterGet("hermite_rhs_storage", "MoL", &parameter_type);
    if (parameter == NULL or parameter_type != PARAMETER_BOOLEAN or
        *static_cast<CCTK_INT const *>(parameter) == 0) {
      CCTK_PARAMWARN("Carpet::hermite_dense_output requires "
                     "MoL::hermite_rhs_storage = yes.");
    }
  }

  if (max_refinement_levels > 1) {

    // InitBase

    enum setup_method_t {
      init_error,
      init_some_levels,
      init_single_level,
      init_two_levels,
      init_all_levels,
    };
    setup_method_t setup_method = init_error;

    if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels")) {
      setup_method = init_some_levels;
    } else if (CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
      setup_method = init_single_level;
    } else if (CCTK_EQUALS(initial_data_setup_method, "init_two_levels")) {
      setup_method = init_two_levels;
    } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {
      setup_method = init_all_levels;
    } else {
      CCTK_PARAMWARN(
          "Unknown value of parameter InitBase::initial_data_setup_method");
    }

    // Carpet

    enum init_method_t {
      error_timelevel,
      each_timelevel,   // Initialise each time level
      fill_timelevels,  // Copy the current to all past timelevels
      three_timelevels, // Carpet's 3 timelevel scheme
      all_timelevels,   // Initial data thorn initialises all timelevels
    };
    init_method_t init_method = error_timelevel;

    if (init_each_timelevel) {
      if (init_fill_timelevels) {
        CCTK_PARAMWARN("Cannot set parameters init_each_timelevel and "
                       "init_fill_timelevels at the same time");
      }
      if (init_3_timelevels) {
        CCTK_PARAMWARN("Cannot set parameters init_each_timelevel and "
                       "init_3_timelevels at the same time");
      }
      init_method = each_timelevel;
    } else if (init_fill_timelevels) {
      if (init_3_timelevels) {
        CCTK_PARAMWARN("Cannot set parameters init_fill_timelevels and "
                       "init_3_timelevels at the same time");
      }
      init_method = fill_timelevels;
    } else if (init_3_timelevels) {
      init_method = three_timelevels;
    } else {
      init_method = all_timelevels;
    }

    switch (init_method) {
    case each_timelevel:
      if (setup_method != init_single_level) {
        CCTK_PARAMWARN(
            "When you set Carpet::init_each_timelevel=yes, then you must also "
            "use InitBase::initial_data_setup_method=\"init_single_level\"");
      }
      break;
    case fill_timelevels:
      // Everything is allowed
      break;
    case three_timelevels:
      // Everything is allowed
      break;
    case all_timelevels:
      if (setup_method != init_all_levels) {
        CCTK_PARAMWARN("When you set neither Carpet::init_each_timelevel=yes "
                       "nor Carpet::init_fill_timelevels=yes nor "
                       "Carpet::init_3_timelevels=yes, then you must also use "
                       "InitBase::initial_data_setup_method=\"init_all_"
                       "levels\"");
      }
      break;
    default:
      assert(0);
    }

  } // if max_refinement_levels > 0
}

} // namespace Carpet
