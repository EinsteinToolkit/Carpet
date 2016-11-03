#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void PeriodicCarpet_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (periodic and (periodic_x or periodic_y or periodic_z))
    CCTK_PARAMWARN("When PeriodicCarpet::periodic is set, all boundaries are "
                   "periodic. It does not make sense to also set any of "
                   "PeriodicCarpet::periodic_[xyz]");
}
