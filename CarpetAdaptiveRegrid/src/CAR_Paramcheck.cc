#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "CAR.hh"
#include "carpet.hh"

namespace CarpetAdaptiveRegrid {

using namespace std;
using namespace Carpet;

void CarpetAdaptiveRegridParamcheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int type;
  const CCTK_INT *const domain_from_coordbase =
      (const CCTK_INT *)CCTK_ParameterGet("domain_from_coordbase", "Carpet",
                                          &type);
  assert(domain_from_coordbase);
  assert(type == PARAMETER_BOOLEAN);
  if (!*domain_from_coordbase) {
    CCTK_PARAMWARN(
        "CarpetAdaptiveRegrid requires that Carpet::domain_from_coordbase=yes");
  }
}
}
