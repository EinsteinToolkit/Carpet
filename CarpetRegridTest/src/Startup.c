#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

int CarpetRegridTest_Startup(void) {
  const char *banner = "CarpetRegridTest: Thoroughly testing PMR";

  CCTK_RegisterBanner(banner);

  return 0;
}
