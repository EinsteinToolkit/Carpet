#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


static const char *rcsid = "$Header:$";

CCTK_FILEVERSION(CarpetExtra_CarpetRegridTest_Startup_c);

void CarpetRegridTest_Startup(void)
{
  const char *banner = "CarpetRegridTest: Thoroughly testing PMR";

  CCTK_RegisterBanner(banner);

}
