#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


static const char *rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/CarpetRegridTest/src/Startup.c,v 1.1 2004/05/23 15:56:14 cott Exp $";

CCTK_FILEVERSION(CarpetExtra_CarpetRegridTest_Startup_c);

void CarpetRegridTest_Startup(void)
{
  const char *banner = "CarpetRegridTest: Thoroughly testing PMR";

  CCTK_RegisterBanner(banner);

}
