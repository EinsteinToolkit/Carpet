#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>


void regridsynctest_output(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  
  CCTK_OutputVarAsByMethod(cctkGH,"regridsynctest::myregridtestgf","IOASCII_1D","myregridtestgf_postregrid");
  //  CCTK_OutputVarAsByMethod(cctkGH,"whisky::rho","IOASCII_1D","rho_postregrid");
    //  CCTK_OutputVarAsByMethod{cctkGH,"whisky::rho","IOASCII_3D","rho_postregrid")



  return;
}
