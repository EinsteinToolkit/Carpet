#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"
#include "regrid.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/paramcheck.cc,v 1.1 2002/05/16 23:25:54 schnetter Exp $";

CCTK_FILEVERSION(CarpetRegrid_paramcheck_cc)

namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  int CarpetRegridParamcheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (refinement_levels > maxreflevels) {
      CCTK_PARAMWARN ("The parameter CarpetRegrid::refinement_levels is larger than Carpet::max_refinement_levels");
    }
    
    return 0;
  }
  
}
