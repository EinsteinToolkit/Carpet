#include <assert.h>
#include <iostream.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/CarpetParamCheck.cc,v 1.2 2001/12/09 16:41:52 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  int CarpetParamCheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    return 0;
  }
  
} // namespace Carpet
