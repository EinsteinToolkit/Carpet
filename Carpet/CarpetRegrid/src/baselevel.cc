#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/baselevel.cc,v 1.1 2004/01/25 14:57:30 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_baselevel_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int BaseLevel (cGH const * const cctkGH,
                 gh<dim> const & hh,
                 int const reflevel,
                 int const map,
                 int const size,
                 jjvect const & nboundaryzones,
                 jjvect const & is_internal,
                 jjvect const & is_staggered,
                 jjvect const & shiftout,
                 gh<dim>::rexts  & bbsss,
                 gh<dim>::rbnds  & obss,
                 gh<dim>::rprocs & pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (reflevel>=0 && reflevel<maxreflevels);
    assert (map>=0 && map<maps);
    
    assert (refinement_levels == 1);
    
    assert (bbsss.size() == 1);
    
    return 0;
  }
  
} // namespace CarpetRegrid
