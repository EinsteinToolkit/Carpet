#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int BaseLevel (cGH const * const cctkGH,
                 gh<dim> const & hh,
                 gh<dim>::rexts  & bbsss,
                 gh<dim>::rbnds  & obss,
                 gh<dim>::rprocs & pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels == 1);
    
    assert (bbsss.size() == 1);
    
    return 0;
  }
  
} // namespace CarpetRegrid
