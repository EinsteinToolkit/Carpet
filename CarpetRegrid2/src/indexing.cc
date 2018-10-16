#include "indexing.hh"

#include <cctk.h>

#include <cassert>

namespace CarpetRegrid2 {

// Get indexing information for a vector grid array
void getvectorindex2(cGH const *const cctkGH, char const *const groupname,
                     int *restrict const lsh, int *restrict const ash) {
  assert(groupname);
  assert(lsh);
  assert(ash);

  int const gi = CCTK_GroupIndex(groupname);
  assert(gi >= 0);

  {
    int const ierr = CCTK_GrouplshGI(cctkGH, 1, lsh, gi);
    assert(not ierr);
  }

  {
    int const ierr = CCTK_GroupashGI(cctkGH, 1, ash, gi);
    assert(not ierr);
  }

  cGroup groupdata;
  {
    int const ierr = CCTK_GroupData(gi, &groupdata);
    assert(not ierr);
  }
  assert(groupdata.vectorgroup);
  assert(groupdata.vectorlength >= 0);
  lsh[1] = groupdata.vectorlength;
  ash[1] = groupdata.vectorlength;
}

} // namespace CarpetRegrid2
