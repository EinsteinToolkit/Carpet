#ifndef INDEXING_HH
#define INDEXING_HH

#include <cctk.h>

#include <cassert>

namespace CarpetRegrid2 {

// Get indexing information for a vector grid array
void getvectorindex2(cGH const *cctkGH, char const *groupname,
                     int *restrict lsh, int *restrict ash);

static inline int index2(int const *const lsh, int const *const ash,
                         int const i, int const j) {
  assert(lsh);
  assert(ash);
  assert(lsh[0] >= 0 && lsh[0] <= ash[0]);
  assert(lsh[1] >= 0 && lsh[1] <= ash[1]);
  assert(i >= 0 and i < lsh[0]);
  assert(j >= 0 and j < lsh[1]);
  return i + ash[0] * j;
}

} // namespace CarpetRegrid2

#endif // #ifndef INDEXING_HH
