#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetMask {

extern "C" {
void CarpetBoxExcludedSetup(CCTK_ARGUMENTS);
}

} // namespace CarpetMask
