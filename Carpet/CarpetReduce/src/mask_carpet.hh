#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetMask {

  extern "C" {
    void
    CarpetMaskSetup (CCTK_ARGUMENTS);
  }
  
} // namespace CarpetMask
