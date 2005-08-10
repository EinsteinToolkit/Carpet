#include "cctk.h"
#include "cctk_Arguments.h"

namespace CarpetEvolutionMask {

  extern "C" {
    void
    CarpetEvolutionMaskSetup (CCTK_ARGUMENTS);
  }
  
} // namespace CarpetEovlutionMask
