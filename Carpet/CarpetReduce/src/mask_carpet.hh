// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/mask_carpet.hh,v 1.1 2004/06/14 07:01:21 schnetter Exp $

#include "cctk.h"
#include "cctk_Arguments.h"

namespace CarpetMask {

  extern "C" {
    void
    CarpetMaskSetup (CCTK_ARGUMENTS);
  }
  
} // namespace CarpetMask
