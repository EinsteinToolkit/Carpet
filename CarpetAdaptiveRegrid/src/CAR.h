#ifndef CARPETADAPTIVEREGRID_H
#define CARPETADAPTIVEREGRID_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetAdaptiveRegrid {
extern "C" {
#endif

/* Scheduled functions */
void CarpetAdaptiveregridParamcheck(CCTK_ARGUMENTS);

#ifdef __cplusplus
} /* extern "C" */
} /* namespace CarpetAdaptiveregrid */
#endif

#endif /* !defined(CARPETADAPTIVEREGRID_H) */
