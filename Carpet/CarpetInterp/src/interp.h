/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.h,v 1.2 2003/11/05 16:18:38 schnetter Exp $ */

#ifndef CARPETINTERP_H
#define CARPETINTERP_H

#ifdef __cplusplus
namespace CarpetInterp {
  extern "C" {
#endif
  
    /* Scheduled functions */
    void CarpetInterpStartup (void);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetInterp */
#endif

#endif /* !defined(CARPETINTERP_H) */
