/* $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/ioflexio.h,v 1.2 2003/09/17 13:47:00 cvs_anon Exp $ */

#ifndef CARPETIOFLEXIO_H
#define CARPETIOFLEXIO_H

#include "cctk_Arguments.h"


    
#ifdef __cplusplus
namespace CarpetIOFlexIO {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetIOFlexIO_Startup (void);
    int CarpetIOFlexIO_ReadData (CCTK_ARGUMENTS);
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOFlexIO */
#endif

#ifdef __cplusplus
namespace CarpetCheckpointRestart {
  extern "C" {
#endif
    
    /* Scheduled functions */
    void CarpetIOFlexIO_EvolutionCheckpoint (const cGH*);
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOFlexIO */
#endif


#endif /* !defined(CARPETIOFLEXIO_H) */





