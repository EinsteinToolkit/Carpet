/* $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/ioflexio.h,v 1.3 2003/12/01 13:15:21 cott Exp $ */

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
    int CarpetIOFlexIO_RecoverParameters (void);
    int CarpetIOFlexIO_Recover (cGH *GH, const char *basefilename, int called_from);
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOFlexIO */
#endif


#endif /* !defined(CARPETIOFLEXIO_H) */





