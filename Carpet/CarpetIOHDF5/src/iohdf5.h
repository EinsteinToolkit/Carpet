/* $Header:$ */

#ifndef CARPETIOHDF5_H
#define CARPETIOHDF5_H

#ifdef __cplusplus
namespace CarpetIOHDF5 {
  extern "C" {
#endif
    
#ifdef USE_CHRISTIANS_ROUTINES
#include "cctk_Arguments.h"

int CarpetIOHDF5Startup (void);
void CarpetIOHDF5Init (CCTK_ARGUMENTS);
void CarpetIOHDF5ReadData (CCTK_ARGUMENTS);
void CarpetIOHDF5_EvolutionCheckpoint (const cGH*);
void CarpetIOHDF5_InitialDataCheckpoint (const cGH*);

#else

/* Scheduled functions */
void CarpetIOHDF5Startup (void);
int CarpetIOHDF5Init (const cGH* const);
int CarpetIOHDF5ReadData (const cGH* const);
int CarpetIOHDF5_InitialDataCheckpoint (const cGH* const);
int CarpetIOHDF5_EvolutionCheckpoint (const cGH* const);
int CarpetIOHDF5_TerminationCheckpoint (const cGH* const);
void CarpetIOHDF5EvolutionCheckpoint (const cGH* const);
void CarpetIOHDF5InitialDataCheckpoint (const cGH* const);
    
#endif

int CarpetIOHDF5_Recover (cGH* cgh, const char *basefilename, int called_from);

int CarpetIOHDF5_RecoverParameters (void);
int CarpetIOHDF5_CloseFile (void);

#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOHDF5 */
#endif

#endif /* !defined(CARPETIOHDF5_H) */
