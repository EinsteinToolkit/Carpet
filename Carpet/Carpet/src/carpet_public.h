/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.13 2004/04/04 19:24:13 schnetter Exp $ */

#ifndef CARPET_PUBLIC_H
#define CARPET_PUBLIC_H

#include <mpi.h>

#include "cctk.h"



/* Tell thorns that the Carpet routines exist */
#define HAVE_CARPET



#ifdef __cplusplus
namespace Carpet {
  extern "C" {
#endif
    
    /* Prolongation management */
    CCTK_INT CarpetEnableProlongating (const CCTK_INT flag);
    
    
    
    /* Call a schedule group */
    int CallScheduleGroup (cGH * const cgh, const char * const group);
    
    /* Call a local function */
    int CallLocalFunction (cGH * const cgh,
                           void (* const function) (cGH * const cgh));
    int CallSinglemapFunction (cGH * const cgh,
                               void (* const function) (cGH * const cgh));
    int CallLevelFunction (cGH * const cgh,
                           void (* const function) (cGH * const cgh));
    int CallGlobalFunction (cGH * const cgh,
                            void (* const function) (cGH * const cgh));
    int CallMetaFunction (cGH * const cgh,
                          void (* const function) (cGH * const cgh));
    
    
    
    /* Helper functions */
    MPI_Comm CarpetMPIComm (void);
    MPI_Datatype CarpetMPIDatatype (int vartype);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace Carpet */
#endif

#endif /* !defined(CARPET_PUBLIC_H) */
