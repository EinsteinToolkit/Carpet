/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.6 2002/09/01 14:52:24 schnetter Exp $ */

#ifndef CARPET_PUBLIC_H
#define CARPET_PUBLIC_H

#include <mpi.h>

#include "cctk_Arguments.h"



#ifdef __cplusplus
namespace Carpet {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetParamCheck (CCTK_ARGUMENTS);
    int CarpetStartup (void);
    
    /* Helper functions */
    MPI_Comm CarpetMPIComm (void);
    MPI_Datatype CarpetMPIDatatype (int vartype);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace Carpet */
#endif

#endif /* !defined(CARPET_PUBLIC_H) */
