/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.5 2001/12/14 16:39:09 schnetter Exp $ */

#ifndef CARPET_PUBLIC_H
#define CARPET_PUBLIC_H

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"



/* Scheduled functions */
int CarpetParamCheck (CCTK_ARGUMENTS);
int CarpetStartup (void);

/* Helper functions */
MPI_Comm CarpetMPIComm (void);
MPI_Datatype CarpetMPIDatatype (int vartype);

#endif /* ! defined(CARPET_PUBLIC_H) */
