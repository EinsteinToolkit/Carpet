/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.4 2001/12/09 16:41:53 schnetter Exp $ */

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"



/* Scheduled functions */
int CarpetParamCheck (CCTK_ARGUMENTS);
int CarpetStartup (void);

/* Helper functions */
MPI_Comm CarpetMPIComm (void);
MPI_Datatype CarpetMPIDatatype (int vartype);
