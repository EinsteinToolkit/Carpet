/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.3 2001/12/07 18:24:17 schnetter Exp $ */

#include <mpi.h>

#include "cctk.h"

/* Scheduled functions */
int CarpetParamCheck (void);
int CarpetStartup (void);

/* Helper functions */
MPI_Comm CarpetMPIComm (void);
MPI_Datatype CarpetMPIDatatype (int vartype);
