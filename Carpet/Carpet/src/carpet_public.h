/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.2 2001/12/05 03:31:56 schnetter Exp $ */

#include <mpi.h>

#include "cctk.h"

/* Scheduled functions */
int CarpetParamCheck (void);
int CarpetStartup (void);

/* Helper functions */
MPI_Comm CarpetMPICommunicator (void);
