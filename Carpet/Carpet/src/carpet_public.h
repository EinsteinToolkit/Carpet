/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.h,v 1.1 2001/07/09 09:00:13 schnetter Exp $ */

#include <mpi.h>

#include "cctk.h"

/* Scheduled functions */
int CarpetStartup (void);

/* Helper functions */
MPI_Comm CarpetMPICommunicator (void);
