// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dist.cc,v 1.6 2003/01/03 15:49:36 schnetter Exp $

#include <assert.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"

#include "dist.hh"

using namespace std;



namespace dist {
  
  MPI_Comm comm;
  
  MPI_Datatype mpi_complex_float;
  MPI_Datatype mpi_complex_double;
  MPI_Datatype mpi_complex_long_double;
  
  void init (int& argc, char**& argv) {
    MPI_Init (&argc, &argv);
    pseudoinit();
  }
  
  void pseudoinit () {
    comm = MPI_COMM_WORLD;
    
    MPI_Type_contiguous (2, MPI_FLOAT, &mpi_complex_float);
    MPI_Type_commit (&mpi_complex_float);
    MPI_Type_contiguous (2, MPI_DOUBLE, &mpi_complex_double);
    MPI_Type_commit (&mpi_complex_double);
    MPI_Type_contiguous (2, MPI_LONG_DOUBLE, &mpi_complex_long_double);
    MPI_Type_commit (&mpi_complex_long_double);
  }
  
  void finalize () {
    MPI_Finalize ();
  }
  
  void checkpoint (const char* file, int line) {
    DECLARE_CCTK_PARAMETERS;
    if (verbose) {
      int rank;
      MPI_Comm_rank (comm, &rank);
      printf ("CHECKPOINT: processor %d, file %s, line %d\n",
	      rank, file, line);
    }
    if (barriers) {
      MPI_Barrier (comm);
    }
  }
  
} // namespace dist
