#include <cassert>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"

#include "dist.hh"

using namespace std;



namespace dist {
  
  MPI_Comm comm = MPI_COMM_NULL;
  
#if 0
  MPI_Datatype mpi_complex_float;
  MPI_Datatype mpi_complex_double;
  MPI_Datatype mpi_complex_long_double;
#else
  MPI_Datatype mpi_complex8;
  MPI_Datatype mpi_complex16;
  MPI_Datatype mpi_complex32;
#endif
  
  void init (int& argc, char**& argv) {
    MPI_Init (&argc, &argv);
    pseudoinit (MPI_COMM_WORLD);
  }
  
  void pseudoinit (MPI_Comm const c) {
    comm = c;
    
#if 0
    MPI_Type_contiguous (2, MPI_FLOAT, &mpi_complex_float);
    MPI_Type_commit (&mpi_complex_float);
    MPI_Type_contiguous (2, MPI_DOUBLE, &mpi_complex_double);
    MPI_Type_commit (&mpi_complex_double);
    MPI_Type_contiguous (2, MPI_LONG_DOUBLE, &mpi_complex_long_double);
    MPI_Type_commit (&mpi_complex_long_double);
#else
#  ifdef CCTK_REAL4
    CCTK_REAL4 dummy4;
    MPI_Type_contiguous (2, datatype(dummy4), &mpi_complex8);
    MPI_Type_commit (&mpi_complex8);
#  endif
#  ifdef CCTK_REAL8
    CCTK_REAL8 dummy8;
    MPI_Type_contiguous (2, datatype(dummy8), &mpi_complex16);
    MPI_Type_commit (&mpi_complex16);
#  endif
#  ifdef CCTK_REAL16
    CCTK_REAL16 dummy16;
    MPI_Type_contiguous (2, datatype(dummy16), &mpi_complex32);
    MPI_Type_commit (&mpi_complex32);
#  endif
#endif
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
