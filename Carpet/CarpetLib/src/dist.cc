#include <cassert>

#include <mpi.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"

#include "dist.hh"

using namespace std;



namespace dist {
  
  MPI_Comm comm_ = MPI_COMM_NULL;
  
  MPI_Datatype mpi_complex8;
  MPI_Datatype mpi_complex16;
  MPI_Datatype mpi_complex32;
  
  void init (int& argc, char**& argv) {
    MPI_Init (&argc, &argv);
    pseudoinit (MPI_COMM_WORLD);
  }
  
  void pseudoinit (MPI_Comm const c) {
    comm_ = c;
    
#ifdef HAVE_CCTK_REAL4
    CCTK_REAL4 dummy4;
    MPI_Type_contiguous (2, datatype(dummy4), &mpi_complex8);
    MPI_Type_commit (&mpi_complex8);
#endif
#ifdef HAVE_CCTK_REAL8
    CCTK_REAL8 dummy8;
    MPI_Type_contiguous (2, datatype(dummy8), &mpi_complex16);
    MPI_Type_commit (&mpi_complex16);
#endif
#ifdef HAVE_CCTK_REAL16
    CCTK_REAL16 dummy16;
    MPI_Type_contiguous (2, datatype(dummy16), &mpi_complex32);
    MPI_Type_commit (&mpi_complex32);
#endif
  }
  
  void finalize () {
    MPI_Finalize ();
  }
  
  
  
  // Create an MPI datatype from a C datatype description
  void create_mpi_datatype (size_t const count,
                            mpi_struct_descr_t const descr[],
                            MPI_Datatype & newtype)
  {
    int blocklengths[count];
    MPI_Aint displacements[count];
    MPI_Datatype types[count];
    for (size_t n=0; n<count; ++n) {
      blocklengths [n] = descr[n].blocklength;
      displacements[n] = descr[n].displacement;
      types        [n] = descr[n].type;
    }
    MPI_Type_struct (count, blocklengths, displacements, types, &newtype);
    MPI_Type_commit (&newtype);
  }
  
  
  
  void checkpoint (const char* file, int line) {
    DECLARE_CCTK_PARAMETERS;
    if (verbose) {
      int rank;
      MPI_Comm_rank (comm(), &rank);
      printf ("CHECKPOINT: processor %d, file %s, line %d\n",
	      rank, file, line);
    }
    if (barriers) {
      MPI_Barrier (comm());
    }
  }
  
  // Set number of threads
  void set_num_threads (int const num_threads)
  {
#ifdef _OPENMP
    if (num_threads > 0) {
      // Set number of threads which should be used
      // TODO: do this at startup, not in this routine
      omp_set_num_threads (num_threads);
    }
#else
    if (num_threads > 0 and num_threads != 1) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "OpenMP is not enabled.  Cannot set the number of threads.");
    }
#endif
  }
  
  // Global number of threads
  int total_num_threads_worker ()
  {
    int total_num_threads_;
    int const mynthreads = num_threads();
    MPI_Allreduce
      (const_cast <int *> (& mynthreads), & total_num_threads_, 1, MPI_INT,
       MPI_SUM, comm());
    assert (total_num_threads_ >= size());
    return total_num_threads_;
  }

} // namespace dist
