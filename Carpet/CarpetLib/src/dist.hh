#ifndef DIST_HH
#define DIST_HH

#include <cassert>
#include <cstdio>
#include <cstdlib>

#if 0
#include <complex>
#endif

#include <mpi.h>

#include "cctk.h"

#include "defs.hh"

using namespace std;



namespace dist {
  
  extern MPI_Comm comm_;
  
#if 0
  extern MPI_Datatype mpi_complex_float;
  extern MPI_Datatype mpi_complex_double;
  extern MPI_Datatype mpi_complex_long_double;
#else
  extern MPI_Datatype mpi_complex8;
  extern MPI_Datatype mpi_complex16;
  extern MPI_Datatype mpi_complex32;
#endif
  
  void init (int& argc, char**& argv);
  void pseudoinit (MPI_Comm const c);
  void finalize ();
  
  // Debugging output
#define CHECKPOINT dist::checkpoint(__FILE__, __LINE__)
  void checkpoint (const char* file, int line);
  
  
  
  // Information about the communicator
  
  // Return the communicator
  inline MPI_Comm comm ()
  {
    return comm_;
  }
  
  // Always return a good communicator
  inline MPI_Comm goodcomm ()
  {
    return comm_ != MPI_COMM_NULL ? comm_ : MPI_COMM_WORLD;
  }
  
  // Rank in the communicator (this processor's number, 0 .. size-1)
  inline int rank ()
  {
    static int rank_ = -1;
    if (rank_ == -1) MPI_Comm_rank (comm(), &rank_);
    return rank_;
  }
  
  // Size of the communicator
  inline int size ()
  {
    static int size_ = -1;
    if (size_ == -1) MPI_Comm_size (comm(), &size_);
    return size_;
  }
  
  
  
  /////////////////////////////////////////////////////////////////////////
  // C Datatype helpers
  // Map a C datatype to a 0-based index running up to c_ndatatypes().
  /////////////////////////////////////////////////////////////////////////
  inline unsigned int c_datatype (const char&)
  { return 0; }
  
  inline unsigned int c_datatype (const signed char&)
  { return 1; }
  
  inline unsigned int c_datatype (const unsigned char&)
  { return 2; }
  
  inline unsigned int c_datatype (const short&)
  { return 3; }
  
  inline unsigned int c_datatype (const unsigned short&)
  { return 4; }
  
  inline unsigned int c_datatype (const int&)
  { return 5; }
  
  inline unsigned int c_datatype (const unsigned int&)
  { return 6; }
  
  inline unsigned int c_datatype (const long&)
  { return 7; }
  
  inline unsigned int c_datatype (const unsigned long&)
  { return 8; }
  
  inline unsigned int c_datatype (const long long&)
  { return 9; }
  
  inline unsigned int c_datatype (const float&)
  { return 10; }
  
  inline unsigned int c_datatype (const double&)
  { return 11; }
  
  inline unsigned int c_datatype (const long double&)
  { return 12; }
  
#if 0
  
  inline unsigned int c_datatype (const complex<float>&)
  { return 13; }
  
  inline unsigned int c_datatype (const complex<double>&)
  { return 14; }
  
  inline unsigned int c_datatype (const complex<long double>&)
  { return 15; }
  
#else
  
#  ifdef CCTK_REAL4
  inline unsigned int c_datatype (const CCTK_COMPLEX8&)
  { return 13; }
#  endif
  
#  ifdef CCTK_REAL8
  inline unsigned int c_datatype (const CCTK_COMPLEX16&)
  { return 14; }
#  endif
  
#  ifdef CCTK_REAL16
  inline unsigned int c_datatype (const CCTK_COMPLEX32&)
  { return 15; }
#  endif
  
#endif
  
  // keep this function's return code consistent with functions above
  inline unsigned int c_ndatatypes ()
  { return 16; }


  /////////////////////////////////////////////////////////////////
  // MPI Datatype helpers
  // Map a C datatype to its corresponding MPI datatype.
  /////////////////////////////////////////////////////////////////
  inline MPI_Datatype datatype (const char&)
  { return MPI_CHAR; }
  
  inline MPI_Datatype datatype (const signed char&)
  { return MPI_UNSIGNED_CHAR; }
  
  inline MPI_Datatype datatype (const unsigned char&)
  { return MPI_BYTE; }
  
  inline MPI_Datatype datatype (const short&)
  { return MPI_SHORT; }
  
  inline MPI_Datatype datatype (const unsigned short&)
  { return MPI_UNSIGNED_SHORT; }
  
  inline MPI_Datatype datatype (const int&)
  { return MPI_INT; }
  
  inline MPI_Datatype datatype (const unsigned int&)
  { return MPI_UNSIGNED; }
  
  inline MPI_Datatype datatype (const long&)
  { return MPI_LONG; }
  
  inline MPI_Datatype datatype (const unsigned long&)
  { return MPI_UNSIGNED_LONG; }
  
  inline MPI_Datatype datatype (const long long&)
  { return MPI_LONG_LONG_INT; }
  
  inline MPI_Datatype datatype (const float&)
  { return MPI_FLOAT; }
  
  inline MPI_Datatype datatype (const double&)
  { return MPI_DOUBLE; }
  
  inline MPI_Datatype datatype (const long double&)
  { return MPI_LONG_DOUBLE; }
  
#if 0
  
  inline MPI_Datatype datatype (const complex<float>&)
  { return mpi_complex_float; }
  
  inline MPI_Datatype datatype (const complex<double>&)
  { return mpi_complex_double; }
  
  inline MPI_Datatype datatype (const complex<long double>&)
  { return mpi_complex_long_double; }
  
#else
  
#  ifdef CCTK_REAL4
  inline MPI_Datatype datatype (const CCTK_COMPLEX8&)
  { return mpi_complex8; }
#  endif
  
#  ifdef CCTK_REAL8
  inline MPI_Datatype datatype (const CCTK_COMPLEX16&)
  { return mpi_complex16; }
#  endif
  
#  ifdef CCTK_REAL16
  inline MPI_Datatype datatype (const CCTK_COMPLEX32&)
  { return mpi_complex32; }
#  endif
  
#endif
  
} // namespace dist



#endif // DIST_HH
