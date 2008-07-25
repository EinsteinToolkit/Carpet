#ifndef DIST_HH
#define DIST_HH

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <mpi.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

#include "cctk.h"

#include "defs.hh"

using namespace std;



namespace dist {
  
  extern MPI_Comm comm_;
  
  extern MPI_Datatype mpi_complex8;
  extern MPI_Datatype mpi_complex16;
  extern MPI_Datatype mpi_complex32;
  
  void init (int& argc, char**& argv);
  void pseudoinit (MPI_Comm const c);
  void finalize ();
  
  // Create MPI datatypes from C structures
  struct mpi_struct_descr_t {
    int          blocklength;
    MPI_Aint     displacement;
    MPI_Datatype type;
  };
  
  void create_mpi_datatype (size_t const count,
                            mpi_struct_descr_t const descr[],
                            MPI_Datatype & newtype);
  
  
  
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
  
  // Set number of threads
  void set_num_threads (int num_threads);
  
  // Local number of threads
  inline int num_threads ()
  {
    static int num_threads_ = -1;
    if (num_threads_ == -1) {
#ifdef _OPENMP
      num_threads_ = omp_get_max_threads();
#else
      num_threads_ = 1;
#endif
      assert (num_threads_ >= 1);
    }
    return num_threads_;
  }
  
  // Global number of threads
  int total_num_threads_worker ();
  inline int total_num_threads ()
  {
    static int total_num_threads_ = -1;
    if (total_num_threads_ == -1) {
      total_num_threads_ = total_num_threads_worker();
    }
    return total_num_threads_;
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
  
#ifdef HAVE_CCTK_COMPLEX8
  inline unsigned int c_datatype (const CCTK_COMPLEX8&)
  { return 13; }
#endif
  
#ifdef HAVE_CCTK_COMPLEX16
  inline unsigned int c_datatype (const CCTK_COMPLEX16&)
  { return 14; }
#endif
  
#ifdef HAVE_CCTK_COMPLEX32
  inline unsigned int c_datatype (const CCTK_COMPLEX32&)
  { return 15; }
#endif
  
  // keep this function's return code consistent with functions above
  inline unsigned int c_ndatatypes ()
  { return 16; }

  template <typename T> unsigned int c_datatype () { abort(); }
  template<> inline unsigned int c_datatype <char>           () { return  0; }
  template<> inline unsigned int c_datatype <signed char>    () { return  1; }
  template<> inline unsigned int c_datatype <unsigned char>  () { return  2; }
  template<> inline unsigned int c_datatype <short>          () { return  3; }
  template<> inline unsigned int c_datatype <unsigned short> () { return  4; }
  template<> inline unsigned int c_datatype <int>            () { return  5; }
  template<> inline unsigned int c_datatype <unsigned int>   () { return  6; }
  template<> inline unsigned int c_datatype <long>           () { return  7; }
  template<> inline unsigned int c_datatype <unsigned long>  () { return  8; }
  template<> inline unsigned int c_datatype <long long>      () { return  9; }
  template<> inline unsigned int c_datatype <float>          () { return 10; }
  template<> inline unsigned int c_datatype <double>         () { return 11; }
  template<> inline unsigned int c_datatype <long double>    () { return 12; }
#ifdef HAVE_CCTK_COMPLEX8
  template<> inline unsigned int c_datatype <CCTK_COMPLEX8>  () { return 13; }
#endif
#ifdef HAVE_CCTK_COMPLEX16
  template<> inline unsigned int c_datatype <CCTK_COMPLEX16> () { return 14; }
#endif
#ifdef HAVE_CCTK_COMPLEX32
  template<> inline unsigned int c_datatype <CCTK_COMPLEX32> () { return 15; }
#endif
  
  /////////////////////////////////////////////////////////////////
  // MPI Datatype helpers
  // Map a C datatype to its corresponding MPI datatype.
  /////////////////////////////////////////////////////////////////
  inline MPI_Datatype datatype (const char&)
  { return MPI_CHAR; }
  
  inline MPI_Datatype datatype (const signed char&)
  { return MPI_CHAR; }
  
  inline MPI_Datatype datatype (const unsigned char&)
  { return MPI_UNSIGNED_CHAR; }
  
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
  
#ifdef HAVE_CCTK_COMPLEX8
  inline MPI_Datatype datatype (const CCTK_COMPLEX8&)
  { return mpi_complex8; }
#endif
  
#ifdef HAVE_CCTK_COMPLEX16
  inline MPI_Datatype datatype (const CCTK_COMPLEX16&)
  { return mpi_complex16; }
#endif
  
#ifdef HAVE_CCTK_COMPLEX32
  inline MPI_Datatype datatype (const CCTK_COMPLEX32&)
  { return mpi_complex32; }
#endif
  
  template <typename T> MPI_Datatype datatype () { abort(); }
  template<> inline MPI_Datatype datatype <char>           () { return MPI_CHAR; }
  template<> inline MPI_Datatype datatype <signed char>    () { return MPI_CHAR; }
  template<> inline MPI_Datatype datatype <unsigned char>  () { return MPI_UNSIGNED_CHAR; }
  template<> inline MPI_Datatype datatype <short>          () { return MPI_SHORT; }
  template<> inline MPI_Datatype datatype <unsigned short> () { return MPI_UNSIGNED_SHORT; }
  template<> inline MPI_Datatype datatype <int>            () { return MPI_INT; }
  template<> inline MPI_Datatype datatype <unsigned int>   () { return MPI_UNSIGNED; }
  template<> inline MPI_Datatype datatype <long>           () { return MPI_LONG; }
  template<> inline MPI_Datatype datatype <unsigned long>  () { return MPI_UNSIGNED_LONG; }
  template<> inline MPI_Datatype datatype <long long>      () { return MPI_LONG_LONG_INT; }
  template<> inline MPI_Datatype datatype <float>          () { return MPI_FLOAT; }
  template<> inline MPI_Datatype datatype <double>         () { return MPI_DOUBLE; }
  template<> inline MPI_Datatype datatype <long double>    () { return MPI_LONG_DOUBLE; }
#ifdef HAVE_CCTK_COMPLEX8
  template<> inline MPI_Datatype datatype <CCTK_COMPLEX8>  () { return mpi_complex8; }
#endif
#ifdef HAVE_CCTK_COMPLEX16
  template<> inline MPI_Datatype datatype <CCTK_COMPLEX16> () { return mpi_complex16; }
#endif
#ifdef HAVE_CCTK_COMPLEX32
  template<> inline MPI_Datatype datatype <CCTK_COMPLEX32> () { return mpi_complex32; }
#endif
  
} // namespace dist



#endif // DIST_HH
