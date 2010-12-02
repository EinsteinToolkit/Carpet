#include <cctk.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <typeinfo>

#ifdef CCTK_MPI
#  include <mpi.h>
#else
#  include "nompi.h"
#endif
#ifdef _OPENMP
#  include <omp.h>
#endif

#include "defs.hh"
#include "limits.hh"
#include "startup_time.hh"

#include "dist.hh"

using namespace std;



namespace dist {
  
  MPI_Comm comm_ = MPI_COMM_NULL;
  
  MPI_Datatype mpi_complex8  = MPI_DATATYPE_NULL;
  MPI_Datatype mpi_complex16 = MPI_DATATYPE_NULL;
  MPI_Datatype mpi_complex32 = MPI_DATATYPE_NULL;
  
  int total_num_threads_ = -1;
  
  void init (int& argc, char**& argv) {
    MPI_Init (&argc, &argv);
    pseudoinit (MPI_COMM_WORLD);
  }
  
  void pseudoinit (MPI_Comm const c) {
    comm_ = c;
    
#ifdef HAVE_CCTK_REAL4
    CCTK_REAL4 dummy4;
    MPI_Type_contiguous (2, mpi_datatype(dummy4), &mpi_complex8);
    MPI_Type_commit (&mpi_complex8);
#endif
#ifdef HAVE_CCTK_REAL8
    CCTK_REAL8 dummy8;
    MPI_Type_contiguous (2, mpi_datatype(dummy8), &mpi_complex16);
    MPI_Type_commit (&mpi_complex16);
#endif
#ifdef HAVE_CCTK_REAL16
    CCTK_REAL16 dummy16;
    MPI_Type_contiguous (2, mpi_datatype(dummy16), &mpi_complex32);
    MPI_Type_commit (&mpi_complex32);
#endif
    
    // Output startup time
    CarpetLib::output_startup_time ();
    
    // Check and/or modify system limits
    CarpetLib::set_system_limits ();
    
    collect_total_num_threads ();
  }
  
  void finalize () {
    MPI_Finalize ();
  }
  
  
  
  // Create an MPI datatype from a C datatype description
  
  ostream& operator<< (ostream& os, mpi_struct_descr_t const& descr)
  {
    int type_size;
    MPI_Type_size (descr.type, &type_size);
    os << "{"
       << "blocklength:" << descr.blocklength << ","
       << "displacement:" << descr.displacement << ","
       << "type:" << descr.type << ","
       << "type_size:" << type_size << ","
       << "field_name:" << descr.field_name << ","
       << "type_name:" << descr.type_name
       << "}";
    return os;
  }
  
  MPI_Datatype create_mpi_datatype (size_t const count,
                                    mpi_struct_descr_t const descr[],
                                    char const * const name, size_t const size)
  {
    DECLARE_CCTK_PARAMETERS;
    int blocklengths[count];
    MPI_Aint displacements[count];
    MPI_Datatype types[count];
    for (size_t n=0; n<count; ++n) {
      blocklengths [n] = descr[n].blocklength;
      displacements[n] = descr[n].displacement;
      types        [n] = descr[n].type;
    }
    MPI_Datatype newtype;
    MPI_Type_struct (count, blocklengths, displacements, types, &newtype);
    MPI_Type_commit (&newtype);
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Creating new MPI type for C type %s:", name);
      cout << "   Type has " << count << " components" << endl;
      for (size_t n=0; n<count; ++n) {
        cout << "   [" << n << "]: " << descr[n] << endl;
      }
      cout << "   New MPI type ID is " << newtype << endl;
      int datatypesize;
      MPI_Type_size (newtype, &datatypesize);
      cout << "   C type size is " << size << endl;
      cout << "   MPI type size is " << datatypesize << endl;
    }
    return newtype;
  }
  
#if 0
  
  ostream&
  generic_mpi_datatype_t::field_t::output (ostream& os) const
  {
    int type_size;
    MPI_Type_size (mpi_datatype, &type_size);
    os << "{"
       << "offset:" << offset << ","
       << "count:" << count << ","
       << "mpi_datatype:" << mpi_datatype << ","
       << "type_size:" << type_size << ","
       << "field_name:" << field_name << ","
       << "type_name:" << type_name
       << "}";
    return os;
  }
  
  generic_mpi_datatype_t::generic_mpi_datatype_t (string const type_name_)
    : type_name (type_name_), type_is_committed (false)
  {
  }
  
  template <typename U>
  void
  generic_mpi_datatype_t::add_field (size_t const offset, size_t const count,
                                     string const field_name)
  {
    assert (not type_is_committed);
    U u;
    entries.push_back (field_t (offset, count, mpi_datatype(u),
                                field_name, typeid(U).name()));
  }
  
  void
  generic_mpi_datatype_t::commit ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Debug output
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Creating new MPI type for C type %s:", type_name.c_str());
      cout << *this;
    }
    
    assert (not type_is_committed);
    type_is_committed = true;
    
    // Out of caution -- this could be allowed
    assert (not entries.empty());
    
    // Create MPI type
    size_t const count = entries.size();
    int          blocklengths [count+1];
    MPI_Aint     displacements[count+1];
    MPI_Datatype types        [count+1];
    {
      size_t n = 0;
      for (list<field_t>::const_iterator ifield =
             entries.begin(); ifield!=entries.end(); ++ifield, ++n)
      {
        blocklengths [n] = ifield->count;
        displacements[n] = ifield->offset;
        types        [n] = ifield->mpi_datatype;
      }
      assert (n == count);
      // Add MPI_UB
      blocklengths [n] = 1;
      displacements[n] = type_size();
      types        [n] = MPI_UB;
    }
    
    MPI_Type_struct
      (count+1, blocklengths, displacements, types, &mpi_datatype);
    MPI_Type_commit (&mpi_datatype);
  }
  
  ostream&
  generic_mpi_datatype_t::output (ostream& os) const
  {
    cout << "Datatype: " << type_name << endl;
    size_t const count = entries.size();
    cout << "   Type has " << count << " components" << endl;
    {
      size_t n = 0;
      for (list<field_t>::const_iterator ifield =
             entries.begin(); ifield!=entries.end(); ++ifield, ++n)
      {
        cout << "   [" << n << "]: " << *ifield << endl;
      }
      assert (n == count);
    }
    cout << "   MPI type ID: " << mpi_datatype << endl;
    int datatypesize;
    MPI_Type_size (mpi_datatype, &datatypesize);
    cout << "   C   type size: " << size << endl;
    cout << "   MPI type size: " << datatypesize << endl;
    return os;
  }
  
#endif
  
  
  
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
  void collect_total_num_threads ()
  {
    int const mynthreads = num_threads();
    // cerr << "QQQ: collect_total_num_threads[1]" << endl;
    MPI_Allreduce
      (const_cast <int *> (& mynthreads), & total_num_threads_, 1, MPI_INT,
       MPI_SUM, comm());
    // cerr << "QQQ: collect_total_num_threads[2]" << endl;
    assert (total_num_threads_ >= size());
  }
  
  
  
  char const * c_datatype_name (unsigned const type)
  {
    switch (type) {
    case  0: return "char";
    case  1: return "signed char";
    case  2: return "unsigned char";
    case  3: return "short";
    case  4: return "unsigned short";
    case  5: return "int";
    case  6: return "unsigned int";
    case  7: return "long";
    case  8: return "unsigned long";
    case  9: return "long long";
    case 10: return "unsigned long long";
    case 11: return "float";
    case 12: return "double";
    case 13: return "long double";
#ifdef HAVE_CCTK_COMPLEX8
    case 14: return "CCTK_COMPLEX8";
#endif
#ifdef HAVE_CCTK_COMPLEX16
    case 15: return "CCTK_COMPLEX16";
#endif
#ifdef HAVE_CCTK_COMPLEX32
    case 16: return "CCTK_COMPLEX32";
#endif
    }
    assert (0); abort();
    return NULL;
  }
  
} // namespace dist
