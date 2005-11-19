#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dist.hh"

#include "mem.hh"



using namespace std;



struct mstat {
  int total_bytes;
  int total_objects;
  int max_bytes;
  int max_objects;
};



// Total number of currently allocated bytes and objects
static size_t total_allocated_bytes   = 0;
static size_t total_allocated_objects = 0;

// Maximum of the above (over time)
static size_t max_allocated_bytes   = 0;
static size_t max_allocated_objects = 0;



template<typename T>
mem<T>::
mem (size_t const vectorlength, size_t const nelems, T * const memptr)
  : storage_ (memptr),
    nelems_ (nelems),
    vectorlength_ (vectorlength),
    owns_storage_ (false),
    clients_ (vectorlength, false),
    num_clients_ (0)
{
  DECLARE_CCTK_PARAMETERS;
  if (memptr == NULL) {
    const size_t nbytes = vectorlength * nelems * sizeof (T);
    if (max_allowed_memory_MB
        and (total_allocated_bytes + nbytes
             > size_t(1000000) * max_allowed_memory_MB))
    {
      T Tdummy;
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Refusing to allocate %.0f bytes (%.3f MB) of memory for type %s.  %.0f bytes (%.3f MB) are currently allocated in %d objects.  The parameter file specifies a maximum of %d MB",
                  double(nbytes), double(nbytes/1.0e6),
                  typestring(Tdummy),
                  double(total_allocated_bytes),
                  double(total_allocated_bytes/1.0e6),
                  int(total_allocated_objects),
                  int(max_allowed_memory_MB));
    }
    try {
      storage_ = new T [vectorlength * nelems];
      owns_storage_ = true;
    } catch (...) {
      T Tdummy;
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Failed to allocate %.0f bytes (%.3f MB) of memory for type %s.  %.0f bytes (%.3f MB) are currently allocated in %d objects",
                  double(nbytes), double(nbytes/1.0e6),
                  typestring(Tdummy),
                  double(total_allocated_bytes),
                  double(total_allocated_bytes/1.0e6),
                  int(total_allocated_objects));
    }
    total_allocated_bytes += nbytes;
    max_allocated_bytes = max (max_allocated_bytes, total_allocated_bytes);
  }
  ++ total_allocated_objects;
  max_allocated_objects = max (max_allocated_objects, total_allocated_objects);
}

template<typename T>
mem<T>::
~mem ()
{
  assert (! has_clients());
  if (owns_storage_) {
    delete [] storage_;
    total_allocated_bytes -= vectorlength_ * nelems_ * sizeof (T);
  }
  -- total_allocated_objects;
}



template<typename T>
void mem<T>::
register_client (size_t const vectorindex)
{
  assert (vectorindex < vectorlength_);
  assert (! clients_.at(vectorindex));
  clients_.at(vectorindex) = true;
  ++ num_clients_;
}

template<typename T>
void mem<T>::
unregister_client (size_t const vectorindex)
{
  assert (vectorindex < vectorlength_);
  assert (clients_.at(vectorindex));
  clients_.at(vectorindex) = false;
  assert (num_clients_ > 0);
  -- num_clients_;
}

template<typename T>
bool mem<T>::
has_clients () const
{
  // return find (clients_.begin(), clients_.end(), true) != clients_.end();
  return num_clients_ > 0;
}



extern "C" void CarpetLib_printmemstats (CCTK_ARGUMENTS);

void CarpetLib_printmemstats (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (print_memstats_every
      and cctk_iteration % print_memstats_every == 0)
  {
    cout << "Memory statistics from CarpetLib:" << endl
         << "   Current number of objects: " << total_allocated_objects << endl
         << "   Current allocated memory:  "
         << setprecision(3) << total_allocated_bytes / 1.0e6 << " MB" << endl
         << "   Maximum number of objects: " << max_allocated_objects << endl
         << "   Maximum allocated memory:  "
         << setprecision(3) << max_allocated_bytes / 1.0e6 << " MB" << endl
         << endl;
    
    if (strcmp (memstat_file, "") != 0) {
      
      mstat mybuf;
      mybuf.total_bytes   = total_allocated_bytes;
      mybuf.total_objects = total_allocated_objects;
      mybuf.max_bytes     = max_allocated_bytes;
      mybuf.max_objects   = max_allocated_objects;
      vector<mstat> allbuf (dist::size());
      MPI_Gather (& mybuf, 4, MPI_INT,
                  & allbuf.front(), 4, MPI_INT,
                  0, dist::comm());
      
      if (dist::rank() == 0) {
        
        int max_max_bytes = 0;
        int avg_max_bytes = 0;
        int cnt_max_bytes = 0;
        for (size_t n=0; n<allbuf.size(); ++n) {
          max_max_bytes = max (max_max_bytes, allbuf[n].max_bytes);
          avg_max_bytes += allbuf[n].max_bytes;
          ++ cnt_max_bytes;
        }
        avg_max_bytes /= cnt_max_bytes;
        
        ostringstream filenamebuf;
        filenamebuf << out_dir << "/" << memstat_file;
        string const filename = filenamebuf.str();
        ofstream file;
        static bool did_truncate = false;
        if (not did_truncate) {
          did_truncate = true;
          file.open (filename.c_str(), ios::out | ios::trunc);
          if (CCTK_IsFunctionAliased ("UniqueBuildID")) {
            char const * const build_id
              = static_cast<char const *> (UniqueBuildID (cctkGH));
            file << "# Build ID: " << build_id << endl;
          }
          if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
            char const * const job_id
              = static_cast<char const *> (UniqueSimulationID (cctkGH));
            file << "# Simulation ID: " << job_id << endl;
          }
          file << "# Running on " << dist::size() << " processors" << endl;
          file << "#" << endl;
          file << "# iteration maxmaxbytes avgmaxbytes" << endl;
        } else {
          file.open (filename.c_str(), ios::out | ios::app);
        }
        
        file << cctk_iteration << " "
             << max_max_bytes << " " << avg_max_bytes << endl;
        
        file.close ();
        
      } // if on root processor
    } // if output to file
    
  }
}



#define INSTANTIATE(T)                          \
  template class mem<T>;

#include "instantiate"

#undef INSTANTIATE
