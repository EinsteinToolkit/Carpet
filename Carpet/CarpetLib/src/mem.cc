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
  double total_bytes;
  double total_objects;
  double max_bytes;
  double max_objects;
};



// Total number of currently allocated bytes and objects
static double total_allocated_bytes   = 0;
static double total_allocated_objects = 0;

// Maximum of the above (over time)
static double max_allocated_bytes   = 0;
static double max_allocated_objects = 0;



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
    const double nbytes = vectorlength * nelems * sizeof (T);
    if (max_allowed_memory_MB > 0
        and (total_allocated_bytes + nbytes > 1.0e6 * max_allowed_memory_MB))
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
      if (poison_new_memory) {
        memset (storage_, poison_value, vectorlength * nelems * sizeof (T));
      }
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
  assert (! clients_.AT(vectorindex));
  clients_.AT(vectorindex) = true;
  ++ num_clients_;
}

template<typename T>
void mem<T>::
unregister_client (size_t const vectorindex)
{
  assert (vectorindex < vectorlength_);
  assert (clients_.AT(vectorindex));
  clients_.AT(vectorindex) = false;
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



size_t const mempool::chunksize;
size_t const mempool::align;

mempool::
mempool ()
  : freeptr (0), freesize (0)
{
}

mempool::
~mempool ()
{
  while (not chunks.empty()) {
    free (chunks.top());
    chunks.pop();
  }
}

void *
mempool::
alloc (size_t nbytes)
{
  // Take a shortcut for silly requests
  if (nbytes == 0) return 0;
  
  // Round up request size
  nbytes = (nbytes + align - 1) / align * align;
  
  // If there is not enough memory left, allocate a new chunk.  Ignore
  // whatever is left in the old chunk.
  if (nbytes > freesize) {
    // Allocate the usual chunk size, or more if more is requested
    freesize = max (chunksize, nbytes);
    freeptr = malloc (freesize);
    if (not freeptr) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Failed to allocate %.3f MB of memory",
                  double(freesize/1.0e6));
    }
    // Remember the pointer so that it can be freed
    chunks.push (freeptr);
  }
  
  // Allocate a piece from the current chunk
  void * const ptr = freeptr;
  assert (freesize >= nbytes);
  freesize -= nbytes;
  assert (freeptr);
  freeptr = static_cast <char *> (freeptr) + nbytes;
  
  return ptr;
}



extern "C" void CarpetLib_printmemstats (CCTK_ARGUMENTS);

void CarpetLib_printmemstats (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (print_memstats_every > 0
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
      MPI_Gather (& mybuf, 4, MPI_DOUBLE,
                  & allbuf.front(), 4, MPI_DOUBLE,
                  0, dist::comm());
      
      if (dist::rank() == 0) {
        
        double max_max_bytes = 0;
        double avg_max_bytes = 0;
        double cnt_max_bytes = 0;
        for (size_t n=0; n<allbuf.size(); ++n) {
          max_max_bytes = max (max_max_bytes, (double) allbuf[n].max_bytes);
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
