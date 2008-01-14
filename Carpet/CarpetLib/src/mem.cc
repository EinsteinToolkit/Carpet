#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#include <sys/resource.h>
#include <sys/time.h>

#include "defs.hh"
#include "dist.hh"

#include "mem.hh"



using namespace std;



struct mstat {
  // Carpet statistics
  double total_bytes;
  double total_objects;
  double max_bytes;
  double max_objects;
  // malloc statistics
  double malloc_used_bytes;
  double malloc_free_bytes;
};
int const mstat_entries = sizeof(mstat) / sizeof(double);



// Total number of currently allocated bytes and objects
static double total_allocated_bytes   = 0;
static double total_allocated_objects = 0;

// Maximum of the above (over time)
static double max_allocated_bytes   = 0;
static double max_allocated_objects = 0;



// TODO: Make this a plain class instead of a template

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



extern "C" void CarpetLib_setmemlimit (CCTK_ARGUMENTS);

void CarpetLib_setmemlimit (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  // Set address space limit
  struct rlimit aslimit;
  {
    int const ierr = getrlimit (RLIMIT_AS, & aslimit);
    assert (not ierr);
  }
  CCTK_VInfo (CCTK_THORNSTRING,
              "Old address space size limit: hard=%lld, soft=%lld",
              (long long) aslimit.rlim_max, (long long) aslimit.rlim_cur);
  if (max_allowed_memory_MB > 0) {
    aslimit.rlim_cur = max_allowed_memory_MB * 1000000LL;
  }
  {
    int const ierr = setrlimit (RLIMIT_AS, & aslimit);
    assert (not ierr);
  }
  {
    int const ierr = getrlimit (RLIMIT_AS, & aslimit);
    assert (not ierr);
  }
  CCTK_VInfo (CCTK_THORNSTRING,
              "Old address space size limit: hard=%lld, soft=%lld",
              (long long) aslimit.rlim_max, (long long) aslimit.rlim_cur);
  CCTK_VInfo (CCTK_THORNSTRING,
              "(Unlimited address space size indicated by %lld)",
              (long long) RLIM_INFINITY);
}



extern "C" void CarpetLib_printmemstats (CCTK_ARGUMENTS);

void CarpetLib_printmemstats (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int const ioproc = 0;
  
  if (print_memstats_every > 0
      and cctk_iteration % print_memstats_every == 0)
  {
    mstat mybuf;
    mybuf.total_bytes   = total_allocated_bytes;
    mybuf.total_objects = total_allocated_objects;
    mybuf.max_bytes     = max_allocated_bytes;
    mybuf.max_objects   = max_allocated_objects;
#ifdef HAVE_MALLINFO
    struct mallinfo minfo = mallinfo ();
    mybuf.malloc_used_bytes = minfo.uordblks;
    mybuf.malloc_free_bytes = minfo.fordblks;
#else
    mybuf.malloc_used_bytes = 0;
    mybuf.malloc_free_bytes = 0;
#endif
    
    cout << "Memory statistics from CarpetLib:" << eol
         << "   Current number of objects: " << total_allocated_objects << eol
         << "   Current allocated memory:  "
         << setprecision(3) << total_allocated_bytes / 1.0e6 << " MB" << eol
         << "   Maximum number of objects: " << max_allocated_objects << eol
         << "   Maximum allocated memory:  "
         << setprecision(3) << max_allocated_bytes / 1.0e6 << " MB" << eol
         << "   Total allocated used system memory: "
         << setprecision(3) << mybuf.malloc_used_bytes / 1.0e6 << " MB" << eol
         << "   Total allocated free system memory: "
         << setprecision(3) << mybuf.malloc_free_bytes / 1.0e6 << " MB" << endl;
    
    if (strcmp (memstat_file, "") != 0) {
      vector<mstat> allbuf (dist::size());
      MPI_Gather (& mybuf, mstat_entries, MPI_DOUBLE,
                  & allbuf.front(), mstat_entries, MPI_DOUBLE,
                  ioproc, dist::comm());
      
      if (dist::rank() == ioproc) {
        
        double max_max_bytes = 0;
        double avg_max_bytes = 0;
        double cnt_max_bytes = 0;
        double max_used_bytes = 0;
        double avg_used_bytes = 0;
        double cnt_used_bytes = 0;
        double max_free_bytes = 0;
        double avg_free_bytes = 0;
        double cnt_free_bytes = 0;
        for (size_t n=0; n<allbuf.size(); ++n) {
          max_max_bytes = max (max_max_bytes, allbuf[n].max_bytes);
          avg_max_bytes += allbuf[n].max_bytes;
          ++ cnt_max_bytes;
          max_used_bytes = max (max_used_bytes, allbuf[n].malloc_used_bytes);
          avg_used_bytes += allbuf[n].malloc_used_bytes;
          ++ cnt_used_bytes;
          max_free_bytes = max (max_free_bytes, allbuf[n].malloc_free_bytes);
          avg_free_bytes += allbuf[n].malloc_free_bytes;
          ++ cnt_free_bytes;
        }
        avg_max_bytes /= cnt_max_bytes;
        avg_used_bytes /= cnt_used_bytes;
        avg_free_bytes /= cnt_free_bytes;
        
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
            file << "# Build ID: " << build_id << eol;
          }
          if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
            char const * const job_id
              = static_cast<char const *> (UniqueSimulationID (cctkGH));
            file << "# Simulation ID: " << job_id << eol;
          }
          file << "# Running on " << dist::size() << " processors" << eol;
          file << "#" << eol;
          file << "# iteration   maxmaxbytes avgmaxbytes   maxusedbytes avgusedbytes   maxfreebytes avgfreebytes" << eol;
        } else {
          file.open (filename.c_str(), ios::out | ios::app);
        }
        
        file << cctk_iteration
              << "\t "<< max_max_bytes << " " << avg_max_bytes
              << "\t "<< max_used_bytes << " " << avg_used_bytes
              << "\t "<< max_free_bytes << " " << avg_free_bytes
             << eol;
        
        file.close ();
        
      } // if on root processor
    } // if output to file
    
  }
}



#define INSTANTIATE(T)                          \
  template class mem<T>;

#include "instantiate"

#undef INSTANTIATE
