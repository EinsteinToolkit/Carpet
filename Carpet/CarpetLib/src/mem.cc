#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "defs.hh"
#include "mem.hh"



using namespace std;



double const gmem::MEGA = 1000*1000;

// Total number of currently allocated bytes and objects
double gmem::total_allocated_bytes   = 0;
double gmem::total_allocated_objects = 0;

// Maximum of the above (over time)
double gmem::max_allocated_bytes   = 0;
double gmem::max_allocated_objects = 0;



// TODO: Make this a plain class instead of a template

template<typename T>
mem<T>::
mem (size_t const vectorlength, size_t const nelems,
     T * const memptr, size_t const memsize)
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
        and (total_allocated_bytes + nbytes > MEGA * max_allowed_memory_MB))
    {
      T Tdummy;
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Refusing to allocate %.0f bytes (%.3f MB) of memory for type %s.  %.0f bytes (%.3f MB) are currently allocated in %d objects.  The parameter file specifies a maximum of %d MB",
                  double(nbytes), double(nbytes/MEGA),
                  typestring(Tdummy),
                  double(total_allocated_bytes),
                  double(total_allocated_bytes/MEGA),
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
                  double(nbytes), double(nbytes/MEGA),
                  typestring(Tdummy),
                  double(total_allocated_bytes),
                  double(total_allocated_bytes/MEGA),
                  int(total_allocated_objects));
    }
    total_allocated_bytes += nbytes;
    max_allocated_bytes = max (max_allocated_bytes, total_allocated_bytes);
    if (poison_new_memory) {
      memset (storage_, poison_value, vectorlength * nelems * sizeof (T));
    }
  } else {
    assert (memsize >= vectorlength * nelems * sizeof (T));
    // Don't poison the memory.  Passing in a pointer allows the
    // pointer to be re-interpreted as a mem object, keeping the
    // previous content.  This is e.g. used to turn communication
    // buffers into mem objects.
  }
  ++ total_allocated_objects;
  max_allocated_objects = max (max_allocated_objects, total_allocated_objects);
}

template<typename T>
mem<T>::
~mem ()
{
  assert (not has_clients());
  if (owns_storage_) {
    delete [] storage_;
    const double nbytes = vectorlength_ * nelems_ * sizeof (T);
    total_allocated_bytes -= nbytes;
  }
  -- total_allocated_objects;
}



template<typename T>
void mem<T>::
register_client (size_t const vectorindex)
{
  assert (vectorindex < vectorlength_);
  assert (not clients_.AT(vectorindex));
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



// Memory usage
template<typename T>
size_t
mem<T>::
memory ()
  const
{
  return
    memoryof (storage_) +
    memoryof (nelems_) +
    memoryof (vectorlength_) +
    memoryof (owns_storage_) +
    memoryof (clients_) +
    memoryof (num_clients_) +
    (owns_storage_ ? sizeof (T) * vectorlength_ * nelems_ : 0);
}



size_t const mempool::chunksize;
size_t const mempool::align;

mempool::
mempool ()
  : allocated (0), freeptr (0), freesize (0)
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
    allocated += freesize;
    if (not freeptr) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Failed to allocate %.3f MB of memory",
                  double(freesize/gmem::MEGA));
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



// Memory usage
size_t
mempool::
memory ()
  const
{
  return
    memoryof (chunks) +
    memoryof (freeptr) +
    memoryof (freesize) +
    memoryof (allocated);
}



#define TYPECASE(N,T)                           \
  template class mem<T>;

#include "typecase.hh"

#undef TYPECASE
