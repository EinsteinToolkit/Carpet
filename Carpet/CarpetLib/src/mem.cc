#include <algorithm>
#include <cassert>

#include "cctk.h"

#include "defs.hh"

#include "mem.hh"



// Total number of currently allocated bytes and objects
static size_t total_allocated_bytes   = 0;
static size_t total_allocated_objects = 0;



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
  if (memptr == NULL) {
    const size_t nbytes = vectorlength * nelems * sizeof (T);
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
  }
  ++ total_allocated_objects;
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



#define INSTANTIATE(T)                          \
  template class mem<T>;

#include "instantiate"

#undef INSTANTIATE
