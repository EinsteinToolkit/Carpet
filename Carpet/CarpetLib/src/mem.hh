#ifndef MEM_HH
#define MEM_HH

#include <cstdlib>
#include <vector>

using namespace std;

// A chunk of memory, possibly shared between some clients
template<typename T>
class mem
{
  T * storage_;
  size_t nelems_;
  size_t vectorlength_;
  bool owns_storage_;
  
  vector<bool> clients_;
  size_t num_clients_;
  
public:
  mem (size_t vectorlength, size_t nelems, T * memptr = NULL);
  ~mem ();
  
  T * storage (size_t vectorindex) const
  {
    assert (vectorindex < vectorlength_);
    assert (clients_.at(vectorindex));
    return & storage_ [vectorindex * nelems_];
  }
  
  void register_client (size_t vectorindex);
  void unregister_client (size_t vectorindex);
  bool has_clients () const;
};

#endif  // ifndef MEM_HH
