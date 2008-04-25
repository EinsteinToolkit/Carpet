#ifndef REGION_HH
#define REGION_HH

#include <iostream>

#include "defs.hh"
#include "bbox.hh"
#include "vect.hh"



// Region description
struct region_t {
  ibbox  extent;                            // extent
  b2vect outer_boundaries;                  // outer boundaries
  int    map;                               // map to which this
                                            // region belongs
  int    processor;                         // processor number
};



bool operator== (region_t const & a, region_t const & b);
inline
bool operator!= (region_t const & a, region_t const & b)
{
  return not (a == b);
}



inline size_t memoryof (region_t const & reg)
{
  return
    memoryof (reg.extent) +
    memoryof (reg.outer_boundaries) +
    memoryof (reg.map) +
    memoryof (reg.processor);
}

istream & operator>> (istream & is, region_t       & reg);
ostream & operator<< (ostream & os, region_t const & reg);



// A pseudoregion is almost a region; it is a bbox that lives on a
// certain processor.  Pseudoregions are a compact way to store
// information about what processors needs to send data to what other
// processors during synchronisation or regridding.
struct pseudoregion_t {
  ibbox extent;
  int processor;
  pseudoregion_t ()
  {
  }
  pseudoregion_t (ibbox const & extent_, int const processor_)
    : extent (extent_), processor (processor_)
  {
  }
};

inline size_t memoryof (pseudoregion_t const & p)
{
  return
    memoryof (p.extent) +
    memoryof (p.processor);
}

ostream & operator<< (ostream & os, pseudoregion_t const & p);

struct sendrecv_pseudoregion_t {
  pseudoregion_t send, recv;
  sendrecv_pseudoregion_t ()
  {
  }
  sendrecv_pseudoregion_t (ibbox const & send_extent, int const send_processor,
                           ibbox const & recv_extent, int const recv_processor)
    : send (pseudoregion_t (send_extent, send_processor)),
      recv (pseudoregion_t (recv_extent, recv_processor))
  {
  }
};

inline size_t memoryof (sendrecv_pseudoregion_t const & srp)
{
  return memoryof (srp.send) + memoryof (srp.recv);
}

ostream & operator<< (ostream & os, sendrecv_pseudoregion_t const & srp);



#endif // #ifndef REGION_HH
