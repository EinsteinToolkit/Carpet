#ifndef REGION_HH
#define REGION_HH

#include <iostream>
#include <vector>

#include "defs.hh"
#include "bbox.hh"
#include "fulltree.hh"
#include "vect.hh"



// Region description
struct region_t {
  ibbox        extent;           // extent
  b2vect       outer_boundaries; // outer boundaries
  int          map;              // map to which this region belongs
  int          processor;        // processor number
  ipfulltree * processors;       // processor decomposition
  
  region_t ();
  region_t (region_t const & a);
  region_t & operator= (region_t const & a);
  ~region_t ();
  
  bool invariant () const;
};



bool operator== (region_t const & a, region_t const & b);
inline
bool operator!= (region_t const & a, region_t const & b)
{
  return not (a == b);
}



void
combine_regions (vector<region_t> const & oldregs,
                 vector<region_t> & newregs);



size_t memoryof (region_t const & reg);

istream & operator>> (istream & is, region_t       & reg);
ostream & operator<< (ostream & os, region_t const & reg);



// A pseudoregion is almost a region; it is a bbox that belongs to a
// certain component.  Pseudoregions are a compact way to store
// information about what components needs to send data to what other
// components during synchronisation or regridding.
struct pseudoregion_t {
  ibbox extent;
  int component;
  pseudoregion_t ()
  {
  }
  pseudoregion_t (ibbox const & extent_, int const component_)
    : extent (extent_), component (component_)
  {
  }
};

bool operator== (pseudoregion_t const & a, pseudoregion_t const & b);
inline
bool operator!= (pseudoregion_t const & a, pseudoregion_t const & b)
{
  return not (a == b);
}

inline size_t memoryof (pseudoregion_t const & p)
{
  return
    memoryof (p.extent) +
    memoryof (p.component);
}

ostream & operator<< (ostream & os, pseudoregion_t const & p);



struct sendrecv_pseudoregion_t {
  pseudoregion_t send, recv;
  sendrecv_pseudoregion_t ()
  {
  }
  sendrecv_pseudoregion_t (ibbox const & send_extent, int const send_component,
                           ibbox const & recv_extent,  int const recv_component)
    : send (pseudoregion_t (send_extent, send_component)),
      recv (pseudoregion_t (recv_extent, recv_component))
  {
  }
};

inline size_t memoryof (sendrecv_pseudoregion_t const & srp)
{
  return memoryof (srp.send) + memoryof (srp.recv);
}

ostream & operator<< (ostream & os, sendrecv_pseudoregion_t const & srp);



#endif // #ifndef REGION_HH
