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
  b2vect refinement_boundaries;             // refinement boundaries
  int    map;                               // map to which this
                                            // region belongs
  int    processor;                         // processor number
};



bool operator== (region_t const & a, region_t const & b);
inline
bool operator!= (region_t const & a, region_t const & b)
{
  return not operator== (a, b);
}



istream & operator>> (istream & is, region_t       & reg);
ostream & operator<< (ostream & os, region_t const & reg);



#endif // #ifndef REGION_HH
