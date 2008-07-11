#include <cassert>
#include <cstdlib>
#include <iostream>

#include "bboxset.hh"
#include "defs.hh"
#include "region.hh"

using namespace std;



region_t::region_t ()
  : processor (-1), processors (NULL)
{
  assert (invariant());
}

region_t::region_t (region_t const & a)
{
  assert (a.invariant());
  extent = a.extent;
  outer_boundaries = a.outer_boundaries;
  map = a.map;
  processor = a.processor;
  if (a.processors == NULL) {
    processors = NULL;
  } else {
    processors = new ipfulltree (*a.processors);
  }
  assert (invariant());
}

region_t &
region_t::operator= (region_t const & a)
{
  assert (invariant());
  if (processors != NULL) {
    delete processors;
  }
  assert (a.invariant());
  extent = a.extent;
  outer_boundaries = a.outer_boundaries;
  map = a.map;
  processor = a.processor;
  if (a.processors == NULL) {
    processors = NULL;
  } else {
    processors = new ipfulltree (*a.processors);
  }
  assert (invariant());
  return *this;
}


region_t::~region_t ()
{
  assert (invariant());
  if (processors != NULL) {
    delete processors;
  }
}



bool
region_t::invariant () const
{
  if (processor >= 0 and processors != NULL) return false;
  return true;
}



// Compare two regions for equality.
bool
operator== (region_t const & a, region_t const & b)
{
  return
    a.extent == b.extent and
    all (all (a.outer_boundaries == b.outer_boundaries)) and
    a.map == b.map and
    a.processor == b.processor and
    ((a.processors == NULL and b.processors == NULL) or
     (a.processors != NULL and b.processors != NULL and
      *a.processors == *b.processors));
}



// Combine a collection of regions.  Regions can be combined if they
// abutt on boundaries which are not outer boundaries, ignoring the
// processor distribution.  This should lead to a canonical
// representations of collections of regions.
//
// We use vectors to represent the collection, but we could also use
// other containers.  oldregs is read, newregs is added-to.  newregs
// is not cleared.
void
combine_regions (vector<region_t> const & oldregs,
                 vector<region_t> & newregs)
{
  // Find the union of all regions' bounding boxes, and the union of
  // all regions' outer boundaries.  Represent the boundaries as the
  // outermost layer of grid points of the corresponding bounding
  // boxes.
  int const m = oldregs.empty() ? -1 : oldregs.AT(0).map;
  ibset comps;
  ibset cobnds[2][dim];
  for (vector<region_t>::const_iterator
         ri = oldregs.begin(); ri != oldregs.end(); ++ ri)
  {
    region_t const & reg = * ri;
    assert (reg.map == m);
    comps += reg.extent;
    for (int f = 0; f < 2; ++ f) {
      for (int d = 0; d < dim; ++ d) {
        if (reg.outer_boundaries[f][d]) {
          ibbox bnd = reg.extent;
          ivect lo = bnd.lower();
          ivect up = bnd.upper();
          if (f==0) {
            up[d] = lo[d];
          } else {
            lo[d] = up[d];
          }
          bnd = ibbox (lo, up, bnd.stride());
          cobnds[f][d] += bnd;
        }
      }
    }
  }
  comps.normalize();
  for (int f = 0; f < 2; ++ f) {
    for (int d = 0; d < dim; ++ d) {
      cobnds[f][d].normalize();
    }
  }
  
  // Reserve (generous) memory for the result
  size_t const needsize = newregs.size() + comps.setsize();
  if (newregs.capacity() < needsize) {
    newregs.reserve (1000 + 2 * needsize);
  }
  
  // Insert the regions
  for (ibset::const_iterator ci = comps.begin(); ci != comps.end(); ++ ci) {
    ibbox const & c = * ci;
    b2vect obnds;
    for (int f = 0; f < 2; ++ f) {
      for (int d = 0; d < dim; ++ d) {
        obnds[f][d] = cobnds[f][d].intersects (c);
        if (obnds[f][d]) {
          ivect lo = c.lower();
          ivect up = c.upper();
          if (f) lo[d]=up[d]; else up[d]=lo[d];
          ibbox const cbnds (lo, up, c.stride());
          if (not ((cobnds[f][d] & ibset(c)) == ibset(cbnds))) {
            cout << "cobnds[f][d] = " << cobnds[f][d] << endl
                 << "ibset(c) = " << ibset(c) << endl
                 << "(cobnds[f][d] & ibset(c)) = " << (cobnds[f][d] & ibset(c)) << endl
                 << "ibset(cbnds) = " << ibset(cbnds) << endl;
          }
          assert ((cobnds[f][d] & ibset(c)) == ibset(cbnds));
        }
      }
    }
    
    region_t reg;
    reg.extent           = c;
    reg.outer_boundaries = obnds;
    reg.map              = m;
    reg.processor        = -1;
    reg.processors       = NULL;
    newregs.push_back (reg);
  }
}



size_t memoryof (region_t const & reg)
{
  return
    memoryof (reg.extent) +
    memoryof (reg.outer_boundaries) +
    memoryof (reg.map) +
    memoryof (reg.processor) +
    memoryof (reg.processors) +
    (reg.processors != NULL ? memoryof (*reg.processors) : 0);
}



istream &
operator>> (istream & is, region_t & reg)
{
  skipws (is);
  consume (is, "region_t");
  skipws (is);
  consume (is, '(');
  
  skipws (is);
  consume (is, "extent");
  skipws (is);
  consume (is, '=');
  is >> reg.extent;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "outer_boundaries");
  skipws (is);
  consume (is, '=');
  is >> reg.outer_boundaries;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "map");
  skipws (is);
  consume (is, '=');
  is >> reg.map;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "processor");
  skipws (is);
  consume (is, '=');
  is >> reg.processor;
  skipws (is);
  consume (is, ')');
  
  reg.processors = NULL;
  
  return is;
}



ostream &
operator<< (ostream & os, region_t const & reg)
{
  os << "region_t("
     << "extent=" << reg.extent << ","
     << "outer_boundaries=" << reg.outer_boundaries << ","
     << "map=" << reg.map << ","
     << "processor=" << reg.processor << ")";
  return os;
}



// Compare two pseudoregions for equality.
bool
operator== (pseudoregion_t const & a, pseudoregion_t const & b)
{
  return
    a.extent == b.extent and
    a.component == b.component;
}



ostream & operator<< (ostream & os, pseudoregion_t const & p)
{
  return os << p.extent << "/c:" << p.component;
}

ostream & operator<< (ostream & os, sendrecv_pseudoregion_t const & srp)
{
  return os << "(send:" << srp.send << ",recv:" << srp.recv << ")";
}
