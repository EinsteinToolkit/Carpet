#ifndef GH_HH
#define GH_HH

#include <cassert>
#include <iostream>
#include <list>
#include <vector>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "dist.hh"
#include "region.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
class dh;
class th;
class gh;



// A refinement hierarchy, where higher levels are finer than the base
// level.  The extents do not include ghost zones.
class gh {
  
public:
  
  // Types
  typedef vector<region_t> cregs; // ... for each component
  typedef vector<cregs> rregs;    // ... for each refinement level
  typedef vector<rregs> mregs;    // ... for each multigrid level
  
public:				// should be readonly
  
  // Fields
  const vector<ivect> reffacts; // refinement factors
  const centering refcent;	// vertex or cell centered
 
  const int mgfact;		// default multigrid factor
  const centering mgcent;	// default (vertex or cell centered)
  
  vector<vector<ibbox> > baseextents; // [ml][rl]
  const i2vect boundary_width;
  
  // Extents of the regions before distributing them over the
  // processors
  rregs superregions;
  
  mregs regions;                // extents and properties of all grids
  mregs oldregions;             // a copy, used during regridding
  
  list<th*> ths;		// list of all time hierarchies
  list<dh*> dhs;		// list of all data hierarchies
  
public:
  
  // Constructors
  gh (vector<ivect> const & reffacts, centering refcent,
      int mgfact, centering mgcent,
      vector<vector<ibbox> > const & baseextents, 
      i2vect const & boundary_width);
  
  // Destructors
  ~gh ();
  
  // Modifiers
  void regrid (rregs const & superregs, mregs const & regs);
  bool recompose (int rl, bool do_prolongate);
  
private:
  
  bool level_did_change (int rl) const;
  
  // Accessors
  
public:
  
  ibbox const & extent (const int ml, const int rl, const int c) const
  {
    return regions.AT(ml).AT(rl).AT(c).extent;
  }
  
  ibbox const & baseextent (const int ml, const int rl) const
  {
    return baseextents.AT(ml).AT(rl);
  }
  
  b2vect const & outer_boundaries (const int rl, const int c) const
  {
    return regions.AT(0).AT(rl).AT(c).outer_boundaries;
  }
  
  int processor (const int rl, const int c) const
  {
    return regions.AT(0).AT(rl).AT(c).processor;
  }

  int old_processor (const int rl, const int c) const
  {
    return oldregions.AT(0).AT(rl).AT(c).processor;
  }

  int mglevels () const
  {
    return (int)regions.size();
  }
  
  int reflevels () const
  {
    if (mglevels() == 0) return 0;
    return (int)regions.AT(0).size();
  }
  
  int components (const int rl) const
  {
    return (int)regions.AT(0).AT(rl).size();
  }
  
  bool is_local (const int rl, const int c) const
  {
    return processor(rl,c) == dist::rank();
  }
  
  int local_components (const int rl) const;
  
  void locate_position (rvect const & rpos,
                        int const ml,
                        int const minrl, int const maxrl,
                        int & rl, int & c, ivect & aligned_ipos) const;
  
  void locate_position (ivect const & ipos,
                        int const ml,
                        int const minrl, int const maxrl,
                        int & rl, int & c, ivect & aligned_ipos) const;
  
  // Time hierarchy management
  void add (th * t);
  void remove (th * t);
  
  // Data hierarchy management
  void add (dh * d);
  void remove (dh * d);
  
  // Output
  size_t memory () const;
  ostream & output (ostream & os) const;

private:
  void do_output_bboxes (ostream & os) const;
  void do_output_bases (ostream & os) const;

};



inline size_t memoryof (gh const & g)
{
  return g.memory ();
}

inline ostream & operator<< (ostream & os, gh const & h)
{
  return h.output(os);
}



#endif // GH_HH
