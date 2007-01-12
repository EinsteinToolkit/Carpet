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

// Output
ostream& operator<< (ostream& os, const gh& h);



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
  
  const ibbox baseextent;
  
private:
  vector<vector<ibbox> > _bases; // [ml][rl]
  
  // Extents and properties of all grids
  mregs _regions;
  // A copy, used during regridding
  mregs _oldregions;
  
  list<th*> ths;		// list of all time hierarchies
  list<dh*> dhs;		// list of all data hierarchies
  
public:
  
  // Constructors
  gh (const vector<ivect> & reffacts, const centering refcent,
      const int mgfact, const centering mgcent,
      const ibbox baseextent);
  
  // Destructors
  ~gh ();
  
  // Modifiers
  void regrid (mregs const & regs);
  bool recompose (const int rl,
                  const bool do_prolongate);
private:
  bool level_did_change (const int rl) const;
  
  // Accessors
  
public:
  mregs const & regions () const
  {
    return _regions;
  }

  const vector<vector<ibbox> > & bases() const
  {
    return _bases;
  }
  
  ibbox extent (const int m, const int rl, const int c) const
  {
    return _regions.at(m).at(rl).at(c).extent;
  }
  
  b2vect outer_boundaries (const int rl, const int c) const
  {
    return _regions.at(0).at(rl).at(c).outer_boundaries;
  }

  b2vect refinement_boundaries (const int rl, const int c) const
  {
    return _regions.at(0).at(rl).at(c).refinement_boundaries;
  }

  int processor (const int rl, const int c) const
  {
    return _regions.at(0).at(rl).at(c).processor;
  }

  int mglevels () const
  {
    return (int)_regions.size();
  }
  
  int reflevels () const
  {
    if (mglevels() == 0) return 0;
    return (int)_regions.at(0).size();
  }
  
  int components (const int rl) const
  {
    return (int)_regions.at(0).at(rl).size();
  }

  bool is_local (const int rl, const int c) const
  {
    return processor(rl,c) == dist::rank();
  }
  
  int local_components (const int rl) const;
  
  // Time hierarchy management
  void add (th* t);
  void remove (th* t);
  
  // Data hierarchy management
  void add (dh* d);
  void remove (dh* d);
  
  // Output
  ostream& output (ostream& os) const;

private:
  void check_multigrid_consistency ();
  void check_component_consistency ();
  void check_base_grid_extent ();
  void check_refinement_levels ();
  void calculate_base_extents_of_all_levels ();
  void do_output_bboxes (ostream& os) const;
  void do_output_bases (ostream& os) const;

};



inline ostream& operator<< (ostream& os, const gh& h) {
  h.output(os);
  return os;
}



#endif // GH_HH
