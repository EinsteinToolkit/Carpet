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
  typedef vector<ibbox> cexts;	// ... for each component
  typedef vector<cexts> rexts;	// ... for each refinement level
  typedef vector<rexts> mexts;	// ... for each multigrid level
  
  typedef vector<bbvect> cbnds;	// ... for each component
  typedef vector<cbnds> rbnds;	// ... for each refinement level
  
  typedef vector<int>    cprocs; // ... for each component
  typedef vector<cprocs> rprocs; // ... for each refinement level
  
public:				// should be readonly
  
  // Fields
  const vector<ivect> reffacts; // refinement factors
  const centering refcent;	// vertex or cell centered
 
  const int mgfact;		// default multigrid factor
  const centering mgcent;	// default (vertex or cell centered)
  
  const ibbox baseextent;
  
  
private:
  vector<vector<ibbox> > _bases; // [ml][rl]
  // TODO: invent structure for this
  mexts _extents;		// extents of all grids
  rbnds _outer_boundaries;	// boundary descriptions of all grids
  rprocs _processors;		// processor numbers of all grids

  list<th*> ths;		// list of all time hierarchies
  list<dh*> dhs;		// list of all data hierarchies
  
public:
  
  // Constructors
  gh (const vector<ivect> & reffacts, const centering refcent,
      const int mgfact, const centering mgcent,
      const ibbox baseextent);
  
  // Destructors
  virtual ~gh ();
  
  // Modifiers
  void recompose (const mexts& exts, const rbnds& outer_bounds,
		  const rprocs& procs,
                  const bool do_prolongate);

  const mexts & extents() const
  {
    return _extents;
  }

  const rbnds & outer_boundaries() const
  {
    return _outer_boundaries;
  }

  const rprocs & processors() const
  {
    return _processors;
  }

  const vector<vector<ibbox> > & bases() const
  {
    return _bases;
  }
  
  // Accessors
  int mglevels () const
  {
    return (int)_extents.size();
  }
  
  int reflevels () const
  {
    if (mglevels() == 0) return 0;
    return (int)_extents.at(0).size();
  }
  
  int components (const int rl) const
  {
    return (int)_extents.at(0).at(rl).size();
  }
  
  bbvect outer_boundary (const int rl, const int c) const
  {
    return _outer_boundaries.at(rl).at(c);
  }
  
  int proc (const int rl, const int c) const
  {
    return _processors.at(rl).at(c);
  }

  bool is_local (const int rl, const int c) const
  {
    return proc(rl,c) == dist::rank();
  }
  
  int local_components (const int rl) const;
  
  // Time hierarchy management
  void add (th* t);
  void remove (th* t);
  
  // Data hierarchy management
  void add (dh* d);
  void remove (dh* d);
  
  // Output
  virtual ostream& output (ostream& os) const;

private:
  void check_processor_number_consistency ();
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
