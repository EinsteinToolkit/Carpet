// $Header:$

#ifndef GH_HH
#define GH_HH

#include <assert.h>

#include <iostream>
#include <list>
#include <vector>

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
template<int D> class dh;
template<int D> class th;
template<int D> class gh;

// Output
template<int D>
ostream& operator<< (ostream& os, const gh<D>& h);



// A refinement hierarchy, where higher levels are finer than the base
// level.  The extents do not include ghost zones.
template<int D>
class gh {
  
public:
  
  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;
  
  typedef vect<vect<bool,2>,D> bvect;
  
  typedef vector<ibbox> mexts;	// ... for each multigrid level
  typedef vector<mexts> cexts;	// ... for each component
  typedef vector<cexts> rexts;	// ... for each refinement level
  
  typedef vector<bvect> cbnds;	// ... for each component
  typedef vector<cbnds> rbnds;	// ... for each refinement level
  
  typedef vector<int>    cprocs; // ... for each component
  typedef vector<cprocs> rprocs; // ... for each refinement level
  
public:				// should be readonly
  
  // Fields
  const int reffact;		// refinement factor
  const centering refcent;	// vertex or cell centered
 
  const int mgfact;		// default multigrid factor
  const centering mgcent;	// default (vertex or cell centered)
  
  const ibbox baseextent;
  
  
private:
  vector<vector<ibbox> > _bases; // [rl][ml]
  // TODO: invent structure for this
  rexts _extents;		// extents of all grids
  rbnds _outer_boundaries;	// boundary descriptions of all grids
  rprocs _processors;		// processor numbers of all grids

  list<th<D>*> ths;		// list of all time hierarchies
  list<dh<D>*> dhs;		// list of all data hierarchies
  
public:
  
  // Constructors
  gh (const int reffact, const centering refcent,
      const int mgfact, const centering mgcent,
      const ibbox baseextent);
  
  // Destructors
  virtual ~gh ();
  
  // Modifiers
  void recompose (const rexts& exts, const rbnds& outer_bounds,
		  const rprocs& procs,
                  const bool do_prolongate);

  const rexts & extents() const {
    return _extents;
  }

  const rbnds & outer_boundaries() const {
    return _outer_boundaries;
  }

  const rprocs & processors() const {
    return _processors;
  }

  const vector<vector<ibbox> > & bases() const {
    return _bases;
  }
  
  // Accessors
  int reflevels () const {
    return (int)_extents.size();
  }
  
  int components (const int rl) const {
    return (int)_extents.at(rl).size();
  }
  
  int mglevels (const int rl, const int c) const {
    return (int)_extents.at(rl).at(c).size();
  }
  
  bvect outer_boundary (const int rl, const int c) const {
    return _outer_boundaries.at(rl).at(c);
  }
  
  int proc (const int rl, const int c) const {
    return _processors.at(rl).at(c);
  }

  bool is_local (const int rl, const int c) const {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    return proc(rl,c) == rank;
  }
  
  int local_components (const int rl) const;
  
  // Time hierarchy management
  void add (th<D>* t);
  void remove (th<D>* t);
  
  // Data hierarchy management
  void add (dh<D>* d);
  void remove (dh<D>* d);
  
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



template<int D>
inline ostream& operator<< (ostream& os, const gh<D>& h) {
  h.output(os);
  return os;
}



#endif // GH_HH
