/***************************************************************************
                          gh.hh  -  Grid Hierarchy
		bounding boxes for each multigrid level of each
		component of each refinement level
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gh.hh,v 1.10 2002/05/05 22:17:02 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GH_HH
#define GH_HH

#include <assert.h>

#include <iostream>
#include <list>
#include <vector>

#include "bbox.hh"
#include "defs.hh"
#include "dggh.hh"
#include "dist.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
template<int D> class dh;
template<int D> class gh;



// A refinement hierarchy, where higher levels are finer than the base
// level.  The extents do not include ghost zones.
template<int D>
class gh: public dimgeneric_gh {
  
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
  
  ibbox baseextent;		// bounds (inclusive) of base level
  vector<vector<ibbox> > bases; // [rl][ml]
  
  // TODO: invent structure for this
  rexts extents;		// extents of all grids
  rbnds outer_boundaries;	// boundary descriptions of all grids
  rprocs processors;		// processor numbers of all grids
  
  list<dh<D>*> dhs;		// list of all data hierarchies
  
public:
  
  // Constructors
  gh (const int reffact, const centering refcent,
      const int mgfact, const centering mgcent,
      const ibbox& baseextent);
  
  // Destructors
  virtual ~gh ();
  
  // Modifiers
  void recompose (const rexts& exts, const rbnds& outer_bounds,
		  const rprocs& procs);
  
  // Helpers
  cexts make_reflevel_multigrid_boxes (const vector<ibbox>& exts,
				       const int mglevels)
    const;
  rexts make_multigrid_boxes (const vector<vector<ibbox> >& exts,
			      const int mglevels)
    const;
  
  // Accessors
  int reflevels () const {
    return (int)extents.size();
  }
  
  int components (const int rl) const {
    assert (rl>=0 && rl<reflevels());
    return (int)extents[rl].size();
  }
  
  int mglevels (const int rl, const int c) const {
    assert (rl>=0 && rl<reflevels());
    assert (c>=0 && c<components(rl));
    return (int)extents[rl][c].size();
  }
  
  bvect outer_boundary (const int rl, const int c) const {
    assert (rl>=0 && rl<reflevels());
    assert (c>=0 && c<components(rl));
    return outer_boundaries[rl][c];
  }
  
  int proc (const int rl, const int c) const {
    assert (rl>=0 && rl<reflevels());
    assert (c>=0 && c<components(rl));
    return processors[rl][c];
  }

  bool is_local (const int rl, const int c) const {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    return proc(rl,c) == rank;
  }
  
  // Data hierarchy management
  void add (dh<D>* d);
  void remove (dh<D>* d);
  
  // Output
  virtual ostream& output (ostream& os) const;
};



#endif // GH_HH
