/***************************************************************************
                          dh.hh  -  Data Hierarchy
			  A grid hierarchy plus ghost zones
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dh.hh,v 1.7 2001/04/23 08:10:15 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DH_HH
#define DH_HH

#include <assert.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "gh.hh"
#include "vect.hh"

using namespace std;



// Forward declaration
template<int D> class dh;
template<int D> class generic_gf;

// Output
template<int D>
ostream& operator<< (ostream& os, const dh<D>& d);



// A data hierarchy (grid hierarchy plus ghost zones)
template<int D>
class dh {
  
  // Types
  typedef vect<int,D>    ivect;
  typedef bbox<int,D>    ibbox;
  typedef bboxset<int,D> ibset;
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect; // vector of lists
  
public:
  struct dboxes {
    ibbox exterior;		// whole region (including boundaries)
    
    ibbox interior;		// interior (without boundaries)
    iblist send_mg_fine;
    iblist send_mg_coarse;
    iblist recv_mg_fine;
    iblist recv_mg_coarse;
    iblistvect send_ref_fine;
    iblistvect send_ref_coarse;
    iblistvect recv_ref_fine;
    iblistvect recv_ref_coarse;
    iblistvect send_sync;	// send while syncing
    iblistvect send_ref_bnd_fine;
    
    ibset boundaries;		// boundaries
    iblistvect recv_sync;	// received while syncing
    iblistvect recv_ref_bnd_coarse; // received from coarser grids
    ibset sync_not;		// not received while syncing (outer boundary of that level)
    ibset recv_not;		// not received while syncing or prolongating (globally outer boundary)
  };
private:
  
  struct dbases {
    ibbox exterior;		// whole region (including boundaries)
    ibbox interior;		// interior (without boundaries)
    ibset boundaries;		// boundaries
  };
  
  typedef vector<dboxes> mboxes; // ... for each multigrid level
  typedef vector<mboxes> cboxes; // ... for each component
  typedef vector<cboxes> rboxes; // ... for each refinement level
  
  typedef vector<dbases> mbases; // ... for each multigrid level
  typedef vector<mbases> rbases; // ... for each refinement level
  
public:				// should be readonly
  
  // Fields
  gh<D> &h;			// hierarchy
  ivect lghosts, ughosts;	// ghost zones
  int prolongation_order;	// order of spatial prolongation operator
  
  rboxes boxes;
  rbases bases;
  
  list<generic_gf<D>*> gfs;
  
public:
  
  // Constructors
  dh (gh<D>& h, const ivect& lghosts, const ivect& ughosts,
      int prolongation_order);
  
  // Destructors
  ~dh ();
  
  // Helpers
  int prolongation_stencil_size () const;
  
  // Modifiers
  void recompose ();
  
  // Grid function management
  void add (generic_gf<D>* f);
  void remove (generic_gf<D>* f);
  
  // Output
  void output (ostream& os) const;
};

template<int D>
inline ostream& operator<< (ostream& os, const dh<D>& d) {
  d.output(os);
  return os;
}



#if defined(TMPL_IMPLICIT)
#  include "dh.cc"
#endif

#endif // DH_HH
