/***************************************************************************
                dggh.hh  -  Dimension Generic Grid Hierarchy
		bounding boxes for each multigrid level of each
		component of each refinement level
                             -------------------
    begin                : Sun Jun 10 2001
    copyright            : (C) 2001 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dggh.hh,v 1.1 2001/06/12 14:56:58 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DGGH_HH
#define DGGH_HH

#include <assert.h>

#include <iostream>
#include <list>

#include "defs.hh"
#include "dist.hh"

using namespace std;



// Forward declaration
class dimgeneric_dh;
class dimgeneric_gh;
class th;

// Output
ostream& operator<< (ostream& os, const dimgeneric_gh& h);



// A refinement hierarchy, where higher levels are finer than the base
// level.  The extents do not include ghost zones.
class dimgeneric_gh {
  
public:				// should be readonly
  
  // Fields
  int reffact;			// refinement factor
  centering refcent;		// vertex or cell centered
  
  int mgfact;			// default multigrid factor
  centering mgcent;		// default (vertex or cell centered)
  
  list<th*> ths;		// list of all time hierarchies
  
public:
  
  // Constructors
  dimgeneric_gh (const int reffact, const centering refcent,
		 const int mgfact, const centering mgcent);
  
  // Destructors
  virtual ~dimgeneric_gh ();
  
  // Accessors
  virtual int reflevels () const = 0;
  virtual int components (const int rl) const = 0;
  virtual int mglevels (const int rl, const int c) const = 0;
  virtual int proc (const int rl, const int c) const = 0;
  virtual bool is_local (const int rl, const int c) const = 0;
  
  // Time hierarchy management
  void add (th* t);
  void remove (th* t);
  
  // Output
  virtual ostream& output (ostream& os) const = 0;
};



inline ostream& operator<< (ostream& os, const dimgeneric_gh& h) {
  return h.output(os);
}



#if defined(TMPL_IMPLICIT)
#  include "dggh.cc"
#endif

#endif // DGGH_HH
