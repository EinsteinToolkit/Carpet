/***************************************************************************
                          dgdh.hh  -  Dimension Generic Data Hierarchy
			  A grid hierarchy plus ghost zones
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dgdh.hh,v 1.1 2001/06/12 14:56:58 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DGDH_HH
#define DGDH_HH

#include <assert.h>

#include <iostream>
#include <list>

#include "defs.hh"

using namespace std;



// Forward declaration
class dimgeneric_dh;
class dimgeneric_gf;

// Output
ostream& operator<< (ostream& os, const dimgeneric_dh& d);



// A data hierarchy (grid hierarchy plus ghost zones)
class dimgeneric_dh {
  
public:				// should be readonly
  
  // Fields
  int prolongation_order;	// order of spatial prolongation operator
  
public:
  
  // Constructors
  dimgeneric_dh (int prolongation_order);
  
  // Destructors
  virtual ~dimgeneric_dh ();
  
  // Helpers
  int prolongation_stencil_size () const;
  
  // Modifiers
  virtual void recompose () = 0;
  
  // Output
  virtual void output (ostream& os) const = 0;
};

inline ostream& operator<< (ostream& os, const dimgeneric_dh& d) {
  d.output(os);
  return os;
}



#if defined(TMPL_IMPLICIT)
#  include "dgdh.cc"
#endif

#endif // DGDH_HH
