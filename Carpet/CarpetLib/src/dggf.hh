/***************************************************************************
                          dggf.hh  -  Dimension Generic Grid Function
                          grid function without type information
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dggf.hh,v 1.3 2001/12/09 16:43:09 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DGGF_HH
#define DGGF_HH

#include <assert.h>

#include <iostream>
#include <string>

#include "defs.hh"
#include "dgdata.hh"
#include "th.hh"

using namespace std;



// Forward declaration
class dimgeneric_gf;

// Output
ostream& operator<< (ostream& os, const dimgeneric_gf& f);



// A generic grid function without type information
class dimgeneric_gf {

public:				// should be readonly

  // Fields
  string name;

  th &t;			// time hierarchy
  int tmin, tmax;		// timelevels
  int prolongation_order_time;	// order of temporal prolongation operator

public:

  // Constructors
  dimgeneric_gf (const string name, th& t, 
		 const int tmin, const int tmax,
		 const int prolongation_order_time);
  
  // Destructors
  virtual ~dimgeneric_gf ();

  // Comparison
  bool operator== (const dimgeneric_gf& f) const;
  
  
  
  // Modifiers
  virtual void recompose () = 0;
  
  // Cycle the time levels by rotating the data sets
  virtual void cycle (int rl, int c, int ml) = 0;
  
  
  
  // TODO:
  
  // are these necessary in dimgeneric_gf, or is it sufficient to have
  // them in generic_gf<D>?
  
  // do i really want all dims, or do i just want to disregard the z
  // dim in carpet?
  
  // i likely also have to make a dimgeneric_data class.
  
  // The grid boundaries have to be updated after calling mg_restrict,
  // mg_prolongate, ref_restrict, or ref_prolongate.

  // "Updating" means here that the boundaries have to be
  // synchronised.  They don't need to be prolongated.

  // Copy a component from the next time level
  virtual void copy (int tl, int rl, int c, int ml) = 0;

  // Synchronise the boundaries of a component
  virtual void sync (int tl, int rl, int c, int ml) = 0;

  // Prolongate the boundaries of a component
  virtual void ref_bnd_prolongate (int tl, int rl, int c, int ml) = 0;
  
  // Restrict a multigrid level
  virtual void mg_restrict (int tl, int rl, int c, int ml) = 0;

  // Prolongate a multigrid level
  virtual void mg_prolongate (int tl, int rl, int c, int ml) = 0;

  // Restrict a refinement level
  virtual void ref_restrict (int tl, int rl, int c, int ml) = 0;

  // Prolongate a refinement level
  virtual void ref_prolongate (int tl, int rl, int c, int ml) = 0;
  
  
  
  // Access to the data
  virtual const dimgeneric_data* operator() (int tl, int rl, int c, int ml)
    const = 0;
  
  virtual dimgeneric_data* operator() (int tl, int rl, int c, int ml) = 0;
  
  
  
  // Output
  virtual ostream& output (ostream& os) const = 0;
};



inline ostream& operator<< (ostream& os, const dimgeneric_gf& f) {
  return f.output(os);
}



#if defined(TMPL_IMPLICIT)
#  include "dggf.cc"
#endif

#endif // DGGF_HH
