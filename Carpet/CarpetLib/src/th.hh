/***************************************************************************
                          th.hh  -  Time Hierarchy
                          information about time levels
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/th.hh,v 1.1 2001/03/01 13:40:10 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TH_HH
#define TH_HH

#include <cassert>
#include <iostream>
#include <vector>

#include "defs.hh"
#include "gh.hh"



// Forward declaration
template<int D> class th;

// Output
template<int D>
ostream& operator<< (ostream& os, const th<D>& t);



// The time hierarchy (information about the current time)
template<int D>
class th {
  
public:				// should be readonly
  
  // Fields
  gh<D> &h;			// hierarchy
  
private:
  
  int delta;			// time step
  vector<vector<int> > times;	// current times
  vector<vector<int> > deltas;	// time steps
  
public:
  
  // Constructors
  th (gh<D>& h, const int basedelta);
  
  // Destructors
  ~th ();
  
  // Modifiers
  void recompose ();
  
  // Time management
  int get_time (const int rl, const int ml) const {
    assert (rl>=0 && rl<h.reflevels());
    assert (ml>=0 && ml<h.mglevels(rl,0));
    return times[rl][ml];
  }
  
  void set_time (const int rl, const int ml, const int t) {
    assert (rl>=0 && rl<h.reflevels());
    assert (ml>=0 && ml<h.mglevels(rl,0));
    times[rl][ml] = t;
  }
  
  void advance_time (const int rl, const int ml) {
    set_time(rl,ml, get_time(rl,ml) + get_delta(rl,ml));
  }
  
  int get_delta (const int rl, const int ml) const {
    assert (rl>=0 && rl<h.reflevels());
    assert (ml>=0 && ml<h.mglevels(rl,0));
    return deltas[rl][ml];
  }
  
  int time (const int tl, const int rl, const int ml) const {
    assert (rl>=0 && rl<h.reflevels());
    assert (ml>=0 && ml<h.mglevels(rl,0));
    return get_time(rl, ml) + tl * get_delta(rl, ml);
  }
  
  // Output
  friend ostream& operator<< <> (ostream& os, const th& d);
};



#if defined(TMPL_IMPLICIT)
#  include "th.cc"
#endif

#endif // TH_HH
