/***************************************************************************
                          th.hh  -  Time Hierarchy
                          information about time levels
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/th.hh,v 1.6 2002/06/07 17:17:19 schnetter Exp $

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

#include <assert.h>

#include <iostream>
#include <vector>

#include "defs.hh"
#include "dggh.hh"

using namespace std;



// Forward declaration
class th;

// Output
ostream& operator<< (ostream& os, const th& t);



// The time hierarchy (information about the current time)
class th {
  
public:				// should be readonly
  
  // Fields
  dimgeneric_gh *h;		// hierarchy
  
private:
  
  int delta;			// time step
  vector<vector<int> > times;	// current times
  vector<vector<int> > deltas;	// time steps
  
public:
  
  // Constructors
  th (dimgeneric_gh* h, const int basedelta);
  
  // Destructors
  ~th ();
  
  // Modifiers
  void recompose ();
  
  // Time management
  int get_time (const int rl, const int ml) const {
    assert (rl>=0 && rl<h->reflevels());
    assert (ml>=0 && ml<h->mglevels(rl,0));
    return times[rl][ml];
  }
  
  void set_time (const int rl, const int ml, const int t) {
    assert (rl>=0 && rl<h->reflevels());
    assert (ml>=0 && ml<h->mglevels(rl,0));
    times[rl][ml] = t;
  }
  
  void advance_time (const int rl, const int ml) {
    set_time(rl,ml, get_time(rl,ml) + get_delta(rl,ml));
  }
  
  int get_delta (const int rl, const int ml) const {
    assert (rl>=0 && rl<h->reflevels());
    assert (ml>=0 && ml<h->mglevels(rl,0));
    return deltas[rl][ml];
  }
  
  int set_delta (const int rl, const int ml, const int dt) {
    assert (rl>=0 && rl<h->reflevels());
    assert (ml>=0 && ml<h->mglevels(rl,0));
    deltas[rl][ml] = dt;
  }
  
  int time (const int tl, const int rl, const int ml) const {
    assert (rl>=0 && rl<h->reflevels());
    assert (ml>=0 && ml<h->mglevels(rl,0));
    return get_time(rl, ml) + tl * get_delta(rl, ml);
  }
  
  // Output
  void output (ostream& os) const;
};



inline ostream& operator<< (ostream& os, const th& t) {
  t.output(os);
  return os;
}



#endif // TH_HH
