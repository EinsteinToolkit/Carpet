/***************************************************************************
                          th.cc  -  Time Hierarchy
                          information about time levels
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/th.cc,v 1.8 2002/09/25 15:49:17 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <assert.h>
#include <math.h>

#include <iostream>

#include "cctk.h"

#include "defs.hh"
#include "dggh.hh"

#include "th.hh"

using namespace std;



// Constructors
th::th (dimgeneric_gh* h, const CCTK_REAL basedelta)
  : h(h), delta(basedelta) {
  h->add(this);
}

// Destructors
th::~th () {
  h->remove(this);
}

// Modifiers
void th::recompose () {
  times.resize(h->reflevels());
  deltas.resize(h->reflevels());
  for (int rl=0; rl<h->reflevels(); ++rl) {
    const int old_mglevels = times[rl].size();
    CCTK_REAL mgtime;
    // Select default time
    if (old_mglevels==0 && rl==0) {
      mgtime = 0;
    } else if (old_mglevels==0) {
      mgtime = times[rl-1][0];
    } else {
      mgtime = times[rl][old_mglevels-1];
    }
    times[rl].resize(h->mglevels(rl,0), mgtime);
    deltas[rl].resize(h->mglevels(rl,0));
    for (int ml=0; ml<h->mglevels(rl,0); ++ml) {
      if (rl==0 && ml==0) {
	deltas[rl][ml] = delta;
      } else if (ml==0) {
// 	assert (deltas[rl-1][ml] % h->reffact == 0);
	assert (fabs(fmod(deltas[rl-1][ml], h->reffact)) < 1e-10);
	deltas[rl][ml] = deltas[rl-1][ml] / h->reffact;
      } else {
	deltas[rl][ml] = deltas[rl][ml-1] * h->mgfact;
      }
    }
  }
}



// Output
void th::output (ostream& os) const {
  os << "th:"
     << "times={";
  for (int rl=0; rl<h->reflevels(); ++rl) {
    for (int ml=0; ml<h->mglevels(rl,0); ++ml) {
      if (!(rl==0 && ml==0)) os << ",";
      os << rl << ":" << ml << ":"
	 << times[rl][ml] << "(" << deltas[rl][ml] << ")";
    }
  }
  os << "}";
}
