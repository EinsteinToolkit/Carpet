// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/th.cc,v 1.12 2003/08/10 21:58:45 schnetter Exp $

#include <assert.h>
#include <math.h>

#include <iostream>

#include "cctk.h"

#include "defs.hh"
#include "gh.hh"

#include "th.hh"

using namespace std;



// Constructors
template<int D>
th<D>::th (gh<D>* h, const CCTK_REAL basedelta)
  : h(h), delta(basedelta) {
  h->add(this);
}

// Destructors
template<int D>
th<D>::~th () {
  h->remove(this);
}

// Modifiers
template<int D>
void th<D>::recompose () {
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
	deltas[rl][ml] = deltas[rl-1][ml] / h->reffact;
      } else {
	deltas[rl][ml] = deltas[rl][ml-1] * h->mgfact;
      }
    }
  }
}



// Output
template<int D>
void th<D>::output (ostream& os) const {
  os << "th<" << D << ">:"
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



template class th<3>;
