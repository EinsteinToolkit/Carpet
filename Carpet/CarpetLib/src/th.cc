// $Header:$

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
th<D>::th (gh<D>& h, const CCTK_REAL basedelta)
  : h(h), delta(basedelta) {
  h.add(this);
}

// Destructors
template<int D>
th<D>::~th () {
  h.remove(this);
}

// Modifiers
template<int D>
void th<D>::recompose () {
  times.resize(h.reflevels());
  deltas.resize(h.reflevels());
  for (int rl=0; rl<h.reflevels(); ++rl) {
    const int old_mglevels = times.at(rl).size();
    CCTK_REAL mgtime;
    // Select default time
    if (old_mglevels==0 && rl==0) {
      mgtime = 0;
    } else if (old_mglevels==0) {
      mgtime = times.at(rl-1).at(0);
    } else {
      mgtime = times.at(rl).at(old_mglevels-1);
    }
    times.at(rl).resize(h.mglevels(rl,0), mgtime);
    deltas.at(rl).resize(h.mglevels(rl,0));
    for (int ml=0; ml<h.mglevels(rl,0); ++ml) {
      if (rl==0 && ml==0) {
	deltas.at(rl).at(ml) = delta;
      } else if (ml==0) {
	deltas.at(rl).at(ml) = deltas.at(rl-1).at(ml) / h.reffact;
      } else {
	deltas.at(rl).at(ml) = deltas.at(rl).at(ml-1) * h.mgfact;
      }
    }
  }
}



// Output
template<int D>
void th<D>::output (ostream& os) const {
  os << "th<" << D << ">:"
     << "times={";
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int ml=0; ml<h.mglevels(rl,0); ++ml) {
      if (!(rl==0 && ml==0)) os << ",";
      os << rl << ":" << ml << ":"
	 << times.at(rl).at(ml) << "(" << deltas.at(rl).at(ml) << ")";
    }
  }
  os << "}";
}



template class th<3>;
