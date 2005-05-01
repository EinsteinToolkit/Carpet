#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "cctk.h"

#include "defs.hh"
#include "gh.hh"

#include "th.hh"

using namespace std;



// Constructors
th::th (gh& h_, const vector<int> & reffacts_, const CCTK_REAL basedelta)
  : h(h_), reffacts(reffacts_), delta(basedelta)
{
  assert (reffacts.size() >= 1);
  assert (reffacts.front() == 1);
  for (size_t n = 1; n < reffacts.size(); ++ n) {
    assert (reffacts.at(n) >= reffacts.at(n-1));
    assert (reffacts.at(n) % reffacts.at(n-1) == 0);
  }
  h.add(this);
}

// Destructors
th::~th ()
{
  h.remove(this);
}

// Modifiers
void th::recompose ()
{
  const int old_mglevels = times.size();
  times.resize(h.mglevels());
  deltas.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    const int old_reflevels = times.at(ml).size();
    times.at(ml).resize(h.reflevels());
    deltas.at(ml).resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      if (old_mglevels==0) {
        times.at(ml).at(rl) = 0;
      } else if (rl < old_reflevels) {
        // do nothing
      } else {
        times.at(ml).at(rl) = times.at(ml).at(rl-1);
      }
      if (ml==0) {
	deltas.at(ml).at(rl) = delta / reffacts.at(rl);
      } else {
	deltas.at(ml).at(rl) = deltas.at(ml-1).at(rl) * h.mgfact;
      }
    }
  }
}



// Output
void th::output (ostream& os) const
{
  os << "th:"
     << "times={";
  const char * sep = "";
  for (int ml=0; ml<h.mglevels(); ++ml) {
    for (int rl=0; rl<h.reflevels(); ++rl) {
      os << sep
         << ml << ":" << rl << ":"
	 << times.at(ml).at(rl) << "(" << deltas.at(ml).at(rl) << ")";
      sep = ",";
    }
  }
  os << "}";
}
