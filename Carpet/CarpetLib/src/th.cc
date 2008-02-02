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
    assert (reffacts.AT(n) >= reffacts.AT(n-1));
    assert (reffacts.AT(n) % reffacts.AT(n-1) == 0);
  }
  h.add(this);
}

// Destructors
th::~th ()
{
  h.remove(this);
}

// Modifiers
void th::regrid ()
{
  const int old_mglevels = times.size();
  times.resize(h.mglevels());
  deltas.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    const int old_reflevels = times.AT(ml).size();
    times.AT(ml).resize(h.reflevels());
    deltas.AT(ml).resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      if (old_mglevels==0) {
        times.AT(ml).AT(rl) = 0;
      } else if (rl < old_reflevels) {
        // do nothing
      } else {
        times.AT(ml).AT(rl) = times.AT(ml).AT(rl-1);
      }
      if (ml==0) {
	deltas.AT(ml).AT(rl) = delta / reffacts.AT(rl);
      } else {
	deltas.AT(ml).AT(rl) = deltas.AT(ml-1).AT(rl) * h.mgfact;
      }
    }
  }
}



// Memory usage
size_t
th::
memory ()
  const
{
  return
    memoryof (reffacts) +
    memoryof (delta) +
    memoryof (times) +
    memoryof (deltas);
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
	 << times.AT(ml).AT(rl) << "(" << deltas.AT(ml).AT(rl) << ")";
      sep = ",";
    }
  }
  os << "}";
}
