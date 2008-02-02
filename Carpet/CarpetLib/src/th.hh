#ifndef TH_HH
#define TH_HH

#include <cassert>
#include <iostream>
#include <vector>

#include "cctk.h"

#include "defs.hh"
#include "gh.hh"

using namespace std;



// Forward declaration
class th;

// Output
ostream& operator<< (ostream& os, const th& t);



// The time hierarchy (information about the current time)
class th {
  
public:				// should be readonly
  
  // Fields
  gh& h;                        // hierarchy
  
private:
  
  const vector<int> reffacts;

  CCTK_REAL delta;		// time step
  vector<vector<CCTK_REAL> > times; // current times [ml][rl]
  vector<vector<CCTK_REAL> > deltas; // time steps [ml][rl]
  
public:
  
  // Constructors
  th (gh& h, const vector<int> & reffacts, const CCTK_REAL basedelta);
  
  // Destructors
  ~th ();
  
  // Modifiers
  void regrid ();
  
  // Time management
  CCTK_REAL get_time (const int rl, const int ml) const
  {
    assert (rl>=0 and rl<h.reflevels());
    assert (ml>=0 and ml<h.mglevels());
    return times.AT(ml).AT(rl);
  }
  
  void set_time (const int rl, const int ml, const CCTK_REAL t)
  {
    assert (rl>=0 and rl<h.reflevels());
    assert (ml>=0 and ml<h.mglevels());
    times.AT(ml).AT(rl) = t;
  }
  
  void advance_time (const int rl, const int ml)
  {
    set_time(rl,ml, get_time(rl,ml) + get_delta(rl,ml));
  }
  
  CCTK_REAL get_delta (const int rl, const int ml) const
  {
    assert (rl>=0 and rl<h.reflevels());
    assert (ml>=0 and ml<h.mglevels());
    return deltas.AT(ml).AT(rl);
  }
  
  void set_delta (const int rl, const int ml, const CCTK_REAL dt)
  {
    assert (rl>=0 and rl<h.reflevels());
    assert (ml>=0 and ml<h.mglevels());
    deltas.AT(ml).AT(rl) = dt;
  }
  
  CCTK_REAL time (const int tl, const int rl, const int ml) const
  {
    assert (rl>=0 and rl<h.reflevels());
    assert (ml>=0 and ml<h.mglevels());
    return get_time(rl, ml) - tl * get_delta(rl, ml);
  }
  
  // Output
  size_t memory () const;
  void output (ostream& os) const;
};



inline size_t memoryof (th const & t)
{
  return t.memory ();
}
inline ostream& operator<< (ostream& os, const th& t) {
  t.output(os);
  return os;
}



#endif // TH_HH
