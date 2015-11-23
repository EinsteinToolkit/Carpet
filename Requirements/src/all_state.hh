#ifndef ALL_STATE_HH
#define ALL_STATE_HH

#include "Requirements.hh"
#include "all_clauses.hh"
#include "gridpoint.hh"

#include <cctk.h>
#include <cctki_Schedule.h>

#include <iostream>
#include <vector>

using namespace std;

namespace Requirements {

// Keep track of which time levels contain good data; modify this
// while time level cycling; routines should specify how many time
// levels they require/provide

// The state (valid/invalid) of parts of the grid for all
// timelevels, maps, refinement levels and variables
class all_state_t {
  typedef vector<gridpoint_t> timelevels_t;
  typedef vector<timelevels_t> maps_t;
  typedef vector<maps_t> reflevels_t;
  typedef vector<reflevels_t> variables_t;
  variables_t vars;
  variables_t old_vars; // for regridding
public:
  void setup(int maps);
  void change_storage(vector<int> const &groups, vector<int> const &timelevels,
                      int reflevel);
  void regrid(int reflevels);
  void recompose(int iteration, int reflevel, valid::valid_t where);
  void regrid_free();
  void cycle(int reflevel);
  void before_routine(cFunctionData const *function_data,
                      all_clauses_t &all_clauses, int iteration, int reflevel,
                      int map, int timelevel, int timelevel_offset) const;
  void after_routine(cFunctionData const *function_data,
                     all_clauses_t &all_clauses, int iteration, int reflevel,
                     int map, int timelevel, int timelevel_offset);
  void sync(cFunctionData const *function_data, vector<int> const &groups,
            int iteration, int reflevel, int timelevel);
  void restrict1(vector<int> const &groups, int iteration, int reflevel);
  void invalidate(vector<int> const &vars, int reflevel, int map,
                  int timelevel);

  // Input/Output helpers
  void input(istream &is);
  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, const all_state_t &a) {
  a.output(os);
  return os;
}
}

#endif
