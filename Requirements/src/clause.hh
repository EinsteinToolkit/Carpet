#ifndef CLAUSE_HH
#define CLAUSE_HH

#include <cctk.h>
#include <cctk_Schedule.h>

#include <vector>

namespace Requirements {

using namespace std;

// Represent scheduled functions and their dependencies
// This reflects exactly what was written in the schedule.ccl file
struct clause_t {
  bool everywhere;          // all grid points (everywhere)
  bool interior;            // all interior points
  bool boundary;            // all boundary points, excluding
                            // ghostzones
  bool boundary_ghostzones; // all boundary ghost points
  bool timelevel0, timelevel1, timelevel2;
  bool all_timelevels; // all time levels
  bool all_maps;       // all maps (i.e. level mode)
  bool all_reflevels;  // all refinement levels (i.e. global mode)
  vector<int> vars;
  clause_t()
      : everywhere(false), interior(false), boundary(false),
        boundary_ghostzones(false), timelevel0(false), timelevel1(false),
        timelevel2(false), all_timelevels(false), all_maps(false),
        all_reflevels(false) {}
  void interpret_options(cFunctionData const *function_data);
  void parse_clause(char const *clause);
  int min_num_timelevels() const;
  bool active_on_timelevel(int tl) const;

  // Input/Output helpers
  void input(istream &is);
  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, const clause_t &a) {
  a.output(os);
  return os;
}
};

#endif
