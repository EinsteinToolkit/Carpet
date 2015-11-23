#ifndef GRIDPOINT_HH
#define GRIDPOINT_HH

#include "clause.hh"
#include "location.hh"

#include <cctk_Schedule.h>

#include <iostream>

using namespace std;

namespace Requirements {

// Represents which have valid information and which do not.
// This will later be indexed by rl, map etc.
// Currently only works with unigrid.
class gridpoint_t {
  bool i_interior, i_boundary, i_ghostzones, i_boundary_ghostzones;

public:
  gridpoint_t()
      : i_interior(false), i_boundary(false), i_ghostzones(false),
        i_boundary_ghostzones(false) {}
  gridpoint_t(bool interior_, bool boundary_, bool ghostzones_,
              bool boundary_ghostzones_)
      : i_interior(interior_), i_boundary(boundary_), i_ghostzones(ghostzones_),
        i_boundary_ghostzones(boundary_ghostzones_) {}

  // Construct an object with information about which points are
  // valid, assuming that a function with the given clause has just
  // been run
  gridpoint_t(clause_t const &clause)
      : i_interior(clause.everywhere or clause.interior),
        i_boundary(clause.everywhere or clause.boundary),
        i_ghostzones(clause.everywhere),
        i_boundary_ghostzones(clause.everywhere or clause.boundary_ghostzones) {
  }

  // Accessors
  bool interior() const { return i_interior; }
  bool boundary() const { return i_boundary; }
  bool ghostzones() const { return i_ghostzones; }
  bool boundary_ghostzones() const { return i_boundary_ghostzones; }
  void set_interior(bool b, const location_t &loc);
  void set_boundary(bool b, const location_t &loc);
  void set_ghostzones(bool b, const location_t &loc);
  void set_boundary_ghostzones(bool b, const location_t &loc);

  void check_state(clause_t const &clause, const location_t &loc) const;
  void update_state(clause_t const &clause, const location_t &loc);

  void report_error(const location_t &loc, char const *where) const;
  void report_warning(const location_t &loc, char const *where) const;

  // Operators
  bool empty() const {
    return i_interior or i_boundary or i_ghostzones or i_boundary_ghostzones;
  }
  gridpoint_t operator^(const gridpoint_t &gp) const {
    return gridpoint_t(i_interior ^ gp.i_interior, i_boundary ^ gp.i_boundary,
                       i_ghostzones ^ gp.i_ghostzones,
                       i_boundary_ghostzones ^ gp.i_boundary_ghostzones);
  }

  // Input/Output helpers
  void input(istream &is);
  void output(ostream &os) const;
  void output_location(const gridpoint_t &oldgp, const location_t &loc) const;

  static bool there_was_an_error;
  static bool there_was_a_warning;
};

inline ostream &operator<<(ostream &os, const gridpoint_t &a) {
  a.output(os);
  return os;
}
}

#endif
