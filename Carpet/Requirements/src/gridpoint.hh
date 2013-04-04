#ifndef GRIDPOINT_HH
#define GRIDPOINT_HH

#include <iostream>

#include <cctk_Schedule.h>
#include <clause.hh>
#include <location.hh>

using namespace std;

namespace Requirements {
  
  extern bool there_was_an_error;
  extern bool there_was_a_warning;

  // Represents which have valid information and which do not.
  // This will later be indexed by rl, map etc.
  // Currently only works with unigrid.
  class gridpoint_t {
    bool i_interior, i_boundary, i_ghostzones, i_boundary_ghostzones;
    public:
    gridpoint_t():
      i_interior(false), i_boundary(false), i_ghostzones(false),
      i_boundary_ghostzones(false)
    {}

    // Construct an object with information about which points are
    // valid, assuming that a function with the given clause has just
    // been run
    gridpoint_t(clause_t const& clause):
      i_interior(clause.everywhere or clause.interior),
      i_boundary(clause.everywhere or clause.boundary),
      i_ghostzones(clause.everywhere),
      i_boundary_ghostzones(clause.everywhere or clause.boundary_ghostzones)
    {}
    // Accessors
    bool interior() const;
    bool boundary() const;
    bool ghostzones() const;
    bool boundary_ghostzones() const;
    void set_interior(bool b, location_t &l);
    void set_boundary(bool b, location_t &l);
    void set_ghostzones(bool b, location_t &l);
    void set_boundary_ghostzones(bool b, location_t &l);

    void check_state(clause_t const& clause,
                     cFunctionData const* function_data,
                     int vi, int rl, int m, int tl) const;
    void report_error(cFunctionData const* function_data,
                      int vi, int rl, int m, int tl,
                      char const* what, char const* where) const;
    void report_warning(cFunctionData const* function_data,
                        int vi, int rl, int m, int tl,
                        char const* what, char const* where) const;
    void update_state(clause_t const& clause, location_t &loc);

    // Input/Output helpers
    void input (istream& is);
    void output (ostream& os) const;
    void output_location (location_t &l, int changed) const;
  };


  inline ostream& operator<< (ostream& os, const gridpoint_t& a) {
    a.output(os);
    return os;
  }
}

#endif
