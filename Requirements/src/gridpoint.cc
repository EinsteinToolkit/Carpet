#include "gridpoint.hh"
#include "util.hh"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cstdlib>

#include <cstdlib>

namespace Requirements {

using namespace std;

bool gridpoint_t::there_was_an_error = false;
bool gridpoint_t::there_was_a_warning = false;

// Accessors
void gridpoint_t::set_interior(bool b, const location_t &loc) {
  const gridpoint_t oldgp = *this;
  i_interior = b;
  output_location(oldgp, loc);
}
void gridpoint_t::set_boundary(bool b, const location_t &loc) {
  const gridpoint_t oldgp = *this;
  i_boundary = b;
  output_location(oldgp, loc);
}
void gridpoint_t::set_ghostzones(bool b, const location_t &loc) {
  const gridpoint_t oldgp = *this;
  i_ghostzones = b;
  output_location(oldgp, loc);
}
void gridpoint_t::set_boundary_ghostzones(bool b, const location_t &loc) {
  const gridpoint_t oldgp = *this;
  i_boundary_ghostzones = b;
  output_location(oldgp, loc);
}

// Check that all the parts of the grid variables read by a function
// are valid.  This will be called before the function is executed.
void gridpoint_t::check_state(clause_t const &clause,
                              const location_t &loc) const {
  if (not i_interior) {
    if (clause.everywhere or clause.interior) {
      report_error(loc, "interior");
    }
  }
  if (not i_boundary) {
    if (clause.everywhere or clause.boundary) {
      report_error(loc, "boundary");
    }
  }
  if (not i_ghostzones) {
    if (clause.everywhere) {
      report_error(loc, "ghostzones");
    }
  }
  if (not i_boundary_ghostzones) {
    if (clause.everywhere or clause.boundary_ghostzones) {
      report_error(loc, "boundary-ghostzones");
    }
  }
}

void gridpoint_t::report_error(const location_t &loc,
                               char const *const where) const {
  char *const fullname = CCTK_FullName(loc.vi);
  ostringstream state;
  state << "Current state: " << *this << std::endl;
  if (loc.fd) {
    // The error is related to a scheduled function
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Schedule READS clause not satisfied: "
               "Function %s::%s in %s: "
               "Variable %s iteration=%d reflevel=%d map=%d timelevel=%d: "
               "%s not valid before %s. %s",
               loc.fd->thorn, loc.fd->routine, loc.fd->where, fullname, loc.it,
               loc.rl, loc.m, loc.tl, where, loc.info.c_str(),
               state.str().c_str());
  } else {
    // The error is not related to a scheduled function
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Schedule READS clause not satisfied: "
               "Variable %s iteration=%d reflevel=%d map=%d timelevel=%d: "
               "%s not valid before %s. %s",
               fullname, loc.it, loc.rl, loc.m, loc.tl, where, loc.info.c_str(),
               state.str().c_str());
  }
  free(fullname);
  there_was_an_error = true;
}

void gridpoint_t::report_warning(const location_t &loc,
                                 char const *const where) const {
  char *const fullname = CCTK_FullName(loc.vi);
  ostringstream state;
  state << "Current state: " << *this << std::endl;
  if (loc.fd) {
    // The error is related to a scheduled function
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Schedule WRITES clause is superfluous: "
               "Function %s::%s in %s: "
               "Variable %s iteration=%d reflevel=%d map=%d timelevel=%d: "
               "%s already valid before %s. %s",
               loc.fd->thorn, loc.fd->routine, loc.fd->where, fullname, loc.it,
               loc.rl, loc.m, loc.tl, where, loc.info.c_str(),
               state.str().c_str());
  } else {
    // The error is not related to a scheduled function
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Schedule WRITES clause already satisfied: "
               "Variable %s iteration=%d reflevel=%d map=%d timelevel=%d: "
               "%s already valid before %s. %s",
               fullname, loc.it, loc.rl, loc.m, loc.tl, where, loc.info.c_str(),
               state.str().c_str());
  }
  free(fullname);
  there_was_a_warning = true;
}

// Update this object to reflect the fact that some parts of some
// variables are now valid after a function has been called
void gridpoint_t::update_state(clause_t const &clause, const location_t &loc) {
  const gridpoint_t oldgp = *this;
  if (clause.everywhere or clause.interior) {
    i_interior = true;
  }
  if (clause.everywhere or clause.boundary) {
    i_boundary = true;
  }
  if (clause.everywhere) {
    i_ghostzones = true;
  }
  if (clause.everywhere or clause.boundary_ghostzones) {
    i_boundary_ghostzones = true;
  }
  output_location(oldgp, loc);
}

void gridpoint_t::output(ostream &os) const {
  os << "(";
  if (i_interior)
    os << "interior;";
  if (i_boundary)
    os << "boundary;";
  if (i_ghostzones)
    os << "ghostzones;";
  if (i_boundary_ghostzones)
    os << "boundary_ghostzones;";
  os << ")";
}

// Some readable and parsable debug output
void gridpoint_t::output_location(const gridpoint_t &oldgp,
                                  const location_t &loc) const {
  DECLARE_CCTK_PARAMETERS;
  if (not output_changes)
    return;

  const gridpoint_t difference = *this ^ oldgp;
  if (difference.empty())
    return;

  cout << loc << " (" << (difference.interior() ? "IN" : "in") << ":"
       << interior() << "," << (difference.boundary() ? "BO" : "bo") << ":"
       << boundary() << "," << (difference.ghostzones() ? "GH" : "gh") << ":"
       << ghostzones() << ","
       << (difference.boundary_ghostzones() ? "BG" : "bg") << ":"
       << boundary_ghostzones() << ")\n";
}
}
