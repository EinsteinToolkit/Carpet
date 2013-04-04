
#include <cctk.h>
#include <cctk_Parameters.h>

#include <Requirements.hh>

#include <gridpoint.hh>
#include <util.hh>

namespace Requirements {

  using namespace std;
  
  bool gridpoint_t::there_was_an_error = false;
  bool gridpoint_t::there_was_a_warning = false;

  // Accessors
  void gridpoint_t::set_interior(bool b, const location_t &l)
  {
    const gridpoint_t oldgp = *this;
    i_interior = b;
    output_location(oldgp, l);
  }
  void gridpoint_t::set_boundary(bool b, const location_t &l)
  {
    const gridpoint_t oldgp = *this;
    i_boundary = b;
    output_location(oldgp, l);
  }
  void gridpoint_t::set_ghostzones(bool b, const location_t &l)
  {
    const gridpoint_t oldgp = *this;
    i_ghostzones = b;
    output_location(oldgp, l);
  }
  void gridpoint_t::set_boundary_ghostzones(bool b, const location_t &l)
  {
    const gridpoint_t oldgp = *this;
    i_boundary_ghostzones = b;
    output_location(oldgp, l);
  }
  
  // Check that all the parts of the grid variables read by a function
  // are valid.  This will be called before the function is executed.
  void gridpoint_t::check_state(clause_t const& clause,
                                cFunctionData const* const function_data,
                                int const vi,
                                int const rl, int const m, int const tl)
    const
  {
    if (not i_interior) {
      if (clause.everywhere or clause.interior) {
        report_error(function_data, vi, rl, m, tl,
                     "calling function", "interior");
      }
    }
    if (not i_boundary) {
      if (clause.everywhere or clause.boundary) {
        report_error(function_data, vi, rl, m, tl,
                     "calling function", "boundary");
      }
    }
    if (not i_ghostzones) {
      if (clause.everywhere) {
        report_error(function_data, vi, rl, m, tl,
                     "calling function", "ghostzones");
      }
    }
    if (not i_boundary_ghostzones) {
      if (clause.everywhere or clause.boundary_ghostzones) {
        report_error(function_data, vi, rl, m, tl,
                     "calling", "boundary-ghostzones");
      }
    }
  }
  
  void gridpoint_t::report_error(cFunctionData const* const function_data,
                                 int const vi,
                                 int const rl, int const m, int const tl,
                                 char const* const what,
                                 char const* const where) const
  {
    char* const fullname = CCTK_FullName(vi);
    ostringstream state;
    state << "current state: " << *this << std::endl;
    if (function_data) {
      // The error is related to a scheduled function
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Schedule READS clause not satisfied: "
                 "Function %s::%s in %s: "
                 "Variable %s reflevel=%d map=%d timelevel=%d: "
                 "%s not valid for %s. %s",
                 function_data->thorn, function_data->routine,
                 function_data->where,
                 fullname, rl, m, tl,
                 where, what, state.str().c_str());
    } else {
      // The error is not related to a scheduled function
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Schedule READS clause not satisfied: "
                 "Variable %s reflevel=%d map=%d timelevel=%d: "
                 "%s not valid for %s. %s",
                 fullname, rl, m, tl,
                 where, what, state.str().c_str());
    }
    free(fullname);
    there_was_an_error = true;
  }
  
  void gridpoint_t::report_warning(cFunctionData const* const function_data,
                                   int const vi,
                                   int const rl, int const m, int const tl,
                                   char const* const what,
                                   char const* const where) const
  {
    char* const fullname = CCTK_FullName(vi);
    ostringstream state;
    state << "current state: " << *this << std::endl;
    if (function_data) {
      // The error is related to a scheduled function
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Schedule WRITES clause is superfluous: "
                 "Function %s::%s in %s: "
                 "Variable %s reflevel=%d map=%d timelevel=%d: "
                 "%s already valid for %s. %s",
                 function_data->thorn, function_data->routine,
                 function_data->where,
                 fullname, rl, m, tl,
                 where, what, state.str().c_str());
    } else {
      // The error is not related to a scheduled function
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Schedule WRITES clause already satisfied: "
                 "Variable %s reflevel=%d map=%d timelevel=%d: "
                 "%s already valid for %s. %s",
                 fullname, rl, m, tl,
                 where, what, state.str().c_str());
    }
    free(fullname);
    there_was_a_warning = true;
  }
  
  // Update this object to reflect the fact that some parts of some
  // variables are now valid after a function has been called
  void gridpoint_t::update_state(clause_t const& clause, const location_t &loc)
  {
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
  
  void gridpoint_t::output(ostream& os) const
  {
    os << "(";
    if (i_interior) os << "interior;";
    if (i_boundary) os << "boundary;";
    if (i_ghostzones) os << "ghostzones;";
    if (i_boundary_ghostzones) os << "boundary_ghostzones;";
    os << ")";
  }
  
  // Some readable and parsable debug output
  void gridpoint_t::output_location(const gridpoint_t& oldgp,
                                    const location_t& l) const
  {
    DECLARE_CCTK_PARAMETERS;
    if (not output_changes) return;
    
    const gridpoint_t difference = *this ^ oldgp;
    if (difference.empty()) return;
    
    cout << "LOC: " << l << " ("
         << (difference.interior()            ? "IN" : "in" ) << ":"
         << interior()            << ","
         << (difference.boundary()            ? "BO" : "bo" ) << ":"
         << boundary()            << ","
         << (difference.ghostzones()          ? "GH" : "gh" ) << ":"
         << ghostzones()          << ","
         << (difference.boundary_ghostzones() ? "BG" : "bg" ) << ":"
         << boundary_ghostzones() << ") "
         << l.info << "\n";
  }

}
