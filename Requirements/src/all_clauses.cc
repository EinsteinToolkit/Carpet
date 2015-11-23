#include "all_clauses.hh"
#include "util.hh"

#include <cctk.h>
#include <cctki_Schedule.h>

#include <iostream>
#include <map>
#include <vector>

using namespace std;

namespace Requirements {

clauses_t const &
all_clauses_t::get_clauses(cFunctionData const *const function_data) {
  clauses_map_t::const_iterator const iclauses =
      clauses_map.find(function_data);
  if (iclauses != clauses_map.end())
    return *iclauses->second;
  clauses_t *const clauses = new clauses_t;
  clauses->setup(function_data);
  pair<clauses_map_t::const_iterator, bool> const ret =
      clauses_map.insert(clauses_map_t::value_type(function_data, clauses));
  assert(ret.second);
  return *ret.first->second;
}

void all_clauses_t::remove_clauses(cFunctionData const *const function_data) {
  clauses_map_t::iterator const iclauses = clauses_map.find(function_data);
  if (iclauses != clauses_map.end()) {
    clauses_map.erase(iclauses);
  }
  return;
}

void all_clauses_t::output(ostream &os) const {
  os << "all_clauses: {" << std::endl;
  for (std::map<cFunctionData const *, clauses_t const *>::const_iterator ti =
           clauses_map.begin();
       ti != clauses_map.end(); ++ti) {
    if (ti != clauses_map.begin())
      os << ",";
    os << ti->first->thorn << "::" << ti->first->routine << " in "
       << ti->first->where << ": " << *ti->second << std::endl;
  }
  os << "}";
}
};
