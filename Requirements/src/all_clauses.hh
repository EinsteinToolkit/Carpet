#ifndef ALL_CLAUSES_HH
#define ALL_CLAUSES_HH

#include "clauses.hh"

#include <cctk.h>
#include <cctki_Schedule.h>

#include <iostream>
#include <map>

namespace Requirements {

using namespace std;

class all_clauses_t {
  // TODO: Represent I/O as well?
  typedef std::map<cFunctionData const *, clauses_t const *> clauses_map_t;
  clauses_map_t clauses_map;
  // Singleton
  all_clauses_t(all_clauses_t const &);
  all_clauses_t &operator=(all_clauses_t const &);

public:
  all_clauses_t() {}
  clauses_t const &get_clauses(cFunctionData const *function_data);
  void remove_clauses(cFunctionData const *function_data);

  // Input/Output helpers
  void input(istream &is);
  void output(ostream &os) const;
};
}

#endif
