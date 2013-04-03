#ifndef CLAUSES_HH
#define CLAUSES_HH

#include <iostream>

#include <cctk.h>
#include <cctki_Schedule.h>

#include <clause.hh>

using namespace std;

namespace Requirements {

  struct clauses_t {
    vector<clause_t> reads, writes;
    clauses_t() {}
    void setup(cFunctionData const* function_data);

    // Input/Output helpers
    void input (istream& is);
    void output (ostream& os) const;
  };

  inline ostream& operator<< (ostream& os, const clauses_t& a) {
    a.output(os);
    return os;
  }

};

#endif
