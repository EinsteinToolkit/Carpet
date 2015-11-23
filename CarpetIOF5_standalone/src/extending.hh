#ifndef EXTENDING_HH
#define EXTENDING_HH

#include <set>
#include <string>
#include <vector>

#include "cctk.h"

namespace CarpetIOF5 {

using std::set;
using std::string;
using std::vector;

class extending_t {

  static char const *const extension_name;

  struct extension_t {
    set<string> did_truncate;
    // [mglevel][reflevel][variable];
    vector<vector<vector<int> > > last_output_iteration;
  };

  extension_t *m_extension;

public:
  static int create(void *(*Setup)(tFleshConfig *config, int convlevel,
                                   cGH *cctkGH));

  static void *setup(cGH *cctkGH, int (*OutputGH)(cGH const *cctkGH),
                     int (*TimeToOutput)(cGH const *cctkGH, int variable),
                     int (*TriggerOutput)(cGH const *cctkGH, int variable),
                     int (*OutputVarAs)(cGH const *cctkGH, char const *varname,
                                        char const *alias));

  extending_t(cGH const *cctkGH);

  bool get_did_truncate(string name) const;

  void set_did_truncate(string name);

  int get_last_output_iteration(int ml, int rl, int vi) const;

  void set_last_output_iteration(int ml, int rl, int vi, int iteration);

private:
  void resize_last_output_iteration(int ml, int rl, int vi) const;
};

} // namespace CarpetIOF5

#endif // #ifndef EXDENDING_HH
