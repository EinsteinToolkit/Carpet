#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <hdf5.h>
#include <F5/F5F.h>
#include <F5/F5R.h>
#include <F5/F5iterate.h>
#include <F5/F5uniform.h>

#include <bbox.hh>
#include <defs.hh>
#include <vect.hh>

#include <carpet.hh>



// This requires defining a boolean local variable "error_flag",
// initialised to false
#define FAILWARN(_expr)                                                 \
  failwarn(error_flag, _expr, __LINE__, __FILE__, CCTK_THORNSTRING, #_expr)

template<typename T>
static
T failwarn (bool& error_flag, T const expr,
            int const line, char const *const file, char const *const thorn,
            char const *const msg)
{
  if (expr < 0) {
    CCTK_VWarn (CCTK_WARN_ALERT, line, file, thorn,
                "Expression \"%s\" return %d", msg, (int)expr);
    error_flag = true;
  }
  return expr;
}

  

namespace CarpetIOF5 {
  class indent_t;
};
ostream& operator<<(ostream& os, CarpetIOF5::indent_t const& indent);



namespace CarpetIOF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Indentation
  class indent_t {
    static int const width = 3;
    static int level;
  public:
    indent_t() { ++level; }
    ~indent_t() { --level; }
    ostream& output(ostream& os) const;
  };
  
  
  
  // File mode for creating directories
  int const mode = 0755;
  
  // A nan
  CCTK_REAL const nan = numeric_limits<CCTK_REAL>::quiet_NaN();
  
  // Special group and attribute names
  char const *const metadata_group = "Parameters and Global Attributes";
  char const *const all_parameters = "All Parameters";
  extern char const *const grid_structure;
  
  
  
  // Data types for HDF5 and F5 types
  typedef vect<hsize_t,dim> hvect;
  
  // Conversion operators for these datatypes
  static inline
  hvect v2h (ivect const& a)
  {
    return hvect(a);
  }
  static inline
  ivect h2v (hvect const& a)
  {
    return ivect(a);
  }
  static inline
  ivect h2v (hsize_t const* const a)
  {
    ivect res;
    for (int d=0; d<dim; ++d) {
      res[d] = a[d];
    }
    return res;
  }
  
  static inline
  F5_vec3_point_t v2p (rvect const& a)
  {
    F5_vec3_point_t res;
    res.x = a[0];
    res.y = a[1];
    res.z = a[2];
    return res;
  }
  
  static inline
  F5_vec3_float_t v2f (rvect const& a)
  {
    F5_vec3_float_t res;
    res.x = a[0];
    res.y = a[1];
    res.z = a[2];
    return res;
  }
  
  static inline
  F5_vec3_double_t v2d (rvect const& a)
  {
    F5_vec3_double_t res;
    res.x = a[0];
    res.y = a[1];
    res.z = a[2];
    return res;
  }
  
  
  
  // Handle HDF5 attributes more comfortably
  bool WriteAttribute (hid_t const group,
                       char const * const name,
                       int const ivalue);
  bool WriteAttribute (hid_t const group,
                       char const * const name,
                       double const dvalue);
  bool WriteAttribute (hid_t const group,
                       char const *const name,
                       char const *const svalue);
  bool WriteAttribute (hid_t const group,
                       char const *const name,
                       int const *const ivalues,
                       int const nvalues);
  bool WriteAttribute (hid_t const group,
                       char const *const name,
                       double const *const dvalues,
                       int const nvalues);
  bool WriteAttribute (hid_t const group,
                       char const *const name,
                       char const *const *const svalues,
                       int const nvalues);
  bool WriteLargeAttribute (hid_t const group,
                            char const *const name,
                            char const *const svalue);
  
  
  
  // Generate a good file name ("alias") for a variable
  string
  generate_basename (cGH const *const cctkGH,
                     int const variable);
  
  // Create the final file name on a particular processor
  string
  create_filename (cGH const *const cctkGH,
                   string const basename,
                   int const proc,
                   bool const create_directories);
  
  // Generate a good grid name (simulation name)
  string
  generate_gridname (cGH const *const cctkGH);
  
  // Generate a good tensor basis name (coordinate system name)
  string
  generate_chartname (cGH const *const cctkGH);
  
  // Generate a good fragment name (map name)
  string
  generate_fragmentname (cGH const *const cctkGH, int const m);
  int
  interpret_fragmentname (cGH const *const cctkGH,
                          char const *const fragmentname);
  
  // Generate a good topology name (map and refinement level name)
  string
  generate_topologyname (cGH const *const cctkGH,
                         int const m, ivect const& reffact);
  
  
  
  void
  write_metadata (cGH *const cctkGH, hid_t const file);
  
  // Handle Carpet's grid structure (this should move to Carpet and/or
  // CarpetLib)
  string
  serialise_grid_structure (cGH const *const cctkGH);
  
} // end namespace CarpetIOF5
