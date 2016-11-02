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
#include <vector>

#include <F5/F5F.h>
#include <F5/F5R.h>
#include <F5/F5iterate.h>
#include <F5/F5uniform.h>
#include <hdf5.h>

#include <bbox.hh>
#include <defs.hh>
#include <vect.hh>

#include <carpet.hh>

#include "distribute.hh"

// This requires defining a boolean local variable "error_flag",
// initialised to false
#define FAILWARN(_expr)                                                        \
  failwarn(error_flag, _expr, __LINE__, __FILE__, CCTK_THORNSTRING, #_expr)
#define FAILWARN0(_expr)                                                       \
  failwarn0(error_flag, _expr, __LINE__, __FILE__, CCTK_THORNSTRING, #_expr)

template <typename T>
static T failwarn(bool &error_flag, T const expr, int const line,
                  char const *const file, char const *const thorn,
                  char const *const msg) {
  static_assert(T(-1) < T(0), "Type T must be signed");
  if (expr < 0) {
    CCTK_VWarn(CCTK_WARN_ALERT, line, file, thorn,
               "Expression \"%s\" return %d", msg, (int)expr);
    error_flag = true;
  }
  return expr;
}

template <typename T>
static T failwarn0(bool &error_flag, T const expr, int const line,
                   char const *const file, char const *const thorn,
                   char const *const msg) {
  static_assert(T(-1) > T(0), "Type T must be unsigned");
  if (expr == 0) {
    CCTK_VWarn(CCTK_WARN_ALERT, line, file, thorn,
               "Expression \"%s\" return %d", msg, (int)expr);
    error_flag = true;
  }
  return expr;
}

namespace CarpetIOF5 {
class indent_t;
}
std::ostream &operator<<(std::ostream &os, CarpetIOF5::indent_t const &indent);

namespace CarpetIOF5 {

using namespace std;
using namespace Carpet;

// Indentation
class indent_t {
  static bool const debug = false;
  static int const width = 3;
  static int level;

private:
  indent_t(indent_t const &);
  indent_t &operator=(indent_t const &);

public:
  indent_t();
  ~indent_t();
  ostream &output(ostream &os) const;
};

// File mode for creating directories
int const mode = 0755;

// A nan
CCTK_REAL const nan = numeric_limits<CCTK_REAL>::quiet_NaN();

// Special group and attribute names
char const *const metadata_group = "Parameters and Global Attributes";
char const *const all_parameters = "All Parameters";
extern char const *const grid_structure;

// Tensor types
enum tensortype_t { tt_scalar, tt_vector, tt_symtensor, tt_tensor, tt_error };

// Data types for HDF5 and F5 types
typedef vect<hsize_t, dim> hvect;

// Conversion operators for these datatypes
static inline hvect v2h(ivect const &a) { return hvect(a); }
static inline ivect h2v(hvect const &a) { return ivect(a); }
static inline ivect h2v(hsize_t const *const a) {
  ivect res;
  for (int d = 0; d < dim; ++d) {
    res[d] = a[d];
  }
  return res;
}

template <typename T, int D>
static inline vect<T, D> filter(vect<T, D> const &a, vect<bool, D> const &c,
                                T const &fill = T(0)) {
  vect<T, D> r(fill);
  int ri = 0;
  for (int ai = 0; ai < dim; ++ai) {
    if (c[ai])
      r[ri++] = a[ai];
  }
  return r;
}

static inline F5_vec3_point_t v2p(rvect const &a) {
  F5_vec3_point_t res;
  res.x = a[0];
  res.y = a[1];
  res.z = a[2];
  return res;
}

static inline F5_vec3_float_t v2f(rvect const &a) {
  F5_vec3_float_t res;
  res.x = a[0];
  res.y = a[1];
  res.z = a[2];
  return res;
}

static inline F5_vec3_double_t v2d(rvect const &a) {
  F5_vec3_double_t res;
  res.x = a[0];
  res.y = a[1];
  res.z = a[2];
  return res;
}

// Handle HDF5 attributes more comfortably
bool WriteAttribute(hid_t const group, char const *const name,
                    int const ivalue);
bool WriteAttribute(hid_t const group, char const *const name,
                    double const dvalue);
bool WriteAttribute(hid_t const group, char const *const name,
                    char const *const svalue);
bool WriteAttribute(hid_t const group, char const *const name,
                    string const svalue);
bool WriteAttribute(hid_t const group, char const *const name,
                    int const *const ivalues, int const nvalues);
bool WriteAttribute(hid_t const group, char const *const name,
                    double const *const dvalues, int const nvalues);
bool WriteAttribute(hid_t const group, char const *const name,
                    char const *const *const svalues, int const nvalues);
bool WriteLargeAttribute(hid_t const group, char const *const name,
                         char const *const svalue);

bool ReadAttribute(hid_t const group, char const *const name, int &ivalue);
bool ReadAttribute(hid_t const group, char const *const name, double &dvalue);
bool ReadAttribute(hid_t const group, char const *const name, string &svalue);
bool ReadLargeAttribute(hid_t const group, char const *const name,
                        string &svalue);

// Generate a good file name ("alias") for a variable
string generate_basename(cGH const *const cctkGH, int const variable);

// Create the final file name on a particular processor
enum io_dir_t {
  io_dir_input,
  io_dir_output,
  io_dir_recover,
  io_dir_checkpoint
};
string create_filename(cGH const *const cctkGH, string const basename,
                       int const iteration, int const proc,
                       io_dir_t const io_dir, bool const create_directories);

// Generate a good grid name (simulation name)
string generate_gridname(cGH const *const cctkGH);

// Generate a good topology name (refinement level name)
string generate_topologyname(cGH const *const cctkGH, int const gi,
                             ivect const &reffact, ivect const &slice_ipos);

// Generate a good chart name (tensor basis name)
string generate_chartname(cGH const *const cctkGH);

// Generate a good fragment name (map and component name)
string generate_fragmentname(cGH const *const cctkGH, int const m, int const c);
void interpret_fragmentname(cGH const *const cctkGH,
                            char const *const fragmentname, int &m, int &c);

// Generate a good field name (group or variable name)
string generate_fieldname(cGH const *const cctkGH, int const vi,
                          tensortype_t const tt);
void interpret_fieldname(cGH const *const cctkGH, string fieldname, int &vi);

// Write/read Cactus metadata to a particular location in an HDF5
// file
void write_metadata(cGH const *const cctkGH, hid_t const group);
void read_metadata(cGH const *const cctkGH, hid_t const group);

// Handle Carpet's grid structure (this should move to Carpet and/or
// CarpetLib)
string serialise_grid_structure(cGH const *const cctkGH);
void deserialise_grid_structure(cGH const *const cctkGH, string const buf);

void output(cGH const *const cctkGH, hid_t const file,
            vector<bool> const &output_var, bool const output_past_timelevels,
            bool const output_metadata);

void input(cGH const *const cctkGH, hid_t const file,
           vector<bool> const &input_var, bool const input_past_timelevels,
           bool const input_metadata, scatter_t &scatter);

// Scheduled routines
extern "C" {
int CarpetIOF5_Startup();
int CarpetIOF5_RecoverParameters();
void CarpetIOF5_Init(CCTK_ARGUMENTS);
void CarpetIOF5_InitialDataCheckpoint(CCTK_ARGUMENTS);
void CarpetIOF5_EvolutionCheckpoint(CCTK_ARGUMENTS);
void CarpetIOF5_TerminationCheckpoint(CCTK_ARGUMENTS);
}

// Registered GH extension setup routine
void *SetupGH(tFleshConfig *const fleshconfig, int const convLevel,
              cGH *const cctkGH);

// Callbacks for CarpetIOF5's I/O method
int Input(cGH *const cctkGH, char const *const basefilename,
          int const called_from);
int OutputGH(cGH const *const cctkGH);
int TimeToOutput(cGH const *const cctkGH, int const vindex);
int TriggerOutput(cGH const *const cctkGH, int const vindex);
int OutputVarAs(cGH const *const cctkGH, const char *const varname,
                const char *const alias);

} // end namespace CarpetIOF5
