#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Version.h"
#include "util_Table.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "CarpetTimers.hh"

#include "CarpetIONirvana.hh"

#include "defs.hh"
#include "gh.hh"

namespace CarpetIONirvana {

using namespace std;
using namespace Carpet;
using namespace Nirvana;

// Variable definitions

// when was the last checkpoint written ?
static int last_checkpoint_iteration = -1;
static CCTK_REAL last_checkpoint_walltime;

// registered GH extension setup routine
static void *SetupGH(tFleshConfig *const fleshconfig, const int convLevel,
                     cGH *const cctkGH);

// callbacks for CarpetIOHDF5's I/O method
static int OutputGH(const cGH *const cctkGH);
static int OutputVarAs(const cGH *const cctkGH, const char *const varname,
                       const char *const alias);
static int TimeToOutput(const cGH *const cctkGH, const int vindex);
static int TriggerOutput(const cGH *const cctkGH, const int vindex);

// general checkpoint routine
static void Checkpoint(const cGH *const cctkGH, int called_from);

// callback for I/O parameter parsing routine
static void GetVarIndex(int vindex, const char *optstring, void *arg);

static void CheckSteerableParameters(const cGH *const cctkGH,
                                     CarpetIONirvanaGH *myGH);

//////////////////////////////////////////////////////////////////////////////
// public routines
//////////////////////////////////////////////////////////////////////////////

int CarpetIONirvana_Startup(void) {
  CCTK_RegisterBanner("Carpet Data Nirvana I/O");

  const int GHExtension = CCTK_RegisterGHExtension(CCTK_THORNSTRING);
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  /*IONirvana<0>::Startup();
  IONirvana<1>::Startup();
  IONirvana<2>::Startup();*/

  return (0);
}

// Called at basegrid during regular startup
void CarpetIONirvana_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = -1;
  *next_output_iteration = 0;
  *next_output_time = cctk_time;

  for (int d = 0; d < 3; ++d) {
    this_iteration_slice[d] = 0;
    last_output_iteration_slice[d] = 0;
    last_output_time_slice[d] = cctk_time;
  }
}

hid_t CCTKtoHDF5_Datatype(const cGH *const cctkGH, int cctk_type,
                          bool single_precision) {
  hid_t retval;

  const CarpetIONirvanaGH *myGH =
      (CarpetIONirvanaGH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);

  /*switch (cctk_type) {

    case CCTK_VARIABLE_CHAR:      retval = HDF5_CHAR; break;

    case CCTK_VARIABLE_INT:       retval = HDF5_INT; break;
#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:      retval = H5T_NATIVE_CHAR; break;
#endif
#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:      retval = H5T_NATIVE_SHORT; break;
#endif
#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:      retval = H5T_NATIVE_INT; break;
#endif
#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:      retval = H5T_NATIVE_LLONG; break;
#endif

    case CCTK_VARIABLE_REAL:      retval = HDF5_REAL;
#ifdef HAVE_CCTK_REAL4
                                  if (single_precision) {
                                    retval = H5T_NATIVE_FLOAT;
                                  }
#endif
                                  break;

    case CCTK_VARIABLE_COMPLEX:   retval = myGH->HDF5_COMPLEX;
#ifdef HAVE_CCTK_REAL4
                                  if (single_precision) {
                                    retval = myGH->HDF5_COMPLEX8;
                                  }
#endif
                                  break;

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:     retval = H5T_NATIVE_FLOAT; break;
    case CCTK_VARIABLE_COMPLEX8:  retval = myGH->HDF5_COMPLEX8; break;
#endif
#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:     retval = H5T_NATIVE_DOUBLE; break;
    case CCTK_VARIABLE_COMPLEX16: retval = myGH->HDF5_COMPLEX16; break;
#endif
#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:    retval = H5T_NATIVE_LDOUBLE; break;
    case CCTK_VARIABLE_COMPLEX32: retval = myGH->HDF5_COMPLEX32; break;
#endif

    default: CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Unsupported CCTK variable datatype %d", cctk_type);
             retval = -1;
  }*/

  return (retval);
}

//////////////////////////////////////////////////////////////////////////////
// private routines
//////////////////////////////////////////////////////////////////////////////
static void *SetupGH(tFleshConfig *const fleshconfig, const int convLevel,
                     cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;

  // register CarpetIONirvana's routines as a new I/O method
  const int IOMethod = CCTK_RegisterIOMethod("IONirvana");
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);

  if (not CCTK_Equals(verbose, "none")) {
    CCTK_INFO("I/O Method 'Nirvana': Sending data to Nirvana ");
  }

  const int numvars = CCTK_NumVars();

  // allocate a new GH extension structure
  CarpetIONirvanaGH *myGH = new CarpetIONirvanaGH;

  myGH->requests.resize(numvars);
  myGH->out_vars = strdup("");
  myGH->out_every_default = out_every - 1;

  // initial I/O parameter check
  myGH->out_dir = 0;
  myGH->stop_on_parse_errors = strict_io_parameter_check;
  CheckSteerableParameters(cctkGH, myGH);
  myGH->stop_on_parse_errors = 0;

  return (myGH);
}

static void CheckSteerableParameters(const cGH *const cctkGH,
                                     CarpetIONirvanaGH *myGH) {
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOHDF5::out_dir' parameter if it has changed
  const char *my_out_dir = *out_dir ? out_dir : io_out_dir;
  char *the_out_dir;
  if (strcmp(my_out_dir, ".")) {
    int i = strlen(my_out_dir);
    if (not strstr(my_out_dir, "%u")) {
      the_out_dir = (char *)malloc(i + 2);
      strcpy(the_out_dir, my_out_dir);
      the_out_dir[i] = '/';
      the_out_dir[i + 1] = 0;
    } else {
      // TODO: ensure that there is exactly one "%u" and no other "%"
      // substrings, except possibly "%%".
      the_out_dir = (char *)malloc(i + 20);
      snprintf(the_out_dir, i + 19, my_out_dir, dist::rank());
      strcat(the_out_dir, "/");
    }
  } else {
    the_out_dir = strdup("");
  }

  if (not myGH->out_dir or strcmp(the_out_dir, myGH->out_dir)) {
    free(myGH->out_dir);
    myGH->out_dir = the_out_dir;

    // create the output directory
    int result =
        IOUtil_CreateDirectory(cctkGH, myGH->out_dir, true, dist::rank());
    if (result < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Problem creating HDF5 output directory '%s'", myGH->out_dir);
    } else if (result > 0 and CCTK_Equals(verbose, "full")) {
      CCTK_VInfo(CCTK_THORNSTRING, "HDF5 output directory '%s' already exists",
                 myGH->out_dir);
    }
  } else {
    free(the_out_dir);
  }

  // re-parse the 'IONirvana::out_vars' parameter if it has changed
  if (strcmp(out_vars, myGH->out_vars)) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING, "IONirvana::out_vars",
                              myGH->stop_on_parse_errors, out_vars, -1, -1.0,
                              &myGH->requests[0]);
#else
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING, "IONirvana::out_vars",
                              myGH->stop_on_parse_errors, out_vars, -1,
                              &myGH->requests[0]);
#endif

    // notify the user about the new setting
    if (not CCTK_Equals(verbose, "none")) {
      int count = 0;
      ostringstream msg;
      msg << "Periodic scalar output requested for:";
      for (int vi = 0; vi < CCTK_NumVars(); ++vi) {
        if (myGH->requests[vi]) {
          ++count;
          char *const fullname = CCTK_FullName(vi);
          msg << eol << "   " << fullname;
          free(fullname);
        }
      }
      if (count > 0) {
        CCTK_INFO(msg.str().c_str());
      }
    }

    // save the last setting of 'IOHDF5::out_vars' parameter
    free(myGH->out_vars);
    myGH->out_vars = strdup(out_vars);
  }
}

static int OutputGH(const cGH *const cctkGH) {
  static Carpet::Timer timer("CarpetIONirvana::OutputGH");
  timer.start();
  for (int vindex = CCTK_NumVars() - 1; vindex >= 0; vindex--) {
    if (TimeToOutput(cctkGH, vindex)) {
      TriggerOutput(cctkGH, vindex);
    }
  }
  timer.stop();

  return (0);
}

static int TimeToOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int numvars = CCTK_NumVars();
  assert(vindex >= 0 and vindex < numvars);

  if (CCTK_GroupTypeFromVarI(vindex) != CCTK_GF and not do_global_mode) {
    return 0;
  }

  CarpetIONirvanaGH *myGH =
      (CarpetIONirvanaGH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);
  CheckSteerableParameters(cctkGH, myGH);

  // check if output for this variable was requested
  if (not myGH->requests[vindex]) {
    return (0);
  }

  // check whether this refinement level should be output
  if (not(myGH->requests[vindex]->refinement_levels & (1 << reflevel))) {
    return (0);
  }

  // check if output for this variable was requested individually
  // by a "<varname>{ out_every = <number> }" option string
  // this will overwrite the output criterion setting
  const char *myoutcriterion =
      CCTK_EQUALS(out_criterion, "default") ? io_out_criterion : out_criterion;
  if (myGH->requests[vindex]->out_every >= 0) {
    myoutcriterion = "divisor";
  }

  if (CCTK_EQUALS(myoutcriterion, "never")) {
    return (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration = false;

  if (CCTK_EQUALS(myoutcriterion, "iteration")) {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myoutevery > 0) {
      if (*this_iteration == cctk_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >= *next_output_iteration) {
        // it is time for the next output
        output_this_iteration = true;
        *this_iteration = cctk_iteration;
        *next_output_iteration = cctk_iteration + myoutevery;
      }
    }
  } else if (CCTK_EQUALS(myoutcriterion, "divisor")) {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myGH->requests[vindex]->out_every >= 0) {
      myoutevery = myGH->requests[vindex]->out_every;
    }
    if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
  } else if (CCTK_EQUALS(myoutcriterion, "time")) {
    CCTK_REAL myoutdt = out_dt == -2 ? io_out_dt : out_dt;
    if (myoutdt == 0 or *this_iteration == cctk_iteration) {
      output_this_iteration = true;
    } else if (myoutdt > 0) {
      int do_output = (cctk_time / cctk_delta_time >=
                       *next_output_time / cctk_delta_time - 1.0e-12);
      MPI_Bcast(&do_output, 1, MPI_INT, 0, dist::comm());
      if (do_output) {
        // it is time for the next output
        output_this_iteration = true;
        *this_iteration = cctk_iteration;
        *next_output_time = cctk_time + myoutdt;
      }
    }
  }

  return output_this_iteration ? 1 : 0;
}

static int TriggerOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_PARAMETERS;
  int retval;

  char *const fullname = CCTK_FullName(vindex);

  const int gindex = CCTK_GroupIndexFromVarI(vindex);
  char *const groupname = CCTK_GroupName(gindex);
  for (char *p = groupname; *p; ++p)
    *p = tolower(*p);
  retval = OutputVarAs(cctkGH, fullname, groupname);
  free(groupname);

  free(fullname);

  return (retval);
}

static void GetVarIndex(int vindex, const char *optstring, void *arg) {
  if (optstring) {
    char *fullname = CCTK_FullName(vindex);
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Option string '%s' will be ignored for Nirvana output of "
               "variable '%s'",
               optstring, fullname);
    free(fullname);
  }

  *((int *)arg) = vindex;
}

static int OutputVarAs(const cGH *const cctkGH, const char *const fullname,
                       const char *const alias) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;
  int vindex = -1;

  if (CCTK_TraverseString(fullname, GetVarIndex, &vindex, CCTK_VAR) < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "error while parsing variable name '%s' (alias name '%s')",
               fullname, alias);
    return (-1);
  }

  if (vindex < 0) {
    return (-1);
  }

  const int group = CCTK_GroupIndexFromVarI(vindex);
  assert(group >= 0);
  cGroup groupdata;
  CCTK_GroupData(group, &groupdata);
  if (groupdata.grouptype == CCTK_SCALAR or groupdata.grouptype == CCTK_ARRAY) {
    assert(do_global_mode);
  }

  // get the default I/O request for this variable
  const CarpetIONirvanaGH *myGH =
      (CarpetIONirvanaGH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);
  ioRequest *request = myGH->requests[vindex];
  if (not request) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    request = IOUtil_DefaultIORequest(cctkGH, vindex, 1, -1.0);
#else
    request = IOUtil_DefaultIORequest(cctkGH, vindex, 1);
#endif
  }

  // Get grid hierarchy extentsion from IOUtil
  const ioGH *const iogh = (const ioGH *)CCTK_GHExtension(cctkGH, "IO");
  assert(iogh);

  // Invent a base file name
  int ioproc = 0, nfiles = 1;
  string filename;
  filename.append(myGH->out_dir);
  // filename.append (alias);
  filename.append("cactus-simulation-output");
  if (not(groupdata.disttype == CCTK_DISTRIB_CONSTANT or dist::size() == 1)) {
    ioproc = dist::rank();
    nfiles = dist::size();
  }

  // check if the file has been created already
  typedef std::map<string, vector<vector<vector<int> > > > filelist;
  static filelist created_files;
  filelist::iterator thisfile = created_files.find(filename);
  bool is_new_file = thisfile == created_files.end();
  if (is_new_file) {
    int const numvars = CCTK_NumVars();
    vector<vector<vector<int> > > last_outputs; // [ml][rl][var]
    last_outputs.resize(mglevels);
    for (int ml = 0; ml < mglevels; ++ml) {
      last_outputs[ml].resize(maxreflevels);
      for (int rl = 0; rl < maxreflevels; ++rl) {
        last_outputs[ml][rl].resize(numvars, cctk_iteration - 1);
      }
    }
    thisfile = created_files.insert(
        thisfile, filelist::value_type(filename, last_outputs));
    assert(thisfile != created_files.end());
  }

  const int firstvar = CCTK_FirstVarIndexI(group);
  const int numvars = CCTK_NumVarsInGroupI(group);

  // check if this variable has been output already during this iteration
  int &last_output = thisfile->second.at(mglevel).at(reflevel).at(vindex);
  if (last_output == cctk_iteration) {
    // Has already been output during this iteration
    if (vindex == firstvar + numvars - 1) {
      char *varname = CCTK_FullName(vindex);
      CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Skipping output for variable \"%s\", because this variable "
                 "has already been output during the current iteration -- "
                 "probably via a trigger during the analysis stage",
                 varname);
      free(varname);
    }
    return (0);
  }
  assert(last_output < cctk_iteration);
  last_output = cctk_iteration;

  // Check for storage
  if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Cannot output variable '%s' because it has no storage",
               fullname);
    return (0);
  }

  // Open the output file if this is a designated I/O processor
  CCTK_REAL io_files = 0;
  CCTK_REAL io_bytes = 0;
  BeginTimingIO(cctkGH);

  if (CCTK_Equals(verbose, "full")) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Writing variable '%s' on mglevel %d reflevel %d", fullname,
               mglevel, reflevel);
  }
  for (int var = firstvar; var < firstvar + numvars; var++) {
    ioRequest *r = myGH->requests[var];
    if (not r) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
      r = IOUtil_DefaultIORequest(cctkGH, var, 1, -1.0);
#else
      r = IOUtil_DefaultIORequest(cctkGH, var, 1);
#endif
    }
    if (groupdata.disttype == CCTK_DISTRIB_CONSTANT) {
      error_count += WriteVar(cctkGH, filename, ioproc, io_bytes, r);
    } else
      error_count += WriteVar(cctkGH, filename, ioproc, io_bytes, r);

    if (r != myGH->requests[var])
      IOUtil_FreeIORequest(&r);

    // mark this variable to have been output at this iteration
    thisfile->second.at(mglevel).at(reflevel).at(var) = cctk_iteration;
  }

  // free I/O request structure
  if (request != myGH->requests[vindex]) {
    IOUtil_FreeIORequest(&request);
  }

  // Close the file

  EndTimingIO(cctkGH, io_files, io_bytes, true);

  if (error_count > 0 and abort_on_io_errors) {
    CCTK_WARN(0, "Aborting simulation due to previous I/O errors");
  }

  return (0);
}

} // namespace CarpetIONirvana
