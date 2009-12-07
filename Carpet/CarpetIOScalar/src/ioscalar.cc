#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Network.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "carpet.hh"
#include "CarpetTimers.hh"



// That's a hack
namespace Carpet {
  void UnsupportedVarType (const int vindex);
}



namespace CarpetIOScalar {

  using namespace std;
  using namespace Carpet;



  // Definition of local types
  struct info {
    string reduction;
    int handle;
  };



  // Begin a new line without flushing the output buffer
  char const * const eol = "\n";



 // Registered functions
  static void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cctkGH);
  static int OutputGH (const cGH* cctkGH);
  static int OutputVarAs (const cGH* cctkGH, const char* varname, const char* alias);
  static int OutputVarAs (const cGH* cctkGH, const char* varname, const char* alias, const char* out_reductions);
  static int TimeToOutput (const cGH* cctkGH, int vindex);
  static int TriggerOutput (const cGH* cctkGH, int vindex);

  // Internal functions
#if 0
  static void SetFlag (int index, const char* optstring, void* arg);
#endif
  static void CheckSteerableParameters (const cGH *const cctkGH);



  // Definition of static members
  vector<bool> do_truncate;
  vector<int> last_output;

  /* CarpetScalar GH extension structure */
  struct
  {
    /* list of variables to output */
    char *out_vars;

    /* stop on I/O parameter parsing errors ? */
    int stop_on_parse_errors;

    /* I/O request description list (for all variables) */
    ioRequest **requests;
  } IOparameters;



  extern "C" int
  CarpetIOScalarStartup ()
  {
    CCTK_RegisterBanner ("AMR scalar I/O provided by CarpetIOScalar");

    int GHExtension = CCTK_RegisterGHExtension("CarpetIOScalar");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);

    int IOMethod = CCTK_RegisterIOMethod ("CarpetIOScalar");
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);

    return 0;
  }



  extern "C" void
  CarpetIOScalarInit (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;

    *this_iteration = 0;
    *last_output_iteration = 0;
    *last_output_time = cctkGH->cctk_time;
  }



  void*
  SetupGH (tFleshConfig* const fc, int const convLevel, cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    const void *dummy;

    dummy = &fc;
    dummy = &convLevel;
    dummy = &cctkGH;
    dummy = &dummy;

    // Truncate all files if this is not a restart
    const int numvars = CCTK_NumVars ();
    do_truncate.resize (numvars, true);

    // No iterations have yet been output
    last_output.resize (numvars, -1);

    IOparameters.requests = (ioRequest **) calloc (numvars, sizeof(ioRequest*));
    IOparameters.out_vars = strdup ("");

    // initial I/O parameter check
    IOparameters.stop_on_parse_errors = strict_io_parameter_check;
    CheckSteerableParameters (cctkGH);
    IOparameters.stop_on_parse_errors = 0;

    // We register only once, ergo we get only one handle.  We store
    // that statically, so there is no need to pass anything to
    // Cactus.
    return NULL;
  }



  int
  OutputGH (const cGH * const cctkGH)
  {
    static Carpet::Timer timer ("CarpetIOScalar::OutputGH");
    timer.start();
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cctkGH, vindex)) {
	TriggerOutput(cctkGH, vindex);
      }
    }
    timer.stop();
    return 0;
  }



  int
  OutputVarAs (const cGH * const cctkGH,
               const char* const varname, const char* const alias)
  {
    DECLARE_CCTK_PARAMETERS;

    int const retval = OutputVarAs (cctkGH, varname, alias,
                                    outScalar_reductions);

    return retval;
  }

  int
  OutputVarAs (const cGH * const cctkGH,
               const char* const varname, const char* const alias,
               const char* out_reductions)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (is_level_mode() or
            (is_singlemap_mode() and maps == 1) or
            (is_local_mode() and maps == 1 and vhh.at(Carpet::map)->local_components(reflevel) == 1));
    BEGIN_LEVEL_MODE (cctkGH) {

    const int n = CCTK_VarIndex(varname);
    if (n<0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Variable \"%s\" does not exist", varname);
      return -1;
    }
    assert (n>=0 and n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 and group<(int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 and n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 and var<CCTK_NumVarsInGroupI(group));
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
    assert (num_tl>=1);

    // Check for storage
    if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot output variable \"%s\" because it has no storage",
		  varname);
      return 0;
    }

    assert (do_global_mode);

    const int vartype = CCTK_VarTypeI(n);
    assert (vartype >= 0);

    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cctkGH, "IO");
    assert (iogh);

    // Create the output directory
    const char* myoutdir = outScalar_dir;
    if (CCTK_EQUALS(myoutdir, "")) {
      myoutdir = out_dir;
    }
    int const iret = IOUtil_CreateDirectory (cctkGH, myoutdir, 0, 0);
    if (iret < 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not create output directory \"%s\"", myoutdir);
    } else if (CCTK_Equals (verbose, "full")) {
      static bool firsttime = true;
      if (firsttime and iret > 0) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Output directory \"%s\" exists already", myoutdir);
      } else if (not firsttime and iret == 0) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Created output directory \"%s\"", myoutdir);
      }
      firsttime = false;
    }

    // Find the set of desired reductions
    list<info> reductions;
    string const redlist (out_reductions);
    string::const_iterator p = redlist.begin();
    while (p!=redlist.end()) {
      while (p!=redlist.end() and isspace(*p)) ++p;
      if (p==redlist.end()) break;
      string::const_iterator const start = p;
      while (p!=redlist.end() and not isspace(*p)) ++p;
      string::const_iterator const end = p;
      string const reduction (start, end);
      int const handle = CCTK_ReductionHandle (reduction.c_str());
      if (handle < 0) {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Reduction operator \"%s\" does not exist (maybe there is no reduction thorn active?)",
                    reduction.c_str());
      } else {
        info i;
        i.reduction = reduction;
        i.handle = handle;
        reductions.push_back (i);
      }
    }

    // Output in global mode
    BEGIN_GLOBAL_MODE(cctkGH) {

      for (list<info>::const_iterator ireduction = reductions.begin();
           ireduction != reductions.end();
           ++ireduction)
      {
        string const reduction = ireduction->reduction;

        fstream file;
        BeginTimingIO (cctkGH);
        CCTK_REAL io_files = 0;
        CCTK_REAL io_bytes_begin = 0, io_bytes_end = 0;
        if (CCTK_MyProc(cctkGH)==0) {

          // Invent a file name
          ostringstream filenamebuf;
          filenamebuf << myoutdir << "/" << alias << "." << reduction
                      << ".asc";
          // we need a persistent temporary here
          string filenamestr = filenamebuf.str();
          const char* const filename = filenamestr.c_str();

          if (do_truncate.at(n) and IO_TruncateOutputFiles (cctkGH)) {
            file.open (filename, ios::out | ios::trunc);
          } else {
            file.open (filename, ios::out | ios::app);
          }
          if (not file.good()) {
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Could not open output file \"%s\" for variable \"%s\"",
                        filename, varname);
          }
          assert (file.is_open());

          io_files += 1;
          io_bytes_begin = file.tellg();

          // If this is the first time, then write a nice header
          if (do_truncate.at(n)) {
            bool want_labels = false;
            bool want_date = false;
            bool want_parfilename = false;
            bool want_other = false;
            if (CCTK_EQUALS (out_fileinfo, "none")) {
              // do nothing
            } else if (CCTK_EQUALS (out_fileinfo, "axis labels")) {
              want_labels = true;
            } else if (CCTK_EQUALS (out_fileinfo, "creation date")) {
              want_date = true;
            } else if (CCTK_EQUALS (out_fileinfo, "parameter filename")) {
              want_parfilename = true;
            } else if (CCTK_EQUALS (out_fileinfo, "all")) {
              want_labels = true;
              want_date = true;
              want_parfilename = true;
              want_other = true;
            } else {
              CCTK_WARN (0, "internal error");
            }
            file << "# Scalar ASCII output created by CarpetIOScalar" << eol;
            if (want_date) {
              char run_host [1000];
              Util_GetHostName (run_host, sizeof run_host);
              char const * const run_user = getenv ("USER");
              char run_date [1000];
              Util_CurrentDate (sizeof run_date, run_date);
              char run_time [1000];
              Util_CurrentTime (sizeof run_time, run_time);
              file << "# created on " << run_host
                   << " by " << run_user
                   << " on " << run_date
                   << " at " << run_time << eol;
            }
            if (want_parfilename) {
              char parameter_filename [10000];
              CCTK_ParameterFilename
                (sizeof parameter_filename, parameter_filename);
              file << "# parameter filename: \"" << parameter_filename << "\"" << eol;
            }
            if (want_other) {
              if (CCTK_IsFunctionAliased ("UniqueBuildID")) {
                char const * const build_id
                  = (char const *) UniqueBuildID (cctkGH);
                file << "# Build ID: " << build_id << eol;
              }
              if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
                char const * const job_id
                  = static_cast<char const *> (UniqueSimulationID (cctkGH));
                file << "# Simulation ID: " << job_id << eol;
              }
              if (CCTK_IsFunctionAliased ("UniqueRunID")) {
                char const * const job_id
                  = static_cast<char const *> (UniqueRunID (cctkGH));
                file << "# Run ID: " << job_id << eol;
              }
            }
            file << "#" << eol;
            if (want_labels) {
              file << "# " << varname << " (" << alias << ")" << eol;
              file << "# 1:iteration 2:time 3:data" << eol;
              int col = 3;
              if (one_file_per_group) {
                file << "# data columns:";
                int const firstvar = CCTK_FirstVarIndexI(group);
                int const numvars = CCTK_NumVarsInGroupI(group);
                for (int n=firstvar; n<firstvar+numvars; ++n) {
                  file << " " << col << ":" << CCTK_VarName(n);
                  col += CarpetSimpleMPIDatatypeLength (vartype);
                }
                file << eol;
              }
            }
          }

          file << setprecision(15);
          assert (file.good());

        } // if on the root processor
        
        if (CCTK_MyProc(cctkGH)==0) {
          file << cctk_iteration << " " << cctk_time;
        }
        
        int const handle = ireduction->handle;

        union {
#define TYPECASE(N,T) T var_##T;
#include "carpet_typecase.hh"
#undef TYPECASE
        } result;

        int const firstvar
          = one_file_per_group ? CCTK_FirstVarIndexI(group) : n;
        int const numvars
          = one_file_per_group ? CCTK_NumVarsInGroupI(group) : 1;
        
        for (int n=firstvar; n<firstvar+numvars; ++n) {
          
          int const ierr
            = CCTK_Reduce (cctkGH, 0, handle, 1, vartype, &result, 1, n);
          if (ierr) {
            char * const fullname = CCTK_FullName (n);
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Error during reduction for variable \"%s\"",
                        fullname);
            free (fullname);
            memset (&result, 0, sizeof result);
          }
          
          if (CCTK_MyProc(cctkGH)==0) {
            file << " ";
            
            switch (vartype) {
#define TYPECASE(N,T)                           \
              case N:                           \
                file << result.var_##T;         \
              break;
#include "carpet_typecase.hh"
#undef TYPECASE
            default:
              UnsupportedVarType (n);
            }
          }
          
        } // for n
        
        if (CCTK_MyProc(cctkGH)==0) {
          file << eol;
          assert (file.good());

          io_bytes_end = file.tellg();
          file.close();
          assert (file.good());
        }

        assert (not file.is_open());

        CCTK_REAL const io_bytes = io_bytes_end - io_bytes_begin;
        EndTimingIO (cctkGH, io_files, io_bytes, false);

      } // for reductions

    } END_GLOBAL_MODE;

    // Don't truncate again
    do_truncate.at(n) = false;

    } END_LEVEL_MODE;

    return 0;
  }



  int
  TimeToOutput (const cGH * const cctkGH, int const vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (vindex>=0 and vindex<CCTK_NumVars());

    if (not do_global_mode) return 0;

    CheckSteerableParameters (cctkGH);

    // check if output for this variable was requested
    if (not IOparameters.requests[vindex])
    {
      return (0);
    }

    // check whether to output at this iteration
    bool output_this_iteration;

    const char* myoutcriterion = outScalar_criterion;
    if (CCTK_EQUALS(myoutcriterion, "default")) {
      myoutcriterion = out_criterion;
    }

    if (CCTK_EQUALS (myoutcriterion, "never")) {

      // Never output
      output_this_iteration = false;

    } else if (CCTK_EQUALS (myoutcriterion, "iteration")) {

      int myoutevery = outScalar_every;
      if (myoutevery == -2) {
        myoutevery = out_every;
      }
      if (myoutevery <= 0) {
        // output is disabled
        output_this_iteration = false;
      } else if (cctk_iteration == *this_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >= *last_output_iteration + myoutevery) {
        // it is time for the next output
        output_this_iteration = true;
        *last_output_iteration = cctk_iteration;
        *this_iteration = cctk_iteration;
      } else {
        // we want no output at this iteration
        output_this_iteration = false;
      }

    } else if (CCTK_EQUALS (myoutcriterion, "divisor")) {

      int myoutevery = outScalar_every;
      if (myoutevery == -2) {
        myoutevery = out_every;
      }
      if (myoutevery <= 0) {
        // output is disabled
        output_this_iteration = false;
      } else if ((cctk_iteration % myoutevery) == 0 ) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else {
        // we want no output at this iteration
        output_this_iteration = false;
      }

    } else if (CCTK_EQUALS (myoutcriterion, "time")) {

      CCTK_REAL myoutdt = outScalar_dt;
      if (myoutdt == -2) {
        myoutdt = out_dt;
      }
      if (myoutdt < 0) {
        // output is disabled
        output_this_iteration = false;
      } else if (myoutdt == 0) {
        // output all iterations
        output_this_iteration = true;
      } else if (cctk_iteration == *this_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else {
        int do_output =
          (cctk_time / cctk_delta_time >=
           (*last_output_time + myoutdt) / cctk_delta_time - 1.0e-12);
        MPI_Bcast (&do_output, 1, MPI_INT, 0, dist::comm());
        if (do_output) {
          // it is time for the next output
          output_this_iteration = true;
          *last_output_time = cctk_time;
          *this_iteration = cctk_iteration;
        } else {
          // we want no output at this iteration
          output_this_iteration = false;
        }
      }

    } else {

      assert (0);

    } // select output criterion

    if (not output_this_iteration) return 0;



#if 0
    // check which variables to output
    static vector<bool> output_variables;
    static int output_variables_iteration = -1;

    if (cctk_iteration > output_variables_iteration) {
      output_variables.resize (CCTK_NumVars());

      const char* const varlist = outScalar_vars;
      if (CCTK_TraverseString (varlist, SetFlag, &output_variables,
                               CCTK_GROUP_OR_VAR) < 0)
      {
        CCTK_WARN (output_variables_iteration < 0 and strict_io_parameter_check ?
                   0 : 1,
                   "error while parsing parameter 'IOScalar::outScalar_vars'");
      }

      output_variables_iteration = cctk_iteration;
    }

    if (not output_variables.at(vindex)) return 0;
#endif


    if (last_output.at(vindex) == cctk_iteration) {
      // Has already been output during this iteration
      char* const varname = CCTK_FullName(vindex);
      CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Skipping output for variable \"%s\", because this variable "
		  "has already been output during the current iteration -- "
		  "probably via a trigger during the analysis stage",
		  varname);
      free (varname);
      return 0;
    }

    assert (last_output.at(vindex) < cctk_iteration);

    // Should be output during this iteration
    return 1;
  }



  int
  TriggerOutput (const cGH * const cctkGH, int const vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (vindex>=0 and vindex<CCTK_NumVars());
    
    // use individual reductions list for this variable when given
    // otherwise IOScalar::outScalar_reductions
    const char* out_reductions = IOparameters.requests[vindex]->reductions;
    if (not out_reductions) out_reductions = outScalar_reductions;

    int retval;
    
    if (one_file_per_group) {
      
      char* const fullname = CCTK_FullName(vindex);
      int const gindex = CCTK_GroupIndexFromVarI(vindex);
      char* const groupname = CCTK_GroupName(gindex);
      for (char* p=groupname; *p; ++p) *p=tolower(*p);
      retval = OutputVarAs (cctkGH, fullname, groupname, out_reductions);
      free (fullname);
      free (groupname);
      
      int const firstvar = CCTK_FirstVarIndexI(gindex);
      int const numvars = CCTK_NumVarsInGroupI(gindex);
      for (int n=firstvar; n<firstvar+numvars; ++n) {
        last_output.at(n) = cctk_iteration;
      }
      
    } else {
      
      char* const fullname = CCTK_FullName(vindex);
      char const* const varname = CCTK_VarName(vindex);
      retval = OutputVarAs (cctkGH, fullname, varname, out_reductions);
      free (fullname);
      
      last_output.at(vindex) = cctk_iteration;
      
    }
    
    return retval;
  }


  static void CheckSteerableParameters (const cGH *const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;

    // re-parse the 'IOScalar::outScalar_vars' parameter if it has changed
    if (strcmp (outScalar_vars, IOparameters.out_vars)) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
      IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                                 "IOScalar::outScalar_vars",
                                 IOparameters.stop_on_parse_errors,
                                 outScalar_vars, -1, -1.0,
                                 IOparameters.requests);
#else
      IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                                 "IOScalar::outScalar_vars",
                                 IOparameters.stop_on_parse_errors,
                                 outScalar_vars, -1,
                                 IOparameters.requests);
#endif

      // notify the user about the new setting
      if (not CCTK_Equals (verbose, "none")) {
        int count = 0;
        string msg ("Periodic scalar output requested for '");
        for (int i = CCTK_NumVars () - 1; i >= 0; i--) {
          if (IOparameters.requests[i]) {
            if (count++) {
              msg += "', '";
            }
            char *fullname = CCTK_FullName (i);
            msg += fullname;
            free (fullname);
          }
        }
        if (count) {
          msg += "'";
          CCTK_INFO (msg.c_str());
        }
      }

      // save the last setting of 'IOScalar::outScalar_vars' parameter
      free (IOparameters.out_vars);
      IOparameters.out_vars = strdup (outScalar_vars);
    }
  }



#if 0
  void
  SetFlag (int const index, const char * const optstring, void * const arg)
  {
    const void *dummy;

    dummy = &optstring;
    dummy = &dummy;

    vector<bool>& flags = *(vector<bool>*)arg;
    flags.at(index) = true;
  }
#endif



} // namespace CarpetIOScalar
