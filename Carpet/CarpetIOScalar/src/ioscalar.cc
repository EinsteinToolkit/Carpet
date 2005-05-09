#include <cassert>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_String.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "carpet.hh"



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


  // Registered functions
  static void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cctkGH);
  static int OutputGH (const cGH* cctkGH);
  static int OutputVarAs (const cGH* cctkGH, const char* varname, const char* alias);
  static int TimeToOutput (const cGH* cctkGH, int vindex);
  static int TriggerOutput (const cGH* cctkGH, int vindex);

  // Internal functions
  static void SetFlag (int index, const char* optstring, void* arg);
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


  // Special output routines for complex numbers

#ifdef CCTK_REAL4
  ostream& operator<< (ostream& os, const CCTK_COMPLEX8& val)
  {
    return os << CCTK_Cmplx8Real(val) << " " << CCTK_Cmplx8Imag(val);
  }
#endif

#ifdef CCTK_REAL8
  ostream& operator<< (ostream& os, const CCTK_COMPLEX16& val)
  {
    return os << CCTK_Cmplx16Real(val) << " " << CCTK_Cmplx16Imag(val);
  }
#endif

#ifdef CCTK_REAL16
  ostream& operator<< (ostream& os, const CCTK_COMPLEX32& val)
  {
    return os << CCTK_Cmplx32Real(val) << " " << CCTK_Cmplx32Imag(val);
  }
#endif



  extern "C" void
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
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cctkGH, vindex)) {
	TriggerOutput(cctkGH, vindex);
      }
    }
    return 0;
  }



  int
  OutputVarAs (const cGH * const cctkGH,
               const char* const varname, const char* const alias)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (is_level_mode());

    const int n = CCTK_VarIndex(varname);
    if (n<0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Variable \"%s\" does not exist", varname);
      return -1;
    }
    assert (n>=0 && n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVarsInGroupI(group));
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
    assert (num_tl>=1);

    // Check for storage
    if (! CCTK_QueryGroupStorageI(cctkGH, group)) {
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
    if (CCTK_MyProc(cctkGH)==0) {
      CCTK_CreateDirectory (0755, myoutdir);
    }

    // Find the set of desired reductions
    list<info> reductions;
    string const redlist (outScalar_reductions);
    string::const_iterator p = redlist.begin();
    while (p!=redlist.end()) {
      while (p!=redlist.end() && isspace(*p)) ++p;
      if (p==redlist.end()) break;
      string::const_iterator const start = p;
      while (p!=redlist.end() && !isspace(*p)) ++p;
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

        ofstream file;
        if (CCTK_MyProc(cctkGH)==0) {

          // Invent a file name
          ostringstream filenamebuf;
          filenamebuf << myoutdir << "/" << alias << "." << reduction
                      << ".asc";
          // we need a persistent temporary here
          string filenamestr = filenamebuf.str();
          const char* const filename = filenamestr.c_str();

          // If this is the first time, then write a nice header
          if (do_truncate.at(n) && IO_TruncateOutputFiles (cctkGH)) {
            file.open (filename, ios::out | ios::trunc);
            file << "# " << varname << " (" << alias << ")" << endl;
            file << "# iteration time data" << endl;
            if (one_file_per_group) {
              file << "# data columns:";
              int const firstvar = CCTK_FirstVarIndexI(group);
              int const numvars = CCTK_NumVarsInGroupI(group);
              for (int n=firstvar; n<firstvar+numvars; ++n) {
                file << " " << CCTK_VarName(n);
              }
              file << endl;
            }
          } else {
            file.open (filename, ios::out | ios::app);
          }
          if (! file.good()) {
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Could not open output file \"%s\" for variable \"%s\"",
                        filename, varname);
          }

          assert (file.is_open());

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
          assert (! ierr);
          
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
          file << endl;
          assert (file.good());
        }
        
        if (CCTK_MyProc(cctkGH)==0) {
          file.close();
          assert (file.good());
        }

        assert (! file.is_open());

      } // for reductions

    } END_GLOBAL_MODE;

    // Don't truncate again
    do_truncate.at(n) = false;

    return 0;
  }



  int
  TimeToOutput (const cGH * const cctkGH, int const vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (vindex>=0 && vindex<CCTK_NumVars());

    if (! do_global_mode) return 0;

    CheckSteerableParameters (cctkGH);

    // check if output for this variable was requested
    if (! IOparameters.requests[vindex])
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
      } else if (cctk_time / cctk_delta_time
                 >= (*last_output_time + myoutdt) / cctk_delta_time - 1.0e-12) {
        // it is time for the next output
        output_this_iteration = true;
        *last_output_time = cctk_time;
        *this_iteration = cctk_iteration;
      } else {
        // we want no output at this iteration
        output_this_iteration = false;
      }

    } else {

      assert (0);

    } // select output criterion

    if (! output_this_iteration) return 0;



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
        CCTK_WARN (output_variables_iteration < 0 && strict_io_parameter_check ?
                   0 : 1,
                   "error while parsing parameter 'IOScalar::outScalar_vars'");
      }

      output_variables_iteration = cctk_iteration;
    }

    if (! output_variables.at(vindex)) return 0;
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

    assert (vindex>=0 && vindex<CCTK_NumVars());
    
    int retval;
    
    if (one_file_per_group) {
      
      char* const fullname = CCTK_FullName(vindex);
      int const gindex = CCTK_GroupIndexFromVarI(vindex);
      char* const groupname = CCTK_GroupName(gindex);
      for (char* p=groupname; *p; ++p) *p=tolower(*p);
      retval = OutputVarAs (cctkGH, fullname, groupname);
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
      const int retval = OutputVarAs (cctkGH, fullname, varname);
      free (fullname);
      
      last_output.at(vindex) = cctk_iteration;
      
    }
    
    return retval;
  }


  static void CheckSteerableParameters (const cGH *const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;

    // re-parse the 'IOScalar::outScalar_vars' parameter if it has changed
    if (strcmp (outScalar_vars, IOparameters.out_vars))
    {
      IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                                 "IOScalar::outScalar_vars",
                                 IOparameters.stop_on_parse_errors,
                                 outScalar_vars, -1, IOparameters.requests);

      // notify the user about the new setting
      if (! CCTK_Equals (verbose, "none"))
      {
        char *msg = NULL;
        for (int i = CCTK_NumVars () - 1; i >= 0; i--)
        {
          if (IOparameters.requests[i])
          {
            char *fullname = CCTK_FullName (i);
            if (! msg)
            {
              Util_asprintf (&msg, "Periodic scalar output requested for '%s'",
                             fullname);
            }
            else
            {
              Util_asprintf (&msg, "%s, '%s'", msg, fullname);
            }
            free (fullname);
          }
        }
        if (msg)
        {
          CCTK_INFO (msg);
          free (msg);
        }
      }

      // save the last setting of 'IOScalar::outScalar_vars' parameter
      free (IOparameters.out_vars);
      IOparameters.out_vars = strdup (outScalar_vars);
    }
  }


  void
  SetFlag (int const index, const char * const optstring, void * const arg)
  {
    const void *dummy;

    dummy = &optstring;
    dummy = &dummy;

    vector<bool>& flags = *(vector<bool>*)arg;
    flags.at(index) = true;
  }



} // namespace CarpetIOScalar
