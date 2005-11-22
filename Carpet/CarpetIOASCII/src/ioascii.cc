#include <cassert>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_Network.h"
#include "util_Table.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "carpet.hh"

#include "ioascii.hh"



// That's a hack
namespace Carpet {
  void UnsupportedVarType (const int vindex);
}



namespace CarpetIOASCII {

  using namespace std;
  using namespace Carpet;

  static void GetVarIndex (int vindex, const char* optstring, void* arg);



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



  template<int outdim>
  void WriteASCII (ostream& os,
		   vector<const gdata*> const gfdatas,
		   const bbox<int,dim>& gfext,
		   const int vi,
		   const int time,
		   const vect<int,dim>& org,
		   const vect<int,outdim>& dirs,
		   const int rl,
		   const int ml,
		   const int m,
		   const int c,
		   const int tl,
		   const CCTK_REAL coord_time,
		   const vect<CCTK_REAL,dim>& coord_lower,
		   const vect<CCTK_REAL,dim>& coord_upper);



  int CarpetIOASCIIStartup()
  {
    IOASCII<0>::Startup();
    IOASCII<1>::Startup();
    IOASCII<2>::Startup();
    IOASCII<3>::Startup();
    return 0;
  }



  void CarpetIOASCIIInit (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;

    for (int d=0; d<4; ++d) {
      this_iteration[d] = 0;
      last_output_iteration[d] = 0;
      last_output_time[d] = cctk_time;
    }
  }



  // Definition of static members
  template<int outdim> const char*        IOASCII<outdim>::my_out_dir;
  template<int outdim> char*              IOASCII<outdim>::my_out_vars;
  template<int outdim> vector<ioRequest*> IOASCII<outdim>::requests;

  static bool stop_on_parse_errors = false;


  template<int outdim>
  int IOASCII<outdim>::Startup()
  {
    ostringstream msg;
    msg << "AMR " << outdim << "D ASCII I/O provided by CarpetIOASCII";
    CCTK_RegisterBanner (msg.str().c_str());

    ostringstream extension_name;
    extension_name << "CarpetIOASCII_" << outdim << "D";
    const int GHExtension = CCTK_RegisterGHExtension(extension_name.str().c_str());
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);

    ostringstream method_name;
    method_name << "IOASCII_" << outdim << "D";
    const int IOMethod = CCTK_RegisterIOMethod (method_name.str().c_str());
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);

    return 0;
  }



  template<int outdim>
  void* IOASCII<outdim>
  ::SetupGH (tFleshConfig* const fc, const int convLevel, cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    const void *dummy;

    dummy = &fc;
    dummy = &convLevel;
    dummy = &cctkGH;
    dummy = &dummy;

    if (not CCTK_Equals (verbose, "none")) {
      CCTK_VInfo (CCTK_THORNSTRING, "I/O Method 'IOASCII_%dD' registered: "
                  "%dD AMR output of grid variables to ASCII files",
                  outdim, outdim);
    }

    // Truncate all files if this is not a restart
    const int numvars = CCTK_NumVars ();
    requests.resize (numvars);

    // create the output directory
    my_out_dir = GetStringParameter("out%dD_dir");
    if (CCTK_EQUALS (my_out_dir, "")) {
      my_out_dir = out_dir;
    }

    const ioGH* const ioUtilGH = (const ioGH*) CCTK_GHExtension (cctkGH, "IO");
    int result = IOUtil_CreateDirectory (cctkGH, my_out_dir, 0, 0);
    if (result < 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Problem creating %dD-output directory '%s'",
                  outdim, my_out_dir);
    } else if (result > 0 and CCTK_Equals (verbose, "full")) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "%dD-output directory '%s' already exists",
                  outdim, my_out_dir);
    }

    // initial I/O parameter check
    my_out_vars = strdup ("");
    stop_on_parse_errors = strict_io_parameter_check != 0;
    CheckSteerableParameters (cctkGH);
    stop_on_parse_errors = false;

    // We register only once, ergo we get only one handle.  We store
    // that statically, so there is no need to pass anything to
    // Cactus.
    return 0;
  }



  template<int outdim>
  void IOASCII<outdim>
  ::CheckSteerableParameters (const cGH *const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;

    // re-parse the 'IOASCII::out%d_vars' parameter if it has changed
    const char* const out_vars = GetStringParameter("out%dD_vars");
    if (strcmp (out_vars, my_out_vars)) {
      ostringstream parameter_name;
      parameter_name << "IOASCII::out" << outdim << "D_vars";
      IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                                 parameter_name.str().c_str(),
                                 stop_on_parse_errors, out_vars,
                                 -1, &requests[0]);

      // notify the user about the new setting
      if (not CCTK_Equals (verbose, "none")) {
        int count = 0;
        ostringstream msg;
        msg << "Periodic " << outdim << "D AMR output requested for '";
        for (int i = CCTK_NumVars () - 1; i >= 0; i--) {
          if (requests[i]) {
            if (count++) {
              msg << "', '";
            }
            char *fullname = CCTK_FullName (i);
            msg << fullname;
            free (fullname);
          }
        }
        if (count) {
          msg << "'";
          CCTK_INFO (msg.str().c_str());
        }
      }

      // save the last setting of 'IOASCII::out%d_vars' parameter
      free (my_out_vars);
      my_out_vars = strdup (out_vars);
    }
  }



  template<int outdim>
  int IOASCII<outdim>
  ::OutputGH (const cGH* const cctkGH)
  {
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cctkGH, vindex)) {
	TriggerOutput(cctkGH, vindex);
      }
    }
    return 0;
  }



  template<int outdim>
  int IOASCII<outdim>
  ::TimeToOutput (const cGH* const cctkGH, const int vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (vindex >= 0 and vindex < CCTK_NumVars ());

    if (CCTK_GroupTypeFromVarI (vindex) != CCTK_GF and not do_global_mode) {
      return 0;
    }

    CheckSteerableParameters (cctkGH);

    // check if output for this variable was requested
    if (not requests[vindex]) {
      return (0);
    }

    // check whether this refinement level should be output
    if (not (requests[vindex]->refinement_levels & (1 << reflevel))) {
      return (0);
    }

    // check if output for this variable was requested individually
    // by a "<varname>{ out_every = <number> }" option string
    // this will overwrite the output criterion setting
    const char* myoutcriterion = GetStringParameter("out%dD_criterion");
    if (CCTK_EQUALS(myoutcriterion, "default")) {
      myoutcriterion = out_criterion;
    }
    if (requests[vindex]->out_every >= 0) {
      myoutcriterion = "divisor";
    }

    if (CCTK_EQUALS (myoutcriterion, "never")) {
      return (0);
    }

    // check whether to output at this iteration
    bool output_this_iteration = false;

    if (CCTK_EQUALS (myoutcriterion, "iteration")) {
      int myoutevery = GetIntParameter("out%dD_every");
      if (myoutevery == -2) {
        myoutevery = out_every;
      }
      if (myoutevery > 0) {
        if (cctk_iteration == this_iteration[outdim]) {
        // we already decided to output this iteration
          output_this_iteration = true;
        } else if (cctk_iteration
                   >= last_output_iteration[outdim] + myoutevery) {
          // it is time for the next output
          output_this_iteration = true;
          last_output_iteration[outdim] = cctk_iteration;
          this_iteration[outdim] = cctk_iteration;
        }
      }
    } else if (CCTK_EQUALS (myoutcriterion, "divisor")) {
      int myoutevery = GetIntParameter("out%dD_every");
      if (myoutevery == -2) {
        myoutevery = out_every;
      }
      if (requests[vindex]->out_every >= 0) {
        myoutevery = requests[vindex]->out_every;
      }
      if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
        // we already decided to output this iteration
        output_this_iteration = true;
      }
    } else if (CCTK_EQUALS (myoutcriterion, "time")) {
      CCTK_REAL myoutdt = GetRealParameter("out%dD_dt");
      if (myoutdt == -2) {
        myoutdt = out_dt;
      }
      if (myoutdt == 0 or cctk_iteration == this_iteration[outdim]) {
        output_this_iteration = true;
      } else if (myoutdt > 0 and (cctk_time / cctk_delta_time
                 >= (last_output_time[outdim] + myoutdt) / cctk_delta_time - 1.0e-12)) {
        // it is time for the next output
        output_this_iteration = true;
        last_output_time[outdim] = cctk_time;
        this_iteration[outdim] = cctk_iteration;
      }
    } // select output criterion

    return output_this_iteration ? 1 : 0;
  }



  template<int outdim>
  int IOASCII<outdim>
  ::OutputVarAs (const cGH* const cctkGH,
		 const char* const varname, const char* const alias)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    int vindex = -1;

    if (CCTK_TraverseString (varname, GetVarIndex, &vindex, CCTK_VAR) < 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "error while parsing variable name '%s' (alias name '%s')",
                  varname, alias);
      return (-1);
    }

    if (vindex < 0) {
      return (-1);
    }

    if (! is_level_mode()) {
      CCTK_WARN (1, "OutputVarAs must be called in level mode");
      return -1;
    }
    assert (is_level_mode());

    const int group = CCTK_GroupIndexFromVarI (vindex);
    assert (group >= 0);
    const int vindex0 = CCTK_FirstVarIndexI (group);
    assert (vindex0 >= 0 and vindex >= vindex0);
    const int var = vindex - vindex0;
    const int num_tl = CCTK_NumTimeLevelsFromVarI (vindex);
    assert (num_tl >= 1);

    const int grouptype = CCTK_GroupTypeI (group);
    if (grouptype != CCTK_GF) {
      assert (do_global_mode);
    }
    const int rl = grouptype == CCTK_GF ? reflevel : 0;

    const int groupdim = CCTK_GroupDimI(group);
    if (outdim > groupdim) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot produce %dD ASCII output file '%s' for variable '%s' "
                  "because it has only %d dimensions",
                  outdim, alias, varname, groupdim);
      return -1;
    }

    // Check for storage
    if (not CCTK_QueryGroupStorageI (cctkGH, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot output variable '%s' because it has no storage",
		  varname);
      return 0;
    }

    // get the default I/O request for this variable
    ioRequest* request = requests[vindex];
    if (not request) {
      request = IOUtil_DefaultIORequest (cctkGH, vindex, 1);
    }

    // check if the file has been created already
    typedef std::map<string, vector<vector<vector<int> > > > filelist;
    static filelist created_files;
    string basefilename (my_out_dir);
    basefilename.append (alias);
    filelist::iterator thisfile = created_files.find (basefilename);
    bool is_new_file = thisfile == created_files.end();
    if (is_new_file) {
      int const numvars = CCTK_NumVars ();
      vector<vector<vector<int> > > last_outputs;   // [ml][rl][var]
      last_outputs.resize (mglevels);
      for (int ml = 0; ml < mglevels; ++ml) {
        last_outputs[ml].resize (maxreflevels);
        for (int rl = 0; rl < maxreflevels; ++rl) {
          last_outputs[ml][rl].resize (numvars, cctk_iteration - 1);
        }
      }
      thisfile = created_files.insert (thisfile,
                                       filelist::value_type (basefilename,
                                                             last_outputs));
      assert (thisfile != created_files.end());
    }
    is_new_file &= IO_TruncateOutputFiles (cctkGH);

    // check if this variable has been output already during this iteration
    int& last_output = thisfile->second.at(mglevel).at(reflevel).at(vindex);
    if (last_output == cctk_iteration) {
      // Has already been output during this iteration
      char* varname = CCTK_FullName (vindex);
      CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Skipping output for variable '%s', because this variable "
                  "has already been output during the current iteration -- "
                  "probably via a trigger during the analysis stage",
                  varname);
      free (varname);
      return (0);
    }
    assert (last_output < cctk_iteration);
    last_output = cctk_iteration;

    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cctkGH, "IO");
    assert (iogh);

    // Loop over all direction combinations
    vect<int,outdim> dirs (0);

    bool done;
    do {

      // Output each combination only once
      bool ascending = true;
      for (int d1=0; d1<outdim; ++d1) {
	for (int d2=d1+1; d2<outdim; ++d2) {
	  ascending = ascending && dirs[d1] < dirs[d2];
	}
      }

      // Skip output if the dimensions are not ascending
      if (ascending) {

	// If this output is desired at all
	bool desired;
	switch (outdim) {
        case 0:
	  // Output is always desired (if switched on)
          desired = true;
          break;
	case 1:
	  switch (dirs[0]) {
	  case 0:
	    desired = out1D_x;
	    break;
	  case 1:
	    desired = out1D_y;
	    break;
	  case 2:
	    desired = out1D_z;
	    break;
	  default:
	    assert (0);
	  }
	  break;
	case 2:
	  if (dirs[0]==0 && dirs[1]==1) {
	    desired = out2D_xy;
	  } else if (dirs[0]==0 && dirs[1]==2) {
	    desired = out2D_xz;
	  } else if (dirs[0]==1 && dirs[1]==2) {
	    desired = out2D_yz;
	  } else {
	    assert (0);
	  }
	  break;
	case 3:
	  // Output is always desired (if switched on)
	  desired = true;
	  break;
	default:
	  assert (0);
	}

	// Skip output if not desired
	if (desired) {

	  // Traverse all maps on this refinement and multigrid level
	  BEGIN_MAP_LOOP(cctkGH, grouptype) {

            // Find the output offset
            ivect offset(0);
            if (grouptype == CCTK_GF) {
              switch (outdim) {
              case 0:
                offset[0] = GetGridOffset
                  (cctkGH, 1,
                   "out%dD_point_xi", /*"out_point_xi"*/ NULL,
                   "out%dD_point_x",  /*"out_point_x"*/  NULL,
                   /*out_point_x*/ 0.0);
                offset[1] = GetGridOffset
                  (cctkGH, 2,
                   "out%dD_point_yi", /*"out_point_yi"*/ NULL,
                   "out%dD_point_y",  /*"out_point_y"*/  NULL,
                   /*out_point_y*/ 0.0);
                offset[2] = GetGridOffset
                  (cctkGH, 3,
                   "out%dD_point_zi", /*"out_point_zi"*/ NULL,
                   "out%dD_point_z",  /*"out_point_z"*/  NULL,
                   /*out_point_z*/ 0.0);
                break;
              case 1:
                switch (dirs[0]) {
                case 0:
                  offset[1] = GetGridOffset (cctkGH, 2,
                                             "out%dD_xline_yi", "out_xline_yi",
                                             "out%dD_xline_y",  "out_xline_y",
                                             out_xline_y);
                  offset[2] = GetGridOffset (cctkGH, 3,
                                             "out%dD_xline_zi", "out_xline_zi",
                                             "out%dD_xline_z",  "out_xline_z",
                                             out_xline_z);
                  break;
                case 1:
                  offset[0] = GetGridOffset (cctkGH, 1,
                                             "out%dD_yline_xi", "out_yline_xi",
                                             "out%dD_yline_x",  "out_yline_x",
                                             out_yline_x);
                  offset[2] = GetGridOffset (cctkGH, 3,
                                             "out%dD_yline_zi", "out_yline_zi",
                                             "out%dD_yline_z",  "out_yline_z",
                                             out_yline_z);
                  break;
                case 2:
                  offset[0] = GetGridOffset (cctkGH, 1,
                                             "out%dD_zline_xi", "out_zline_xi",
                                             "out%dD_zline_x",  "out_zline_x",
                                             out_zline_x);
                  offset[1] = GetGridOffset (cctkGH, 2,
                                             "out%dD_zline_yi", "out_zline_yi",
                                             "out%dD_zline_y",  "out_zline_y",
                                             out_zline_y);
                  break;
                default:
                  assert (0);
                }
                break;
              case 2:
                if (dirs[0]==0 && dirs[1]==1) {
                  offset[2] = GetGridOffset
                    (cctkGH, 3,
                     "out%dD_xyplane_zi", "out_xyplane_zi",
                     "out%dD_xyplane_z",  "out_xyplane_z",
                     out_xyplane_z);
                } else if (dirs[0]==0 && dirs[1]==2) {
                  offset[1] = GetGridOffset
                    (cctkGH, 2,
                     "out%dD_xzplane_yi", "out_xzplane_yi",
                     "out%dD_xzplane_y",  "out_xzplane_y",
                     out_xzplane_y);
                } else if (dirs[0]==1 && dirs[1]==2) {
                  offset[0] = GetGridOffset
                    (cctkGH, 1,
                     "out%dD_yzplane_xi", "out_yzplane_xi",
                     "out%dD_yzplane_x",  "out_yzplane_x",
                     out_yzplane_x);
                } else {
                  assert (0);
                }
                break;
              case 3:
                // The offset doesn't matter in this case
                break;
              default:
                assert (0);
              }
            } // if grouptype is GF

            ofstream file;
            if (CCTK_MyProc(cctkGH)==0) {

              // Invent a file name
              ostringstream filenamebuf;
              filenamebuf << my_out_dir << "/" << alias << ".";
              if (maps > 1) {
                filenamebuf << Carpet::map << ".";
              }
              if (new_filename_scheme) {
                for (int d=0; d<outdim; ++d) {
                  const char* const coords = "xyz";
                  filenamebuf << coords[dirs[d]];
                }
// The offsets differ per level
//                 for (int dd=0; dd<groupdim; ++dd) {
//                   bool print_dir = true;
//                   for (int d=0; d<outdim; ++d) {
//                     print_dir = print_dir && dirs[d] != dd;
//                   }
//                   if (print_dir) {
//                     filenamebuf << "." << offset[dd];
//                   }
//                 }
                filenamebuf << ".asc";
              } else {
                for (int d=0; d<outdim; ++d) {
                  assert (dirs[d]>=0 && dirs[d]<3);
                  const char* const coords = "xyz";
                  filenamebuf << coords[dirs[d]];
                }
                const char* const suffixes = "plpv";
                filenamebuf << suffixes[outdim];
              }
              // we need a persistent temporary here
              string filenamestr = filenamebuf.str();
              const char* const filename = filenamestr.c_str();


              // Open the file
              file.open (filename, ios::out |
                         (is_new_file ? ios::trunc : ios::app));
              if (! file.good()) {
                CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                            "Could not open output file '%s' for variable '%s'",
                            filename, varname);
              }
              // If this is the first time, then write a nice header
              if (is_new_file) {
                bool want_date = false;
                bool want_parfilename = false;
                bool want_other = false;
                if (CCTK_EQUALS (out_fileinfo, "none")) {
                  // do nothing
                } else if (CCTK_EQUALS (out_fileinfo, "axis labels")) {
                  // do nothing
                } else if (CCTK_EQUALS (out_fileinfo, "creation date")) {
                  want_date = true;
                } else if (CCTK_EQUALS (out_fileinfo, "parameter filename")) {
                  want_parfilename = true;
                } else if (CCTK_EQUALS (out_fileinfo, "all")) {
                  want_date = true;
                  want_parfilename = true;
                  want_other = true;
                } else {
                  CCTK_WARN (0, "internal error");
                }
                file << "# "<< outdim << "D ASCII output created by CarpetIOASCII" << endl;
                if (want_date) {
                  char run_host [1000];
                  Util_GetHostName (run_host, sizeof run_host);
#if 0
                  char const * const run_user = CCTK_RunUser();
#else
                  char const * const run_user = getenv ("USER");
#endif
                  char run_date [1000];
                  Util_CurrentDate (sizeof run_date, run_date);
                  char run_time [1000];
                  Util_CurrentTime (sizeof run_time, run_time);
                  file << "# created on " << run_host
                       << " by " << run_user
                       << " on " << run_date
                       << " at " << run_time << endl;
                }
                if (want_parfilename) {
                  char parameter_filename [10000];
                  CCTK_ParameterFilename (sizeof parameter_filename, parameter_filename);
                  file << "# parameter filename: \"" << parameter_filename << "\"" << endl;
                }
                if (want_other) {
                  if (CCTK_IsFunctionAliased ("UniqueBuildID")) {
                    char const * const build_id
                      = (char const *) UniqueBuildID (cctkGH);
                    file << "# Build ID: " << build_id << endl;
                  }
                  if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
                    char const * const job_id
                      = (char const *) UniqueSimulationID (cctkGH);
                    file << "# Simulation ID: " << job_id << endl;
                  }
                }
                file << "#" << endl;
                if (one_file_per_group) {
                  char* groupname = CCTK_GroupNameFromVarI(vindex);
                  file << "# " << groupname;
                  free (groupname);
                } else {
                  file << "# " << varname;
                }
                for (int d=0; d<outdim; ++d) {
                  file << " " << "xyz"[dirs[d]];
                }
                file << " (" << alias << ")" << endl;
                file << "#" << endl;
              } // if is_new_file

              file << setprecision(out_precision);

            } // if on the root processor

            // Traverse and components on this multigrid and
            // refinement level and map
            BEGIN_COMPONENT_LOOP(cctkGH, grouptype) {

              cGroup groupdata;
              int const ierr = CCTK_GroupData (group, & groupdata);
              assert (! ierr);
              if (groupdata.disttype != CCTK_DISTRIB_CONSTANT
                  or component == 0) {

                const ggf* const ff
                  = arrdata.at(group).at(Carpet::map).data.at(var);

                const int maxtl = output_all_timelevels ? num_tl : 1;
                for (int tl=0; tl<maxtl; ++tl) {

                  const gdata* const data
                    = (*ff) (tl, rl, component, mglevel);
                  ibbox ext = data->extent();

                  ivect lo = ext.lower();
                  ivect hi = ext.upper();
                  ivect str = ext.stride();

                  // Ignore symmetry boundaries if desired
                  if (! output_symmetry_points) {
                    if (grouptype == CCTK_GF) {
                      CCTK_INT const symtable
                        = SymmetryTableHandleForGrid (cctkGH);
                      if (symtable < 0) CCTK_WARN (0, "internal error");
                      CCTK_INT symbnd[2*dim];
                      int const ierr = Util_TableGetIntArray
                        (symtable, 2*dim, symbnd, "symmetry_handle");
                      if (ierr != 2*dim) CCTK_WARN (0, "internal error");
                      for (int d=0; d<dim; ++d) {
                        if (symbnd[2*d] < 0) {
                          // lower boundary is a symmetry boundary
                          lo[d] += cctkGH->cctk_nghostzones[d] * str[d];
                        }
                        if (symbnd[2*d+1] < 0) {
                          // upper boundary is a symmetry boundary
                          hi[d] -= cctkGH->cctk_nghostzones[d] * str[d];
                        }
                      }
                    }
                  }

                  // Ignore ghost zones if desired
                  for (int d=0; d<dim; ++d) {
                    bool output_lower_ghosts
                      = cctkGH->cctk_bbox[2*d] ? out3D_outer_ghosts : out3D_ghosts;
                    bool output_upper_ghosts
                      = cctkGH->cctk_bbox[2*d+1] ? out3D_outer_ghosts : out3D_ghosts;

                    if (! output_lower_ghosts) {
                      lo[d] += cctkGH->cctk_nghostzones[d] * str[d];
                    }
                    if (! output_upper_ghosts) {
                      hi[d] -= cctkGH->cctk_nghostzones[d] * str[d];
                    }
                  }
                  ext = ibbox(lo,hi,str);

                  // coordinates
                  const CCTK_REAL coord_time = cctkGH->cctk_time;
                  rvect global_lower;
                  rvect coord_delta;
                  if (grouptype == CCTK_GF) {
                    for (int d=0; d<dim; ++d) {
                      global_lower[d] = cctkGH->cctk_origin_space[d];
                      coord_delta[d]
                        = cctkGH->cctk_delta_space[d] / maxspacereflevelfact[d];
                    }
                  } else {
                    for (int d=0; d<dim; ++d) {
                      global_lower[d] = 0.0;
                      coord_delta[d] = 1.0;
                    }
                  }
                  const rvect coord_lower
                    = global_lower + coord_delta * rvect(lo);
                  const rvect coord_upper
                    = global_lower + coord_delta * rvect(hi);

                  ivect offset1;
                  if (grouptype == CCTK_GF) {
                    const ibbox& baseext
                      = vdd.at(Carpet::map)->bases.at(mglevel).at(0).exterior;
                    offset1 = baseext.lower() + offset * ext.stride();
                  } else {
                    offset1 = offset * ext.stride();
                  }
                  for (int d=0; d<outdim; ++d) {
                    offset1[dirs[d]] = ext.lower()[dirs[d]];
                  }

                  vector<const gdata*> datas;
                  if (one_file_per_group) {
                    int const numvars = CCTK_NumVarsInGroupI(group);
                    datas.resize (numvars);
                    for (int n=0; n<numvars; ++n) {
                      const ggf* const ff1
                        = arrdata.at(group).at(Carpet::map).data.at(n);
                      datas.at(n) = (*ff1) (tl, rl, component, mglevel);
                    }
                  } else {
                    datas.resize (1);
                    datas.at(0) = data;
                  }
                  WriteASCII (file, datas, ext, vindex, cctkGH->cctk_iteration,
                              offset1, dirs,
                              rl, mglevel, Carpet::map, component, tl,
                              coord_time, coord_lower, coord_upper);

                  // Append EOL after every component
                  if (CCTK_MyProc(cctkGH)==0) {
                    if (separate_components) {
                      assert (file.good());
                      file << endl;
                    }
                  }
                  assert (file.good());

                } // for tl

              } // if distrib!=CONST or component==0

            } END_COMPONENT_LOOP;

            // Append EOL after every complete set of components
            if (CCTK_MyProc(cctkGH)==0) {
              if (separate_grids) {
                assert (file.good());
                file << endl;
              }
              file.close();
              assert (file.good());
            }

            assert (! file.is_open());

          } END_MAP_LOOP;

	} // if (desired)

      } // if (ascending)

      // Next direction combination
      done = true;
      for (int d=0; d<outdim; ++d) {
	++dirs[d];
	if (dirs[d]<groupdim) {
	  done = false;
	  break;
	}
	dirs[d] = 0;
      }

    } while (! done);		// all directions

    return 0;
  }



  template<int outdim>
  int IOASCII<outdim>
  ::TriggerOutput (const cGH* const cctkGH, const int vindex)
  {
    DECLARE_CCTK_PARAMETERS;

    assert (vindex >= 0 and vindex < CCTK_NumVars ());

    char* const fullname = CCTK_FullName(vindex);

    int retval;
    if (one_file_per_group) {
      char* alias = CCTK_GroupNameFromVarI (vindex);
      for (char* p = alias; *p; ++p) *p = (char) tolower (*p);
      retval = OutputVarAs (cctkGH, fullname, alias);
      free (alias);
    } else {
      const char* const alias = CCTK_VarName (vindex);
      retval = OutputVarAs (cctkGH, fullname, alias);
    }

    free (fullname);

    return retval;
  }



  template<int outdim>
  int IOASCII<outdim>
  ::GetGridOffset (const cGH* const cctkGH, const int dir,
		   const char* const itempl, const char* const iglobal,
		   const char* const ctempl, const char* const cglobal,
		   const CCTK_REAL cfallback)
  {
    // First choice: explicit coordinate
    char cparam[1000];
    snprintf (cparam, sizeof cparam, ctempl, outdim);
    const int ncparam = CCTK_ParameterQueryTimesSet (cparam, CCTK_THORNSTRING);
    assert (ncparam >= 0);
    if (ncparam > 0) {
      int ptype;
      const CCTK_REAL* const pcoord
	= ((const CCTK_REAL*)CCTK_ParameterGet
           (cparam, CCTK_THORNSTRING, &ptype));
      assert (pcoord);
      const CCTK_REAL coord = *pcoord;
      assert (ptype == PARAMETER_REAL);
      return CoordToOffset (cctkGH, dir, coord, 0);
    }

    // Second choice: explicit index
    char iparam[1000];
    snprintf (iparam, sizeof iparam, itempl, outdim);
    const int niparam = CCTK_ParameterQueryTimesSet (iparam, CCTK_THORNSTRING);
    assert (niparam >= 0);
    if (niparam > 0) {
      int ptype;
      const int* const pindex
	= (const int*)CCTK_ParameterGet (iparam, CCTK_THORNSTRING, &ptype);
      assert (pindex);
      const int index = *pindex;
      assert (ptype == PARAMETER_INT);
      return index;
    }

    // Third choice: explicit global coordinate
    const char* iothorn = CCTK_ImplementationThorn ("IO");
    assert (iothorn);
    if (cglobal) {
      const int ncglobal = CCTK_ParameterQueryTimesSet (cglobal, iothorn);
      assert (ncglobal >= 0);
      if (ncglobal > 0) {
        int ptype;
        const CCTK_REAL* const pcoord
          = (const CCTK_REAL*)CCTK_ParameterGet (cglobal, iothorn, &ptype);
        assert (pcoord);
        const CCTK_REAL coord = *pcoord;
        assert (ptype == PARAMETER_REAL);
        return CoordToOffset (cctkGH, dir, coord, 0);
      }
    }

    // Fourth choice: explicit global index
    if (iglobal) {
      const int niglobal = CCTK_ParameterQueryTimesSet (iglobal, iothorn);
      assert (niglobal >= 0);
      if (niglobal > 0) {
        int ptype;
        const int* const pindex
          = (const int*)CCTK_ParameterGet (iglobal, iothorn, &ptype);
        assert (pindex);
        const int index = *pindex;
        assert (ptype == PARAMETER_INT);
        return index;
      }
    }

    // Fifth choice: default coordinate
    return CoordToOffset (cctkGH, dir, cfallback, 0);
  }



  template<int outdim>
  int IOASCII<outdim>
  ::CoordToOffset (const cGH* cctkGH, const int dir, const CCTK_REAL coord,
		   const int ifallback)
  {
    assert (dir>=1 && dir<=dim);

    assert (mglevel!=-1 && reflevel!=-1 && Carpet::map!=-1);

    const CCTK_REAL delta = cctkGH->cctk_delta_space[dir-1] / cctkGH->cctk_levfac[dir-1];
    const CCTK_REAL lower = cctkGH->cctk_origin_space[dir-1];
#if 0
    const int npoints = cctkGH->cctk_gsh[dir-1];
    const CCTK_REAL upper = lower + (npoints-1) * delta;
#endif

    const CCTK_REAL rindex = (coord - lower) / delta;
    int cindex = (int)floor(rindex + 0.75);

#if 0
    if (cindex<0 || cindex>=npoints) {
      cindex = ifallback;

      assert (dir>=1 && dir<=3);
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The specified coordinate value %g for the %c-direction is not within the grid range [%g,%g] on convergence level %d, refinement level %d, map %d; using %g instead",
                  coord, "xyz"[dir-1], lower, upper,
                  mglevel, reflevel, Carpet::map, lower + delta * cindex);
    }

    assert (cindex>=0 && cindex<npoints);
#else
    const void *dummy;
    dummy = &ifallback;
    dummy = &dummy;
#endif

    return cindex;
  }



  template<int outdim>
  const char* IOASCII<outdim>
  ::GetStringParameter (const char* const parametertemplate)
  {
    char parametername[1000];
    snprintf (parametername, sizeof parametername, parametertemplate, outdim);
    int ptype;
    const char* const* const ppval = (const char* const*)CCTK_ParameterGet
      (parametername, CCTK_THORNSTRING, &ptype);
    assert (ppval);
    const char* const pval = *ppval;
    assert (ptype == PARAMETER_STRING || ptype == PARAMETER_KEYWORD);
    return pval;
  }



  template<int outdim>
  CCTK_INT IOASCII<outdim>
  ::GetIntParameter (const char* const parametertemplate)
  {
    char parametername[1000];
    snprintf (parametername, sizeof parametername, parametertemplate, outdim);
    int ptype;
    const CCTK_INT* const ppval
      = (const CCTK_INT*)CCTK_ParameterGet
      (parametername, CCTK_THORNSTRING, &ptype);
    assert (ppval);
    assert (ptype == PARAMETER_INT || ptype == PARAMETER_BOOLEAN);
    const CCTK_INT pval = *ppval;
    return pval;
  }



  template<int outdim>
  CCTK_REAL IOASCII<outdim>
  ::GetRealParameter (const char* const parametertemplate)
  {
    char parametername[1000];
    snprintf (parametername, sizeof parametername, parametertemplate, outdim);
    int ptype;
    const CCTK_REAL* const ppval
      = (const CCTK_REAL*)CCTK_ParameterGet
      (parametername, CCTK_THORNSTRING, &ptype);
    assert (ppval);
    assert (ptype == PARAMETER_REAL);
    const CCTK_REAL pval = *ppval;
    return pval;
  }



  static void GetVarIndex (int vindex, const char* optstring, void* arg)
  {
    if (optstring) {
      char *fullname = CCTK_FullName (vindex);
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Option string '%s' will be ignored for ASCII output of "
                  "variable '%s'", optstring, fullname);
      free (fullname);
    }

    *((int *) arg) = vindex;
  }



  static CCTK_REAL
  nicelooking (CCTK_REAL const val,
               CCTK_REAL const base)
  {
    return floor(val / base + 0.5) * base;
  }



  // Output
  template<int outdim>
  void WriteASCII (ostream& os,
		   vector<const gdata*> const gfdatas,
		   const bbox<int,dim>& gfext,
		   const int vi,
		   const int time,
		   const vect<int,dim>& org,
		   const vect<int,outdim>& dirs,
		   const int rl,
		   const int ml,
                   const int m,
		   const int c,
		   const int tl,
		   const CCTK_REAL coord_time,
		   const vect<CCTK_REAL,dim>& coord_lower,
		   const vect<CCTK_REAL,dim>& coord_upper)
  {
    DECLARE_CCTK_PARAMETERS;

    assert (outdim<=dim);
    const int vartype = CCTK_VarTypeI(vi);

    bool all_on_root = true;
    for (size_t n=0; n<gfdatas.size(); ++n) {
      all_on_root &= gfdatas.at(n)->proc() == 0;
    }
    if (all_on_root) {
      // output on processor 0

      if (dist::rank() == 0) {

	assert (os.good());

	os << "# iteration " << time << endl
	   << "# refinement level " << rl
           << "   multigrid level " << ml
	   << "   map " << m
	   << "   component " << c
           << "   time level " << tl
           << endl
	   << "# column format: 1:it\t2:tl 3:rl 4:c 5:ml";
        int col=6;
	assert (dim>=1 && dim<=3);
	const char* const coords = "xyz";
	for (int d=0; d<dim; ++d) {
          os << (d==0 ? "\t" : " ") << col++ << ":i" << coords[d];
        }
	os << "\t" << col++ << ":time";
	for (int d=0; d<dim; ++d) {
          os << (d==0 ? "\t" : " ") << col++ << ":" << coords[d];
        }
	os << "\t" << col << ":data" << endl;
        if (one_file_per_group) {
          os << "# data columns:";
          int const gindex = CCTK_GroupIndexFromVarI(vi);
          int const firstvar = CCTK_FirstVarIndexI(gindex);
          int const numvars = CCTK_NumVarsInGroupI(gindex);
          for (int n=firstvar; n<firstvar+numvars; ++n) {
            os << " " << col << ":" << CCTK_VarName(n);
            col += CarpetSimpleMPIDatatypeLength (vartype);
          }
          os << endl;
        }

	const vect<int,outdim> lo = gfext.lower()[dirs];
	const vect<int,outdim> up = gfext.upper()[dirs];
	const vect<int,outdim> str = gfext.stride()[dirs];
	const bbox<int,outdim> ext(lo,up,str);

	// Check whether the output origin is contained in the extent
	// of the data that should be output
	ivect org1(org);
	for (int d=0; d<outdim; ++d) org1[dirs[d]] = ext.lower()[d];
	if (gfext.contains(org1)) {

          typename bbox<int,outdim>::iterator it=ext.begin();
          do {

	    ivect index(org);
	    for (int d=0; d<outdim; ++d) index[dirs[d]] = (*it)[d];
	    os << time << "\t" << tl << " " << rl << " " << c << " " << ml
               << "\t";
	    for (int d=0; d<dim-1; ++d) os << index[d] << " "; os << index[dim-1];
	    os << "\t" << coord_time << "\t";
	    for (int d=0; d<dim; ++d) {
              if (d > 0) os << " ";
	      assert (gfext.upper()[d] - gfext.lower()[d] >= 0);
	      if (gfext.upper()[d] - gfext.lower()[d] == 0) {
                os << coord_lower[d];
              } else {
                CCTK_REAL const dx = ((coord_upper[d] - coord_lower[d])
                                      / (gfext.upper()[d] - gfext.lower()[d]));
                os << (nicelooking
                       (coord_lower[d] + (index[d] - gfext.lower()[d]) * dx,
                        dx * 1.0e-8));
              }
	    }
            os << "\t";
            for (size_t n=0; n<gfdatas.size(); ++n) {
              const gdata* gfdata = gfdatas.at(n);
              if (n > 0) os << " ";
              switch (vartype) {
#define TYPECASE(N,T)                                           \
                case N:                                         \
                  os << (*(const data<T>*)gfdata)[index];	\
                break;
#include "carpet_typecase.hh"
#undef TYPECASE
              default:
                UnsupportedVarType(vi);
              }
            } // for n
	    os << endl;

            ++it;

	    for (int d=0; d<outdim; ++d) {
	      if ((*it)[d]!=(*ext.end())[d]) break;
	      os << endl;
	    }

          } while (it!=ext.end());

	} else {

          os << "#" << endl;

	} // if ! ext contains org

	assert (os.good());

      }

    } else {
      // copy to processor 0 and output there

      vector<const gdata*> tmps (gfdatas.size());
      for (size_t n=0; n<gfdatas.size(); ++n) {
        gdata * const tmp = gfdatas.at(n)->make_typed(vi);
        tmp->allocate(gfdatas.at(n)->extent(), 0);
        for (comm_state state; !state.done(); state.step()) {
          tmp->copy_from (state, gfdatas.at(n), gfdatas.at(n)->extent());
        }
        tmps.at(n) = tmp;
      }
      WriteASCII (os, tmps, gfext, vi, time, org, dirs, rl, ml, m, c, tl,
		  coord_time, coord_lower, coord_upper);
      for (size_t n=0; n<gfdatas.size(); ++n) {
        delete tmps.at(n);
      }

    }
  }





  // Explicit instantiation for all output dimensions
  template class IOASCII<0>;
  template class IOASCII<1>;
  template class IOASCII<2>;
  template class IOASCII<3>;

  template
  void WriteASCII (ostream& os,
		   vector<const gdata*> const gfdatas,
		   const bbox<int,dim>& gfext,
		   const int vi,
		   const int time,
		   const vect<int,dim>& org,
		   const vect<int,0>& dirs,
		   const int rl,
		   const int ml,
		   const int m,
		   const int c,
		   const int tl,
		   const CCTK_REAL coord_time,
		   const vect<CCTK_REAL,dim>& coord_lower,
		   const vect<CCTK_REAL,dim>& coord_upper);

  template
  void WriteASCII (ostream& os,
		   vector<const gdata*> const gfdatas,
		   const bbox<int,dim>& gfext,
		   const int vi,
		   const int time,
		   const vect<int,dim>& org,
		   const vect<int,1>& dirs,
		   const int rl,
		   const int ml,
		   const int m,
		   const int c,
		   const int tl,
		   const CCTK_REAL coord_time,
		   const vect<CCTK_REAL,dim>& coord_lower,
		   const vect<CCTK_REAL,dim>& coord_upper);

  template
  void WriteASCII (ostream& os,
		   vector<const gdata*> const gfdatas,
		   const bbox<int,dim>& gfext,
		   const int vi,
		   const int time,
		   const vect<int,dim>& org,
		   const vect<int,2>& dirs,
		   const int rl,
		   const int ml,
		   const int m,
		   const int c,
		   const int tl,
		   const CCTK_REAL coord_time,
		   const vect<CCTK_REAL,dim>& coord_lower,
		   const vect<CCTK_REAL,dim>& coord_upper);

  template
  void WriteASCII (ostream& os,
		   vector<const gdata*> const gfdatas,
		   const bbox<int,dim>& gfext,
		   const int vi,
		   const int time,
		   const vect<int,dim>& org,
		   const vect<int,3>& dirs,
		   const int rl,
		   const int ml,
		   const int m,
		   const int c,
		   const int tl,
		   const CCTK_REAL coord_time,
		   const vect<CCTK_REAL,dim>& coord_lower,
		   const vect<CCTK_REAL,dim>& coord_upper);

} // namespace CarpetIOASCII
