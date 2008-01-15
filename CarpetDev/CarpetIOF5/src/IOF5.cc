#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "extending.hh"
#include "file.hh"
#include "utils.hh"
#include "writer.hh"



namespace CarpetIOF5 {
  
  
  
  int const Error_none                        =  0; // no error
  
  int const Error_illegal_varname             = -1;
  int const Error_group_has_no_storage        = -2;
  
  int const Error_lonely_option_string        = -3;
  int const Error_unterminated_option_string  = -4;
  int const Error_garbage_after_option_string = -5;
  int const Error_invalid_variable_name       = -6;
  
  
  
  extern "C" int
  CarpetIOF5_Startup ();
  
  static void *
  Setup (tFleshConfig * const fleshconfig,
         int const convlevel,
         cGH * const cctkGH);
  
  extern "C" void
  CarpetIOF5_Init (CCTK_ARGUMENTS);
  
  
  
  static int
  OutputGH (cGH const * cctkGH);
  
  static void
  mark_variables (int variable,
                  char const * options,
                  void * ptr);
  
  
  
  static int
  TimeToOutput (cGH const * cctkGH,
                int variable);
  
  
  
  static int
  TriggerOutput (cGH const * cctkGH,
                 int variable);
  
  
  
  static int
  OutputVarAs (cGH const * cctkGH,
               char const * varname,
               char const * alias);
  
  
  
  static void
  WriteParameters (F5::file_t & file);
  
  
  
  int
  CarpetIOF5_Startup ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose) CCTK_INFO ("Startup");
    return extending_t::create (Setup);
  }
  
  
  
  void *
  Setup (tFleshConfig * const fleshconfig,
         int const convlevel,
         cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (fleshconfig != 0);
    if (verbose) CCTK_INFO ("Setup");
    return extending_t::setup
      (cctkGH, OutputGH, TimeToOutput, TriggerOutput, OutputVarAs);
  }
  
  
  
  void
  CarpetIOF5_Init (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose) CCTK_INFO ("Init");

    * next_output_iteration = 0;
    * next_output_time = cctk_time;
    * this_iteration = -1;
  }
  
  
  
  int
  OutputGH (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    
    if (verbose) CCTK_INFO ("OutputGH");
    
    int ierr;
    
    vector<bool> want_variables (CCTK_NumVars());
    ierr
      = CCTK_TraverseString (out_vars,
                             mark_variables,
                             static_cast<void *> (& want_variables),
                             CCTK_GROUP_OR_VAR);
    switch (ierr) {
    case -2: return Error_lonely_option_string;
    case -3: return Error_unterminated_option_string;
    case -4: return Error_garbage_after_option_string;
    case -5: return Error_invalid_variable_name;
    }
    assert (ierr >= 0);
    
    ierr = Error_none;
    for (int variable = 0; variable < CCTK_NumVars(); ++ variable)
    {
      if (want_variables.at(variable))
      {
        if (TimeToOutput (cctkGH, variable))
        {
          ierr = TriggerOutput (cctkGH, variable);
        }
      }
    }
    
    return ierr;
  }
  
  
  
  void
  mark_variables (int const variable,
                  char const * const options,
                  void * const ptr)
  {
    vector<bool> & want_variables = * static_cast<vector<bool> *> (ptr);
    want_variables.at (variable) = true;
  }
  
  
  
  int
  TimeToOutput (cGH const * const cctkGH,
                int const variable)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    assert (variable >= 0 and variable < CCTK_NumVars());
    
    assert (Carpet::is_level_mode());
    
    if (verbose)
    {
      char * const fullname = CCTK_FullName(variable);
      CCTK_VInfo (CCTK_THORNSTRING, "TimeToOutput \"%s\"", fullname);
      free (fullname);
    }
    
    bool should_output;
    
    char const * const my_out_criterion
      = (CCTK_EQUALS (out_criterion, "default")
         ? IO_out_criterion
         : out_criterion);
    
    if (CCTK_EQUALS (my_out_criterion, "always"))
    {
      should_output = true;
    }
    else if (CCTK_EQUALS (my_out_criterion, "never"))
    {
      should_output = false;
    }
    else if (CCTK_EQUALS (my_out_criterion, "iteration"))
    {
      int const my_out_every = out_every == -2 ? IO_out_every : out_every;
      switch (my_out_every)
      {
      case 0:
        should_output = true;
        break;
      case -1:
        should_output = false;
        break;
      default:
        if (* this_iteration == cctk_iteration)
        {
          // we already decided to output this iteration
          should_output = true;
        }
        else if (cctk_iteration >= * next_output_iteration)
        {
          // it is time for the next output
          should_output = true;
          * this_iteration = cctk_iteration;
          * next_output_iteration = cctk_iteration + my_out_every;
        }
        else
        {
          should_output = false;
        }
        break;
      }
    }
    else if (CCTK_EQUALS (my_out_criterion, "time"))
    {
      CCTK_REAL const my_out_dt = out_dt == -2 ? IO_out_dt : out_dt;
      if (out_dt == 0)
      {
        should_output = true;
      }
      else if (out_dt == -1)
      {
        should_output = false;
      }
      else
      {
        if (* this_iteration == cctk_iteration)
        {
          // we already decided to output this iteration
          should_output = true;
        }
        else if (cctk_time / cctk_delta_time
                 >= * next_output_time / cctk_delta_time - dt_fudge)
        {
          // it is time for the next output
          should_output = true;
          * this_iteration = cctk_iteration;
          * next_output_time = cctk_time + my_out_dt;
        }
        else
        {
          should_output = false;
        }
      }
    }
    else
    {
      CCTK_WARN (1, "internal error");
      should_output = false;
    }
    
    if (should_output)
    {
      extending_t extending (cctkGH);
      int const last_output_iteration
        = (extending.get_last_output_iteration
           (Carpet::mglevel, Carpet::reflevel, variable));
      assert (last_output_iteration <= cctk_iteration);
      if (last_output_iteration == cctk_iteration)
      {
        // Skipping output for variable, because this variable has
        // already been output during the current iteration --
        // probably via a trigger during the analysis stage
        should_output = false;
      }
      else
      {
        extending.set_last_output_iteration
          (Carpet::mglevel, Carpet::reflevel, variable, cctk_iteration);
      }
    }

    return should_output;
  }
  
  
  
  int
  TriggerOutput (cGH const * const cctkGH,
                 int const variable)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    assert (variable >= 0 and variable < CCTK_NumVars());
    
    if (verbose)
    {
      char * const fullname = CCTK_FullName(variable);
      CCTK_VInfo (CCTK_THORNSTRING, "TriggerOutput \"%s\"", fullname);
      free (fullname);
    }
    
    char * fullname = CCTK_FullName (variable);
    assert (fullname);
    
    int const ierr = OutputVarAs (cctkGH, fullname, out_filename);
    
    free (fullname);
    
    return ierr;
  }
  
  
  
  int
  OutputVarAs (cGH const * const cctkGH,
               char const * const varname,
               char const * const alias)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    assert (varname != 0);
    assert (alias != 0);
    
    if (verbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "OutputVarAs \"%s\" \"%s\"",
                  varname, alias);
    }
    
    int const variable = CCTK_VarIndex (varname);
    if (variable < 0)
    {
      return Error_illegal_varname;
    }
    assert (variable >= 0 and variable < CCTK_NumVars());
    
    int const group = CCTK_GroupIndexFromVarI (variable);
    assert (group >= 0 and group < CCTK_NumGroups());
    
    if (CCTK_ActiveTimeLevelsGI (cctkGH, group) == 0)
    {
      return Error_group_has_no_storage;
    }
    
    extending_t extending (cctkGH);
    
    bool const use_IO_out_dir = std::strcmp (out_dir, "") == 0;
    string const basename
      = (use_IO_out_dir ? IO_out_dir : out_dir) + string ("/") + alias;
    bool const want_metafile = CCTK_MyProc (cctkGH) == 0;
    
    bool const did_truncate = extending.get_did_truncate (basename);
    bool const do_truncate
      = not did_truncate and IO_TruncateOutputFiles (cctkGH);
    extending.set_did_truncate (basename);
    
    F5::file_t file (cctkGH, basename, string (out_extension),
                     want_metafile, do_truncate);
    
    writer_t writer (cctkGH, variable);
    writer.write (file);
    
    return Error_none;
  }
  
  
  
  void
  WriteParameters (F5::file_t & file)
  {
    hid_t const parameter_group
      = F5::open_or_create_group (file.get_hdf5_file(), "Cactus parameters");
    
    int first = 1;
    for (;;)
    {
      cParamData const * parameter_data;
      char * parameter_fullname;
      
      int const ierr
        = CCTK_ParameterWalk (first, 0,
                              & parameter_fullname, & parameter_data);
      if (ierr > 0) break;
      assert (ierr >= 0);
      
      int type;
      void const * const parameter_value
        = CCTK_ParameterGet (parameter_data->name, parameter_data->thorn,
                             & type);
      assert (type == parameter_data->type);
      assert (parameter_value != 0);
      
      switch (type)
      {
      case PARAMETER_BOOLEAN:
      case PARAMETER_INT:
        {
          CCTK_INT const value
            = * static_cast<CCTK_INT const *> (parameter_value);
          F5::write_or_check_attribute
            (parameter_group, parameter_fullname, value);
        }
        break;
      case PARAMETER_REAL:
        {
          CCTK_REAL const value
            = * static_cast<CCTK_REAL const *> (parameter_value);
          F5::write_or_check_attribute
            (parameter_group, parameter_fullname, value);
        }
        break;
      case PARAMETER_KEYWORD:
      case PARAMETER_STRING:
        {
          char const * const value
            = * static_cast<char const * const *> (parameter_value);
          F5::write_or_check_attribute
            (parameter_group, parameter_fullname, value);
        }
        break;
      default:
        assert (0);
      }
      
      free (parameter_fullname);
      
      first = 0;
    }
    
    herr_t const herr = H5Gclose (parameter_group);
    assert (not herr);
  }
  
} // namespace CarpetIOF5
