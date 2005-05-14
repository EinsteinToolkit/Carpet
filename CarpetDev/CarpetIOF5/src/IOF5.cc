#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "coordinate_system.hh"
#include "data_region.hh"
#include "extending.hh"
#include "file.hh"
#include "physical_quantity.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"
#include "utils.hh"



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
  write_one_mglevel (cGH const * cctkGH,
                     F5::timestep_t & timestep,
                     int group,
                     int variable);
  
  static void
  write_global (cGH const * cctkGH,
                F5::simulation_t & simulation,
                int group,
                int variable);

  static void
  write_one_reflevel (cGH const * cctkGH,
                      F5::simulation_t & simulation,
                      int group,
                      int variable);
  
  static void
  write_one_map (cGH const * cctkGH,
                 F5::tensor_component_t & tensor_component,
                 int group);
  
  static void
  write_one_component (cGH const * cctkGH,
                       F5::tensor_component_t & tensor_component);
  
  
  
  static void
  WriteParameters (F5::file_t & file);
  
  
  
  int
  CarpetIOF5_Startup ()
  {
    CCTK_INFO ("Startup");
    return extending_t::create (Setup);
  }
  
  
  
  void *
  Setup (tFleshConfig * const fleshconfig,
         int const convlevel,
         cGH * const cctkGH)
  {
    assert (fleshconfig != 0);
    CCTK_INFO ("Setup");
    return extending_t::setup
      (cctkGH, OutputGH, TimeToOutput, TriggerOutput, OutputVarAs);
  }
  
  
  
  int
  OutputGH (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    
    CCTK_INFO ("OutputGH");
    
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
    
    CCTK_VInfo (CCTK_THORNSTRING, "OutputVarAs \"%s\" \"%s\"",
                varname, alias);
    
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
    
    ostringstream filenamebuf;
    bool const use_IO_out_dir = strcmp (out_dir, "") == 0;
    filenamebuf << (use_IO_out_dir ? IO_out_dir : out_dir)
                << "/"
                << alias
                << out_extension;
    string const filename = filenamebuf.str();
    
    bool const did_truncate = extending.get_did_truncate (filename);
    bool const do_truncate
      = ! did_truncate and IO_TruncateOutputFiles (cctkGH);
    extending.set_did_truncate (filename);
    
    F5::file_t file (cctkGH, filename, do_truncate);
    
    F5::timestep_t timestep (file, cctkGH->cctk_time);
    
    if (Carpet::is_meta_mode())
    {
      for (Carpet::mglevel_iterator mglevel_iter (cctkGH);
           ! mglevel_iter.done();
           mglevel_iter.step())
      {
        write_one_mglevel (cctkGH, timestep, group, variable);
      }
    }
    else
    {
      write_one_mglevel (cctkGH, timestep, group, variable);
    }
    
    return Error_none;
  }
  
  
  
  void
  write_one_mglevel (cGH const * const cctkGH,
                     F5::timestep_t & timestep,
                     int const group,
                     int const variable)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "OutputVarAs/write_one_mglevel mglevel=%d", Carpet::mglevel);
    
    ostringstream namebuf;
    namebuf << "convlevel=" << cctkGH->cctk_convlevel;
    string const namestr = namebuf.str();
    char const * const name = namestr.c_str();
    F5::simulation_t simulation (timestep, name);
    
    int const grouptype = CCTK_GroupTypeI (group);
    assert (grouptype >= 0);
    switch (grouptype)
    {
    case CCTK_ARRAY:
    case CCTK_SCALAR:
      {
        if (Carpet::do_global_mode)
        {
          write_global (cctkGH, simulation, group, variable);
        }
      }
      break;
    case CCTK_GF:
      {
        if (Carpet::is_global_mode())
        {
          for (Carpet::reflevel_iterator reflevel_iter (cctkGH);
               ! reflevel_iter.done();
               reflevel_iter.step())
          {
            write_one_reflevel (cctkGH, simulation, group, variable);
          }
        }
        else
        {
          write_one_reflevel (cctkGH, simulation, group, variable);
        }
      }
      break;
    default:
      assert (0);
    }
  }
  
  
  
  void
  write_global (cGH const * const cctkGH,
                F5::simulation_t & simulation,
                int const group,
                int const variable)
  {
    CCTK_INFO ("OutputVarAs/write_global");
    
    F5::unigrid_topology_t topology (simulation);
    
    int const grouptype = CCTK_GroupTypeI (group);
    assert (grouptype >= 0);
    assert (grouptype == CCTK_SCALAR or grouptype == CCTK_ARRAY);
    
    vect<CCTK_REAL, dim> level_origin, level_delta;
    for (int d=0; d<dim; ++d)
    {
      level_origin[d] = 0.0;
      level_delta[d]  = 1.0;
    }
    F5::Cartesian_coordinate_system_t coordinate_system
      (topology, level_origin, level_delta);
    
    F5::physical_quantity_t physical_quantity (coordinate_system, group);
    
    F5::tensor_component_t tensor_component (physical_quantity, variable);
    
    int const myproc = CCTK_MyProc (cctkGH);
    dh * const dd = Carpet::arrdata.at(group).at(0).dd;
    bbox<int, dim> const & region
      = dd->boxes.at(Carpet::mglevel).at(0).at(myproc).exterior;
    
    F5::data_region_t data_region (tensor_component, region);
    
    void const * const varptr = CCTK_VarDataPtrI (cctkGH, 0, variable);
    assert (varptr != 0);
    int const vartype = CCTK_VarTypeI (variable);
    assert (vartype >= 0);
    data_region.write (varptr, vartype);
  }
  
  
  
  void
  write_one_reflevel (cGH const * const cctkGH,
                      F5::simulation_t & simulation,
                      int const group,
                      int const variable)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "OutputVarAs/write_one_reflevel reflevel=%d",
                Carpet::reflevel);
    
    int const grouptype = CCTK_GroupTypeI (group);
    assert (grouptype >= 0);
    assert (grouptype == CCTK_GF);
    
    F5::mesh_refinement_topology_t topology
      (simulation, Carpet::reflevel, Carpet::maxreflevels,
       Carpet::spacereflevelfact, Carpet::maxspacereflevelfact);
    
    vect<CCTK_REAL, dim> level_origin, level_delta;
    for (int d=0; d<dim; ++d)
    {
      DECLARE_CCTK_ARGUMENTS;
      level_origin[d] = CCTK_ORIGIN_SPACE(d);
      level_delta[d]  = CCTK_DELTA_SPACE(d);
    }
    F5::Cartesian_coordinate_system_t coordinate_system
      (topology, level_origin, level_delta);
    
    F5::physical_quantity_t physical_quantity (coordinate_system, group);
    
    F5::tensor_component_t tensor_component (physical_quantity, variable);
    
    if (Carpet::is_level_mode())
    {
      for (Carpet::map_iterator map_iter (cctkGH, grouptype);
           ! map_iter.done();
           map_iter.step())
      {
        write_one_map (cctkGH, tensor_component, group);
      }
    }
    else
    {
      write_one_map (cctkGH, tensor_component, group);
    }
  }
  
  
  
  void
  write_one_map (cGH const * const cctkGH,
                 F5::tensor_component_t & tensor_component,
                 int const group)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "OutputVarAs/write_one_map map=%d", Carpet::map);
    
    if (Carpet::is_singlemap_mode())
    {
      int const grouptype = CCTK_GroupTypeI (group);
      assert (grouptype >= 0);
      
      for (Carpet::component_iterator component_iter (cctkGH, grouptype);
           ! component_iter.done();
           component_iter.step())
      {
        write_one_component (cctkGH, tensor_component);
      }
    }
    else
    {
      write_one_component (cctkGH, tensor_component);
    }
  }
  
  
  
  void
  write_one_component (cGH const * const cctkGH,
                       F5::tensor_component_t & tensor_component)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "OutputVarAs/write_one_component component=%d",
                Carpet::component);
    
    gh * const hh = Carpet::vhh.at(Carpet::map);
    if (hh->is_local (Carpet::reflevel, Carpet::component))
    {
      dh * const dd = Carpet::vdd.at(Carpet::map);
      bbox<int, dim> const & region
        = (dd->boxes.at(Carpet::mglevel).at(Carpet::reflevel)
           .at(Carpet::component).exterior);
      
      F5::data_region_t data_region (tensor_component, region);
      
      int const variable = tensor_component.get_variable();
      void const * const varptr = CCTK_VarDataPtrI (cctkGH, 0, variable);
      assert (varptr != 0);
      int const vartype = CCTK_VarTypeI (variable);
      assert (vartype >= 0);
      data_region.write (varptr, vartype);
    }
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
    assert (! herr);
  }
  
} // namespace CarpetIOF5
