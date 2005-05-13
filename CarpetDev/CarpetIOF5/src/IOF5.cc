#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sstream>
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
  SetupGH (tFleshConfig * fleshconfig,
           int convlevel,
           cGH * cctkGH);
  
  
  
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
  
  
  
  extern "C" int
  CarpetIOF5_Startup ()
  {
    int const ierr
      = CCTK_RegisterBanner ("AMR HDF5 I/O provided by CarpetIOF5");
    assert (! ierr);
    
    extending_t::create();
    
    return 0;                   // no error
  }
  
  
  
  void *
  SetupGH (tFleshConfig * const fleshconfig,
           int const convlevel,
           cGH * const cctkGH)
  {
    int ierr;
    int const io_method = CCTK_RegisterIOMethod ("IOF5");
    ierr = CCTK_RegisterIOMethodOutputGH (io_method, OutputGH);
    assert (! ierr);
    ierr = CCTK_RegisterIOMethodOutputVarAs (io_method, OutputVarAs);
    assert (! ierr);
    ierr = CCTK_RegisterIOMethodTimeToOutput (io_method, TimeToOutput);
    assert (! ierr);
    ierr = CCTK_RegisterIOMethodTriggerOutput (io_method, TriggerOutput);
    assert (! ierr);
    
    return extending_t::setup (fleshconfig, convlevel, cctkGH);
  }
  
  
  
  int
  OutputGH (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    
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
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH != 0);
    assert (variable >= 0 and variable < CCTK_NumVars());
    
    assert (Carpet::is_level_mode());
    extending_t extending (cctkGH);
    
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
      int const last_output_iteration
        = (extending.get_last_output_iteration
           (Carpet::mglevel, Carpet::reflevel, variable));
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
        should_output
          = cctkGH->cctk_iteration >= last_output_iteration + my_out_every;
        break;
      }
    }
    else if (CCTK_EQUALS (my_out_criterion, "time"))
    {
      CCTK_REAL const last_output_time
        = (extending.get_last_output_time
           (Carpet::mglevel, Carpet::reflevel, variable));
      CCTK_REAL const my_out_dt
        = abs (out_dt - (-2)) <= dt_fudge ? IO_out_dt : out_dt;
      if (out_dt == 0)
      {
        should_output = true;
      }
      else if (abs (out_dt - (-1)) <= dt_fudge)
      {
        should_output = false;
      }
      else
      {
        should_output
          = cctkGH->cctk_time - dt_fudge >= last_output_time + my_out_dt;
      }
    }
    else
    {
      CCTK_WARN (1, "internal error");
      should_output = false;
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
    
    char * fullname = CCTK_FullName (variable);
    assert (fullname);
    
    int const ierr = OutputVarAs (cctkGH, fullname, out_filename);
    
    free (fullname);
    
    assert (Carpet::is_level_mode());
    extending_t extending (cctkGH);
    extending.set_last_output_iteration
      (Carpet::mglevel, Carpet::reflevel, variable, cctkGH->cctk_iteration);
    extending.set_last_output_time
      (Carpet::mglevel, Carpet::reflevel, variable, cctkGH->cctk_time);
    
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
    
    std::ostringstream filenamebuf;
    bool const use_IO_out_dir = strcmp (out_dir, "") == 0;
    filenamebuf << (use_IO_out_dir ? IO_out_dir : out_dir)
                << "/"
                << alias
                << out_extension;
    char const * const filename = filenamebuf.str().c_str();
    
    bool const did_truncate = extending.get_did_truncate (filename);
    bool const do_truncate
      = ! did_truncate and IO_TruncateOutputFiles (cctkGH);
    extending.set_did_truncate (filename);
    
    F5::file_t file (cctkGH, filenamebuf.str().c_str(), do_truncate);
    
    F5::timestep_t timestep (file, cctkGH->cctk_time);
    
    if (Carpet::is_meta_mode())
    {
      int const grouptype = CCTK_GroupTypeI (group);
      assert (grouptype >= 0);
      switch (grouptype)
      {
      case CCTK_GF:
        {
          for (Carpet::mglevel_iterator mglevel_iter (cctkGH);
               ! mglevel_iter.done();
               mglevel_iter.step())
          {
            write_one_mglevel (cctkGH, timestep, group, variable);
          }
        }
        break;
      case CCTK_ARRAY:
      case CCTK_SCALAR:
        {
          Carpet::enter_global_mode (cctkGH, 0);
          write_one_mglevel (cctkGH, timestep, group, variable);
          Carpet::leave_global_mode (cctkGH);
        }
        break;
      default:
        assert (0);
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
    std::ostringstream namebuf;
    namebuf << "convlevel " << cctkGH->cctk_convlevel;
    F5::simulation_t simulation (timestep, namebuf.str().c_str());
    
    if (grouptype == CCTK_GF and Carpet::is_global_mode())
    {
      int const grouptype = CCTK_GroupTypeI (group);
      assert (grouptype >= 0);
      switch (grouptype)
      {
      case CCTK_GF:
        {
          for (Carpet::reflevel_iterator reflevel_iter (cctkGH);
               ! reflevel_iter.done();
               reflevel_iter.step())
          {
            write_one_reflevel (cctkGH, simulation, group, variable);
          }
        }
        break;
      case CCTK_ARRAY:
      case CCTK_SCALAR:
        {
          Carpet::enter_level_mode (cctkGH, 0);
          write_one_reflevel (cctkGH, timestep, group, variable);
          Carpet::leave_level_mode (cctkGH);
        }
        break;
      default:
        assert (0);
      }
    }
    else
    {
      write_one_reflevel (cctkGH, simulation, group, variable);
    }
  }
  
  
  
  void
  write_one_reflevel (cGH const * const cctkGH,
                      F5::simulation_t & simulation,
                      int const group,
                      int const variable)
  {
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
      int const grouptype = CCTK_GroupTypeI (group);
      assert (grouptype >= 0);
      
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
    int const variable = tensor_component.get_variable();
    int const group = tensor_component.get_physical_quantity().get_group();
    
    dh * const dd = Carpet::arrdata.at(group).at(Carpet::map).dd;
    bbox<int, dim> const & region
      = (dd->boxes.at(Carpet::mglevel).at(Carpet::reflevel)
         .at(Carpet::component).exterior);
    
    F5::data_region_t data_region (tensor_component, region);
    
    void const * const varptr = CCTK_VarDataPtrI (cctkGH, 0, variable);
    assert (varptr != 0);
    int const vartype = CCTK_VarTypeI (variable);
    assert (vartype >= 0);
    data_region.write (varptr, vartype);
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
