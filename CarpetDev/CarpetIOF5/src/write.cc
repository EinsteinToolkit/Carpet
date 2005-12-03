#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"
  
#include "data_region.hh"
#include "file.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"
#include "write.hh"  


  
namespace CarpetIOF5 {
  
  namespace write {
    
    void
    write_meta (cGH const * const cctkGH,
                F5::file_t & file,
                int const group,
                int const variable)
    {
      DECLARE_CCTK_PARAMETERS;
      
      if (verbose or veryverbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "OutputVarAs/write_meta");
      }
      
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
    }
    
    
    
    void
    write_one_mglevel (cGH const * const cctkGH,
                       F5::timestep_t & timestep,
                       int const group,
                       int const variable)
    {
      DECLARE_CCTK_PARAMETERS;
      
      if (verbose or veryverbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "OutputVarAs/write_one_mglevel mglevel=%d",
                    Carpet::mglevel);
      }
      
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
      DECLARE_CCTK_PARAMETERS;
      
      if (verbose or veryverbose)
      {
        CCTK_INFO ("OutputVarAs/write_global");
      }
      
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
      DECLARE_CCTK_PARAMETERS;
      
      if (verbose or veryverbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "OutputVarAs/write_one_reflevel reflevel=%d",
                    Carpet::reflevel);
      }
      
      int const grouptype = CCTK_GroupTypeI (group);
      assert (grouptype >= 0);
      assert (grouptype == CCTK_GF);
      
      F5::mesh_refinement_topology_t topology
        (simulation, Carpet::reflevel, Carpet::maxreflevels,
         Carpet::spacereflevelfact, Carpet::maxspacereflevelfact);
      
      if (Carpet::is_level_mode())
      {
        for (Carpet::map_iterator map_iter (cctkGH, grouptype);
             ! map_iter.done();
             map_iter.step())
        {
          write_one_map (cctkGH, topology, group, variable);
        }
      }
      else
      {
        write_one_map (cctkGH, topology, group, variable);
      }
    }
    
    
    
    void
    write_one_map (cGH const * const cctkGH,
                   F5::topology_t & topology,
                   int const group,
                   int const variable)
    {
      DECLARE_CCTK_PARAMETERS;
      
      if (verbose or veryverbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "OutputVarAs/write_one_map map=%d", Carpet::map);
      }
      
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
      DECLARE_CCTK_PARAMETERS;
      
      if (verbose or veryverbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "OutputVarAs/write_one_component component=%d",
                    Carpet::component);
      }
      
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
    
  } // namespace write
  
} // namespace CarpetIOF5
