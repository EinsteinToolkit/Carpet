#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "data_region.hh"
#include "file.hh"
#include "meta_data_region.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"
#include "writer.hh"


  
namespace CarpetIOF5 {
  
  writer_t::
  writer_t (cGH const * const cctkGH,
            int const variable)
    : m_cctkGH (cctkGH),
      m_variable (variable)
  {
  }
  
  
  
  void writer_t::
  write (F5::file_t & file)
    const
  {
    write_meta (file, file.get_have_metafile());
  }
  
  
  
  void writer_t::
  write_meta (F5::file_t & file,
              bool const have_metafile)
    const
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "OutputVarAs/write_meta");
    }
    
    F5::timestep_t timestep (file, m_cctkGH->cctk_time);
    
    if (Carpet::is_meta_mode())
    {
      for (Carpet::mglevel_iterator mglevel_iter (m_cctkGH);
           not mglevel_iter.done();
           mglevel_iter.step())
      {
        write_one_mglevel (timestep, have_metafile);
      }
    }
    else
    {
      write_one_mglevel (timestep, have_metafile);
    }
  }
  
  
  
  void writer_t::
  write_one_mglevel (F5::timestep_t & timestep,
                     bool const have_metafile)
    const
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "OutputVarAs/write_one_mglevel mglevel=%d",
                  Carpet::mglevel);
    }
    
    ostringstream namebuf;
    namebuf << "convlevel=" << m_cctkGH->cctk_convlevel;
    string const namestr = namebuf.str();
    char const * const name = namestr.c_str();
    F5::simulation_t simulation (timestep, name);
    
    int const grouptype = CCTK_GroupTypeFromVarI (m_variable);
    assert (grouptype >= 0);
    switch (grouptype)
    {
    case CCTK_ARRAY:
    case CCTK_SCALAR:
      {
        if (Carpet::do_global_mode)
        {
          write_global (simulation, have_metafile);
        }
      }
      break;
    case CCTK_GF:
      {
        if (Carpet::is_global_mode())
        {
          for (Carpet::reflevel_iterator reflevel_iter (m_cctkGH);
               not reflevel_iter.done();
               reflevel_iter.step())
          {
            write_one_reflevel (simulation, have_metafile);
          }
        }
        else
        {
          write_one_reflevel (simulation, have_metafile);
        }
      }
      break;
    default:
      assert (0);
    }
  }
  
  
  
  void writer_t::
  write_global (F5::simulation_t & simulation,
                bool const have_metafile)
    const
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose)
    {
      CCTK_INFO ("OutputVarAs/write_global");
    }
    
    F5::unigrid_topology_t topology (simulation);
    
    int const grouptype = CCTK_GroupTypeFromVarI (m_variable);
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
    
    int const group = CCTK_GroupIndexFromVarI (m_variable);
    assert (group >= 0 and group < CCTK_NumGroups());
    F5::physical_quantity_t physical_quantity (coordinate_system, group);
    
    F5::tensor_component_t tensor_component (physical_quantity, m_variable);
    
    int const map = 0;
    int const reflevel = 0;
    int const myproc = CCTK_MyProc (m_cctkGH);
    dh * const dd = Carpet::arrdata.at(group).at(map).dd;
    bbox<int, dim> const & region
      = dd->boxes.at(Carpet::mglevel).at(reflevel).at(myproc).exterior;
    
    if (have_metafile)
    {
      F5::meta_data_region_t meta_data_region (tensor_component, region);
      gh * const hh = Carpet::vhh.at(Carpet::map);
      int const proc = hh->processor (Carpet::reflevel, Carpet::component);
      meta_data_region.write (proc);
    }
    
    F5::data_region_t data_region (tensor_component, region);
    int const timelevel = 0;
    void const * const varptr
      = CCTK_VarDataPtrI (m_cctkGH, timelevel, m_variable);
    assert (varptr != 0);
    int const vartype = CCTK_VarTypeI (m_variable);
    assert (vartype >= 0);
    data_region.write (varptr, vartype);
  }
  
  
  
  void writer_t::
  write_one_reflevel (F5::simulation_t & simulation,
                      bool const have_metafile)
    const
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "OutputVarAs/write_one_reflevel reflevel=%d",
                  Carpet::reflevel);
    }
    
    int const grouptype = CCTK_GroupTypeFromVarI (m_variable);
    assert (grouptype >= 0);
    assert (grouptype == CCTK_GF);
    
    if (Carpet::is_level_mode())
    {
      for (Carpet::map_iterator map_iter (m_cctkGH, grouptype);
           not map_iter.done();
           map_iter.step())
      {
        write_one_map (simulation, have_metafile);
      }
    }
    else
    {
      write_one_map (simulation, have_metafile);
    }
  }
  
  
  
  void writer_t::
  write_one_map (F5::simulation_t & simulation,
                 bool const have_metafile)
    const
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "OutputVarAs/write_one_map map=%d", Carpet::map);
    }
    
    F5::mesh_refinement_topology_t topology
      (simulation, Carpet::map, Carpet::reflevel, Carpet::maxreflevels,
       Carpet::spacereflevelfact, Carpet::maxspacereflevelfact);
    
    vect<CCTK_REAL, dim> level_origin, level_delta;
    for (int d=0; d<dim; ++d)
    {
      cGH const * const cctkGH = m_cctkGH;
      DECLARE_CCTK_ARGUMENTS;
      level_origin[d] = CCTK_ORIGIN_SPACE(d);
      level_delta[d]  = CCTK_DELTA_SPACE(d);
    }
    F5::Cartesian_coordinate_system_t coordinate_system
      (topology, level_origin, level_delta);
    
    int const group = CCTK_GroupIndexFromVarI (m_variable);
    assert (group >= 0 and group < CCTK_NumGroups());
    F5::physical_quantity_t physical_quantity (coordinate_system, group);
    
    F5::tensor_component_t tensor_component (physical_quantity, m_variable);
    
    if (Carpet::is_singlemap_mode())
    {
      int const grouptype = CCTK_GroupTypeI (group);
      assert (grouptype >= 0);
      
      for (Carpet::component_iterator component_iter (m_cctkGH, grouptype);
           not component_iter.done();
           component_iter.step())
      {
        write_one_component (tensor_component, have_metafile);
      }
    }
    else
    {
      write_one_component (tensor_component, have_metafile);
    }
  }
  
  
  
  void writer_t::
  write_one_component (F5::tensor_component_t & tensor_component,
                       bool const have_metafile)
    const
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "OutputVarAs/write_one_component component=%d",
                  Carpet::component);
    }
    
    gh * const hh = Carpet::vhh.at(Carpet::map);
    bool const is_local = hh->is_local (Carpet::reflevel, Carpet::component);
    if (have_metafile or is_local)
    {
      dh * const dd = Carpet::vdd.at(Carpet::map);
      bbox<int, dim> const & region
        = (dd->boxes.at(Carpet::mglevel).at(Carpet::reflevel)
           .at(Carpet::component).exterior);
      
      if (have_metafile)
      {
        F5::meta_data_region_t meta_data_region (tensor_component, region);
        int const proc = hh->processor (Carpet::reflevel, Carpet::component);
        meta_data_region.write (proc);
      }
      
      if (is_local)
      {
        F5::data_region_t data_region (tensor_component, region);
        int const timelevel = 0;
        void const * const varptr
          = CCTK_VarDataPtrI (m_cctkGH, timelevel, m_variable);
        assert (varptr != 0);
        int const vartype = CCTK_VarTypeI (m_variable);
        assert (vartype >= 0);
        data_region.write (varptr, vartype);
      }
    }
  }
  
} // namespace CarpetIOF5
