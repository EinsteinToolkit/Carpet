#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "coordinate_system.hh"
#include "utils.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    coordinate_system_t::
    coordinate_system_t (topology_t & topology)
      : m_topology (topology)
    {
    }
    
    
    
    coordinate_system_t::
    ~ coordinate_system_t ()
    {
    }
    
    
    
    topology_t & coordinate_system_t::
    get_topology ()
      const
    {
      return m_topology;
    }
    
    
    
    hid_t coordinate_system_t::
    get_hdf5_coordinate_system ()
      const
    {
      return m_hdf5_coordinate_system;
    }
    
    
    
    bool coordinate_system_t::
    invariant ()
      const
    {
      return true;
    }
    
    
    
    Cartesian_coordinate_system_t::
    Cartesian_coordinate_system_t (topology_t & topology,
                                   vect<CCTK_REAL, dim> const & level_origin,
                                   vect<CCTK_REAL, dim> const & level_delta)
      : coordinate_system_t (topology),
        m_level_origin (level_origin),
        m_level_delta (level_delta)
    {
      assert (all (m_level_delta > 0));
      
      init();
      
      assert (invariant());
    }
    
    
    
    Cartesian_coordinate_system_t::
    Cartesian_coordinate_system_t (topology_t & topology,
                                   vect<CCTK_REAL, dim> const & coarse_origin,
                                   vect<CCTK_REAL, dim> const & coarse_delta,
                                   vect<int, dim> const & level_offset,
                                   vect<int, dim> const & level_offset_denominator)
      : coordinate_system_t (topology)
    {
      assert (all (coarse_delta > 0));
      assert (all (level_offset_denominator > 0));
      
      mesh_refinement_topology_t * mesh_refinement_topology
        = dynamic_cast<mesh_refinement_topology_t *> (& topology);
      assert (mesh_refinement_topology != 0);
      mesh_refinement_topology
        ->calculate_level_origin_delta (coarse_origin, coarse_delta,
                                        level_offset, level_offset_denominator,
                                        m_level_origin, m_level_delta);
      
      init();
      
      assert (invariant());
    }
    
    
    
    Cartesian_coordinate_system_t::
    ~ Cartesian_coordinate_system_t ()
    {
      herr_t const herr = H5Gclose (m_hdf5_coordinate_system);
      assert (! herr);
    }
    
    
    
    void Cartesian_coordinate_system_t::
    init ()
    {
      assert (all (m_level_delta > 0));
      
      ostringstream namebuf;
      namebuf << "Cartesian 3D, x0=" << m_level_origin
              << ", dx=" << m_level_delta;
      string const namestr = namebuf.str();
      char const * const name = namestr.c_str();
      
      m_hdf5_coordinate_system
        = open_or_create_group (m_topology.get_hdf5_topology(), name);
      assert (m_hdf5_coordinate_system >= 0);
      
      write_or_check_attribute
        (m_hdf5_coordinate_system, "origin", m_level_origin);
      write_or_check_attribute
        (m_hdf5_coordinate_system, "delta", m_level_delta);
    }
    
    
    
    bool Cartesian_coordinate_system_t::
    invariant ()
      const
    {
      return (coordinate_system_t::invariant()
              and all (m_level_delta > 0));
    }
    
  } // namespace F5
  
} // namespace CarpetIOF5
