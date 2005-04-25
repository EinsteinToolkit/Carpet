#ifndef TOPOLOGY_HH
#define TOPOLOGY_HH

#include <hdf5.h>

#include "cctk.h"

#include "defs.hh"
#include "vect.hh"

#include "simulation.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    class topology_t {
      
    protected:
      
      simulation_t & m_simulation;
      
      hid_t m_hdf5_topology;
      
      topology_t (simulation_t & simulation);
      
    public:
      
      virtual
      ~ topology_t ();
      
    public:
      
      hid_t
      get_hdf5_topology ()
        const;
      
      virtual bool
      invariant ()
        const;
    };
    
    
    
    class unigrid_topology_t : public topology_t {
      
    public:
      
      // Create unigrid topology
      unigrid_topology_t (simulation_t & simulation);
      
      virtual
      ~ unigrid_topology_t ();
      
      virtual bool
      invariant ()
        const;
    };
  
  
  
    class mesh_refinement_topology_t : public topology_t {
      
      int const m_refinement_level;
      int const m_max_refinement_levels;
      int const m_level_refinement_factor;
      int const m_max_refinement_factor;
      
    public:
      
      mesh_refinement_topology_t (simulation_t & simulation,
                                  int refinement_level,
                                  int max_refinement_levels,
                                  int level_refinement_factor,
                                  int max_refinement_factor);
      
      void
      calculate_level_origin_delta (vect<CCTK_REAL, dim> const & coarse_origin,
                                    vect<CCTK_REAL, dim> const & coarse_delta,
                                    vect<int, dim> const & level_offset,
                                    vect<int, dim> const & level_offset_denominator,
                                    vect<CCTK_REAL, dim> & level_origin,
                                    vect<CCTK_REAL, dim> & level_delta)
        const;
      
      virtual
      ~ mesh_refinement_topology_t ();
      
      virtual bool
      invariant ()
        const;
    };
    
  } // namespace F5
  
} // namespace CarpetIOF5

#endif // #ifndef TOPOLOGY_HH
