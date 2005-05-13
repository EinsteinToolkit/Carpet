#ifndef COORDINATE_SYSTEM_HH
#define COORDINATE_SYSTEM_HH

#include <hdf5.h>

#include "cctk.h"

#include "defs.hh"
#include "vect.hh"

#include "topology.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    class coordinate_system_t {
      
    protected:
      
      topology_t & m_topology;
      
      hid_t m_hdf5_coordinate_system;
      
      coordinate_system_t (topology_t & topology);
      
    public:
      
      virtual
      ~ coordinate_system_t ();
      
      topology_t &
      get_topology ()
        const;
      
      hid_t
      get_hdf5_coordinate_system ()
        const;
      
      virtual bool
      invariant ()
        const;
      
    };
    
    
    
    class Cartesian_coordinate_system_t : public coordinate_system_t {
      
      vect<CCTK_REAL, dim> m_level_origin;
      vect<CCTK_REAL, dim> m_level_delta;
      
    public:
      
      Cartesian_coordinate_system_t (topology_t & topology,
                                     vect<CCTK_REAL, dim> const & level_origin,
                                     vect<CCTK_REAL, dim> const & level_delta);
      
      Cartesian_coordinate_system_t (topology_t & topology,
                                     vect<CCTK_REAL, dim> const & coarse_origin,
                                     vect<CCTK_REAL, dim> const & coarse_delta,
                                     vect<int, dim> const & level_offset,
                                     vect<int, dim> const & level_offset_denominator);
      
    private:
      
      void
      init ();
      
    public:
      
      virtual
      ~ Cartesian_coordinate_system_t ();
      
      virtual bool
      invariant ()
        const;
      
    };

  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef COORDINATE_SYSTEM_HH
