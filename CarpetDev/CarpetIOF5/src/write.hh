#ifndef WRITE_HH
#define WRITE_HH

#include "cctk.h"

#include "file.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"



namespace CarpetIOF5 {
  
  namespace write {
    
    void
    write_meta (cGH const * cctkGH,
                F5::file_t & file,
                int group,
                int variable);
    
    void
    write_one_mglevel (cGH const * cctkGH,
                       F5::timestep_t & timestep,
                       int group,
                       int variable);
    
    void
    write_global (cGH const * cctkGH,
                  F5::simulation_t & simulation,
                  int group,
                  int variable);
    
    void
    write_one_reflevel (cGH const * cctkGH,
                        F5::simulation_t & simulation,
                        int group,
                        int variable);
    
    void
    write_one_map (cGH const * cctkGH,
                   F5::topology_t & topology,
                   int group,
                   int variable);
    
    void
    write_one_component (cGH const * cctkGH,
                         F5::tensor_component_t & tensor_component);
    
  } // namespace write
  
} // namespace CarpetIOF5

#endif  // #ifndef WRITE_HH
