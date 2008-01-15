#ifndef WRITER_HH
#define WRITER_HH

#include "cctk.h"

#include "file.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"



namespace CarpetIOF5 {
  
  class writer_t {
    
    cGH const * const m_cctkGH;
    int const m_variable;
    
  public:
    
    writer_t (cGH const * cctkGH,
              int variable);
    
    void
    write (F5::file_t & file)
      const;
    
  private:
    
    void
    write_meta (F5::file_t & file,
                bool have_metafile)
      const;
    
    void
    write_one_mglevel (F5::timestep_t & timestep,
                       bool have_metafile)
      const;
    
    void
    write_global (F5::simulation_t & simulation,
                  bool have_metafile)
      const;
    
    void
    write_one_reflevel (F5::simulation_t & simulation,
                        bool have_metafile)
      const;
    
    void
    write_one_map (F5::simulation_t & simulation,
                   bool have_metafile)
      const;
    
    void
    write_one_component (F5::tensor_component_t & tensor_component,
                         bool have_metafile)
      const;
    
  };
  
} // namespace CarpetIOF5

#endif  // #ifndef WRITER_HH
