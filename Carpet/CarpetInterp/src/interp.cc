// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.1 2003/04/30 12:37:31 schnetter Exp $

#include <assert.h>

#include <complex>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/vect.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "interp.hh"

extern "C" {
  static char const * const rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.1 2003/04/30 12:37:31 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetInterp_interp_cc);
}



namespace CarpetInterp {
  
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_REAL,dim> rvect;
  
  
  
  void CarpetInterpStartup ()
  {
    CCTK_OverloadInterpGridArrays (InterpGridArrays);
  }
  
  
  
  int InterpGridArrays (cGH const * const cGH,
                        int const N_dims,
                        int const local_interp_handle,
                        int const param_table_handle,
                        int const coord_system_handle,
                        int const N_interp_points,
                        int const interp_coords_type_code,
                        void const * const interp_coords [],
                        int const N_input_arrays,
                        CCTK_INT const input_array_variable_indices [],
                        int const N_output_arrays,
                        CCTK_INT const output_array_type_codes [],
                        void * const output_arrays [])
  {
    int ierr;
    
    assert (cGH);
    assert (N_dims == dim);
    
    const char * coord_system_name
      = CCTK_CoordSystemName (coord_system_handle);
    assert (coord_system_name);
    rvect lower, upper;
    for (int d=0; d<dim; ++d) {
      ierr = CCTK_CoordRange
        (cGH, &lower[d], &upper[d], d+1, 0, coord_system_name);
      assert (!ierr);
    }
    
    const ivect lbnd (cGH->cctk_lbnd);
    const ivect lsh  (cGH->cctk_lsh );
    const ivect gsh  (cGH->cctk_gsh );
    
    const rvect coord_delta
      = (upper - lower) / rvect((gsh - 1) * reflevelfact);
    const rvect coord_origin = lower + rvect(lbnd) * coord_delta;
    
    vector<CCTK_INT> input_array_type_codes(N_input_arrays);
    vector<const void *> input_arrays(N_input_arrays);
    for (int n=0; n<N_input_arrays; ++n) {
      if (input_array_variable_indices[n] == -1) {
        
        // Ignore this entry
        input_array_type_codes[n] = -1;
        input_arrays[n] = 0;
        
      } else {
        
        int const vi = input_array_variable_indices[n];
        assert (vi>=0 && vi<CCTK_NumVars());
        
        int const gi = CCTK_GroupIndexFromVarI (vi);
        assert (gi>=0 && gi<CCTK_NumGroups());
        
        cGroup group;
        ierr = CCTK_GroupData (gi, &group);
        assert (!ierr);
        assert (group.grouptype == CCTK_GF);
        assert (group.dim == dim);
        assert (group.disttype == CCTK_DISTRIB_DEFAULT);
        assert (group.stagtype == 0); // not staggered
        
        input_array_type_codes[n] = group.vartype;
        input_arrays[n] = CCTK_VarDataPtrI (cGH, 0, vi);
        
      }
    }
    
    ierr = CCTK_InterpLocalUniform (N_dims,
                                    local_interp_handle,
                                    param_table_handle,
                                    &coord_origin[0],
                                    &coord_delta[0],
                                    N_interp_points,
                                    interp_coords_type_code,
                                    interp_coords,
                                    N_input_arrays,
                                    &lsh[0],
                                    &input_array_type_codes[0],
                                    &input_arrays[0],
                                    N_output_arrays,
                                    output_array_type_codes,
                                    output_arrays);
    
    return ierr;
  }
  
} // namespace CarpetInterp
