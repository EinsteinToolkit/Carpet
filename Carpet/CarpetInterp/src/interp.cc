// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.3 2003/05/08 15:35:49 schnetter Exp $

#include <assert.h>
#include <math.h>

#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/data.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "interp.hh"

extern "C" {
  static char const * const rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.3 2003/05/08 15:35:49 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetInterp_interp_cc);
}



namespace CarpetInterp {
  
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_REAL,dim> rvect;
  
  typedef bbox<int,dim> ibbox;
  
  
  
  ivect rvect2ivect (rvect const & pos,
                     rvect const & lower,
                     rvect const finest_delta)
  {
    static CCTK_REAL (* const rfloor) (CCTK_REAL const) = floor;
    int const stride = maxreflevelfact / reflevelfact;
    return (ivect(map(rfloor, (pos - lower) / finest_delta / stride + 0.5))
            * stride);
  }
  
  
  
  void CarpetInterpStartup ()
  {
    CCTK_OverloadInterpGridArrays (InterpGridArrays);
  }
  
  
  
  int InterpGridArrays (cGH const * const cgh,
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
    
    assert (cgh);
    assert (N_dims == dim);
    
    
    
    // We want to be in level mode
    assert (reflevel != -1);
    if (hh->local_components(reflevel) != 1 && component != -1) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot interpolate in local mode");
    }
    if (hh->local_components(reflevel) != 1) assert (component == -1);
    
    
    
    // Find out about the coordinates
    const char * coord_system_name
      = CCTK_CoordSystemName (coord_system_handle);
    assert (coord_system_name);
    rvect lower, upper;
    for (int d=0; d<dim; ++d) {
//       ierr = CCTK_CoordRange
//         (cgh, &lower[d], &upper[d], d+1, 0, coord_system_name);
//       assert (!ierr);
      lower[d] = cgh->cctk_origin_space[d];
      upper[d] = cgh->cctk_origin_space[d] + (cgh->cctk_gsh[d] - 1) * cgh->cctk_delta_space[d] / cgh->cctk_levfac[d];
    }
    
    const ivect gsh   (cgh->cctk_gsh );
    const rvect delta ((upper - lower) / rvect((gsh - 1) * maxreflevelfact));
    
    assert (interp_coords);
    for (int d=0; d<dim; ++d) {
      assert (interp_coords[d]);
    }
    
    
    
    // Get some information
    MPI_Comm const comm = CarpetMPIComm ();
    int const myproc = CCTK_MyProc (cgh);
    int const nprocs = CCTK_nProcs (cgh);
    assert (myproc>=0 && myproc<nprocs);
    
    int const ncomps = hh->components(reflevel);
    
    
    
    // Assign interpolation points to components
    vector<int> home (N_interp_points); // component of point n
    vector<int> homecnts (ncomps); // number of points in component c
    
    for (int n=0; n<N_interp_points; ++n) {
      
      // Find the grid point closest to the interpolation point
      assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
      rvect pos;
      for (int d=0; d<dim; ++d) {
        pos[d] = ((CCTK_REAL const *)interp_coords[d])[n];
      }
      ivect const ipos = rvect2ivect (pos, lower, delta);
      
      // Find the component that this grid point belongs to
      home[n] = -1;
      // TODO: use something faster than a linear search
      for (int c=0; c<ncomps; ++c) {
        if (hh->extents[reflevel][c][mglevel].contains(ipos)) {
          home[n] = c;
          break;
        }
      }
      if (! (home[n]>=0 && home[n]<ncomps)) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Interpolation point #%d at real [%g,%g,%g] is at integer [%d,%d,%d], which is not on any grid patch",
                    n, pos[0], pos[1], pos[2], ipos[0], ipos[1], ipos[2]);
      }
      assert (home[n]>=0 && home[n]<ncomps);
      ++ homecnts [home[n]];
      
    } // for n
    
    // Communicate counts
    vector<int> allhomecnts(nprocs * ncomps);
    MPI_Allgather (&homecnts   [0], ncomps, MPI_INT,
                   &allhomecnts[0], ncomps, MPI_INT, comm);
    
    
    
    // Create coordinate patches
    vector<data<CCTK_REAL,dim> > allcoords (nprocs * ncomps);
    for (int p=0; p<nprocs; ++p) {
      for (int c=0; c<ncomps; ++c) {
        ivect lo (0);
        ivect up (1);
        up[0] = allhomecnts[c+ncomps*p];
        up[1] = dim;
        ivect str (1);
        ibbox extent (lo, up-str, str);
        allcoords [c+ncomps*p].allocate (extent, p);
      }
    }
    
    // Fill in local coordinate patches
    {
      vector<int> tmpcnts (ncomps);
      for (int n=0; n<N_interp_points; ++n) {
        int const c = home[n];
        assert (c>=0 && c<ncomps);
        assert (tmpcnts[c]>=0 && tmpcnts[c]<homecnts[c]);
        assert (dim==3);
        for (int d=0; d<dim; ++d) {
          allcoords[c+ncomps*myproc][ivect(tmpcnts[c],d,0)]
            = ((CCTK_REAL const *)interp_coords[d])[n];
        }
        ++ tmpcnts[c];
      }
      for (int c=0; c<ncomps; ++c) {
        assert (tmpcnts[c] == homecnts[c]);
      }
    }
    
    // Transfer coordinate patches
    for (int p=0; p<nprocs; ++p) {
      for (int c=0; c<ncomps; ++c) {
        allcoords[c+ncomps*p].change_processor (hh->processors[reflevel][c]);
      }
    }
    
    
    
    // Create output patches
    vector<data<CCTK_REAL,dim> > alloutputs (nprocs * ncomps);
    for (int p=0; p<nprocs; ++p) {
      for (int c=0; c<ncomps; ++c) {
        ivect lo (0);
        ivect up (1);
        up[0] = allhomecnts[c+ncomps*p];
        up[1] = N_output_arrays;
        ivect str (1);
        ibbox extent (lo, up-str, str);
        alloutputs[c+ncomps*p].allocate (extent, hh->processors[reflevel][c]);
      }
    }
    
    
    
    //
    // Do the local interpolation
    //
    BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
      
      // Find out about the local geometry
      const ivect lbnd (cgh->cctk_lbnd);
      const ivect lsh  (cgh->cctk_lsh );
      
      const rvect coord_delta  = delta * (maxreflevelfact / reflevelfact);
      const rvect coord_origin = lower + rvect(lbnd) * coord_delta;
      
      
      
      // Find out about the grid functions
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
          input_arrays[n] = CCTK_VarDataPtrI (cgh, 0, vi);
          
        }
      }
      
      
      
      // Work on the data from all processors
      for (int p=0; p<nprocs; ++p) {
        assert (allcoords[component+ncomps*p].owns_storage());
        assert (allhomecnts[component+ncomps*p]
                == allcoords[component+ncomps*p].shape()[0]);
        assert (allhomecnts[component+ncomps*p]
                == alloutputs[component+ncomps*p].shape()[0]);
        
        // Do the processor-local interpolation
        vector<const void *> tmp_interp_coords (dim);
        for (int d=0; d<dim; ++d) {
          tmp_interp_coords[d]
            = &allcoords[component+ncomps*p][ivect(0,d,0)];
        }
        vector<void *> tmp_output_arrays (N_output_arrays);
        for (int m=0; m<N_output_arrays; ++m) {
          assert (output_array_type_codes[m] == CCTK_VARIABLE_REAL);
          tmp_output_arrays[m]
            = &alloutputs[component+ncomps*p][ivect(0,m,0)];
        }
        ierr = CCTK_InterpLocalUniform (N_dims,
                                        local_interp_handle,
                                        param_table_handle,
                                        &coord_origin[0],
                                        &coord_delta[0],
                                        allhomecnts[component+ncomps*p],
                                        interp_coords_type_code,
                                        &tmp_interp_coords[0],
                                        N_input_arrays,
                                        &lsh[0],
                                        &input_array_type_codes[0],
                                        &input_arrays[0],
                                        N_output_arrays,
                                        output_array_type_codes,
                                        &tmp_output_arrays[0]);
        
      } // for p
      
    } END_LOCAL_COMPONENT_LOOP(cgh);
    
    
    
    // Transfer output patches back
    for (int p=0; p<nprocs; ++p) {
      for (int c=0; c<ncomps; ++c) {
        alloutputs[c+ncomps*p].change_processor (p);
      }
    }
    
    // Read out local output patches
    {
      vector<int> tmpcnts (ncomps);
      for (int n=0; n<N_interp_points; ++n) {
        int const c = home[n];
        for (int m=0; m<N_output_arrays; ++m) {
          assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
          assert (dim==3);
          ((CCTK_REAL *)output_arrays[m])[n] =
            alloutputs[c+ncomps*myproc][ivect(tmpcnts[c],m,0)];
        }
        ++ tmpcnts[c];
      }
      for (int c=0; c<ncomps; ++c) {
        assert (tmpcnts[c] == homecnts[c]);
      }
    }
    
    
    
    // Done.
    return 0;
  }
  
} // namespace CarpetInterp
