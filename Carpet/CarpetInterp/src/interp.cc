// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.5 2003/05/13 12:14:00 schnetter Exp $

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/data.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "interp.hh"

extern "C" {
  static char const * const rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.5 2003/05/13 12:14:00 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetInterp_interp_cc);
}



namespace CarpetInterp {
  
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_REAL,dim> rvect;
  
  typedef bbox<int,dim> ibbox;
  
  
  
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
    
    
    
    // Find out about the coordinates
    const char * coord_system_name
      = CCTK_CoordSystemName (coord_system_handle);
    assert (coord_system_name);
    rvect clower, cupper, cdelta;
    for (int d=0; d<dim; ++d) {
      ierr = CCTK_CoordRange
        (cgh, &clower[d], &cupper[d], d+1, 0, coord_system_name);
      assert (!ierr);
      cdelta[d] = cgh->cctk_delta_space[d];
    }
    
    assert (interp_coords);
    for (int d=0; d<dim; ++d) {
      assert (interp_coords[d]);
    }
    
    
    
    // Get some information
    MPI_Comm const comm = CarpetMPIComm ();
    int const myproc = CCTK_MyProc (cgh);
    int const nprocs = CCTK_nProcs (cgh);
    assert (myproc>=0 && myproc<nprocs);
    
    int const minrl = reflevel==-1 ? 0               : reflevel;
    int const maxrl = reflevel==-1 ? hh->reflevels() : reflevel+1;
    int maxncomps = 0;
    for (int rl=minrl; rl<maxrl; ++rl) {
      maxncomps = max(maxncomps, hh->components(rl));
    }
    
    
    
    // Assign interpolation points to components
    vector<int> rlev (N_interp_points); // refinement level of point n
    vector<int> home (N_interp_points); // component of point n
    vector<int> homecnts ((maxrl-minrl) * maxncomps); // number of points in component c
    
    for (int n=0; n<N_interp_points; ++n) {
      
      // Find the grid point closest to the interpolation point
      assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
      rvect pos;
      for (int d=0; d<dim; ++d) {
        pos[d] = ((CCTK_REAL const *)interp_coords[d])[n];
      }
      
      // Find the component that this grid point belongs to
      rlev[n] = -1;
      home[n] = -1;
      for (int rl=maxrl-1; rl>=minrl; --rl) {
        
        bbox<int,dim> const& baseext = dd->bases[rl][0].exterior;
        rvect const lower = clower + cdelta / maxreflevelfact * rvect(baseext.lower());
        rvect const delta = cdelta / maxreflevelfact;
        
        CCTK_REAL (* const rfloor) (CCTK_REAL const) = floor;
        int const stride = maxreflevelfact / reflevelfact;
        ivect const ipos = ivect(map(rfloor, (pos - lower) / delta / stride + 0.5)) * stride;
        
        // TODO: use something faster than a linear search
        for (int c=0; c<hh->components(rl); ++c) {
          if (hh->extents[reflevel][c][mglevel].contains(ipos)) {
            rlev[n] = rl;
            home[n] = c;
            goto found;
          }
        }
      }
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Interpolation point #%d at [%g,%g,%g] is not on any grid patch",
                  n, pos[0], pos[1], pos[2]);
    found:
      assert (rlev[n]>=minrl && rlev[n]<maxrl);
      assert (home[n]>=0 && home[n]<hh->components(rlev[n]));
      ++ homecnts [home[n] + (rlev[n]-minrl) * maxncomps];
      
    } // for n
    
    // Communicate counts
    vector<int> allhomecnts(nprocs * (maxrl-minrl) * maxncomps);
    MPI_Allgather (&homecnts   [0], (maxrl-minrl) * maxncomps, MPI_INT,
                   &allhomecnts[0], (maxrl-minrl) * maxncomps, MPI_INT, comm);
    
    
    
    // Create coordinate patches
    vector<data<CCTK_REAL,dim> > allcoords (nprocs * (maxrl-minrl) * maxncomps);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          ivect lo (0);
          ivect up (1);
          up[0] = allhomecnts[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p];
          up[1] = dim;
          ivect str (1);
          ibbox extent (lo, up-str, str);
          allcoords [c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].allocate (extent, p);
        }
      }
    }
    
    // Fill in local coordinate patches
    {
      vector<int> tmpcnts ((maxrl-minrl) * maxncomps);
      for (int n=0; n<N_interp_points; ++n) {
        int const rl = rlev[n];
        int const c = home[n];
        assert (rl>=minrl && rl<maxrl);
        assert (c>=0 && c<hh->components(rl));
        assert (tmpcnts[c + (rl-minrl)*maxncomps] >= 0
                && tmpcnts[c + (rl-minrl)*maxncomps] < homecnts[c + (rl-minrl)*maxncomps]);
        assert (dim==3);
        for (int d=0; d<dim; ++d) {
          allcoords[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*myproc][ivect(tmpcnts[c + (rl-minrl)*maxncomps],d,0)]
            = ((CCTK_REAL const *)interp_coords[d])[n];
        }
        ++ tmpcnts[c + (rl-minrl)*maxncomps];
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          assert (tmpcnts[c + (rl-minrl)*maxncomps] == homecnts[c + (rl-minrl)*maxncomps]);
        }
      }
    }
    
    // Transfer coordinate patches
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          allcoords[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].change_processor (hh->processors[rl][c]);
        }
      }
    }
    
    
    
    // Create output patches
    vector<data<CCTK_REAL,dim> > alloutputs (nprocs * (maxrl-minrl) * maxncomps);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          ivect lo (0);
          ivect up (1);
          up[0] = allhomecnts[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p];
          up[1] = N_output_arrays;
          ivect str (1);
          ibbox extent (lo, up-str, str);
          alloutputs[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].allocate (extent, hh->processors[rl][c]);
        }
      }
    }
    
    
    
    //
    // Do the local interpolation
    //
    int const saved_reflevel = reflevel;
    int const saved_mglevel = mglevel;
    int const saved_component = component;
    if (component!=-1) {
      set_component (cgh, -1);
    }
    if (mglevel!=-1) {
      set_mglevel (cgh, -1);
    }
    if (reflevel!=-1) {
      set_reflevel (cgh, -1);
    }
    BEGIN_REFLEVEL_LOOP(cgh) {
      if (reflevel>=minrl && reflevel<maxrl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          
          BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
            
            // Find out about the local geometry
            ivect lsh;
            rvect coord_origin, coord_delta;
            for (int d=0; d<dim; ++d) {
              lsh[d] = cgh->cctk_lsh[d];
              coord_origin[d] = cgh->cctk_origin_space[d];
              coord_delta[d] = cgh->cctk_delta_space[d] / cgh->cctk_levfac[d];
            }
            
            
            
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
            } // for input arrays
            
            
            
            // Work on the data from all processors
            for (int p=0; p<nprocs; ++p) {
              assert (allcoords[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].owns_storage());
              assert (allhomecnts[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p]
                      == allcoords[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].shape()[0]);
              assert (allhomecnts[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p]
                      == alloutputs[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].shape()[0]);
              
              // Do the processor-local interpolation
              vector<const void *> tmp_interp_coords (dim);
              for (int d=0; d<dim; ++d) {
                tmp_interp_coords[d]
                  = &allcoords[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p][ivect(0,d,0)];
              }
              vector<void *> tmp_output_arrays (N_output_arrays);
              for (int m=0; m<N_output_arrays; ++m) {
                assert (output_array_type_codes[m] == CCTK_VARIABLE_REAL);
                tmp_output_arrays[m]
                  = &alloutputs[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p][ivect(0,m,0)];
              }
              ierr = CCTK_InterpLocalUniform (N_dims,
                                              local_interp_handle,
                                              param_table_handle,
                                              &coord_origin[0],
                                              &coord_delta[0],
                                              allhomecnts[component + (reflevel-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p],
                                              interp_coords_type_code,
                                              &tmp_interp_coords[0],
                                              N_input_arrays,
                                              &lsh[0],
                                              &input_array_type_codes[0],
                                              &input_arrays[0],
                                              N_output_arrays,
                                              output_array_type_codes,
                                              &tmp_output_arrays[0]);
              
            } // for processors
            
          } END_LOCAL_COMPONENT_LOOP(cgh);
        } END_MGLEVEL_LOOP(cgh);
      }
    } END_REFLEVEL_LOOP(cgh);
    if (saved_reflevel!=-1) {
      set_reflevel (cgh, saved_reflevel);
    }
    if (saved_mglevel!=-1) {
      set_mglevel (cgh, saved_mglevel);
    }
    if (saved_component!=-1) {
      set_component (cgh, saved_component);
    }
    
    
    
    // Transfer output patches back
    for (int p=0; p<nprocs; ++p) {
      for (int rl=0; rl<minrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          alloutputs[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*p].change_processor (p);
        }
      }
    }
    
    // Read out local output patches
    {
      vector<int> tmpcnts ((maxrl-minrl) * maxncomps);
      for (int n=0; n<N_interp_points; ++n) {
        int const rl = rlev[n];
        int const c = home[n];
        for (int m=0; m<N_output_arrays; ++m) {
          assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
          assert (dim==3);
          ((CCTK_REAL *)output_arrays[m])[n] =
            alloutputs[c + (rl-minrl)*maxncomps + (maxrl-minrl)*maxncomps*myproc][ivect(tmpcnts[c + (rl-minrl)*maxncomps],m,0)];
        }
        ++ tmpcnts[c + (rl-minrl)*maxncomps];
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          assert (tmpcnts[c + (rl-minrl)*maxncomps] == homecnts[c + (rl-minrl)*maxncomps]);
        }
      }
    }
    
    
    
    // Done.
    return 0;
  }
  
} // namespace CarpetInterp
