// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.13 2003/10/14 16:39:16 schnetter Exp $

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "bbox.hh"
#include "data.hh"
#include "vect.hh"

#include "carpet.hh"

#include "interp.hh"

extern "C" {
  static char const * const rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.13 2003/10/14 16:39:16 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetInterp_interp_cc);
}



namespace CarpetInterp {
  
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_REAL,dim> rvect;
  
  typedef bbox<int,dim> ibbox;
  
  
  
#define ind_rc(rl,c) ind_rc_(rl,minrl,maxrl,c,maxncomps,hh)
  static inline int ind_rc_(int const rl, int const minrl, int const maxrl,
                            int const c, int const maxncomps,
                            gh<dim> const * const hh)
  {
    assert (rl>=minrl && rl<maxrl);
    assert (minrl>=0 && maxrl<=hh->reflevels());
    assert (c>=0 && c<maxncomps);
    assert (maxncomps>=0 && maxncomps<=hh->components(rl));
    int const ind = rl * maxncomps + c;
    assert (ind>=0 && ind < (maxrl-minrl) * maxncomps);
    return ind;
  }
  
#define ind_prc(p,rl,c) ind_prc_(p,nprocs,rl,minrl,maxrl,c,maxncomps,cgh,hh)
  static inline int ind_prc_(int const p, int const nprocs,
                             int const rl, int const minrl, int const maxrl,
                             int const c, int const maxncomps,
                             cGH const * const cgh,
                             gh<dim> const * const hh)
  {
    assert (p>=0 && p<nprocs);
    assert (nprocs==CCTK_nProcs(cgh));
    assert (rl>=minrl && rl<maxrl);
    assert (minrl>=0 && maxrl<=hh->reflevels());
    assert (c>=0 && c<maxncomps);
    assert (maxncomps>=0 && maxncomps<=hh->components(rl));
    int const ind = (p * (maxrl-minrl) + rl) * maxncomps + c;
    assert (ind>=0 && ind < nprocs * (maxrl-minrl) * maxncomps);
    return ind;
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
    
    
    
    // Find out about the coordinates
    const char * coord_system_name
      = CCTK_CoordSystemName (coord_system_handle);
    assert (coord_system_name);
    rvect lower, upper, delta;
    for (int d=0; d<dim; ++d) {
      ierr = CCTK_CoordRange
        (cgh, &lower[d], &upper[d], d+1, 0, coord_system_name);
      assert (!ierr);
      delta[d] = (upper[d] - lower[d]) / (hh->baseextent.shape()[d] - hh->baseextent.stride()[d]);
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
    int const ml = 0;
    int maxncomps = 0;
    for (int rl=minrl; rl<maxrl; ++rl) {
      maxncomps = max(maxncomps, hh->components(rl));
    }
    
    
    
    // Assert that all refinement levels have the same time
    // TODO: interpolate in time
    {
      bool can_interp = true;
      for (int rl=minrl; rl<maxrl; ++rl) {
        int const tl = 0;
        CCTK_REAL const time1 = tt->time (tl, rl, ml);
        CCTK_REAL const time2 = cgh->cctk_time / cgh->cctk_delta_time;
//         assert (fabs((time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time))) < 1e-12);
        can_interp = can_interp && (fabs((time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time))) < 1e-12);
      }
      if (! can_interp) {
        if (reflevel == -1) {
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Cannot interpolate in global mode at iteration %d (time %g), because not all refinement levels exist at this time.  Interpolation in time is not yet implemented.",
                      cgh->cctk_iteration, (double)cgh->cctk_time);
        } else {
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Cannot interpolate in level mode at iteration %d (time %g), because refinement level %d does not exist at this time.  Interpolation in time is not yet implemented.",
                      cgh->cctk_iteration, (double)cgh->cctk_time, reflevel);
        }
        assert (0);
      }
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
        pos[d] = static_cast<CCTK_REAL const *>(interp_coords[d])[n];
      }
      
      // Find the component that this grid point belongs to
      rlev[n] = -1;
      home[n] = -1;
      for (int rl=maxrl-1; rl>=minrl; --rl) {
        
        CCTK_REAL (* const rfloor) (CCTK_REAL const) = floor;
        ivect const ipos = ivect(map(rfloor, (pos - lower) / delta + 0.5));
        
        // TODO: use something faster than a linear search
        for (int c=0; c<hh->components(rl); ++c) {
          if (hh->extents[rl][c][ml].contains(ipos)) {
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
      ++ homecnts [ind_rc(rlev[n], home[n])];
      
    } // for n
    
    // Communicate counts
    vector<int> allhomecnts(nprocs * (maxrl-minrl) * maxncomps);
    MPI_Allgather (&homecnts   [0], (maxrl-minrl) * maxncomps, MPI_INT,
                   &allhomecnts[0], (maxrl-minrl) * maxncomps, MPI_INT, comm);
    
    
    
    // Create coordinate patches
    vector<data<CCTK_REAL,dim> > allcoords (nprocs * (maxrl-minrl) * maxncomps, -1);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          ivect lo (0);
          ivect up (1);
          up[0] = allhomecnts[ind_prc(p,rl,c)];
          up[1] = dim;
          ivect str (1);
          ibbox extent (lo, up-str, str);
          allcoords [ind_prc(p,rl,c)].allocate (extent, p);
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
        assert (tmpcnts[ind_rc(rl,c)] >= 0);
        assert (tmpcnts[ind_rc(rl,c)] < homecnts[ind_rc(rl,c)]);
        assert (dim==3);
        for (int d=0; d<dim; ++d) {
          allcoords[ind_prc(myproc,rl,c)][ivect(tmpcnts[ind_rc(rl,c)],d,0)]
            = static_cast<CCTK_REAL const *>(interp_coords[d])[n];
        }
        ++ tmpcnts[c + (rl-minrl)*maxncomps];
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          assert (tmpcnts[ind_rc(rl,c)] == homecnts[ind_rc(rl,c)]);
        }
      }
    }
    
    // Transfer coordinate patches
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          allcoords[ind_prc(p,rl,c)].change_processor (hh->processors[rl][c]);
        }
      }
    }
    
    
    
    // Create output patches
    vector<data<CCTK_REAL,dim> > alloutputs (nprocs * (maxrl-minrl) * maxncomps, -1);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          ivect lo (0);
          ivect up (1);
          up[0] = allhomecnts[ind_prc(p,rl,c)];
          up[1] = N_output_arrays;
          ivect str (1);
          ibbox extent (lo, up-str, str);
          alloutputs[ind_prc(p,rl,c)].allocate (extent, hh->processors[rl][c]);
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
      set_component ((cGH*)cgh, -1);
    }
    if (mglevel!=-1) {
      set_mglevel ((cGH*)cgh, -1);
    }
    if (reflevel!=-1) {
      set_reflevel ((cGH*)cgh, -1);
    }
    BEGIN_REFLEVEL_LOOP(cgh) {
      if (reflevel>=minrl && reflevel<maxrl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          
          BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
            
            // Find out about the local geometry
            ivect lsh;
            rvect coord_origin, coord_delta;
            for (int d=0; d<dim; ++d) {
              lsh[d] = cgh->cctk_lsh[d];
              coord_delta[d] = cgh->cctk_delta_space[d] / cgh->cctk_levfac[d];
              coord_origin[d] = cgh->cctk_origin_space[d] + (1.0 * cgh->cctk_levoff[d] / cgh->cctk_levoffdenom[d] + cgh->cctk_lbnd[d]) * coord_delta[d];
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
              assert (allcoords[ind_prc(p,reflevel,component)].owns_storage());
              assert (allhomecnts[ind_prc(p,reflevel,component)]
                      == allcoords[ind_prc(p,reflevel,component)].shape()[0]);
              assert (allhomecnts[ind_prc(p,reflevel,component)]
                      == alloutputs[ind_prc(p,reflevel,component)].shape()[0]);
              
              // Do the processor-local interpolation
              vector<const void *> tmp_interp_coords (dim);
              for (int d=0; d<dim; ++d) {
                tmp_interp_coords[d]
                  = &allcoords[ind_prc(p,reflevel,component)][ivect(0,d,0)];
              }
              vector<void *> tmp_output_arrays (N_output_arrays);
              for (int m=0; m<N_output_arrays; ++m) {
                assert (output_array_type_codes[m] == CCTK_VARIABLE_REAL);
                tmp_output_arrays[m]
                  = &alloutputs[ind_prc(p,reflevel,component)][ivect(0,m,0)];
              }
              
              ierr = CCTK_InterpLocalUniform
                (N_dims, local_interp_handle, param_table_handle,
                 &coord_origin[0], &coord_delta[0],
                 allhomecnts[ind_prc(p,reflevel,component)],
                 interp_coords_type_code, &tmp_interp_coords[0],
                 N_input_arrays, &lsh[0],
                 &input_array_type_codes[0], &input_arrays[0],
                 N_output_arrays,
                 output_array_type_codes, &tmp_output_arrays[0]);
              assert (!ierr);
              
            } // for processors
            
          } END_LOCAL_COMPONENT_LOOP;
        } END_MGLEVEL_LOOP;
      } // if reflevel active
    } END_REFLEVEL_LOOP;
    if (saved_reflevel!=-1) {
      set_reflevel ((cGH*)cgh, saved_reflevel);
    }
    if (saved_mglevel!=-1) {
      set_mglevel ((cGH*)cgh, saved_mglevel);
    }
    if (saved_component!=-1) {
      set_component ((cGH*)cgh, saved_component);
    }
    
    
    
    // Transfer output patches back
    for (int p=0; p<nprocs; ++p) {
      for (int rl=0; rl<minrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          alloutputs[ind_prc(p,rl,c)].change_processor (p);
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
          static_cast<CCTK_REAL *>(output_arrays[m])[n] =
            alloutputs[ind_prc(myproc,rl,c)][ivect(tmpcnts[ind_rc(rl,c)],m,0)];
        }
        ++ tmpcnts[ind_rc(rl,c)];
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<hh->components(rl); ++c) {
          assert (tmpcnts[ind_rc(rl,c)] == homecnts[ind_rc(rl,c)]);
        }
      }
    }
    
    
    
    // Done.
    return 0;
  }
  
} // namespace CarpetInterp
