// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.21 2004/02/15 10:34:14 schnetter Exp $

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "bbox.hh"
#include "data.hh"
#include "defs.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "interp.hh"

extern "C" {
  static char const * const rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetInterp/src/interp.cc,v 1.21 2004/02/15 10:34:14 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetInterp_interp_cc);
}



namespace CarpetInterp {
  
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_REAL,dim> rvect;
  
  typedef bbox<int,dim> ibbox;
  
  
  
#define ind_rc(m,rl,c) ind_rc_(m,rl,minrl,maxrl,c,maxncomps,vhh)
  static inline int ind_rc_(int const m,
                            int const rl, int const minrl, int const maxrl,
                            int const c, int const maxncomps,
                            vector<gh<dim> *> const hh)
  {
    assert (m>=0 && m<maps);
    assert (rl>=minrl && rl<maxrl);
    assert (minrl>=0 && maxrl<=hh.at(m)->reflevels());
    assert (c>=0 && c<maxncomps);
    assert (maxncomps>=0 && maxncomps<=hh.at(m)->components(rl));
    int const ind = (rl-minrl) * maxncomps + c;
    assert (ind>=0 && ind < (maxrl-minrl) * maxncomps);
    return ind;
  }
  
#define ind_prc(p,m,rl,c) \
  ind_prc_(p,nprocs,m,rl,minrl,maxrl,c,maxncomps,cgh,vhh)
  static inline int ind_prc_(int const p, int const nprocs,
                             int const m,
                             int const rl, int const minrl, int const maxrl,
                             int const c, int const maxncomps,
                             cGH const * const cgh,
                             vector<gh<dim> *> const hh)
  {
    assert (p>=0 && p<nprocs);
    assert (nprocs==CCTK_nProcs(cgh));
    assert (m>=0 && m<maps);
    assert (rl>=minrl && rl<maxrl);
    assert (minrl>=0 && maxrl<=hh.at(m)->reflevels());
    assert (c>=0 && c<maxncomps);
    assert (maxncomps>=0 && maxncomps<=hh.at(m)->components(rl));
    int const ind = (p * (maxrl-minrl) + (rl-minrl)) * maxncomps + c;
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
    
    
    
    if (is_meta_mode()) {
      CCTK_WARN (0, "It is not possible to interpolate in meta mode");
    }
    
    
    
    // Get some information
    MPI_Comm const comm = CarpetMPIComm ();
    int const myproc = CCTK_MyProc (cgh);
    int const nprocs = CCTK_nProcs (cgh);
    assert (myproc>=0 && myproc<nprocs);
    
    // Multiple maps are not supported
    // (because we don't know how to select a map)
    assert (maps == 1);
    const int m = 0;
    
    int const minrl = reflevel==-1 ? 0                      : reflevel;
    int const maxrl = reflevel==-1 ? vhh.at(m)->reflevels() : reflevel+1;
    
    // Multiple convergence levels are not supported
    assert (mglevels == 1);
    int const ml = 0;
    
    int maxncomps = 0;
    for (int rl=minrl; rl<maxrl; ++rl) {
      maxncomps = max(maxncomps, vhh.at(m)->components(rl));
    }
    
    
    
    // Find the time interpolation order
    int partype;
    void const * const parptr
      = CCTK_ParameterGet ("prolongation_order_time", "Carpet", &partype);
    assert (parptr);
    assert (partype == PARAMETER_INTEGER);
    int const prolongation_order_time = * (CCTK_INT const *) parptr;
    
    
    
    // Find out about the coordinates
    const char * coord_system_name
      = CCTK_CoordSystemName (coord_system_handle);
    assert (coord_system_name);
    rvect lower, upper, delta;
    for (int d=0; d<dim; ++d) {
      ierr = CCTK_CoordRange
        (cgh, &lower[d], &upper[d], d+1, 0, coord_system_name);
      assert (!ierr);
      const ibbox & baseext = vhh.at(m)->baseextent;
      delta[d]
        = (upper[d] - lower[d]) / (baseext.upper()[d] - baseext.lower()[d]);
    }
    
    assert (N_interp_points >= 0);
    assert (interp_coords);
    for (int d=0; d<dim; ++d) {
      assert (N_interp_points==0 || interp_coords[d]);
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
      rlev.at(n) = -1;
      home.at(n) = -1;
      for (int rl=maxrl-1; rl>=minrl; --rl) {
        
        const int fact = maxreflevelfact * ipow(mgfact, basemglevel + mglevel) / ipow(reffact, rl);
        CCTK_REAL (* const rfloor) (CCTK_REAL const) = floor;
        ivect const ipos = ivect(::map(rfloor, (pos - lower) / (delta * fact) + 0.5)) * fact;
        assert (all(ipos % vhh.at(m)->bases.at(rl).at(ml).stride() == 0));
        
        // TODO: use something faster than a linear search
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          if (vhh.at(m)->extents.at(rl).at(c).at(ml).contains(ipos)) {
            rlev.at(n) = rl;
            home.at(n) = c;
            goto found;
          }
        }
      }
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Interpolation point #%d at [%g,%g,%g] is not on any grid patch",
                  n, pos[0], pos[1], pos[2]);
    found:
      assert (rlev.at(n)>=minrl && rlev.at(n)<maxrl);
      assert (home.at(n)>=0 && home.at(n)<vhh.at(m)->components(rlev.at(n)));
      ++ homecnts.at(ind_rc(m, rlev.at(n), home.at(n)));
      
    } // for n
    
    // Communicate counts
    vector<int> allhomecnts(nprocs * (maxrl-minrl) * maxncomps);
    MPI_Allgather
      (&homecnts   .at(0), (maxrl-minrl) * maxncomps, MPI_INT,
       &allhomecnts.at(0), (maxrl-minrl) * maxncomps, MPI_INT, comm);
    
    
    
    // Create coordinate patches
    vector<data<CCTK_REAL,dim> > allcoords
      (nprocs * (maxrl-minrl) * maxncomps);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          ivect lo (0);
          ivect up (1);
          up[0] = allhomecnts.at(ind_prc(p,m,rl,c));
          up[1] = dim;
          ivect str (1);
          ibbox extent (lo, up-str, str);
          allcoords.at(ind_prc(p,m,rl,c)).allocate (extent, p);
        }
      }
    }
    
    // Fill in local coordinate patches
    {
      vector<int> tmpcnts ((maxrl-minrl) * maxncomps);
      for (int n=0; n<N_interp_points; ++n) {
        int const rl = rlev.at(n);
        int const c = home.at(n);
        assert (rl>=minrl && rl<maxrl);
        assert (c>=0 && c<vhh.at(m)->components(rl));
        assert (tmpcnts.at(ind_rc(m,rl,c)) >= 0);
        assert (tmpcnts.at(ind_rc(m,rl,c)) < homecnts.at(ind_rc(m,rl,c)));
        assert (dim==3);
        for (int d=0; d<dim; ++d) {
          allcoords.at(ind_prc(myproc,m,rl,c))[ivect(tmpcnts.at(ind_rc(m,rl,c)),d,0)]
            = static_cast<CCTK_REAL const *>(interp_coords[d])[n];
        }
        ++ tmpcnts.at(c + (rl-minrl)*maxncomps);
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          assert (tmpcnts.at(ind_rc(m,rl,c)) == homecnts.at(ind_rc(m,rl,c)));
        }
      }
    }
    
    // Transfer coordinate patches
    for (comm_state<dim> state; !state.done(); state.step()) {
      for (int p=0; p<nprocs; ++p) {
        for (int rl=minrl; rl<maxrl; ++rl) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            allcoords.at(ind_prc(p,m,rl,c)).change_processor
              (state, vhh.at(m)->processors.at(rl).at(c));
          }
        }
      }
    }
    
    
    
    // Create output patches
    vector<data<CCTK_REAL,dim> > alloutputs
      (nprocs * (maxrl-minrl) * maxncomps, -1);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          ivect lo (0);
          ivect up (1);
          up[0] = allhomecnts.at(ind_prc(p,m,rl,c));
          up[1] = N_output_arrays;
          ivect str (1);
          ibbox extent (lo, up-str, str);
          alloutputs.at(ind_prc(p,m,rl,c)).allocate
            (extent, vhh.at(m)->processors.at(rl).at(c));
        }
      }
    }
    
    
    
    //
    // Do the local interpolation
    //
    int overall_ierr = 0;
    BEGIN_GLOBAL_MODE(cgh) {

      BEGIN_REFLEVEL_LOOP(cgh) {
        if (reflevel>=minrl && reflevel<maxrl) {
        
          // Number of necessary time levels
          int const tl = 0;
          CCTK_REAL const time1
            = vtt.at(m)->time (tl, reflevel, mglevel);
          CCTK_REAL const time2 = cgh->cctk_time / cgh->cctk_delta_time;
          bool const need_time_interp
            = fabs((time1 - time2) / (fabs(time1) + fabs(time2)
                                      + fabs(cgh->cctk_delta_time))) > 1e-12;
          int const num_tl
            = need_time_interp ? prolongation_order_time + 1 : 1;
        
          BEGIN_MAP_LOOP(cgh, CCTK_GF) {
            BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
            
              // Find out about the local geometry
              ivect lsh;
              rvect coord_origin, coord_delta;
              for (int d=0; d<dim; ++d) {
                lsh[d] = cgh->cctk_lsh[d];
                coord_delta[d] = cgh->cctk_delta_space[d] / cgh->cctk_levfac[d];
                coord_origin[d] = cgh->cctk_origin_space[d] + (cgh->cctk_lbnd[d] + 1.0 * cgh->cctk_levoff[d] / cgh->cctk_levoffdenom[d]) * coord_delta[d];
                
              }
            
            
            
              // Find out about the grid functions
              vector<CCTK_INT> input_array_type_codes(N_input_arrays * num_tl);
              vector<const void *> input_arrays(N_input_arrays * num_tl);
              for (int n=0; n<N_input_arrays; ++n) {
                if (input_array_variable_indices[n] == -1) {
                
                  // Ignore this entry
                  input_array_type_codes.at(tl+n*num_tl) = -1;
                  input_arrays.at(tl+n*num_tl) = 0;
                
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
                
                  assert (group.numtimelevels >= num_tl);
                
                  for (int tl=0; tl<num_tl; ++tl) {
                    input_array_type_codes.at(tl+n*num_tl) = group.vartype;
                    input_arrays.at(tl+n*num_tl) = CCTK_VarDataPtrI (cgh, tl, vi);
                  }
                
                }
              } // for input arrays
            
            
            
              // Work on the data from all processors
              for (int p=0; p<nprocs; ++p) {
                assert (allcoords.at(ind_prc(p,m,reflevel,component)).owns_storage());
                assert (allhomecnts.at(ind_prc(p,m,reflevel,component)) == allcoords.at(ind_prc(p,m,reflevel,component)).shape()[0]);
                assert (allhomecnts.at(ind_prc(p,m,reflevel,component)) == alloutputs.at(ind_prc(p,m,reflevel,component)).shape()[0]);
                
                int const npoints = allhomecnts.at(ind_prc(p,m,reflevel,component));
                
                // Do the processor-local interpolation
                vector<const void *> tmp_interp_coords (dim);
                for (int d=0; d<dim; ++d) {
                  tmp_interp_coords.at(d) = &allcoords.at(ind_prc(p,m,reflevel,component))[ivect(0,d,0)];
                }
                vector<void *> tmp_output_arrays (N_output_arrays * num_tl);
                if (need_time_interp) {
                  for (int j=0; j<N_output_arrays; ++j) {
                    assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
                    for (int tl=0; tl<num_tl; ++tl) {
                      tmp_output_arrays.at(tl+j*num_tl)
                        = new CCTK_REAL [npoints];
                    }
                  }
                } else {
                  for (int j=0; j<N_output_arrays; ++j) {
                    assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
                    tmp_output_arrays.at(j) = &alloutputs.at(ind_prc(p,m,reflevel,component))[ivect(0,j,0)];
                  }
                }
              
                ierr = CCTK_InterpLocalUniform
                  (N_dims, local_interp_handle, param_table_handle,
                   &coord_origin[0], &coord_delta[0],
                   npoints,
                   interp_coords_type_code, &tmp_interp_coords[0],
                   N_input_arrays * num_tl, &lsh[0],
                   &input_array_type_codes[0], &input_arrays[0],
                   N_output_arrays * num_tl,
                   output_array_type_codes, &tmp_output_arrays[0]);
                if (ierr) {
                  CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "The local interpolator returned the error code %d", ierr);
                }
                overall_ierr = min(overall_ierr, ierr);
              
                // Interpolate in time, if necessary
                if (need_time_interp) {
                  CCTK_REAL const time = cgh->cctk_time / cgh->cctk_delta_time;
                  vector<CCTK_REAL> times(num_tl);
                  for (int tl=0; tl<num_tl; ++tl) {
                    times.at(tl) = vtt.at(m)->time (tl, reflevel, mglevel);
                  }
                  vector<CCTK_REAL> tfacs(num_tl);
                  switch (num_tl) {
                  case 1:
                    // no interpolation
                    assert (fabs((time - times.at(0)) / fabs(time + times.at(0) + cgh->cctk_delta_time)) < 1e-12);
                    tfacs.at(0) = 1.0;
                    break;
                  case 2:
                    // linear (2-point) interpolation
                    tfacs.at(0) = (time - times.at(1)) / (times.at(0) - times.at(1));
                    tfacs.at(1) = (time - times.at(0)) / (times.at(1) - times.at(0));
                    break;
                  case 3:
                    // quadratic (3-point) interpolation
                    tfacs.at(0) = (time - times.at(1)) * (time - times.at(2)) / ((times.at(0) - times.at(1)) * (times.at(0) - times.at(2)));
                    tfacs.at(1) = (time - times.at(0)) * (time - times.at(2)) / ((times.at(1) - times.at(0)) * (times.at(1) - times.at(2)));
                    tfacs.at(2) = (time - times.at(0)) * (time - times.at(1)) / ((times.at(2) - times.at(0)) * (times.at(2) - times.at(1)));
                    break;
                  default:
                    assert (0);
                  }
                  for (int j=0; j<N_output_arrays; ++j) {
                    assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
                    for (int k=0; k<npoints; ++k) {
                      CCTK_REAL & dest = alloutputs.at(ind_prc(p,m,reflevel,component))[ivect(k,j,0)];
                      dest = 0;
                      for (int tl=0; tl<num_tl; ++tl) {
                        CCTK_REAL const src = ((CCTK_REAL const *)tmp_output_arrays.at(tl+j*num_tl))[k];
                        dest += tfacs[tl] * src;
                      }
                    }
                    for (int tl=0; tl<num_tl; ++tl) {
                      delete [] (CCTK_REAL *)tmp_output_arrays.at(tl+j*num_tl);
                    }
                  }
                }
              
              } // for processors
            
            } END_LOCAL_COMPONENT_LOOP;
          } END_MAP_LOOP;
        
        } // if reflevel active
      } END_REFLEVEL_LOOP;
      
    } END_GLOBAL_MODE;
    
    
    
    // Transfer output patches back
    for (comm_state<dim> state; !state.done(); state.step()) {
      for (int p=0; p<nprocs; ++p) {
        for (int rl=minrl; rl<maxrl; ++rl) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            alloutputs.at(ind_prc(p,m,rl,c)).change_processor (state, p);
          }
        }
      }
    }
    
    
    
    // Read out local output patches
    {
      vector<int> tmpcnts ((maxrl-minrl) * maxncomps);
      for (int n=0; n<N_interp_points; ++n) {
        int const rl = rlev.at(n);
        int const c = home.at(n);
        for (int j=0; j<N_output_arrays; ++j) {
          assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
          assert (alloutputs.at(ind_prc(myproc,m,rl,c)).owns_storage());
          static_cast<CCTK_REAL *>(output_arrays[j])[n] = alloutputs.at(ind_prc(myproc,m,rl,c))[ivect(tmpcnts.at(ind_rc(m,rl,c)),j,0)];
        }
        ++ tmpcnts.at(ind_rc(m,rl,c));
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          assert (tmpcnts.at(ind_rc(m,rl,c)) == homecnts.at(ind_rc(m,rl,c)));
        }
      }
    }
    
    
    
    int global_overall_ierr;
    MPI_Allreduce
      (&overall_ierr, &global_overall_ierr, 1, MPI_INT, MPI_MIN, comm);
    
    
    
    // Done.
    return global_overall_ierr;
  }
  
} // namespace CarpetInterp
