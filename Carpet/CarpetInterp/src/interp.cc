#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "bbox.hh"
#include "data.hh"
#include "defs.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "interp.hh"



namespace CarpetInterp {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_REAL,dim> rvect;
  
  typedef bbox<int,dim> ibbox;
  
  
  
#define ind_rc(m,rl,c) ind_rc_(m,rl,minrl,maxrl,c,maxncomps,vhh)
  static inline int ind_rc_(int const m,
                            int const rl, int const minrl, int const maxrl,
                            int const c, int const maxncomps,
                            vector<gh*> const hh)
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
                             vector<gh*> const hh)
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
    if (CCTK_IsFunctionAliased ("SymmetryInterpolate")) {
      return SymmetryInterpolate
        (cgh, N_dims,
         local_interp_handle, param_table_handle, coord_system_handle,
         N_interp_points, interp_coords_type_code, interp_coords,
         N_input_arrays, input_array_variable_indices,
         N_output_arrays, output_array_type_codes, output_arrays);
    } else {
      return Carpet_DriverInterpolate
        (cgh, N_dims,
         local_interp_handle, param_table_handle, coord_system_handle,
         N_interp_points, interp_coords_type_code, interp_coords,
         N_input_arrays, input_array_variable_indices,
         N_output_arrays, output_array_type_codes, output_arrays);
    }
  }
  
  
  
  extern "C" CCTK_INT
  Carpet_DriverInterpolate (CCTK_POINTER_TO_CONST const cgh_,
                            CCTK_INT const N_dims,
                            CCTK_INT const local_interp_handle,
                            CCTK_INT const param_table_handle,
                            CCTK_INT const coord_system_handle,
                            CCTK_INT const N_interp_points,
                            CCTK_INT const interp_coords_type_code,
                            CCTK_POINTER_TO_CONST const interp_coords [],
                            CCTK_INT const N_input_arrays,
                            CCTK_INT const input_array_variable_indices [],
                            CCTK_INT const N_output_arrays,
                            CCTK_INT const output_array_type_codes [],
                            CCTK_POINTER const output_arrays [])
  {
    cGH const * const cgh = static_cast<cGH const *> (cgh_);
    int ierr;
    
    assert (cgh);
    assert (N_dims == dim);
    
    
    
    if (is_meta_mode()) {
      CCTK_WARN (0, "It is not possible to interpolate in meta mode");
    }
    
    bool const want_global_mode = is_global_mode();
    
    
    
    // Get some information
    MPI_Comm const comm = CarpetMPIComm ();
    int const myproc = CCTK_MyProc (cgh);
    int const nprocs = CCTK_nProcs (cgh);
    assert (myproc>=0 && myproc<nprocs);
    
    // Multiple maps are not supported
    // (because we don't know how to select a map)
    assert (maps == 1);
    const int m = 0;
    
    int const minrl = want_global_mode ? 0                      : reflevel;
    int const maxrl = want_global_mode ? vhh.at(m)->reflevels() : reflevel+1;
    
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
    
    CCTK_REAL const current_time = cgh->cctk_time / cgh->cctk_delta_time;
    
    vector<CCTK_REAL> interpolation_times (N_interp_points);
    bool interpolation_times_differ;
    ierr = Util_TableGetRealArray
      (param_table_handle, N_interp_points,
       &interpolation_times.front(), "interpolation_times");
    if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      for (int n=0; n<N_interp_points; ++n) {
        interpolation_times.at(n) = current_time;
      }
      interpolation_times_differ = false;
    } else {
      assert (ierr == N_interp_points);
      interpolation_times_differ = true;
    }
    
    vector<CCTK_INT> time_deriv_order (N_interp_points);
    bool have_time_derivs;
    ierr = Util_TableGetIntArray
      (param_table_handle, N_interp_points,
       &time_deriv_order.front(), "time_deriv_order");
    if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      for (int n=0; n<N_interp_points; ++n) {
        time_deriv_order.at(n) = 0;
      }
      have_time_derivs = false;
    } else {
      assert (ierr == N_interp_points);
      have_time_derivs = true;
    }
    
    
    
    // Find out about the coordinates: origin and delta for the Carpet
    // grid indices
#if 0
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
#else
    rvect const lower = rvect::ref(cgh->cctk_origin_space);
    rvect const delta = rvect::ref(cgh->cctk_delta_space) / maxreflevelfact;
    rvect const upper = lower + delta * (vhh.at(m)->baseextent.upper() - vhh.at(m)->baseextent.lower());
#endif
    
    assert (N_interp_points >= 0);
    assert (interp_coords);
    for (int d=0; d<dim; ++d) {
      assert (N_interp_points==0 || interp_coords[d]);
    }
    
    assert (N_output_arrays >= 0);
    if (N_interp_points > 0) {
      assert (output_arrays);
      for (int j=0; j<N_output_arrays; ++j) {
        assert (output_arrays[j]);
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
      rlev.at(n) = -1;
      home.at(n) = -1;
      if (all(pos>=lower && pos<=upper)) {
        for (int rl=maxrl-1; rl>=minrl; --rl) {
          
          int const fact = maxreflevelfact / ipow(reffact, rl) * ipow(mgfact, mglevel);
          ivect const ipos = ivect(floor((pos - lower) / (delta * fact) + 0.5)) * fact;
          assert (all(ipos % vhh.at(m)->bases().at(rl).at(ml).stride() == 0));
          
          // TODO: use something faster than a linear search
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            if (vhh.at(m)->extents().at(rl).at(c).at(ml).contains(ipos)) {
              rlev.at(n) = rl;
              home.at(n) = c;
              goto found;
            }
          }
        }
      }
      CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Interpolation point #%d at [%g,%g,%g] is not on any grid patch",
                  n, pos[0], pos[1], pos[2]);
      // Map the point (arbitrarily) to the first component of the
      // coarsest grid
      rlev.at(n) = minrl;
      home.at(n) = 0;
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
    vector<data<CCTK_REAL> > allcoords (nprocs * (maxrl-minrl) * maxncomps);
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
    for (comm_state state; !state.done(); state.step()) {
      for (int p=0; p<nprocs; ++p) {
        for (int rl=minrl; rl<maxrl; ++rl) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            allcoords.at(ind_prc(p,m,rl,c)).change_processor
              (state, vhh.at(m)->processors().at(rl).at(c));
          }
        }
      }
    }
    
    
    
    // Create output patches
    vector<data<CCTK_REAL> > alloutputs
      (nprocs * (maxrl-minrl) * maxncomps, -1);
    vector<data<CCTK_INT> > allstatuses
      (nprocs * (maxrl-minrl) * maxncomps, -1);
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          ivect const lo (0);
          ivect up (1);
          up[0] = allhomecnts.at(ind_prc(p,m,rl,c));
          up[1] = N_output_arrays;
          ivect const str (1);
          ibbox const extent (lo, up-str, str);
          alloutputs.at(ind_prc(p,m,rl,c)).allocate
            (extent, vhh.at(m)->processors().at(rl).at(c));
          
          ivect const slo (0);
          ivect sup (1);
          sup[0] = allhomecnts.at(ind_prc(p,m,rl,c));
          ivect const sstr (1);
          ibbox const sextent (lo, up-str, str);
          allstatuses.at(ind_prc(p,m,rl,c)).allocate
            (sextent, vhh.at(m)->processors().at(rl).at(c));
        }
      }
    }
    
    
    
    //
    // Do the local interpolation
    //
    int overall_ierr = 0;
    BEGIN_GLOBAL_MODE(cgh) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        enter_level_mode (const_cast<cGH*>(cgh), rl);
        
        // Number of necessary time levels
        CCTK_REAL const level_time = cgh->cctk_time / cgh->cctk_delta_time;
        bool const need_time_interp
          = (interpolation_times_differ || have_time_derivs
             || (fabs(current_time - level_time)
                 > 1e-12 * (fabs(level_time) + fabs(current_time)
                            + fabs(cgh->cctk_delta_time))));
        assert (! (! want_global_mode
                   && ! interpolation_times_differ && ! have_time_derivs
                   && need_time_interp));
        int const num_tl = need_time_interp ? prolongation_order_time + 1 : 1;
        
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
            vector<CCTK_INT> input_array_type_codes(N_input_arrays);
            vector<const void *> input_arrays(N_input_arrays);
            for (int n=0; n<N_input_arrays; ++n) {
              if (input_array_variable_indices[n] == -1) {
                
                // Ignore this entry
                input_array_type_codes.at(n) = -1;
                
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
                
                input_array_type_codes.at(n) = group.vartype;
                
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
              
              
              
              vector<vector<void *> > tmp_output_arrays (num_tl);
              vector<void *> tmp_status_array (num_tl);
              
              for (int tl=0; tl<num_tl; ++tl) {
                
                for (int n=0; n<N_input_arrays; ++n) {
                  
                  int const vi = input_array_variable_indices[n];
                  assert (vi>=0 && vi<CCTK_NumVars());
                  
                  int const gi = CCTK_GroupIndexFromVarI (vi);
                  assert (gi>=0 && gi<CCTK_NumGroups());
                  
                  int const table = CCTK_GroupTagsTableI (gi);
                  assert (table>=0);
                  
                  int my_num_tl;
                  CCTK_INT interp_num_time_levels;
                  int const ilen = Util_TableGetInt
                    (table, &interp_num_time_levels, "InterpNumTimelevels");
                  if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
                    my_num_tl = num_tl;
                  } else if (ilen >= 0) {
                    my_num_tl = interp_num_time_levels;
                    assert (my_num_tl>0 && my_num_tl<=num_tl);
                  } else {
                    assert (0);
                  }
                  
                  if (input_array_variable_indices[n] == -1) {
                    input_arrays.at(n) = 0;
                  } else if (tl >= my_num_tl) {
                    // Do a dummy interpolation from a later timelevel
                    int const vi = input_array_variable_indices[n];
                    assert (vi>=0 && vi<CCTK_NumVars());
                    input_arrays.at(n) = CCTK_VarDataPtrI (cgh, 0, vi);
                  } else {
                    int const vi = input_array_variable_indices[n];
                    assert (vi>=0 && vi<CCTK_NumVars());
                    input_arrays.at(n) = CCTK_VarDataPtrI (cgh, tl, vi);
                  }
                } // for input arrays
                
                tmp_output_arrays.at(tl).resize (N_output_arrays);
                for (int j=0; j<N_output_arrays; ++j) {
                  assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
                  if (need_time_interp) {
                    tmp_output_arrays.at(tl).at(j) = new CCTK_REAL [npoints];
                  } else {
                    tmp_output_arrays.at(tl).at(j) = &alloutputs.at(ind_prc(p,m,reflevel,component))[ivect(0,j,0)];
                  }
                }
                tmp_status_array.at(tl) = &allstatuses.at(ind_prc(p,m,reflevel,component))[ivect(0,0,0)];
                
                ierr = Util_TableSetPointer
                  (param_table_handle,
                   tmp_status_array.at(tl), "per_point_status");
                assert (ierr >= 0);
                
                ierr = CCTK_InterpLocalUniform
                  (N_dims, local_interp_handle, param_table_handle,
                   &coord_origin[0], &coord_delta[0],
                   npoints,
                   interp_coords_type_code, &tmp_interp_coords[0],
                   N_input_arrays, &lsh[0],
                   &input_array_type_codes[0], &input_arrays[0],
                   N_output_arrays,
                   output_array_type_codes, &tmp_output_arrays.at(tl)[0]);
                if (ierr) {
                  CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "The local interpolator returned the error code %d", ierr);
                }
                overall_ierr = min(overall_ierr, ierr);
                
                ierr = Util_TableDeleteKey
                  (param_table_handle, "per_point_status");
                assert (! ierr);
                
              } // for tl
              
              // Interpolate in time, if necessary
              if (need_time_interp) {
                
                for (int j=0; j<N_output_arrays; ++j) {
                  
                  // Find output variable indices
                  vector<CCTK_INT> operand_indices (N_output_arrays);
                  ierr = Util_TableGetIntArray
                    (param_table_handle, N_output_arrays,
                     &operand_indices.front(), "operand_indices");
                  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
                    assert (N_output_arrays == N_input_arrays);
                    for (int m=0; m<N_output_arrays; ++m) {
                      operand_indices.at(m) = m;
                    }
                  } else {
                    assert (ierr == N_output_arrays);
                  }
                  
                  // Find input array for this output array
                  int const n = operand_indices.at(j);
                  
                  int const vi = input_array_variable_indices[n];
                  assert (vi>=0 && vi<CCTK_NumVars());
                  
                  int const gi = CCTK_GroupIndexFromVarI (vi);
                  assert (gi>=0 && gi<CCTK_NumGroups());
                  
                  int const table = CCTK_GroupTagsTableI (gi);
                  assert (table>=0);
                  
                  int my_num_tl;
                  CCTK_INT interp_num_time_levels;
                  int const ilen = Util_TableGetInt
                    (table, &interp_num_time_levels, "InterpNumTimelevels");
                  if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
                    my_num_tl = num_tl;
                  } else if (ilen >= 0) {
                    my_num_tl = interp_num_time_levels;
                    assert (my_num_tl>0 && my_num_tl<=num_tl);
                  } else {
                    assert (0);
                  }
                  
                  if (! interpolation_times_differ && ! have_time_derivs) {
                    
                    // Get interpolation times
                    CCTK_REAL const time = current_time;
                    vector<CCTK_REAL> times(my_num_tl);
                    for (int tl=0; tl<my_num_tl; ++tl) {
                      times.at(tl) = vtt.at(m)->time (-tl, reflevel, mglevel);
                    }
                    
                    // Calculate interpolation weights
                    vector<CCTK_REAL> tfacs(my_num_tl);
                    switch (my_num_tl) {
                    case 1:
                      // no interpolation
                      // We have to assume that any GF with one timelevel
                      // is constant in time!!!
                      // assert (fabs((time - times.at(0)) / fabs(time + times.at(0) + cgh->cctk_delta_time)) < 1e-12);
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
                    
                    // Interpolate
                    assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
                    for (int k=0; k<npoints; ++k) {
                      CCTK_REAL & dest = alloutputs.at(ind_prc(p,m,reflevel,component))[ivect(k,j,0)];
                      dest = 0;
                      for (int tl=0; tl<my_num_tl; ++tl) {
                        CCTK_REAL const src = ((CCTK_REAL const *)tmp_output_arrays.at(tl).at(j))[k];
                        dest += tfacs[tl] * src;
                      }
                    }
                    
                  } else {
                    
                    // Interpolate
                    assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
                    for (int k=0; k<npoints; ++k) {
                      
                      CCTK_REAL const time = interpolation_times.at(k);
                      CCTK_INT const deriv_order = time_deriv_order.at(k);
                      
                      // Get interpolation times
                      vector<CCTK_REAL> times(my_num_tl);
                      for (int tl=0; tl<my_num_tl; ++tl) {
                        times.at(tl) = vtt.at(m)->time (-tl, reflevel, mglevel);
                      }
                      
                      // Calculate interpolation weights
                      vector<CCTK_REAL> tfacs(my_num_tl);
                      switch (deriv_order) {
                      case 0:
                        switch (my_num_tl) {
                        case 1:
                          // no interpolation
                          // We have to assume that any GF with one timelevel
                          // is constant in time!!!
                          // assert (fabs((time - times.at(0)) / fabs(time + times.at(0) + cgh->cctk_delta_time)) < 1e-12);
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
                        break;
                        
                      case 1:
                        switch (my_num_tl) {
                        case 2:
                          // linear (2-point) interpolation
                          // first order accurate
                          tfacs.at(0) = (1 - times.at(1)) / (times.at(0) - times.at(1));
                          tfacs.at(1) = (1 - times.at(0)) / (times.at(1) - times.at(0));
                          break;
                        case 3:
                          // quadratic (3-point) interpolation
                          // second order accurate
                          tfacs.at(0) = ((1 - times.at(1)) * (time - times.at(2)) + (time - times.at(1)) * (1 - times.at(2))) / ((times.at(0) - times.at(1)) * (times.at(0) - times.at(2)));
                          tfacs.at(1) = ((1 - times.at(0)) * (time - times.at(2)) + (time - times.at(0)) * (1 - times.at(2))) / ((times.at(1) - times.at(0)) * (times.at(1) - times.at(2)));
                          tfacs.at(2) = ((1 - times.at(0)) * (time - times.at(1)) + (time - times.at(0)) * (1 - times.at(1))) / ((times.at(2) - times.at(0)) * (times.at(2) - times.at(1)));
                          break;
                        default:
                          assert (0);
                        }
                        break;
                        
                      case 2:
                        switch (my_num_tl) {
                        case 3:
                          // quadratic (3-point) interpolation
                          // second order accurate
                          tfacs.at(0) = 2 * (1 - times.at(1)) * (1 - times.at(2)) / ((times.at(0) - times.at(1)) * (times.at(0) - times.at(2)));
                          tfacs.at(1) = 2 * (1 - times.at(0)) * (1 - times.at(2)) / ((times.at(1) - times.at(0)) * (times.at(1) - times.at(2)));
                          tfacs.at(2) = 2 * (1 - times.at(0)) * (1 - times.at(1)) / ((times.at(2) - times.at(0)) * (times.at(2) - times.at(1)));
                          break;
                        default:
                          assert (0);
                        }
                        break;
                        
                      default:
                        assert (0);
                      } // switch time_deriv_order
                      
                      CCTK_REAL & dest = alloutputs.at(ind_prc(p,m,reflevel,component))[ivect(k,j,0)];
                      dest = 0;
                      for (int tl=0; tl<my_num_tl; ++tl) {
                        CCTK_REAL const src = ((CCTK_REAL const *)tmp_output_arrays.at(tl).at(j))[k];
                        dest += tfacs[tl] * src;
                      }
                    }
                    
                  } // if interpolation_times_differ
                  
                } // for j
                  
                for (int tl=0; tl<num_tl; ++tl) {
                  for (int j=0; j<N_output_arrays; ++j) {
                    delete [] (CCTK_REAL *)tmp_output_arrays.at(tl).at(j);
                  }
                }
                
              } // if need_time_interp
              
            } // for processors
            
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
        
        leave_level_mode (const_cast<cGH*>(cgh));
      } // for rl
      
    } END_GLOBAL_MODE;
    
    
    
    // Transfer output patches back
    for (comm_state state; !state.done(); state.step()) {
      for (int p=0; p<nprocs; ++p) {
        for (int rl=minrl; rl<maxrl; ++rl) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            alloutputs.at(ind_prc(p,m,rl,c)).change_processor (state, p);
            allstatuses.at(ind_prc(p,m,rl,c)).change_processor (state, p);
          }
        }
      }
    }
    
    
    
    // Read out local output patches
    {
      vector<int> tmpcnts ((maxrl-minrl) * maxncomps);
      CCTK_INT local_interpolator_status = 0;
      for (int n=0; n<N_interp_points; ++n) {
        int const rl = rlev.at(n);
        int const c = home.at(n);
        for (int j=0; j<N_output_arrays; ++j) {
          assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
          assert (alloutputs.at(ind_prc(myproc,m,rl,c)).owns_storage());
          assert (output_arrays);
          assert (output_arrays[j]);
          static_cast<CCTK_REAL *>(output_arrays[j])[n] = alloutputs.at(ind_prc(myproc,m,rl,c))[ivect(tmpcnts.at(ind_rc(m,rl,c)),j,0)];
        }
        local_interpolator_status
          = min (local_interpolator_status,
                 allstatuses.at(ind_prc(myproc,m,rl,c))[ivect(tmpcnts.at(ind_rc(m,rl,c)),0,0)]);
        ++ tmpcnts.at(ind_rc(m,rl,c));
      }
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          assert (tmpcnts.at(ind_rc(m,rl,c)) == homecnts.at(ind_rc(m,rl,c)));
        }
      }
      ierr = Util_TableSetInt
        (param_table_handle,
         local_interpolator_status, "local_interpolator_status");
      assert (ierr >= 0);
    }
    
    
    
    int global_overall_ierr;
    MPI_Allreduce
      (&overall_ierr, &global_overall_ierr, 1, MPI_INT, MPI_MIN, comm);
    
    
    
    // Done.
    return global_overall_ierr;
  }
  
} // namespace CarpetInterp
