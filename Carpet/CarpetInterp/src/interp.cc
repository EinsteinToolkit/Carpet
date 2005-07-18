#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "bbox.hh"
#include "data.hh"
#include "defs.hh"
#include "dist.hh"
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
    int const ind = ((rl-minrl) * maps + m) * maxncomps + c;
    assert (ind>=0 && ind < (maxrl-minrl) * maps * maxncomps);
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
    int const ind
      = ((p * (maxrl-minrl) + (rl-minrl)) * maps + m) * maxncomps + c;
    assert (ind>=0 && ind < nprocs * (maxrl-minrl) * maps * maxncomps);
    return ind;
  }
  
  
  
  static void
  find_time_interpolation_order (cGH const * const cgh,
                                 int const param_table_handle,
                                 int const N_interp_points,
                                 int & prolongation_order_time,
                                 CCTK_REAL & current_time,
                                 vector<CCTK_INT> & time_deriv_order,
                                 vector<CCTK_REAL> & interpolation_times,
                                 bool & interpolation_times_differ,
                                 bool & have_time_derivs);
  
  static void
  assign_interpolation_points_to_components (vector<int> & rlev,
                                             vector<int> & home,
                                             vector<int> & homecnts,
                                             int const ml,
                                             vector<CCTK_INT> const & source_map,
                                             int const interp_coords_type_code,
                                             void const * const interp_coords [],
                                             vector<rvect> const & lower,
                                             vector<rvect> const & upper,
                                             vector<rvect> const & delta,
                                             int const minrl,
                                             int const maxrl,
                                             int const maxncomps,
                                             int const N_interp_points);
  
  static void
  create_coordinate_patches (cGH const * const cgh,
                             int const nprocs,
                             int const minrl,
                             int const maxrl,
                             int const maxncomps,
                             vector<int> const & allhomecnts,
                             vector<data<CCTK_INT> *> & allmaps,
                             vector<data<CCTK_REAL> *> & allcoords);
  
  static void
  fill_local_coordinate_patches (cGH const * const cgh,
                                 int const nprocs,
                                 int const minrl,
                                 int const maxrl,
                                 int const maxncomps,
                                 CCTK_INT const N_interp_points,
                                 vector<int> const & rlev,
                                 vector<CCTK_INT> const & source_map,
                                 vector<int> const & home,
                                 vector<int> const & homecnts,
                                 int const myproc,
                                 vector<int> const & allhomecnts,
                                 vector<data<CCTK_INT> *> & allmaps,
                                 vector<data<CCTK_REAL> *> & allcoords,
                                 CCTK_POINTER_TO_CONST const interp_coords []);
  
  static void
  transfer_coordinate_patches (cGH const * const cgh,
                               int const nprocs,
                               int const minrl,
                               int const maxrl,
                               int const maxncomps,
                               vector<data<CCTK_INT> *> & allmaps,
                               vector<data<CCTK_REAL> *> & allcoords);
  
  static void
  create_output_patches (cGH const * const cgh,
                         int const nprocs,
                         int const minrl,
                         int const maxrl,
                         int const maxncomps,
                         CCTK_INT const N_output_arrays,
                         vector<int> const & allhomecnts,
                         vector<data<CCTK_REAL> *> & allcoords,
                         vector<data<CCTK_REAL> *> & alloutputs,
                         vector<data<CCTK_INT> *> & allstatuses);
  
  static int
  do_local_interpolation (cGH const * const cgh,
                          int const minrl,
                          int const maxrl,
                          int const maxncomps,
                          bool const want_global_mode,
                          int const prolongation_order_time,
                          int const N_dims,
                          vector<data<CCTK_REAL> *> & allcoords,
                          vector<data<CCTK_REAL> *> & alloutputs,
                          vector<data<CCTK_INT> *> & allstatuses,
                          vector<CCTK_INT> & time_deriv_order,
                          vector<int> & allhomecnts,
                          vector<CCTK_REAL> & interpolation_times,
                          CCTK_INT const local_interp_handle,
                          CCTK_INT const param_table_handle,
                          CCTK_REAL const current_time,
                          int const interp_coords_type_code,
                          bool const interpolation_times_differ,
                          bool const have_time_derivs,
                          int const N_input_arrays,
                          int const N_output_arrays,
                          CCTK_INT const output_array_type_codes [],
                          CCTK_INT const input_array_variable_indices []);
  
  static void
  transfer_output_patches_back (cGH const * const cgh,
                                int const nprocs,
                                int const minrl,
                                int const maxrl,
                                int const maxncomps,
                                vector<data<CCTK_REAL> *> & alloutputs,
                                vector<data<CCTK_INT> *> & allstatuses);
  
  static void
  read_out_local_output_patches (cGH const * const cgh,
                                 int const nprocs,
                                 int const minrl,
                                 int const maxrl,
                                 int const maxncomps,
                                 int const myproc,
                                 CCTK_INT const N_interp_points,
                                 vector<int> const & rlev,
                                 vector<CCTK_INT> const & source_map,
                                 vector<int> const & home,
                                 CCTK_INT const N_output_arrays,
                                 CCTK_INT const interp_coords_type_code,
                                 vector<data<CCTK_REAL> *> const & alloutputs,
                                 CCTK_POINTER const output_arrays [],
                                 vector<data<CCTK_INT> *> const & allstatuses,
                                 vector<int> const & allhomecnts,
                                 vector<int> const & homecnts,
                                 CCTK_INT const param_table_handle);
  
  static int
  interpolate_within_component (cGH const * const cgh,
                                int const minrl,
                                int const maxrl,
                                int const maxncomps,
                                int const N_dims,
                                vector<data<CCTK_REAL> *> & allcoords,
                                vector<data<CCTK_REAL> *> & alloutputs,
                                vector<data<CCTK_INT> *> & allstatuses,
                                vector<CCTK_INT> & time_deriv_order,
                                vector<int> & allhomecnts,
                                vector<CCTK_REAL> &interpolation_times,
                                CCTK_INT const local_interp_handle,
                                CCTK_INT const param_table_handle,
                                CCTK_REAL const current_time,
                                int const interp_coords_type_code,
                                int const num_tl,
                                bool const need_time_interp,
                                bool const interpolation_times_differ,
                                bool const have_time_derivs,
                                int const N_input_arrays,
                                int const N_output_arrays,
                                CCTK_INT const output_array_type_codes [],
                                CCTK_INT const input_array_variable_indices []);
  
  
  
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
    
    // Multiple convergence levels are not supported
    assert (mglevels == 1);
    int const ml = 0;
    
    // Find range of refinement levels
    assert (maps > 0);
    for (int m=1; m<maps; ++m) {
      assert (vhh.at(0)->reflevels() == vhh.at(m)->reflevels());
    }
    int const minrl = want_global_mode ? 0                      : reflevel;
    int const maxrl = want_global_mode ? vhh.at(0)->reflevels() : reflevel+1;
    
    // Find source maps
    vector<CCTK_INT> source_map (N_interp_points);
    if (Carpet::map != -1) {
      // Interpolate from the current map
      for (int n=0; n<N_interp_points; ++n) {
        source_map.at(n) = Carpet::map;
      }
    } else if (maps == 1) {
      // Interpolate from map 0 if this is the only one
      // (for backwards compatibility)
      for (int n=0; n<N_interp_points; ++n) {
        source_map.at(n) = 0;
      }
    } else {
      source_map.resize (N_interp_points);
      ierr = Util_TableGetIntArray
        (param_table_handle, N_interp_points,
         &source_map.front(), "source_map");
      if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        CCTK_WARN (1, "No source map specified");
        return -1;
      } else if (ierr < 0) {
        CCTK_WARN (1, "internal error");
        return -1;
      } else if (ierr != N_interp_points) {
        CCTK_WARN (1, "Source map array has wrong size");
        return -1;
      }
      for (int n=0; n<N_interp_points; ++n) {
        assert (source_map.at(n) >=0 and source_map.at(n) < maps);
      }
    }
    
    int maxncomps = 0;
    for (int rl=minrl; rl<maxrl; ++rl) {
      for (int m=0; m<maps; ++m) {
        maxncomps = max(maxncomps, vhh.at(m)->components(rl));
      }
    }
    
    
    
    // Find the time interpolation order
    vector<CCTK_INT> time_deriv_order (N_interp_points);
    vector<CCTK_REAL> interpolation_times (N_interp_points);
    bool interpolation_times_differ;
    bool have_time_derivs;
    CCTK_REAL current_time;
    int prolongation_order_time;
    
    find_time_interpolation_order
      (cgh,
       param_table_handle,
       N_interp_points,
       prolongation_order_time, current_time,
       time_deriv_order, interpolation_times,
       interpolation_times_differ, have_time_derivs);
    
    
    
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
    vector<rvect> lower (maps);
    vector<rvect> delta (maps);
    vector<rvect> upper (maps);
    for (int m=0; m<maps; ++m) {
      lower.at(m) = rvect::ref(cgh->cctk_origin_space);
      delta.at(m) = rvect::ref(cgh->cctk_delta_space) / maxspacereflevelfact;
      upper.at(m) = lower.at(m) + delta.at(m) * (vhh.at(m)->baseextent.upper() - vhh.at(m)->baseextent.lower());
    }
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
    vector<int> homecnts ((maxrl-minrl) * maps * maxncomps); // number of points in component c
    
    assign_interpolation_points_to_components
      (rlev, home, homecnts,
       ml,
       source_map,
       interp_coords_type_code, interp_coords,
       lower, upper, delta,
       minrl, maxrl,
       maxncomps,
       N_interp_points);
    
    // Communicate counts
    vector<int> allhomecnts(nprocs * (maxrl-minrl) * maps * maxncomps);
    MPI_Allgather
      (&homecnts   .at(0), (maxrl-minrl) * maps * maxncomps, MPI_INT,
       &allhomecnts.at(0), (maxrl-minrl) * maps * maxncomps, MPI_INT, comm);
    
    
    
    // Create coordinate patches
    vector<data<CCTK_INT> *>  allmaps   (nprocs * (maxrl-minrl) * maps * maxncomps);
    vector<data<CCTK_REAL> *> allcoords (nprocs * (maxrl-minrl) * maps * maxncomps);
    
    create_coordinate_patches
      (cgh, nprocs,
       minrl, maxrl,
       maxncomps,
       allhomecnts,
       allmaps, allcoords);
    
    fill_local_coordinate_patches
      (cgh, nprocs,
       minrl, maxrl,
       maxncomps,
       N_interp_points,
       rlev,
       source_map,
       home,
       homecnts,
       myproc,
       allhomecnts,
       allmaps,
       allcoords,
       interp_coords);
    
    transfer_coordinate_patches
      (cgh, nprocs,
       minrl, maxrl,
       maxncomps,
       allmaps, allcoords);
    
    
    
    // Create output patches
    vector<data<CCTK_REAL> *> alloutputs
      (nprocs * (maxrl-minrl) * maps * maxncomps, static_cast<data<CCTK_REAL> *> (0));
    vector<data<CCTK_INT> *> allstatuses
      (nprocs * (maxrl-minrl) * maps * maxncomps, static_cast<data<CCTK_INT> *> (0));
    
    create_output_patches
      (cgh, nprocs,
       minrl, maxrl,
       maxncomps,
       N_output_arrays,
       allhomecnts,
       allcoords,
       alloutputs,
       allstatuses);
    
    
    
    //
    // Do the local interpolation
    //
    int const overall_ierr = do_local_interpolation
      (cgh,
       minrl, maxrl,
       maxncomps,
       want_global_mode,
       prolongation_order_time,
       N_dims,
       allcoords, alloutputs, allstatuses,
       time_deriv_order,
       allhomecnts,
       interpolation_times,
       local_interp_handle, param_table_handle,
       current_time,
       interp_coords_type_code,
       interpolation_times_differ,
       have_time_derivs,
       N_input_arrays, N_output_arrays,
       output_array_type_codes, input_array_variable_indices);
    
    
    
    // Free coordinate patches
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int m=0; m<maps; ++m) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            delete allmaps.at(ind_prc(p,m,rl,c));
            allmaps.at(ind_prc(p,m,rl,c)) = NULL;
            delete allcoords.at(ind_prc(p,m,rl,c));
            allcoords.at(ind_prc(p,m,rl,c)) = NULL;
          }
        }
      }
    }
    
    
    
    transfer_output_patches_back
      (cgh, nprocs,
       minrl, maxrl,
       maxncomps,
       alloutputs,
       allstatuses);
    
    
    
    read_out_local_output_patches
      (cgh, nprocs,
       minrl, maxrl,
       maxncomps,
       myproc,
       N_interp_points,
       rlev,
       source_map,
       home,
       N_output_arrays,
       interp_coords_type_code,
       alloutputs,
       output_arrays,
       allstatuses,
       allhomecnts,
       homecnts,
       param_table_handle);
    
    
    
    // Free local output patches
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int m=0; m<maps; ++m) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            delete alloutputs.at(ind_prc(p,m,rl,c));
            alloutputs.at(ind_prc(p,m,rl,c)) = NULL;
            delete allstatuses.at(ind_prc(p,m,rl,c));
            allstatuses.at(ind_prc(p,m,rl,c)) = NULL;
          }
        }
      }
    }
    
    
    
    int global_overall_ierr;
    MPI_Allreduce
      (const_cast<int *> (& overall_ierr), & global_overall_ierr,
       1, MPI_INT, MPI_MIN, comm);
    
    
    
    // Done.
    return global_overall_ierr;
  }
  
  
  
  static void
  find_time_interpolation_order (cGH const * const cgh,
                                 int const param_table_handle,
                                 int const N_interp_points,
                                 int & prolongation_order_time,
                                 CCTK_REAL & current_time,
                                 vector<CCTK_INT> & time_deriv_order,
                                 vector<CCTK_REAL> & interpolation_times,
                                 bool & interpolation_times_differ,
                                 bool & have_time_derivs)
  {
    // Find the time interpolation order
    int ierr;
    
    int partype;
    void const * const parptr
      = CCTK_ParameterGet ("prolongation_order_time", "Carpet", &partype);
    assert (parptr);
    assert (partype == PARAMETER_INTEGER);
    prolongation_order_time = * (CCTK_INT const *) parptr;
    
    current_time = cgh->cctk_time / cgh->cctk_delta_time;
    
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
  }
  
  
  
  static void
  assign_interpolation_points_to_components (vector<int> & rlev,
                                             vector<int> & home,
                                             vector<int> & homecnts,
                                             int const ml,
                                             vector<CCTK_INT> const & source_map,
                                             int const interp_coords_type_code,
                                             void const * const interp_coords [],
                                             vector<rvect> const & lower,
                                             vector<rvect> const & upper,
                                             vector<rvect> const & delta,
                                             int const minrl,
                                             int const maxrl,
                                             int const maxncomps,
                                             int const N_interp_points)
  {
    for (int n=0; n<N_interp_points; ++n) {
      
      int const m = source_map.at(n);
      
      // Find the grid point closest to the interpolation point
      assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
      rvect pos;
      for (int d=0; d<dim; ++d) {
        pos[d] = static_cast<CCTK_REAL const *>(interp_coords[d])[n];
      }
      
      // Find the component that this grid point belongs to
      rlev.at(n) = -1;
      home.at(n) = -1;
      if (all(pos>=lower.at(m) && pos<=upper.at(m))) {
        for (int rl=maxrl-1; rl>=minrl; --rl) {
          
          ivect const fact = maxspacereflevelfact / spacereffacts.at(rl) * ipow(mgfact, mglevel);
          ivect const ipos = ivect(floor((pos - lower.at(m)) / (delta.at(m) * fact) + 0.5)) * fact;
          assert (all(ipos % vhh.at(m)->bases().at(ml).at(rl).stride() == 0));
          
          // TODO: use something faster than a linear search
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            if (vhh.at(m)->extents().at(ml).at(rl).at(c).contains(ipos)) {
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
  }
  
  
  
  static void
  create_coordinate_patches (cGH const * const cgh,
                             int const nprocs,
                             int const minrl,
                             int const maxrl,
                             int const maxncomps,
                             vector<int> const & allhomecnts,
                             vector<data<CCTK_INT> *> & allmaps,
                             vector<data<CCTK_REAL> *> & allcoords)
  {
    // Create coordinate patches
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int m=0; m<maps; ++m) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            {
              ivect lo (0);
              ivect up (1);
              up[0] = allhomecnts.at(ind_prc(p,m,rl,c));
              up[1] = 1;
              ivect str (1);
              ibbox extent (lo, up-str, str);
              allmaps.at(ind_prc(p,m,rl,c)) = new data<CCTK_INT>;
              allmaps.at(ind_prc(p,m,rl,c))->allocate (extent, p);
            }
            {
              ivect lo (0);
              ivect up (1);
              up[0] = allhomecnts.at(ind_prc(p,m,rl,c));
              up[1] = dim;
              ivect str (1);
              ibbox extent (lo, up-str, str);
              allcoords.at(ind_prc(p,m,rl,c)) = new data<CCTK_REAL>;
              allcoords.at(ind_prc(p,m,rl,c))->allocate (extent, p);
            }
          }
        }
      }
    }
  }
  
  
  
  static void
  fill_local_coordinate_patches (cGH const * const cgh,
                                 int const nprocs,
                                 int const minrl,
                                 int const maxrl,
                                 int const maxncomps,
                                 CCTK_INT const N_interp_points,
                                 vector<int> const & rlev,
                                 vector<CCTK_INT> const & source_map,
                                 vector<int> const & home,
                                 vector<int> const & homecnts,
                                 int const myproc,
                                 vector<int> const & allhomecnts,
                                 vector<data<CCTK_INT> *> & allmaps,
                                 vector<data<CCTK_REAL> *> & allcoords,
                                 CCTK_POINTER_TO_CONST const interp_coords [])
  {
    // Fill in local coordinate patches
    vector<int> tmpcnts ((maxrl-minrl) * maps * maxncomps);
    for (int n=0; n<N_interp_points; ++n) {
      int const rl = rlev.at(n);
      int const m = source_map.at(n);
      int const c = home.at(n);
      assert (rl>=minrl && rl<maxrl);
      assert (m>=0 && m<maps);
      assert (c>=0 && c<vhh.at(m)->components(rl));
      assert (tmpcnts.at(ind_rc(m,rl,c)) >= 0);
      assert (tmpcnts.at(ind_rc(m,rl,c)) < homecnts.at(ind_rc(m,rl,c)));
      assert (dim==3);
      (*allmaps.at(ind_prc(myproc,m,rl,c)))[ivect(tmpcnts.at(ind_rc(m,rl,c)),0,0)]
        = source_map.at(n);
      for (int d=0; d<dim; ++d) {
        (*allcoords.at(ind_prc(myproc,m,rl,c)))[ivect(tmpcnts.at(ind_rc(m,rl,c)),d,0)]
          = static_cast<CCTK_REAL const *>(interp_coords[d])[n];
      }
      ++ tmpcnts.at(c + maxncomps * (m + maps * (rl-minrl)));
    }
    for (int rl=minrl; rl<maxrl; ++rl) {
      for (int m=0; m<maps; ++m) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          assert (tmpcnts.at(ind_rc(m,rl,c)) == homecnts.at(ind_rc(m,rl,c)));
        }
      }
    }
  }
  
  
  
  static void
  transfer_coordinate_patches (cGH const * const cgh,
                               int const nprocs,
                               int const minrl,
                               int const maxrl,
                               int const maxncomps,
                               vector<data<CCTK_INT> *> & allmaps,
                               vector<data<CCTK_REAL> *> & allcoords)
  {
    // Transfer coordinate patches
    for (comm_state state; !state.done(); state.step()) {
      for (int p=0; p<nprocs; ++p) {
        for (int rl=minrl; rl<maxrl; ++rl) {
          for (int m=0; m<maps; ++m) {
            for (int c=0; c<vhh.at(m)->components(rl); ++c) {
              allmaps.at(ind_prc(p,m,rl,c))->change_processor
                (state, vhh.at(m)->processors().at(rl).at(c));
              allcoords.at(ind_prc(p,m,rl,c))->change_processor
                (state, vhh.at(m)->processors().at(rl).at(c));
            }
          }
        }
      }
    }
  }
  
  
  
  static void
  create_output_patches (cGH const * const cgh,
                         int const nprocs,
                         int const minrl,
                         int const maxrl,
                         int const maxncomps,
                         CCTK_INT const N_output_arrays,
                         vector<int> const & allhomecnts,
                         vector<data<CCTK_REAL> *> & allcoords,
                         vector<data<CCTK_REAL> *> & alloutputs,
                         vector<data<CCTK_INT> *> & allstatuses)
  {
    for (int p=0; p<nprocs; ++p) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        for (int m=0; m<maps; ++m) {
          for (int c=0; c<vhh.at(m)->components(rl); ++c) {
            ivect const lo (0);
            ivect up (1);
            up[0] = allhomecnts.at(ind_prc(p,m,rl,c));
            up[1] = N_output_arrays;
            ivect const str (1);
            ibbox const extent (lo, up-str, str);
            alloutputs.at(ind_prc(p,m,rl,c)) = new data<CCTK_REAL>;
            alloutputs.at(ind_prc(p,m,rl,c))->allocate
              (extent, vhh.at(m)->processors().at(rl).at(c));
            
            ivect const slo (0);
            ivect sup (1);
            sup[0] = allhomecnts.at(ind_prc(p,m,rl,c));
            ivect const sstr (1);
            ibbox const sextent (lo, up-str, str);
            allstatuses.at(ind_prc(p,m,rl,c)) = new data<CCTK_INT>;
            allstatuses.at(ind_prc(p,m,rl,c))->allocate
              (sextent, vhh.at(m)->processors().at(rl).at(c));
          }
        }
      }
    }
  }
  
  
  
  static int
  do_local_interpolation (cGH const * const cgh,
                          int const minrl,
                          int const maxrl,
                          int const maxncomps,
                          bool const want_global_mode,
                          int const prolongation_order_time,
                          int const N_dims,
                          vector<data<CCTK_REAL> *> & allcoords,
                          vector<data<CCTK_REAL> *> & alloutputs,
                          vector<data<CCTK_INT> *> & allstatuses,
                          vector<CCTK_INT> & time_deriv_order,
                          vector<int> & allhomecnts,
                          vector<CCTK_REAL> & interpolation_times,
                          CCTK_INT const local_interp_handle,
                          CCTK_INT const param_table_handle,
                          CCTK_REAL const current_time,
                          int const interp_coords_type_code,
                          bool const interpolation_times_differ,
                          bool const have_time_derivs,
                          int const N_input_arrays,
                          int const N_output_arrays,
                          CCTK_INT const output_array_type_codes [],
                          CCTK_INT const input_array_variable_indices [])
  {
    // Do the local interpolation
    int ierr;
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
            
            ierr = interpolate_within_component
              (cgh,
               minrl, maxrl, maxncomps,
               N_dims,
               allcoords, alloutputs, allstatuses,
               time_deriv_order,
               allhomecnts,
               interpolation_times,
               local_interp_handle, param_table_handle,
               current_time,
               interp_coords_type_code,
               num_tl,
               need_time_interp,
               interpolation_times_differ,
               have_time_derivs,
               N_input_arrays, N_output_arrays,
               output_array_type_codes,
               input_array_variable_indices);
            overall_ierr = min(overall_ierr, ierr);
            
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
        
        leave_level_mode (const_cast<cGH*>(cgh));
      } // for rl
      
    } END_GLOBAL_MODE;
    
    return overall_ierr;
  }
  
  
  
  static void
  transfer_output_patches_back (cGH const * const cgh,
                                int const nprocs,
                                int const minrl,
                                int const maxrl,
                                int const maxncomps,
                                vector<data<CCTK_REAL> *> & alloutputs,
                                vector<data<CCTK_INT> *> & allstatuses)
  {
    // Transfer output patches back
    for (comm_state state; !state.done(); state.step()) {
      for (int p=0; p<nprocs; ++p) {
        for (int rl=minrl; rl<maxrl; ++rl) {
          for (int m=0; m<maps; ++m) {
            for (int c=0; c<vhh.at(m)->components(rl); ++c) {
              alloutputs.at(ind_prc(p,m,rl,c))->change_processor (state, p);
              allstatuses.at(ind_prc(p,m,rl,c))->change_processor (state, p);
            }
          }
        }
      }
    }
  }
  
  
  
  static void
  read_out_local_output_patches (cGH const * const cgh,
                                 int const nprocs,
                                 int const minrl,
                                 int const maxrl,
                                 int const maxncomps,
                                 int const myproc,
                                 CCTK_INT const N_interp_points,
                                 vector<int> const & rlev,
                                 vector<CCTK_INT> const & source_map,
                                 vector<int> const & home,
                                 CCTK_INT const N_output_arrays,
                                 CCTK_INT const interp_coords_type_code,
                                 vector<data<CCTK_REAL> *> const & alloutputs,
                                 CCTK_POINTER const output_arrays [],
                                 vector<data<CCTK_INT> *> const & allstatuses,
                                 vector<int> const & allhomecnts,
                                 vector<int> const & homecnts,
                                 CCTK_INT const param_table_handle)
  {
    // Read out local output patches
    vector<int> tmpcnts ((maxrl-minrl) * maps * maxncomps);
    CCTK_INT local_interpolator_status = 0;
    for (int n=0; n<N_interp_points; ++n) {
      int const rl = rlev.at(n);
      int const m = source_map.at(n);
      int const c = home.at(n);
      for (int j=0; j<N_output_arrays; ++j) {
        assert (interp_coords_type_code == CCTK_VARIABLE_REAL);
        assert (alloutputs.at(ind_prc(myproc,m,rl,c))->proc() == dist::rank());
        assert (output_arrays);
        assert (output_arrays[j]);
        static_cast<CCTK_REAL *>(output_arrays[j])[n]
          = (*alloutputs.at(ind_prc(myproc,m,rl,c)))[ivect(tmpcnts.at(ind_rc(m,rl,c)),j,0)];
      }
      local_interpolator_status
        = min (local_interpolator_status,
               (*allstatuses.at(ind_prc(myproc,m,rl,c)))[ivect(tmpcnts.at(ind_rc(m,rl,c)),0,0)]);
      ++ tmpcnts.at(ind_rc(m,rl,c));
    }
    for (int rl=minrl; rl<maxrl; ++rl) {
      for (int m=0; m<maps; ++m) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          assert (tmpcnts.at(ind_rc(m,rl,c)) == homecnts.at(ind_rc(m,rl,c)));
        }
      }
    }
    int ierr = Util_TableSetInt
      (param_table_handle,
       local_interpolator_status, "local_interpolator_status");
    assert (ierr >= 0);
  }
  
 
  /** A list of time values corresoponding to the last few timelevels
   * on the given patch.
   *
   * Values are determined by secret magic.
   * */
  class InterpolationTimes : private vector<CCTK_REAL>
  {
  public:
    InterpolationTimes (CCTK_INT num_timelevels )
      : vector<CCTK_REAL> (num_timelevels )
    {
      for (int tl=0; tl<num_timelevels; ++tl) {
        // SW: vtt, Carpet::map, reflevel and mglevel are all globals :<
 //       at(tl) = vtt.at(map)->time (-tl, reflevel, mglevel);
        // SW oh dear, this formula changed under my feet!
        at(tl) = vtt.at(Carpet::map)->time (tl, reflevel, mglevel);
      }
    }

    CCTK_REAL const & operator [] (size_type i ) const {
      return at (i );
    }

    int num_timelevels () const {
      return (int)size ();
    }
  };

  /** Represents interpolstion coefficients, or weights, for a given
   *  derivative order and set  of time levels.
   *  
   * */
  class InterpolationWeights : private vector<CCTK_REAL>
  {
  public:
    InterpolationWeights (CCTK_INT deriv_order, 
		    const InterpolationTimes & times, CCTK_REAL current_time )
      : vector<CCTK_REAL> (times.num_timelevels () )
    {
      CCTK_INT num_timelevels = times.num_timelevels ();
      // Initialise weights
      switch (deriv_order) {
      case 0:
        switch (num_timelevels) {
        case 1:
          for_no_interp (times, current_time );
          break;
        case 2:
          for_linear_2_pt_interp (times, current_time );
          break;
        case 3:
          for_quadratic_3_pt_interp (times, current_time );
          break;
        default:
          assert (0);
        }
      break;
                
      case 1:
        switch (num_timelevels) {
        case 2:
          for_1st_deriv_linear_2_pt_interp_1st_order (times, current_time );
          break;
        case 3:
          for_1st_deriv_quadratic_3_pt_interp_2nd_order (times, current_time );
          break;
        default:
          assert (0);
        }
        break;
        
      case 2:
        switch (num_timelevels) {
        case 3:
          for_2nd_deriv_quadratic_3_pt_interp_2nd_order (times, current_time );
          break;
        default:
        assert (0);
      }
      break;
        
      default:
        assert (0);
      } // switch time_deriv_order
    }

    CCTK_REAL const & operator [] (size_type i ) const {
      return at (i );
    }

  private:
    void for_no_interp (const InterpolationTimes & t, CCTK_REAL time )
    {
      // We have to assume that any GF with one timelevel is constant in time!!!
      // assert (fabs((time - t.at(0)) / fabs(time + t.at(0) + cgh->cctk_delta_time)) < 1e-12);
      at(0) = 1.0;
    }

    void for_linear_2_pt_interp (const InterpolationTimes & t, CCTK_REAL t0 )
    {
      at(0) = (t0 - t[1]) / (t[0] - t[1]);
      at(1) = (t0 - t[0]) / (t[1] - t[0]);
    }

    void for_quadratic_3_pt_interp (const InterpolationTimes & t, CCTK_REAL t0 )
    {
      at(0) = (t0 - t[1]) * (t0 - t[2]) / ((t[0] - t[1]) * (t[0] - t[2]));
      at(1) = (t0 - t[0]) * (t0 - t[2]) / ((t[1] - t[0]) * (t[1] - t[2]));
      at(2) = (t0 - t[0]) * (t0 - t[1]) / ((t[2] - t[0]) * (t[2] - t[1]));
    }

    void for_1st_deriv_linear_2_pt_interp_1st_order (
                                               const InterpolationTimes & t,
                                               CCTK_REAL t0 )
    {
      at(0) = (1 - t[1]) / (t[0] - t[1]);
      at(1) = (1 - t[0]) / (t[1] - t[0]);
    }

    void for_1st_deriv_quadratic_3_pt_interp_2nd_order (
                                               const InterpolationTimes & t,
                                               CCTK_REAL t0 )
    {
      at(0) = ((1 - t[1]) * (t0 - t[2]) + (t0 - t[1]) * (1 - t[2]))
                  / ((t[0] - t[1]) * (t[0] - t[2]));
      at(1) = ((1 - t[0]) * (t0 - t[2]) + (t0 - t[0]) * (1 - t[2]))
                  / ((t[1] - t[0]) * (t[1] - t[2]));
      at(2) = ((1 - t[0]) * (t0 - t[1]) + (t0 - t[0]) * (1 - t[1]))
                  / ((t[2] - t[0]) * (t[2] - t[1]));
    }

    void for_2nd_deriv_quadratic_3_pt_interp_2nd_order (
                                               const InterpolationTimes & t,
                                               CCTK_REAL t0 )
    {
      at(0) = 2 * (1 - t[1]) * (1 - t[2]) / ((t[0] - t[1]) * (t[0] - t[2]));
      at(1) = 2 * (1 - t[0]) * (1 - t[2]) / ((t[1] - t[0]) * (t[1] - t[2]));
      at(2) = 2 * (1 - t[0]) * (1 - t[1]) / ((t[2] - t[0]) * (t[2] - t[1]));
    }
  };
  
  
  static int
  interpolate_within_component (cGH const * const cgh,
                                int const minrl,
                                int const maxrl,
                                int const maxncomps,
                                int const N_dims,
                                vector<data<CCTK_REAL> *> & allcoords,
                                vector<data<CCTK_REAL> *> & alloutputs,
                                vector<data<CCTK_INT> *> & allstatuses,
                                vector<CCTK_INT> & time_deriv_order,
                                vector<int> & allhomecnts,
                                vector<CCTK_REAL> &interpolation_times,
                                CCTK_INT const local_interp_handle,
                                CCTK_INT const param_table_handle,
                                CCTK_REAL const current_time,
                                int const interp_coords_type_code,
                                int const num_tl,
                                bool const need_time_interp,
                                bool const interpolation_times_differ,
                                bool const have_time_derivs,
                                int const N_input_arrays,
                                int const N_output_arrays,
                                CCTK_INT const output_array_type_codes [],
                                CCTK_INT const input_array_variable_indices [])
  {
    int ierr;
    int overall_ierr = 0;
    int const nprocs = CCTK_nProcs (cgh);
    
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
      assert (allcoords.at(ind_prc(p,Carpet::map,reflevel,component))->proc() == dist::rank());
      assert (allhomecnts.at(ind_prc(p,Carpet::map,reflevel,component)) == allcoords.at(ind_prc(p,Carpet::map,reflevel,component))->shape()[0]);
      assert (allhomecnts.at(ind_prc(p,Carpet::map,reflevel,component)) == alloutputs.at(ind_prc(p,Carpet::map,reflevel,component))->shape()[0]);
      
      int const npoints = allhomecnts.at(ind_prc(p,Carpet::map,reflevel,component));
      
      // Do the processor-local interpolation
      const void * tmp_interp_coords[dim];
      for (int d=0; d<dim; ++d) {
        tmp_interp_coords[d] = &(*allcoords.at(ind_prc(p,Carpet::map,reflevel,component)))[ivect(0,d,0)];
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
          
          if (vi == -1) {
            input_arrays.at(n) = NULL;
          } else {
            // Do a dummy interpolation from a later timelevel
            // if the desired timelevel does not exist
            int const my_tl = tl >= my_num_tl ? 0 : tl;
            assert (vi>=0 && vi<CCTK_NumVars());
#if 0
            input_arrays.at(n) = CCTK_VarDataPtrI (cgh, my_tl, vi);
#else
            int const vi0 = CCTK_FirstVarIndexI (gi);
            input_arrays.at(n)
              = ((*arrdata.at(gi).at(Carpet::map).data.at(vi-vi0))
                 (my_tl, reflevel, component, mglevel)->storage());
#endif
          }
        } // for input arrays
        
        tmp_output_arrays.at(tl).resize (N_output_arrays);
        for (int j=0; j<N_output_arrays; ++j) {
          assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
          if (need_time_interp) {
            tmp_output_arrays.at(tl).at(j) = new CCTK_REAL [npoints];
          } else {
            tmp_output_arrays.at(tl).at(j) = &(*alloutputs.at(ind_prc(p,Carpet::map,reflevel,component)))[ivect(0,j,0)];
          }
        }
        tmp_status_array.at(tl) = &(*allstatuses.at(ind_prc(p,Carpet::map,reflevel,component)))[ivect(0,0,0)];
        
        ierr = Util_TableSetPointer
          (param_table_handle, tmp_status_array.at(tl), "per_point_status");
        assert (ierr >= 0);
        
        ierr = CCTK_InterpLocalUniform
          (N_dims, local_interp_handle, param_table_handle,
           &coord_origin[0], &coord_delta[0],
           npoints,
           interp_coords_type_code, tmp_interp_coords,
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
            const InterpolationTimes times (my_num_tl );
            const InterpolationWeights tfacs (0, times, current_time );

            // Interpolate
            assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
            for (int k=0; k<npoints; ++k) {
              CCTK_REAL & dest = (*alloutputs.at(ind_prc(p,Carpet::map,reflevel,component)))[ivect(k,j,0)];
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
              
              const InterpolationTimes times (my_num_tl );
              const InterpolationWeights tfacs (deriv_order, times,
			      current_time );
             
              CCTK_REAL & dest = (*alloutputs.at(ind_prc(p,Carpet::map,reflevel,component)))[ivect(k,j,0)];
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
    
    return overall_ierr;
  }
  
} // namespace CarpetInterp
