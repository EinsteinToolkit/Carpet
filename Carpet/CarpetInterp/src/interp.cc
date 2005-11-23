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


//////////////////////////////////////////////////////////////////////////
// For a general overview on the implementation
// please see CarpetInterp's thorn documentation.
//////////////////////////////////////////////////////////////////////////


namespace CarpetInterp {

  using namespace std;
  using namespace Carpet;


// Return a unique index for component c on reflevel rl, map m and processor p
#define component_idx(p,m,rl,c)    \
        component_idx_ (p, m, rl, c, minrl, maxrl, maxncomps)
  static inline int component_idx_ (int const p,
                                    int const m,
                                    int const rl,
                                    int const c,
                                    int const minrl, int const maxrl,
                                    int const maxncomps)
  {
    assert (p  >= 0     and p  < dist::size());
    assert (m  >= 0     and m  < maps);
    assert (rl >= minrl and rl < maxrl);
    assert (c  >= 0     and c  < maxncomps);
    int const local_idx  = ((rl-minrl) * maps + m) * maxncomps + c;
    int const global_idx = p * (maxrl-minrl) * maps * maxncomps + local_idx;
    return global_idx;
  }


  static int
  extract_parameter_table_options (cGH const * const cctkGH,
                                   int const param_table_handle,
                                   int const N_interp_points,
                                   int const N_input_arrays,
                                   int const N_output_arrays,
                                   bool& have_source_map,
                                   bool& have_time_derivs,
                                   int & prolongation_order_time,
                                   CCTK_REAL & current_time,
                                   CCTK_REAL & delta_time,
                                   vector<CCTK_INT> & source_map,
                                   vector<CCTK_INT> & operand_indices,
                                   vector<CCTK_INT> & time_deriv_order);

  static void
  map_points (cGH const * const cctkGH,
              int const coord_system_handle,
              int const coord_group,
              int const ml,
              int const minrl,
              int const maxrl,
              int const maxncomps,
              int const N_dims,
              int const npoints,
              vector<CCTK_INT> const& source_map,
              void const * const coords_list[],
              vector<CCTK_REAL> const& coords,
              vector<int>& procs,
              vector<int>& N_points_to,
              vector<int>& rlev,
              vector<int>& home,
              vector<int>& homecnts);

  static void
  interpolate_components (cGH const * const cctkGH,
                          int const coord_system_handle,
                          int const coord_group,
                          int const minrl,
                          int const maxrl,
                          int const maxncomps,
                          bool const want_global_mode,
                          int const prolongation_order_time,
                          int const N_dims,
                          vector<int> & homecnts,
                          vector<CCTK_REAL*> & coords,
                          vector<char*> & outputs,
                          CCTK_INT* const per_proc_statuses,
                          CCTK_INT* const per_proc_retvals,
                          vector<CCTK_INT> & operand_indices,
                          vector<CCTK_INT> & time_deriv_order,
                          bool const have_time_derivs,
                          CCTK_INT const local_interp_handle,
                          CCTK_INT const param_table_handle,
                          CCTK_REAL const current_time,
                          CCTK_REAL const delta_time,
                          int const N_input_arrays,
                          int const N_output_arrays,
                          CCTK_INT const output_array_type_codes [],
                          CCTK_INT const input_array_variable_indices []);

  static void
  interpolate_single_component (cGH const * const cctkGH,
                                int const coord_system_handle, 
                                int const coord_group,
                                int const minrl,
                                int const maxrl,
                                int const maxncomps,
                                int const N_dims,
                                int const npoints,
                                CCTK_REAL const* const coords,
                                char* const outputs,
                                CCTK_INT& overall_status,
                                CCTK_INT& overall_retval,
                                vector<CCTK_INT> & operand_indices,
                                vector<CCTK_INT> & time_deriv_order,
                                vector<CCTK_INT> & interp_num_time_levels,
                                CCTK_INT const local_interp_handle,
                                CCTK_INT const param_table_handle,
                                int const num_tl,
                                bool const need_time_interp,
                                CCTK_REAL const current_time,
                                CCTK_REAL const delta_time,
                                int const N_input_arrays,
                                int const N_output_arrays,
                                CCTK_INT const output_array_type_codes [],
                                CCTK_INT const input_array_variable_indices []);



  void CarpetInterpStartup ()
  {
    CCTK_OverloadInterpGridArrays (InterpGridArrays);
  }



  int InterpGridArrays (cGH const * const cctkGH,
                        int const N_dims,
                        int const local_interp_handle,
                        int const param_table_handle,
                        int const coord_system_handle,
                        int const N_interp_points,
                        int const interp_coords_type_code,
                        void const * const coords [],
                        int const N_input_arrays,
                        CCTK_INT const input_array_variable_indices [],
                        int const N_output_arrays,
                        CCTK_INT const output_array_type_codes [],
                        void * const output_arrays [])
  {
    if (CCTK_IsFunctionAliased ("SymmetryInterpolate")) {
      return SymmetryInterpolate
        (cctkGH, N_dims,
         local_interp_handle, param_table_handle, coord_system_handle,
         N_interp_points, interp_coords_type_code, coords,
         N_input_arrays, input_array_variable_indices,
         N_output_arrays, output_array_type_codes, output_arrays);
    } else {
      return Carpet_DriverInterpolate
        (cctkGH, N_dims,
         local_interp_handle, param_table_handle, coord_system_handle,
         N_interp_points, interp_coords_type_code, coords,
         N_input_arrays, input_array_variable_indices,
         N_output_arrays, output_array_type_codes, output_arrays);
    }
  }



  extern "C" CCTK_INT
  Carpet_DriverInterpolate (CCTK_POINTER_TO_CONST const cctkGH_,
                            CCTK_INT const N_dims,
                            CCTK_INT const local_interp_handle,
                            CCTK_INT const param_table_handle,
                            CCTK_INT const coord_system_handle,
                            CCTK_INT const N_interp_points,
                            CCTK_INT const interp_coords_type_code,
                            CCTK_POINTER_TO_CONST const coords_list [],
                            CCTK_INT const N_input_arrays,
                            CCTK_INT const input_array_variable_indices [],
                            CCTK_INT const N_output_arrays,
                            CCTK_INT const output_array_type_codes [],
                            CCTK_POINTER const output_arrays [])
  {
    cGH const * const cctkGH = static_cast<cGH const *> (cctkGH_);
    assert (cctkGH);
    assert (0 <= N_dims and N_dims <= dim);

    // Check input arrays
    int coord_group = -1;
    cGroupDynamicData coord_group_data;
    for (int n = 0; n < N_input_arrays; n++) {

      // Negative indices are ignored
      const int vindex = input_array_variable_indices[n];
      if (vindex < 0) {
        continue;
      }

      const int group = CCTK_GroupIndexFromVarI (vindex);
      if (group < 0) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "input array variable %d = %d is not a valid "
                   "variable index", n, vindex);
      }
      const int gtype = CCTK_GroupTypeI (group);
      if (gtype != CCTK_GF and gtype != CCTK_ARRAY) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "input array variable %d is not of type CCTK_GF or "
                    "CCTK_ARRY", n);
      }
      if (CCTK_GroupStaggerIndexGI (group) != 0) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "interpolation of staggered input array variable %d "
                    "is not supported", n);
      }
      if (coord_group < 0) {
        coord_group = group;
        CCTK_GroupDynamicData (cctkGH, coord_group, &coord_group_data);
      } else {
        // Check that group and coord_group have the same layout
        cGroupDynamicData gdata;
        CCTK_GroupDynamicData (cctkGH, group, &gdata);
        const int size = gdata.dim * sizeof (int);
        if (gdata.dim != coord_group_data.dim or
            memcmp (gdata.lsh,         coord_group_data.lsh,  size) or
            memcmp (gdata.lbnd,        coord_group_data.lbnd, size) or
            memcmp (gdata.ubnd,        coord_group_data.ubnd, size) or
            memcmp (gdata.bbox,        coord_group_data.bbox, 2*size) or
            memcmp (gdata.nghostzones, coord_group_data.nghostzones, size)) {
          CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "input array variable %d has different layout than "
                      "the underlying coordinate system", n);
        }
      }
    }
    assert (coord_group >= 0);

    // Check output arrays
    assert (N_output_arrays > 0);
    const int output_array_type = output_array_type_codes[0];
    for (int n = 1; n < N_output_arrays; n++) {
      if (output_array_type != output_array_type_codes[n]) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Currently all output arrays have have the same datatype. "
                    "Array 0 has type '%s' but array %d has type '%s'",
                    CCTK_VarTypeName (output_array_type), n,
                    CCTK_VarTypeName (output_array_type_codes[n]));
      }
    }


    if (is_meta_mode()) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "It is not possible to interpolate in meta mode");
    }

    bool const want_global_mode = is_global_mode();


    // Multiple convergence levels are not supported
    assert (mglevels == 1);
    int const ml = 0;

    assert (N_interp_points >= 0);
    assert (coords_list);
    for (int d = 0; d < N_dims; ++d) {
      assert (N_interp_points == 0 or coords_list[d]);
    }

    if (interp_coords_type_code != CCTK_VARIABLE_REAL) {
      CCTK_WARN (CCTK_WARN_ABORT, "CarpetInterp does not support interpolation "
                 "coordinates other than datatype CCTK_VARIABLE_REAL");
    }

    assert (N_output_arrays >= 0);
    if (N_interp_points > 0) {
      assert (output_arrays);
      for (int j=0; j<N_output_arrays; ++j) {
        assert (output_arrays[j]);
        for (int jj=0; jj<j; ++jj) {
          assert (output_arrays[j] != output_arrays[jj]);
        }
      }
    }


    // Find range of refinement levels
    assert (maps > 0);
    for (int m=1; m<maps; ++m) {
      assert (arrdata.at(coord_group).at(0).hh->reflevels() ==
              arrdata.at(coord_group).at(m).hh->reflevels());
    }
    int const minrl = want_global_mode ? 0 : reflevel;
    int const maxrl = want_global_mode ?
                        arrdata[coord_group][0].hh->reflevels() : reflevel+1;

    // Find maximum number of components over all levels and maps
    int maxncomps = 0;
    for (int rl=minrl; rl<maxrl; ++rl) {
      for (int m=0; m<maps; ++m) {
        maxncomps = max(maxncomps, arrdata[coord_group][m].hh->components(rl));
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Extract parameter table options:
    //   - source map
    //   - output array operand indices
    //   - time interpolation order
    //////////////////////////////////////////////////////////////////////
    vector<CCTK_INT> source_map (N_interp_points);
    vector<CCTK_INT> operand_indices (N_output_arrays);
    vector<CCTK_INT> time_deriv_order (N_output_arrays);
    bool have_source_map, have_time_derivs;
    CCTK_REAL current_time, delta_time;
    int prolongation_order_time;

    int iret = extract_parameter_table_options
      (cctkGH,
       param_table_handle,
       N_interp_points, N_input_arrays, N_output_arrays,
       have_source_map, have_time_derivs,
       prolongation_order_time, current_time, delta_time,
       source_map, operand_indices, time_deriv_order);
    if (iret < 0) {
      return iret;
    }

    //////////////////////////////////////////////////////////////////////
    // Map interpolation points to processors
    //////////////////////////////////////////////////////////////////////
    vector<int> N_points_to (dist::size());
    vector<int> rlev (N_interp_points); // refinement level of point n
    vector<int> home (N_interp_points); // component of point n
    vector<int> dstprocs (N_interp_points); // which processor owns point n
    vector<int> allhomecnts ((maxrl-minrl) * maps * maxncomps * dist::size());
                                            // number of points in component c
    vector<CCTK_REAL> coords_buffer (N_dims * N_interp_points);

    // Each point from coord_list is mapped onto the processor
    // that owns it (dstprocs)
    // Also accumulate the number of points per processor (N_points_to)
    map_points (cctkGH, coord_system_handle, coord_group,
                ml, minrl, maxrl, maxncomps, N_dims, N_interp_points,
                source_map,
                coords_list, coords_buffer,
                dstprocs, N_points_to,
                rlev, home, allhomecnts);

    //////////////////////////////////////////////////////////////////////
    // Communicate the number of points each processor is going to communicate
    //////////////////////////////////////////////////////////////////////

    // N_points_from denotes the number of points
    // that this processor is to receive from others
    vector<int> N_points_from (dist::size());
    MPI_Alltoall (&N_points_to[0],   1, dist::datatype (N_points_to[0]),
                  &N_points_from[0], 1, dist::datatype (N_points_from[0]),
                  dist::comm());

    //////////////////////////////////////////////////////////////////////
    // Communicate the interpolation coordinates
    //////////////////////////////////////////////////////////////////////
    vector<int> sendcnt (dist::size());
    vector<int> recvcnt (dist::size());
    vector<int> senddispl (dist::size());
    vector<int> recvdispl (dist::size());

    // Set up counts and displacements for MPI_Alltoallv()
    sendcnt[0] = N_dims * N_points_to[0];
    recvcnt[0] = N_dims * N_points_from[0];
    senddispl[0] = recvdispl[0] = 0;
    int N_points_local = recvcnt[0];
    for (int p = 1; p < dist::size(); p++) {
      sendcnt[p] = N_dims * N_points_to[p];
      recvcnt[p] = N_dims * N_points_from[p];
      N_points_local += recvcnt[p];
      senddispl[p] = senddispl[p-1] + sendcnt[p-1];
      recvdispl[p] = recvdispl[p-1] + recvcnt[p-1];
    }
    if (N_dims > 0) assert (N_points_local % N_dims == 0);
    // N_points_local is the total number of points to receive
    // and thus the total number of points to interpolate on this processor
    if (N_dims > 0) N_points_local /= N_dims;

    // Set up the per-component coordinates
    // as offset into the single communication send buffer
    vector<CCTK_REAL*> coords (allhomecnts.size());
    int offset = 0;
    for (int c = 0; c < allhomecnts.size(); c++) {
      coords[c] = &coords_buffer.front() + N_dims*offset;
      offset += allhomecnts[c];
    }
    assert (offset == N_interp_points);

    // Copy the input coordinates into the communication send buffer
    // also remember the position of each point in the original input arrays
    vector<int> indices (N_interp_points);
    {
      // totalhomecnts is the accumulated number of points over all components
      vector<int> totalhomecnts (allhomecnts.size());
      for (int idx = 0; idx < allhomecnts.size() - 1; idx++) {
        totalhomecnts[idx + 1] = totalhomecnts[idx] + allhomecnts[idx];
      }

      vector<int> tmpcnts (allhomecnts.size());
      for (int n = 0; n < N_interp_points; n++) {
        int const idx = component_idx (dstprocs[n], source_map[n], rlev[n],
                                       home[n]);
        for (int d = 0; d < N_dims; d++) {
          coords[idx][d + N_dims*tmpcnts[idx]] =
            static_cast<CCTK_REAL const *>(coords_list[d])[n];
        }
        indices[n] = totalhomecnts[idx] + tmpcnts[idx];
        tmpcnts[idx]++;
      }
      assert (tmpcnts == allhomecnts);
    }

    // Send this processor's points to owning processors,
    // receive other processors' points for interpolation on this processor
    {
      vector<CCTK_REAL> tmp (N_dims * N_points_local);
      MPI_Datatype const datatype = dist::datatype (tmp[0]);
      MPI_Alltoallv (&coords_buffer[0], &sendcnt[0], &senddispl[0], datatype,
                     &tmp[0],           &recvcnt[0], &recvdispl[0], datatype,
                     dist::comm());
      coords_buffer.swap (tmp);
    }

    //////////////////////////////////////////////////////////////////////
    // Communicate the source map (if neccessary)
    //////////////////////////////////////////////////////////////////////
    if (have_source_map) {
      vector<CCTK_INT> tmp (N_interp_points);
      vector<int> tmpcnts (N_points_to.size());
      for (int n = 0; n < N_interp_points; n++) {
        int const p = dstprocs[n];
        int const offset = senddispl[p] + tmpcnts[p];
        tmp[offset] = source_map[n];
        tmpcnts[p]++;
      }
      assert (tmpcnts == N_points_to);

      sendcnt[0] = N_points_to[0];
      recvcnt[0] = N_points_from[0];
      senddispl[0] = recvdispl[0] = 0;
      for (int p = 1; p < dist::size(); p++) {
        sendcnt[p] = N_points_to[p];
        recvcnt[p] = N_points_from[p];
        senddispl[p] = senddispl[p-1] + sendcnt[p-1];
        recvdispl[p] = recvdispl[p-1] + recvcnt[p-1];
      }

      source_map.resize (N_points_local);
      MPI_Datatype const datatype = dist::datatype (tmp[0]);
      MPI_Alltoallv (&tmp[0],        &sendcnt[0], &senddispl[0], datatype,
                     &source_map[0], &recvcnt[0], &recvdispl[0], datatype,
                     dist::comm());
    } else {
      // No explicit source map specified
      if (Carpet::map != -1) {
        // Interpolate from the current map
        source_map.resize (N_points_local, Carpet::map);
      } else {
        // Interpolate from map 0 if this is the only one
        // (for backwards compatibility)
        assert (maps == 1);
        source_map.resize (N_points_local, 0);
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Map (local) interpolation points to components
    //////////////////////////////////////////////////////////////////////
    rlev.resize (N_points_local);          // reflevel of (local) point n
    home.resize (N_points_local);          // component of (local) point n
    vector<int> srcprocs (N_points_local); // which processor requested point n
    vector<int> homecnts (allhomecnts.size());  // points per components

    // Remember from where point n came from
    offset = 0;
    for (int p = 0; p < N_points_from.size(); p++) {
      for (int n = 0; n < N_points_from[p]; n++) {
        srcprocs[offset++] = p;
      }
    }
    assert (offset == N_points_local);

    // map each point in coords_buffer onto its component (srcproc, rlev, home)
    // also accumulate the number of points in each component (homecnts)
    //
    // In order to be somewhat more efficient, we could also map
    // all processors' (rlev, home) points onto the same component
    // and call the local interpolator on that (thus saving potentially
    // nprocs-1 calls).
    // The reason why this hasn't been implemented here is because
    // CCTK_InterpGridArrays() is supposed to return a value which is
    // the minimum over all local interpolator invocations' return values
    // for the points requested by this processor. This overall minimum
    // isn't defined anymore if we call the local interpolator on a component
    // now containing points from different processors.
    // (One could argue though that the per-point status codes as returned
    // by the local interpolator could be used to determine the global
    // interpolator return value instead.)
    map_points (cctkGH, coord_system_handle, coord_group,
                ml, minrl, maxrl, maxncomps, N_dims, N_points_local,
                source_map,
                NULL, coords_buffer,
                srcprocs, N_points_to,
                rlev, home, homecnts);

    // Free some memory
    source_map.clear();
    srcprocs.clear();

    // Reorder the coordinates from <N_dims>-tupels into <N_dims> vectors
    // as expected by CCTK_InterpLocalUniform()
    {
      int offset = 0;
      vector<CCTK_REAL> tmp (coords_buffer.size());
      for (int c = 0; c < homecnts.size(); c++) {
        for (int n = 0; n < homecnts[c]; n++) {
          for (int d = 0; d < N_dims; d++) {
            tmp[d*homecnts[c] + n + offset] = coords_buffer[n*N_dims + d + offset];
          }
        }
        offset += N_dims * homecnts[c];
      }
      assert (offset == N_dims * N_points_local);
      coords_buffer.swap (tmp);
    }

    //////////////////////////////////////////////////////////////////////
    // Do the local interpolation on individual components
    //////////////////////////////////////////////////////////////////////
    const int vtype = output_array_type_codes[0];
    const int vtypesize = CCTK_VarTypeSize (vtype);
    assert (vtypesize > 0);
    vector<char>  outputs_buffer (N_points_local * N_output_arrays * vtypesize);
    vector<char*> outputs (homecnts.size());
    vector<CCTK_INT> status_and_retval_buffer (2 * dist::size());
    CCTK_INT* per_proc_statuses = &status_and_retval_buffer.front();
    CCTK_INT* per_proc_retvals  = per_proc_statuses + dist::size();

    // Set up the per-component coordinates and output arrays as offsets
    // into the single communication buffers
    offset = 0;
    for (int c = 0; c < homecnts.size(); c++) {
      coords[c]   = &coords_buffer.front()  + N_dims*offset;
      outputs[c]  = &outputs_buffer.front() + N_output_arrays*offset*vtypesize;
      offset += homecnts[c];
    }
    assert (offset == N_points_local);

    interpolate_components
      (cctkGH, coord_system_handle, coord_group,
       minrl, maxrl,
       maxncomps,
       want_global_mode,
       prolongation_order_time,
       N_dims,
       homecnts,
       coords, outputs, per_proc_statuses, per_proc_retvals,
       operand_indices, time_deriv_order, have_time_derivs,
       local_interp_handle, param_table_handle,
       current_time, delta_time,
       N_input_arrays, N_output_arrays,
       output_array_type_codes, input_array_variable_indices);

    // Free some memory
    coords_buffer.clear();
    coords.clear();
    homecnts.clear();
    home.clear();
    rlev.clear();

    //////////////////////////////////////////////////////////////////////
    // Communicate interpolation results
    //////////////////////////////////////////////////////////////////////
    sendcnt[0] = N_output_arrays * N_points_from[0];
    recvcnt[0] = N_output_arrays * N_points_to[0];
    senddispl[0] = recvdispl[0] = 0;
    for (int p = 1; p < dist::size(); p++) {
      sendcnt[p] = N_output_arrays * N_points_from[p];
      recvcnt[p] = N_output_arrays * N_points_to[p];
      senddispl[p] = senddispl[p-1] + sendcnt[p-1];
      recvdispl[p] = recvdispl[p-1] + recvcnt[p-1];
    }

    {
      MPI_Datatype datatype;
      switch (vtype) {
#define TYPECASE(N,T)                                                     \
        case  N: { T dummy; datatype = dist::datatype(dummy); } break;
#include "carpet_typecase.hh"
#undef TYPECASE
        default: assert (0 and "invalid datatype");
      }

      vector<char> tmp (N_interp_points * N_output_arrays * vtypesize);
      MPI_Alltoallv (&outputs_buffer[0], &sendcnt[0], &senddispl[0], datatype,
                     &tmp[0],            &recvcnt[0], &recvdispl[0], datatype,
                     dist::comm());
      outputs_buffer.swap (tmp);
    }

    //////////////////////////////////////////////////////////////////////
    // Communicate interpolation status codes and return values
    //////////////////////////////////////////////////////////////////////
    {
      // A processor's overall status and return code
      // is defined as the minimum over all local interpolator status and
      // return codes across all processors for that processor
      vector<CCTK_INT> tmp (status_and_retval_buffer.size());
      MPI_Allreduce (&status_and_retval_buffer[0], &tmp[0], tmp.size(),
                     dist::datatype (tmp[0]), MPI_MIN, dist::comm());
      status_and_retval_buffer.swap (tmp);
      per_proc_statuses = &status_and_retval_buffer.front();
      per_proc_retvals  = per_proc_statuses + dist::size();
    }

    //////////////////////////////////////////////////////////////////////
    // Finally, sort the received outputs back into the caller's arrays
    //////////////////////////////////////////////////////////////////////

    // Sorting is done with the help of the inverse indices vector
    vector<int> reverse_indices(indices.size());
    for (int i = 0; i < indices.size(); i++) {
      reverse_indices.at(indices[i]) = i;
    }

    for (int d = 0; d < N_output_arrays; d++) {
      char* output_array = static_cast<char*>(output_arrays[d]);
      for (int c = 0, i = 0, offset = 0; c < allhomecnts.size(); c++) {
        assert (allhomecnts[c]*(d+1) + offset <=
                N_output_arrays*N_interp_points);
        for (int n = 0; n < allhomecnts[c]; n++, i++) {
          memcpy (output_array + reverse_indices[i]*vtypesize,
                  &outputs_buffer.front() +
                    (allhomecnts[c]*d + offset + n) * vtypesize,
                  vtypesize);
        }
        offset += N_output_arrays * allhomecnts[c];
      }
    }

    // set this processor's overall local interpolator status code
    int ierr = Util_TableSetInt (param_table_handle,
                                 per_proc_statuses[dist::rank()],
                                 "local_interpolator_status");
    assert (ierr >= 0);

    // Done.
    return per_proc_retvals[dist::rank()];
  }



  static int
  extract_parameter_table_options (cGH const* const cctkGH,
                                   int const param_table_handle,
                                   int const N_interp_points,
                                   int const N_input_arrays,
                                   int const N_output_arrays,
                                   bool& have_source_map,
                                   bool& have_time_derivs,
                                   int& prolongation_order_time,
                                   CCTK_REAL& current_time,
                                   CCTK_REAL& delta_time,
                                   vector<CCTK_INT>& source_map,
                                   vector<CCTK_INT>& operand_indices,
                                   vector<CCTK_INT>& time_deriv_order)
  {
    // Find source map
    assert (source_map.size() == N_interp_points);
    int iret = Util_TableGetIntArray (param_table_handle, N_interp_points,
                                      &source_map.front(), "source_map");
    have_source_map = not (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    if (not have_source_map) {
      // No explicit source map specified
      if (Carpet::map != -1) {
        // Interpolate from the current map
        source_map.assign (source_map.size(), Carpet::map);
      } else if (maps == 1) {
        // Interpolate from map 0 if this is the only one
        // (for backwards compatibility)
        source_map.assign (source_map.size(), 0);
      } else {
        CCTK_WARN (CCTK_WARN_ALERT, "No source map specified");
        return -1;
      }
    } else if (iret < 0) {
      CCTK_WARN (CCTK_WARN_ALERT, "internal error");
      return -1;
    } else if (iret != N_interp_points) {
      CCTK_WARN (CCTK_WARN_ALERT, "Source map array has wrong size");
      return -1;
    } else {
      iret = Util_TableGetIntArray (param_table_handle, source_map.size(),
                                    &source_map.front(), "source_map");
      assert (iret == source_map.size());
      // Check source map
      for (int n = 0; n < source_map.size(); ++n) {
        assert (source_map[n] >= 0 and source_map[n] < maps);
      }
    }

    // Find the time interpolation order
    int partype;
    void const* const parptr = CCTK_ParameterGet ("prolongation_order_time",
                                                  "Carpet", &partype);
    assert (parptr);
    assert (partype == PARAMETER_INTEGER);
    prolongation_order_time = *(CCTK_INT const*) parptr;

    current_time = cctkGH->cctk_time / cctkGH->cctk_delta_time;
    delta_time = cctkGH->cctk_delta_time;

    iret = Util_TableGetIntArray (param_table_handle, N_output_arrays,
                                 &time_deriv_order.front(), "time_deriv_order");
    have_time_derivs = not (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    if (not have_time_derivs) {
      time_deriv_order.assign (time_deriv_order.size(), 0);
    } else {
      assert (iret == N_output_arrays);
    }

    // Find output variable indices
    iret = Util_TableGetIntArray (param_table_handle, N_output_arrays,
                                  &operand_indices.front(), "operand_indices");
    if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      assert (N_output_arrays == N_input_arrays);
      for (int m = 0; m < N_output_arrays; ++m) {
        operand_indices[m] = m;
      }
    } else {
      assert (iret == N_output_arrays);
    }

    return 0;
  }


  static void
  map_points (cGH const* const cctkGH,
              int const coord_system_handle,
              int const coord_group,
              int const ml,
              int const minrl,
              int const maxrl,
              int const maxncomps,
              int const N_dims,
              int const npoints,
              vector<CCTK_INT> const& source_map,
              void const* const coords_list[],
              vector<CCTK_REAL> const& coords,
              vector<int>& procs,
              vector<int>& N_points_to,
              vector<int>& rlev,
              vector<int>& home,
              vector<int>& homecnts)
  {
    bool const map_onto_processors = coords_list != NULL;

    if (not map_onto_processors) {
      assert (coords.size() == N_dims * npoints);
    }
    assert (procs.size() == npoints);
    assert (N_points_to.size() == dist::size());
    assert (rlev.size() == npoints);
    assert (home.size() == npoints);
    assert (source_map.size() == npoints);


    // Find out about the coordinates: origin and delta
    // for the Carpet grid indices
    vector<rvect> lower (maps);
    vector<rvect> delta (maps);
    vector<rvect> upper (maps);
    rvect coord_lower, coord_upper;
    const char* coord_system_name = CCTK_CoordSystemName (coord_system_handle);

    for (int d = 0; d < N_dims; d++) {
      const int iret = CCTK_CoordRange (cctkGH, &coord_lower[d],
                                        &coord_upper[d], d+1,
                                        NULL, coord_system_name);
      assert (iret == 0);
    }

    for (int m = 0; m < maps; ++m) {
      const ibbox& baseextent = arrdata.at(coord_group).at(m).hh->baseextent;
      lower.at(m) = coord_lower;
      delta.at(m) = (coord_upper - coord_lower) /
                    (baseextent.upper() * maxspacereflevelfact);
      upper.at(m) = lower.at(m) +
                    delta.at(m) * (baseextent.upper() - baseextent.lower());
    }

    // Assign interpolation points to processors/components
    for (int n = 0; n < npoints; ++n) {

      int const m = source_map[n];

      // Find the grid point closest to the interpolation point
      rvect pos;
      for (int d = 0; d < N_dims; ++d) {
        if (map_onto_processors) {
          pos[d] = static_cast<CCTK_REAL const *>(coords_list[d])[n];
        } else {
          pos[d] = coords[N_dims*n + d];
        }
      }

      // Find the component that this grid point belongs to
      int rl = -1, c = -1;
      if (all (pos >= lower.at(m) and pos <= upper.at(m))) {
        for (rl = maxrl-1; rl >= minrl; --rl) {

          ivect const fact = maxspacereflevelfact / spacereffacts.at(rl) *
                             ipow(mgfact, mglevel);
          ivect const ipos = ivect(floor((pos - lower.at(m)) / (delta.at(m) *
                             fact) + 0.5)) * fact;
          const gh* hh = arrdata[coord_group][m].hh;
          assert (all (ipos % hh->bases().at(ml).at(rl).stride() == 0));

          // TODO: use something faster than a linear search
          for (c = 0; c < hh->components(rl); ++c) {
            if (hh->extents().at(ml).at(rl).at(c).contains(ipos)) {
              goto found;
            }
          }
        }
      }
      // Point could not be mapped onto any processor's domain.
      // Map the point (arbitrarily) to the first component of the
      // coarsest grid
      rl = minrl;
      c = 0;
      // Only warn once
      if (map_onto_processors) {
        CCTK_VWarn (CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Interpolation point #%d at [%g,%g,%g] is not on "
                    "any grid patch", n, pos[0], pos[1], pos[2]);
      }
    found:
      assert (rl >= minrl and rl < maxrl);
      assert (c >= 0 and c < arrdata[coord_group][m].hh->components(rl));

      if (map_onto_processors) {
        procs[n] = arrdata[coord_group][m].hh->proc(rl, c);
        ++ N_points_to[procs[n]];
      }
      rlev[n] = rl;
      home[n] = c;
      ++ homecnts.at(component_idx (procs[n], m, rl, c));

    } // for n
  }


  static void
  interpolate_components (cGH const* const cctkGH,
                          int const coord_system_handle,
                          int const coord_group,
                          int const minrl,
                          int const maxrl,
                          int const maxncomps,
                          bool const want_global_mode,
                          int const prolongation_order_time,
                          int const N_dims,
                          vector<int> & homecnts,
                          vector<CCTK_REAL*>& coords,
                          vector<char*>& outputs,
                          CCTK_INT* const per_proc_statuses,
                          CCTK_INT* const per_proc_retvals,
                          vector<CCTK_INT>& operand_indices,
                          vector<CCTK_INT>& time_deriv_order,
                          bool const have_time_derivs,
                          CCTK_INT const local_interp_handle,
                          CCTK_INT const param_table_handle,
                          CCTK_REAL const current_time,
                          CCTK_REAL const delta_time,
                          int const N_input_arrays,
                          int const N_output_arrays,
                          CCTK_INT const output_array_type_codes [],
                          CCTK_INT const input_array_variable_indices [])
  {
    // Find out about the number of time levels for interpolation
    // for all input arrays
    vector<CCTK_INT> interp_num_time_levels (N_input_arrays, 0);
    for (int n = 0; n < N_input_arrays; n++) {
      int const vi = input_array_variable_indices[n];
      if (vi >= 0) {
        int const gi = CCTK_GroupIndexFromVarI (vi);
        assert (gi >= 0 and gi < CCTK_NumGroups ());
        int const table = CCTK_GroupTagsTableI (gi);
        assert (table >= 0);

        Util_TableGetInt (table, &interp_num_time_levels[n],
                          "InterpNumTimelevels");
      }
    }

    BEGIN_GLOBAL_MODE(cctkGH) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        enter_level_mode (const_cast<cGH*>(cctkGH), rl);

        // Number of neccessary time levels
        CCTK_REAL const level_time = cctkGH->cctk_time / cctkGH->cctk_delta_time;
        bool const need_time_interp
          = (have_time_derivs
             or (fabs(current_time - level_time)
                 > 1e-12 * (fabs(level_time) + fabs(current_time)
                            + fabs(cctkGH->cctk_delta_time))));
        assert (! (! want_global_mode
                   and ! have_time_derivs
                   and need_time_interp));
        int const num_tl = need_time_interp ? prolongation_order_time + 1 : 1;

        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {

            for (int p = 0; p < dist::size(); p++) {
              int const idx = component_idx (p, Carpet::map,reflevel,component);
              if (homecnts[idx] > 0) {
                interpolate_single_component
                  (cctkGH, coord_system_handle, coord_group,
                   minrl, maxrl, maxncomps, N_dims,
                   homecnts[idx], coords[idx], outputs[idx],
                   per_proc_statuses[p], per_proc_retvals[p],
                   operand_indices, time_deriv_order, interp_num_time_levels,
                   local_interp_handle, param_table_handle,
                   num_tl, need_time_interp, current_time, delta_time,
                   N_input_arrays, N_output_arrays,
                   output_array_type_codes,
                   input_array_variable_indices);
              }
            }

          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;

        leave_level_mode (const_cast<cGH*>(cctkGH));
      } // for rl

    } END_GLOBAL_MODE;
  }

  ///////////////////////////////////////////////////////////////////////////


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

  /** Represents interpolation coefficients, or weights, for a given
   *  derivative order and set  of time levels.
   *
   * */
  class InterpolationWeights : private vector<CCTK_REAL>
  {
  public:
    InterpolationWeights (CCTK_INT deriv_order,
                          const InterpolationTimes & times,
                          const CCTK_REAL current_time,
                          const CCTK_REAL delta_time)
      : vector<CCTK_REAL> (times.num_timelevels () )
    {
      CCTK_INT num_timelevels = times.num_timelevels ();
      // Initialise weights
      switch (deriv_order) {
      case 0:
        switch (num_timelevels) {
        case 1:
          for_no_interp (times, current_time, delta_time);
          break;
        case 2:
          for_linear_2_pt_interp (times, current_time, delta_time);
          break;
        case 3:
          for_quadratic_3_pt_interp (times, current_time, delta_time);
          break;
        default:
          assert (0);
        }
        break;

      case 1:
        switch (num_timelevels) {
        case 2:
          for_1st_deriv_linear_2_pt_interp_1st_order (times, current_time, delta_time);
          break;
        case 3:
          for_1st_deriv_quadratic_3_pt_interp_2nd_order (times, current_time, delta_time);
          break;
        default:
          assert (0);
        }
        break;

      case 2:
        switch (num_timelevels) {
        case 3:
          for_2nd_deriv_quadratic_3_pt_interp_2nd_order (times, current_time, delta_time);
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
    void for_no_interp (const InterpolationTimes & t,
                        const CCTK_REAL t0, const CCTK_REAL dt)
    {
      // We have to assume that any GF with one timelevel is constant in time
      at(0) = 1.0;
    }

    void for_linear_2_pt_interp (const InterpolationTimes & t,
                                 const CCTK_REAL t0, const CCTK_REAL dt)
    {
      at(0) = (t0 - t[1]) / (t[0] - t[1]);
      at(1) = (t0 - t[0]) / (t[1] - t[0]);
    }

    void for_quadratic_3_pt_interp (const InterpolationTimes & t,
                                    const CCTK_REAL t0, const CCTK_REAL dt)
    {
      at(0) = (t0 - t[1]) * (t0 - t[2]) / ((t[0] - t[1]) * (t[0] - t[2]));
      at(1) = (t0 - t[0]) * (t0 - t[2]) / ((t[1] - t[0]) * (t[1] - t[2]));
      at(2) = (t0 - t[0]) * (t0 - t[1]) / ((t[2] - t[0]) * (t[2] - t[1]));
    }

    void for_1st_deriv_linear_2_pt_interp_1st_order (
                                               const InterpolationTimes & t,
                                               const CCTK_REAL t0, const CCTK_REAL dt)
    {
      at(0) = 1 / (dt * (t[0] - t[1]));
      at(1) = 1 / (dt * (t[1] - t[0]));
    }

    void for_1st_deriv_quadratic_3_pt_interp_2nd_order (
                                               const InterpolationTimes & t,
                                               const CCTK_REAL t0, const CCTK_REAL dt )
    {
      at(0) = ((t0 - t[2]) + (t0 - t[1]))
                  / (dt * (t[0] - t[1]) * (t[0] - t[2]));
      at(1) = ((t0 - t[2]) + (t0 - t[0]))
                  / (dt * (t[1] - t[0]) * (t[1] - t[2]));
      at(2) = ((t0 - t[1]) + (t0 - t[0]))
                  / (dt * (t[2] - t[0]) * (t[2] - t[1]));
    }

    void for_2nd_deriv_quadratic_3_pt_interp_2nd_order (
                                               const InterpolationTimes & t,
                                               const CCTK_REAL t0, const CCTK_REAL dt )
    {
      at(0) = 2 / (dt * dt * (t[0] - t[1]) * (t[0] - t[2]));
      at(1) = 2 / (dt * dt * (t[1] - t[0]) * (t[1] - t[2]));
      at(2) = 2 / (dt * dt * (t[2] - t[0]) * (t[2] - t[1]));
    }
  };


  static void
  interpolate_single_component (cGH const* const cctkGH,
                                int const coord_system_handle,
                                int const coord_group,
                                int const minrl,
                                int const maxrl,
                                int const maxncomps,
                                int const N_dims,
                                int const npoints,
                                CCTK_REAL const* const coords,
                                char* const outputs,
                                CCTK_INT& overall_status,
                                CCTK_INT& overall_retval,
                                vector<CCTK_INT>& operand_indices,
                                vector<CCTK_INT>& time_deriv_order,
                                vector<CCTK_INT>& interp_num_time_levels,
                                CCTK_INT const local_interp_handle,
                                CCTK_INT const param_table_handle,
                                int const num_tl,
                                bool const need_time_interp,
                                CCTK_REAL const current_time,
                                CCTK_REAL const delta_time,
                                int const N_input_arrays,
                                int const N_output_arrays,
                                CCTK_INT const output_array_type_codes [],
                                CCTK_INT const input_array_variable_indices [])
  {
    // Find out the datatypes of all input grid arrays
    vector<CCTK_INT> input_array_type_codes(N_input_arrays, -1);
    vector<const void *> input_arrays(N_input_arrays);
    for (int n = 0; n < N_input_arrays; ++n) {
      if (input_array_variable_indices[n] >= 0) {
        input_array_type_codes[n] =
          CCTK_VarTypeI (input_array_variable_indices[n]);
      }
    }

    // Do the processor-local interpolation
    // Find out about the local geometry
    rvect coord_origin, coord_delta;

    // Get global origin and spacing of the underlying coordinate system
#ifdef NEW_COORD_API
    const iret1 = Util_TableGetRealArray (coord_system_handle, N_dims,
                                          &coord_origin[0], "COMPMIN");
    const iret2 = Util_TableGetRealArray (coord_system_handle, N_dims,
                                          &coord_delta[0], "DELTA");
    assert (iret1 == N_dims and iret2 == N_dims);
#else
    const char* coord_system_name = CCTK_CoordSystemName(coord_system_handle);
    assert (CCTK_CoordSystemDim (coord_system_name) >= N_dims);

    rvect lower, upper;
    for (int d = 0; d < N_dims; d++) {
      const int iret = CCTK_CoordRange (cctkGH, &lower[d], &upper[d], d+1,
                                        NULL, coord_system_name);
      assert (iret == 0);
    }
    const ibbox& baseextent =
      arrdata.at(coord_group).at(Carpet::map).hh->baseextent;
    coord_origin = lower;
    coord_delta  = (upper - lower) / baseextent.upper();
#endif

    // Get processor-local origin and spacing
    cGroupDynamicData coord_group_data;
    CCTK_GroupDynamicData (cctkGH, coord_group, &coord_group_data);
    for (int d = 0; d < N_dims; ++d) {
      coord_delta[d]  /= cctkGH->cctk_levfac[d];
      coord_origin[d] += (coord_group_data.lbnd[d] +
                          cctkGH->cctk_levoff[d] / cctkGH->cctk_levoffdenom[d])
                       * coord_delta[d];
    }

    void const* tmp_coords[dim];
    for (int d = 0; d < N_dims; ++d) {
      tmp_coords[d] = coords + d*npoints;
    }

    vector<vector<void *> > tmp_output_arrays (num_tl);

    for (int tl=0; tl<num_tl; ++tl) {

      for (int n=0; n<N_input_arrays; ++n) {

        int const vi = input_array_variable_indices[n];
        assert (vi >= 0 and vi < CCTK_NumVars());

        if (vi == -1) {
          input_arrays[n] = NULL;
        } else {
          int const interp_num_tl = interp_num_time_levels[n] > 0 ?
                                    interp_num_time_levels[n] : num_tl;

          // Do a dummy interpolation from a later timelevel
          // if the desired timelevel does not exist
          int const my_tl = tl >= interp_num_tl ? 0 : tl;
#if 0
          input_arrays[n] = CCTK_VarDataPtrI (cctkGH, my_tl, vi);
#else
          int const gi = CCTK_GroupIndexFromVarI (vi);
          int const vi0 = CCTK_FirstVarIndexI (gi);
          input_arrays[n]
            = ((*arrdata.at(gi).at(Carpet::map).data.at(vi-vi0))
               (my_tl, reflevel, component, mglevel)->storage());
#endif
        }
      } // for input arrays

      tmp_output_arrays[tl].resize (N_output_arrays);
      for (int j=0; j<N_output_arrays; ++j) {
        if (need_time_interp) {
          if (output_array_type_codes[j] != CCTK_VARIABLE_REAL) {
            CCTK_WARN (CCTK_WARN_ABORT,
                       "time interpolation into output arrays of datatype "
                       "other than CCTK_VARIABLE_REAL is not supported");
          }
          tmp_output_arrays[tl][j] = new CCTK_REAL [npoints];
        } else {
          const int vartypesize = CCTK_VarTypeSize (output_array_type_codes[j]);
          tmp_output_arrays[tl][j] = outputs + j*npoints*vartypesize;
        }
      }

      vector<CCTK_INT> per_point_status (npoints);
      int ierr = Util_TableSetPointer
        (param_table_handle, &per_point_status.front(), "per_point_status");
      assert (ierr >= 0);

      const int retval = CCTK_InterpLocalUniform
        (N_dims, local_interp_handle, param_table_handle,
         &coord_origin[0], &coord_delta[0],
         npoints,
         CCTK_VARIABLE_REAL, tmp_coords,
         N_input_arrays, coord_group_data.lsh,
         &input_array_type_codes.front(), &input_arrays.front(),
         N_output_arrays,
         output_array_type_codes, &tmp_output_arrays[tl].front());
      if (retval) {
        CCTK_VWarn (CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The local interpolator returned the error code %d",retval);
      }

      overall_retval = min (overall_retval, retval);
      for (int n = 0; n < per_point_status.size(); n++) {
        overall_status = min (overall_status, per_point_status[n]);
      }
      ierr = Util_TableDeleteKey (param_table_handle, "per_point_status");
      assert (! ierr);

    } // for tl

    // Interpolate in time, if neccessary
    if (need_time_interp) {

      for (int j=0; j<N_output_arrays; ++j) {

        // Find input array for this output array
        int const n = operand_indices.at(j);
        CCTK_INT const deriv_order = time_deriv_order.at(j);

        int const interp_num_tl = interp_num_time_levels.at(n) > 0 ?
                                  interp_num_time_levels.at(n) : num_tl;
        const InterpolationTimes times (interp_num_tl);
        const InterpolationWeights tfacs (deriv_order, times, current_time,
                                          delta_time);

        // Interpolate
        assert (output_array_type_codes[j] == CCTK_VARIABLE_REAL);
        for (int k=0; k<npoints; ++k) {
          CCTK_REAL& dest =
            static_cast<CCTK_REAL*>(static_cast<void*>(outputs))[k + j*npoints];
          dest = 0;
          for (int tl = 0; tl < interp_num_tl; ++tl) {
            CCTK_REAL const src =
              ((CCTK_REAL const *)tmp_output_arrays[tl][j])[k];
            dest += tfacs[tl] * src;
          }
        }

      } // for j

      for (int tl=0; tl<num_tl; ++tl) {
        for (int j=0; j<N_output_arrays; ++j) {
          delete [] (CCTK_REAL *)tmp_output_arrays[tl][j];
        }
      }

    } // if need_time_interp
  }

} // namespace CarpetInterp
