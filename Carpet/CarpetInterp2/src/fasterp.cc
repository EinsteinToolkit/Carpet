#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>

#include <cctk.h>

#include <mpi.h>

#include <carpet.hh>

#include "fasterp.hh"



namespace CarpetInterp2 {
  
  
  
  // Create an MPI datatype for location information
  MPI_Datatype
  fasterp_iloc_t::mpi_datatype ()
  {
    static bool initialised = false;
    static MPI_Datatype newtype;
    if (not initialised) {
      static fasterp_iloc_t s;
#define ARRSIZE(m) (sizeof(s.m) / sizeof(s.m[0]))
#define OFFSET(m)  ((char*)&(s.m) - (char*)&(s)) // offsetof doesn't work (why?)
#define SIZE       (sizeof(s))
      CCTK_REAL rdummy;
      dist::mpi_struct_descr_t const descr[] = {
        {              1, OFFSET(m     ), MPI_INT               },
        {              1, OFFSET(rl    ), MPI_INT               },
        {              1, OFFSET(c     ), MPI_INT               },
        {              1, OFFSET(ind3d ), MPI_INT               },
        {ARRSIZE(offset), OFFSET(offset), dist::datatype(rdummy)},
        {              1, SIZE          , MPI_UB                }
      };
#undef ARRSIZE
#undef OFFSET
#undef SIZE
      dist::create_mpi_datatype
        (sizeof(descr) / sizeof(descr[0]), descr, newtype);
      initialised = true;
    }
    return newtype;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  // Calculate the coefficients for one interpolation stencil
  void
  fasterp_src_loc_t::
  calc_stencil (fasterp_iloc_t const & iloc,
                int const order)
  {
    assert (order <= max_order);
    CCTK_REAL const rone = 1.0;
    CCTK_REAL const eps = 1.0e-12;
    int const n0 = (order+1) / 2;
    CCTK_REAL const offset = order % 2 == 0 ? 0.5 : 0.0;
    for (int d=0; d<dim; ++d) {
      // C_n = Pi_m[m!=n] [(x_i - x_m) / (x_n - x_m)]
      CCTK_REAL const xi = iloc.offset[d] - offset;
      if (order % 2 == 0) {
        assert (xi >= -rone/2 and xi < +rone/2);
      } else {
        assert (xi >= 0 and xi < rone);
      }
      if (fabs(xi) >= eps) {
        assert (fabs(fmod(xi, rone)) <= eps/2);
        for (int n=0; n<=order; ++n) {
          CCTK_REAL const xn = n - n0;
          CCTK_REAL c = 1.0;
          for (int m=0; m<=order; ++m) {
            if (m != n) {
              CCTK_REAL const xm = m - n0;
              c *= (xi - xm) / (xn - xm);
            }
          }
          coeffs[d][n] = c;
        }
        exact[d] = false;
        // The sum of the coefficients must be one
        CCTK_REAL s = 0.0;
        for (int n=0; n<=order; ++n) {
          s += coeffs[d][n];
        }
        assert (fabs(s - rone) <= eps);
      } else {
        exact[d] = true;
      }
    }
    ind3d = iloc.ind3d;
  }
  
  
  
  // Interpolate a set of variables at the same location, with a given
  // set of interpolation coefficients.  The template mechanism allows
  // the order in each direction to be different; in particular, the
  // order can be 0 if the interpolation location coincides with a
  // source grid point.  For interpolation to order O, this function
  // should be instantiated eight times, with On=O and On=0, for
  // n=[0,1,2].  We hope that the compiler optimises the for loops
  // away if On=0.
  template <int O0, int O1, int O2>
  void
  fasterp_src_loc_t::
  interpolate (ivect const & lsh,
               vector<CCTK_REAL const *> const & varptrs,
               CCTK_REAL * restrict const vals)
    const
  {
    size_t const di = 1;
    size_t const dj = di * lsh[0];
    size_t const dk = dj * lsh[1];
    
    for (size_t v=0; v<varptrs.size(); ++v) {
      vals[v] = 0.0;
    }
    
    for (size_t k=0; k<=O2; ++k) {
      CCTK_REAL const coeff_k = coeffs[2][k];
      for (size_t j=0; j<=O1; ++j) {
        CCTK_REAL const coeff_jk = coeff_k * coeffs[1][j];
        for (size_t i=0; i<=O0; ++i) {
          CCTK_REAL const coeff_ijk = coeff_jk * coeffs[0][i];
          
          for (size_t v=0; v<varptrs.size(); ++v) {
            vals[v] += coeff_ijk * varptrs.AT(v)[ind3d + i*di + j*dj + k*dk];
          } // for v
          
        }
      }
    }
  }
  
  
  
  // Interpolate a set of variables at the same location, with a given
  // set of interpolation coefficients.  This calls the specialised
  // interpolation function, depending on whether interpolation is
  // exact in each of the directions.
  template <int O>
  void
  fasterp_src_loc_t::
  interpolate (ivect const & lsh,
               vector<CCTK_REAL const *> const & varptrs,
               CCTK_REAL * restrict const vals)
    const
  {
    int const Z = 0;
    if (exact[2]) {
      if (exact[1]) {
        if (exact[0]) {
          interpolate<Z,Z,Z> (lsh, varptrs, vals);
        } else {
          interpolate<Z,Z,O> (lsh, varptrs, vals);
        }
      } else {
        if (exact[0]) {
          interpolate<Z,O,Z> (lsh, varptrs, vals);
        } else {
          interpolate<Z,O,O> (lsh, varptrs, vals);
        }
      }
    } else {
      if (exact[1]) {
        if (exact[0]) {
          interpolate<O,Z,Z> (lsh, varptrs, vals);
        } else {
          interpolate<O,Z,O> (lsh, varptrs, vals);
        }
      } else {
        if (exact[0]) {
          interpolate<O,O,Z> (lsh, varptrs, vals);
        } else {
          interpolate<O,O,O> (lsh, varptrs, vals);
        }
      }
    }
  }
  
  
  
  // Interpolate a set of variables at the same location, with a given
  // set of interpolation coefficients.  This calls a specialised
  // interpolation function, depending on whether interpolation is
  // exact in each of the directions.
  void
  fasterp_src_loc_t::
  interpolate (ivect const & lsh,
               int const order,
               vector<CCTK_REAL const *> const & varptrs,
               CCTK_REAL * restrict const vals)
    const
  {
    switch (order) {
    case  0: interpolate< 0> (lsh, varptrs, vals); break;
    case  1: interpolate< 1> (lsh, varptrs, vals); break;
    case  2: interpolate< 2> (lsh, varptrs, vals); break;
    case  3: interpolate< 3> (lsh, varptrs, vals); break;
    case  4: interpolate< 4> (lsh, varptrs, vals); break;
    case  5: interpolate< 5> (lsh, varptrs, vals); break;
    case  6: interpolate< 6> (lsh, varptrs, vals); break;
    case  7: interpolate< 7> (lsh, varptrs, vals); break;
    case  8: interpolate< 8> (lsh, varptrs, vals); break;
    case  9: interpolate< 9> (lsh, varptrs, vals); break;
    case 10: interpolate<10> (lsh, varptrs, vals); break;
    case 11: interpolate<11> (lsh, varptrs, vals); break;
    default:
      // Add higher orders here as desired
      CCTK_WARN (CCTK_WARN_ABORT, "Interpolation orders larger than 11 are not yet implemented");
      assert (0);
    }
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  // Set up an interpolation
  fasterp_setup_t::
  fasterp_setup_t (cGH const * restrict const cctkGH,
                   fasterp_glocs_t const & locations,
                   int const order_)
    : order (order_)
  {
    // Some global properties
    int const npoints = locations.size();
    int const nprocs = CCTK_nProcs (cctkGH);
    // int const myproc = CCTK_MyProc (cctkGH);
    
    
    
    // Calculate patch numbers and local coordinates
    fasterp_llocs_t local_locations (npoints);
    
    if (CCTK_IsFunctionAliased ("MultiPatch_GlobalToLocal")) {
      // This is a muulti-patch simulation: convert global to local
      // coordinates
      
      CCTK_REAL const * coords[dim];
      CCTK_REAL       * local_coords[dim];
      for (int d=0; d<dim; ++d) {
        coords[d]       = &locations.coords[d].front();
        local_coords[d] = &local_locations.coords[d].front();
      }
      
      int const ierr = MultiPatch_GlobalToLocal
        (cctkGH, dim, npoints,
         coords,
         &local_locations.maps.front(), local_coords, NULL, NULL);
      assert (not ierr);
      
    } else {
      // This is a single-patch simulation
      
      // TODO: Optimise this: don't copy coordinates, and don't store
      // map numbers if there is only one map.
      for (int n=0; n<npoints; ++n) {
        int const m = 0;
        local_locations.maps.AT(n) = m;
        for (int d=0; d<dim; ++d) {
          local_locations.coords[d].AT(n) = locations.coords[d].AT(n);
        }
      }
      
    } // if not multi-patch
    
    // Obtain the coordinate ranges for all patches
    vector<rvect> lower (Carpet::maps);
    vector<rvect> upper (Carpet::maps);
    vector<rvect> delta (Carpet::maps); // spacing on finest possible grid
    for (int m=0; m<Carpet::maps; ++m) {
      jvect gsh;
      int const ierr =
        GetCoordRange (cctkGH, m, Carpet::mglevel, dim,
                       & gsh[0],
                       & lower.AT(m)[0], & upper.AT(m)[0], & delta.AT(m)[0]);
      assert (not ierr);
      delta.AT(m) /= Carpet::maxspacereflevelfact;
    }
    
    // Calculate refinement levels, components, and integer grid point
    // indices
    vector<fasterp_iloc_t> ilocs (npoints);
    vector<int> proc (npoints);
    vector<int> nlocs (nprocs, 0);
    int n_nz_nlocs = 0;
    assert (Carpet::is_level_mode());
    for (int n=0; n<npoints; ++n) {
      int const m = local_locations.maps.AT(n);
      rvect const pos (local_locations.coords[0].AT(n),
                       local_locations.coords[1].AT(n),
                       local_locations.coords[2].AT(n));
      
      gh const * const hh = Carpet::vhh.AT(m);
      ibbox const & baseext = hh->baseextent(Carpet::mglevel, 0);
      rvect const rpos =
        (pos - lower.AT(m)) / (upper.AT(m) - lower.AT(m)) *
        rvect (baseext.upper() - baseext.lower());
      
      // Find refinement level and component
      int rl, c;
      ivect ipos;
      hh->locate_position (rpos,
                           Carpet::mglevel,
                           Carpet::reflevel, Carpet::reflevel+1,
                           rl, c, ipos);
      
      ibbox const & ext =
        Carpet::vdd.AT(m)->boxes.AT(Carpet::mglevel).AT(rl).AT(c).exterior;
      
      rvect dpos = rpos - rvect(ipos);
      if (order % 2 != 0) {
        // Potentially shift the stencil anchor for odd interpolation
        // orders (i.e., for even numbers of stencil points)
        ipos -= either (dpos > rvect(0), ext.stride(), ivect(0));
        dpos = rpos - rvect(ipos);
      }
      
      ivect const ind = ipos - ext.lower() / ext.stride();
      ivect const lsh = ext.shape() / ext.stride();
      int const ind3d = ind[0] + lsh[0] * (ind[1] + lsh[1] * ind[2]);
      
      // Store result
      fasterp_iloc_t & iloc = ilocs.AT(n);
      iloc.m      = m;
      iloc.rl     = rl;
      iloc.c      = c;
      iloc.ind3d  = ind3d;
      iloc.offset = dpos;
      
      // Find source processor
      int const p = Carpet::vhh.AT(m)->processor(rl,c);
      proc.AT(n) = p;
      if (nlocs.AT(p) == 0) ++n_nz_nlocs;
      ++ nlocs.AT(p);
    }
    
    
    
    // Find mapping from processors to "processor indices": It may be
    // too expensive to store data for all processors, so we store
    // only data about those processors with which we actually
    // communicate.
    {
      recv_descr.procs.resize (n_nz_nlocs);
      recv_descr.procinds.resize (nprocs, -1);
      int pp = 0;
      int offset = 0;
      for (int p=0; p<nprocs; ++p) {
        if (nlocs.AT(p) > 0) {
          recv_descr.procs.AT(pp).p       = p;
          recv_descr.procs.AT(pp).offset  = offset;
          recv_descr.procs.AT(pp).npoints = nlocs.AT(p);
          recv_descr.procinds.AT(p) = pp;
          ++ pp;
          offset += nlocs.AT(p);
        }
      }
      assert (pp == n_nz_nlocs);
      assert (offset == npoints);
      recv_descr.npoints = npoints;
    }
    
    // Create a mapping "index" from location index, as specified by
    // the user, to the index order in which the data are received
    // from the other processors.
    {
      vector<int> index (recv_descr.procs.size());
      for (int pp=0; pp<int(recv_descr.procs.size()); ++pp) {
        index.AT(pp) = recv_descr.procs.AT(pp).offset;
      }
      recv_descr.index.resize (npoints);
      for (int n=0; n<npoints; ++n) {
        int const p = proc.AT(n);
        int const pp = recv_descr.procinds.AT(p);
        assert (pp >= 0);
        recv_descr.index.AT(n) = index.AT(pp);
        ++index.AT(pp);
      }
      for (int pp=0; pp<int(recv_descr.procs.size()); ++pp) {
        int const recv_npoints = index.AT(pp) - recv_descr.procs.AT(pp).offset;
        assert (recv_npoints == recv_descr.procs.AT(pp).npoints);
      }
    }
    
    
    
    // Count the number of points which have to be sent to other
    // processors, and exchange this information with MPI
    vector<int> send_npoints (nprocs, 0), send_offsets (nprocs);
    {
      int offset = 0;
      for (int pp=0; pp<int(recv_descr.procs.size()); ++pp) {
        int const p = recv_descr.procs.AT(pp).p;
        send_npoints.AT(p) = recv_descr.procs.AT(pp).npoints;
        send_offsets.AT(p) = offset;
        offset += send_npoints.AT(p);
      }
      assert (offset == npoints);
    }
    vector<int> recv_npoints (nprocs), recv_offsets (nprocs);
    MPI_Alltoall (&send_npoints.front(), 1, MPI_INT,
                  &recv_npoints.front(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    int npoints_source = 0;
    for (int p=0; p<nprocs; ++p) {
      recv_offsets.AT(p) = npoints_source;
      npoints_source += recv_npoints.AT(p);
    }
    
    // Scatter the location information into a send-array, and
    // exchange it with MPI
    vector<fasterp_iloc_t> scattered_ilocs(npoints);
    for (int n=0; n<npoints; ++n) {
      scattered_ilocs.AT(recv_descr.index.AT(n)) = ilocs.AT(n);
    }
    vector<fasterp_iloc_t> gathered_ilocs(npoints_source);
    MPI_Comm & comm_world = * (MPI_Comm *) GetMPICommWorld (cctkGH);
    MPI_Alltoallv
      (&scattered_ilocs.front(), &send_npoints.front(), &send_offsets.front(),
       fasterp_iloc_t::mpi_datatype(),
       &gathered_ilocs.front(),  &recv_npoints.front(), &recv_offsets.front(),
       fasterp_iloc_t::mpi_datatype(),
       comm_world);
    
    
    
    // Fill in send descriptors
    send_descr.npoints = npoints_source;
    
    {
      int n_nz_recv_npoints = 0;
      for (int p=0; p<nprocs; ++p) {
        if (recv_npoints.AT(p) > 0) ++n_nz_recv_npoints;
      }
      send_descr.procs.resize (n_nz_recv_npoints);
      int pp = 0;
      for (int p=0; p<nprocs; ++p) {
        if (recv_npoints.AT(p) > 0) {
          send_proc_t & send_proc = send_descr.procs.AT(pp);
          send_proc.p       = p;
          send_proc.offset  = recv_offsets.AT(p);
          send_proc.npoints = recv_npoints.AT(p);
          ++pp;
        }
      }
      assert (pp == n_nz_recv_npoints);
    }
    
    for (size_t pp=0; pp<send_descr.procs.size(); ++pp) {
      send_proc_t & send_proc = send_descr.procs.AT(pp);
      
      send_proc.maps.resize (Carpet::maps);
      for (int m=0; m<Carpet::maps; ++m) {
        send_proc.maps.AT(m).rls.resize (Carpet::reflevels);
        for (int rl=0; rl<Carpet::reflevels; ++rl) {
          int const ncomps = Carpet::vhh.AT(m)->components (rl);
          send_proc.maps.AT(m).rls.AT(rl).compinds.resize (ncomps, -1);
        }
      }
      
      vector<vector<vector<int> > > npoints_comp (Carpet::maps);
      for (int m=0; m<Carpet::maps; ++m) {
        npoints_comp.AT(m).resize(Carpet::reflevels);
        for (int rl=0; rl<Carpet::reflevels; ++rl) {
          int const ncomps = Carpet::vhh.AT(m)->components (rl);
          npoints_comp.AT(m).AT(rl).resize(ncomps, 0);
        }
      }
      
      vector<vector<int> > n_nz_npoints_comp (Carpet::maps);
      for (int m=0; m<Carpet::maps; ++m) {
        n_nz_npoints_comp.AT(m).resize(Carpet::reflevels, 0);
      }
      
      for (int n=0; n<send_proc.npoints; ++n) {
        fasterp_iloc_t const & iloc = gathered_ilocs.AT(send_proc.offset + pp);
        int const m  = iloc.m;
        int const rl = iloc.rl;
        int const c  = iloc.c;
        if (npoints_comp.AT(m).AT(rl).AT(c) == 0) {
          ++ n_nz_npoints_comp.AT(m).AT(rl);
        }
        ++ npoints_comp.AT(m).AT(rl).AT(c);
      }
      
      for (int m=0; m<Carpet::maps; ++m) {
        for (int rl=0; rl<Carpet::reflevels; ++rl) {
          send_proc.maps.AT(m).rls.AT(rl).comps.resize
            (n_nz_npoints_comp.AT(m).AT(rl));
          int cc = 0;
          int const ncomps = Carpet::vhh.AT(m)->components (rl);
          for (int c=0; c<ncomps; ++c) {
            if (npoints_comp.AT(m).AT(rl).AT(c) > 0) {
              send_comp_t & comp = send_proc.maps.AT(m).rls.AT(rl).comps.AT(cc);
              comp.locs.clear();
              comp.locs.reserve (npoints_comp.AT(m).AT(rl).AT(c));
              comp.c = c;
              send_proc.maps.AT(m).rls.AT(rl).compinds.AT(c) = cc;
              ++cc;
            }
          }
        }
      }
      
    } // for pp
    
    
    
    // Calculate stencil coefficients
    for (size_t pp=0; pp<send_descr.procs.size(); ++pp) {
      send_proc_t & send_proc = send_descr.procs.AT(pp);
      for (int n=0; n<send_proc.npoints; ++n) {
        fasterp_iloc_t const & iloc = gathered_ilocs.AT(send_proc.offset+pp);
        int const m  = iloc.m;
        int const rl = iloc.rl;
        int const c  = iloc.c;
        
        fasterp_src_loc_t sloc;
        sloc.calc_stencil (iloc, order);
        
        int const cc = send_proc.maps.AT(m).rls.AT(rl).compinds.AT(c);
        send_comp_t & comp = send_proc.maps.AT(m).rls.AT(rl).comps.AT(cc);
        comp.locs.push_back (sloc);
      }
    } // for pp
  }
  
  
  
  // Free the setup for one interpolation
  fasterp_setup_t::
  ~fasterp_setup_t ()
  {
    // do nothing -- C++ calls destructors automatically
  }
  
  
  
  // Interpolate
  void 
  fasterp_setup_t::
  interpolate (cGH const * restrict const cctkGH,
               vector<int> const & varinds,
               vector<CCTK_REAL *> & values)
    const
  {
    size_t const nvars = varinds.size();
    assert (values.size() == nvars);
    for (size_t v=0; v<values.size(); ++v) {
      assert (varinds.AT(v) >= 0 and varinds.AT(v) <= CCTK_NumVars());
      // Ensure that variable is GF
      // Ensure that variable has "good" type
      // Ensure that variable has storage
      assert (values.AT(v) != NULL);
    }
    
    MPI_Comm & comm_world = * (MPI_Comm *) GetMPICommWorld (cctkGH);
    
    // Post Irecvs
    vector<CCTK_REAL> recv_points (recv_descr.npoints * nvars);
    vector<MPI_Request> recv_reqs (recv_descr.procs.size());
    for (size_t pp=0; pp<recv_descr.procs.size(); ++pp) {
      recv_proc_t const & recv_proc = recv_descr.procs.AT(pp);
      
      CCTK_REAL rdummy;
      int const tag = 0;
      MPI_Irecv (& recv_points.AT(recv_proc.offset * nvars),
                 recv_proc.npoints * nvars,
                 dist::datatype (rdummy), recv_proc.p, tag,
                 comm_world, & recv_reqs.AT(pp));
    }
    
    // Interpolate data and post Isends
    vector<CCTK_REAL> send_points (send_descr.npoints);
    vector<CCTK_REAL>::iterator destit = send_points.begin();
    vector<MPI_Request> send_reqs (send_descr.procs.size());
    for (size_t pp=0; pp<send_descr.procs.size(); ++pp) {
      send_proc_t const & send_proc = send_descr.procs.AT(pp);
      
      for (size_t m=0; m<send_proc.maps.size(); ++m) {
        send_map_t const & send_map = send_proc.maps.AT(m);
        for (size_t rl=0; rl<send_map.rls.size(); ++rl) {
          send_rl_t const & send_rl = send_map.rls.AT(rl);
          for (size_t cc=0; cc<send_rl.comps.size(); ++cc) {
            send_comp_t const & send_comp = send_rl.comps.AT(cc);
            int const c = send_comp.c;
            
            dh const & dd = * Carpet::vdd.AT(m);
            ibbox const &
              ext = dd.boxes.AT(Carpet::mglevel).AT(rl).AT(c).exterior;
            ivect const lsh = ext.shape() / ext.stride();
            
            // Get pointers to 3D data
            vector<CCTK_REAL const *> varptrs (nvars);
            for (size_t v=0; v<nvars; ++v) {
              int const gi = CCTK_GroupIndexFromVarI (varinds.AT(v));
              assert (gi >= 0);
              int const vi = varinds.AT(v) - CCTK_FirstVarIndexI (gi);
              assert (vi >= 0);
              int const tl = 0;
              varptrs.AT(v) =
                (CCTK_REAL const *)
                (* Carpet::arrdata.AT(gi).AT(m).data.AT(vi))
                (tl, rl, c, Carpet::mglevel)->storage();
              assert (varptrs.AT(v));
            }
            
            for (size_t n=0; n<send_comp.locs.size(); ++n) {
              assert (destit + nvars <= send_points.end());
              send_comp.locs.AT(n).interpolate
                (lsh, order, varptrs, & * destit);
              destit += nvars;
            }
            
          } // for cc
        }   // for rl
      }     // for m
      
      CCTK_REAL rdummy;
      int const tag = 0;
      MPI_Isend (& send_points.AT(send_proc.offset * nvars),
                 send_proc.npoints * nvars,
                 dist::datatype (rdummy), send_proc.p, tag,
                 comm_world, & send_reqs.AT(pp));
    } // for pp
    assert (destit == send_points.end());
    
    // Wait for Irecvs to complete
    MPI_Waitall (recv_reqs.size(), & recv_reqs.front(), MPI_STATUSES_IGNORE);
    
    // Gather data
    for (int n=0; n<recv_descr.npoints; ++n) {
      int const nn = recv_descr.index.AT(n);
      for (size_t v=0; v<nvars; ++v) {
        values.AT(v)[n] = recv_points.AT(nn);
      }
    }
    
    // Wait for Isends to complete
    MPI_Waitall (send_reqs.size(), & send_reqs.front(), MPI_STATUSES_IGNORE);
  }
  
  
  
} // namespace CarpetInterp2
