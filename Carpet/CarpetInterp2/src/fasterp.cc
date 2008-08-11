#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <mpi.h>

#include <carpet.hh>
#include <vect.hh>

#include "fasterp.hh"



namespace CarpetInterp2 {
  
  
  
  // Erik's gdb cannot print local variables in functions where an
  // assert fails.  Hence the calls to assert are temporarily moved
  // into a function of their own.
  void ASSERT (bool const cond)
  {
    if (cond) return;
    abort();
  }
  
  
  
  int const ipoison = -1234567890;
  CCTK_REAL const poison = -1.0e+12;
  
  template <typename T> T get_poison () { abort(); }
  template<> int get_poison () { return ipoison; }
  template<> CCTK_REAL get_poison () { return poison; }
  
  template <typename T>
  void fill (vector<T> & v, T const & val)
  {
    // fill (v.begin(), v.end(), val);
#pragma omp parallel for
    for (int n=0; n<int(v.size()); ++n) {
      v.AT(n) = val;
    }
  }
  
  template <typename T>
  void fill_with_poison (vector<T> & v)
  {
#ifndef NDEBUG
    fill (v, get_poison<T>());
    // fill (v.begin(), v.end(), get_poison<T>());
#endif
  }
  
  
  
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
      dist::mpi_struct_descr_t const descr[] = {
        {              3, OFFSET(mrc   ), MPI_INT                    },
        {              1, OFFSET(ind3d ), MPI_INT                    },
        {ARRSIZE(offset), OFFSET(offset), dist::datatype<CCTK_REAL>()},
        {              1, SIZE          , MPI_UB                     }
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
  
  void
  fasterp_iloc_t::
  output (ostream& os)
    const
  {
    os << "fasterp_iloc_t{"
       << "mrc=" << mrc << ","
       << "ind3d=" << ind3d << ","
       << "offset=" << offset << "}";
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  int mrc_t::maps       = -1;
  int mrc_t::reflevels  = -1;
  int mrc_t::components = -1;
  
  void
  mrc_t::
  determine_mrc_info ()
  {
    maps = Carpet::maps;
    reflevels = Carpet::reflevels;
    components = 0;
    for (int m=0; m<maps; ++m) {
      for (int rl=0; rl<reflevels; ++rl) {
        int const ncomps = Carpet::vhh.AT(m)->components (rl);
        components = max (components, ncomps);
      }
    }
  }
  
  mrc_t::
  mrc_t (int ind)
  {
    ASSERT (ind >= 0);
    int const ind0 = ind;
    c = ind % components;
    ind /= components;
    rl = ind % reflevels;
    ind /= reflevels;
    m = ind % maps;
    ind /= maps;
    ASSERT (ind == 0);
    ASSERT (get_ind() == ind0);
  }
  
  void
  mrc_t::
  output (ostream& os)
    const
  {
    os << "mrc{m=" << m << ",rl=" << rl << ",c=" << c << "}";
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  // Calculate the coefficients for one interpolation stencil
  // TODO: Could templatify this function on the order to improve
  // efficiency
  void
  fasterp_src_loc_t::
  calc_stencil (fasterp_iloc_t const & iloc,
                ivect const & lsh,
                int const order)
  {
    ASSERT (order <= max_order);
    CCTK_REAL const eps = 1.0e-12;
    
#ifndef NDEBUG
    // Poison all coefficients
    for (int d=0; d<dim; ++d) {
      for (int n=0; n<=max_order; ++n) {
        coeffs[d][n] = poison;
      }
    }
#endif
    
    if (order == 0) {
      for (int d=0; d<dim; ++d) {
        exact[d] = true;
      }
      ind3d = iloc.ind3d;
      return;
    }
    
    // Choose stencil anchor (first stencil point)
    ivect iorigin;
    if (order % 2 == 0) {
      iorigin = - ivect(order/2);
    } else {
      // Potentially shift the stencil anchor for odd interpolation
      // orders (i.e., for even numbers of stencil points)
      ivect const ioffset (iloc.offset < 0.0);
      iorigin = - ivect((order-1)/2) - ioffset;
    }
    rvect const offset = iloc.offset - rvect(iorigin);
    // Ensure that the stencil is centred
    ASSERT (all (offset >= 0.5*(order-1) and offset < 0.5*(order+1)));
    
    for (int d=0; d<dim; ++d) {
      // C_n = PRODUCT_m,m!=n [(x - x_m) / (x_n - x_m)]
      CCTK_REAL const x = offset[d];
      // round is not available with PGI compilers
      // CCTK_REAL const rx = round(x);
      CCTK_REAL const rx = floor(x+0.5);
      if (abs(x - rx) < eps) {
        // The interpolation point coincides with a grid point; no
        // interpolation is necessary (this is a special case)
        iorigin[d] += int(rx);
        exact[d] = true;
      } else {
        for (int n=0; n<=order; ++n) {
          CCTK_REAL const xn = n;
          CCTK_REAL c = 1.0;
          for (int m=0; m<=order; ++m) {
            if (m != n) {
              CCTK_REAL const xm = m;
              c *= (x - xm) / (xn - xm);
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
        ASSERT (abs(s - 1.0) <= eps);
      }
    }
    
    // Set 3D location of stencil anchor
    ind3d = iloc.ind3d + index(lsh, iorigin);
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
    ASSERT (O0 == 0 or not exact[0]);
    ASSERT (O1 == 0 or not exact[1]);
    ASSERT (O2 == 0 or not exact[2]);
    
    size_t const di = 1;
    size_t const dj = di * lsh[0];
    size_t const dk = dj * lsh[1];
    
    for (size_t v=0; v<varptrs.size(); ++v) {
      vals[v] = 0.0;
    }
    
    for (size_t k=0; k<=O2; ++k) {
      ASSERT (O2 == 0 or coeffs[2][k] != poison);
      CCTK_REAL const coeff_k = (O2==0 ? 1.0 : coeffs[2][k]);
      for (size_t j=0; j<=O1; ++j) {
        ASSERT (O1 == 0 or coeffs[1][j] != poison);
        CCTK_REAL const coeff_jk = coeff_k * (O1==0 ? 1.0 : coeffs[1][j]);
        for (size_t i=0; i<=O0; ++i) {
          ASSERT (O0 == 0 or coeffs[0][i] != poison);
          CCTK_REAL const coeff_ijk = coeff_jk * (O0==0 ? 1.0 : coeffs[0][i]);
          
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
          interpolate<O,Z,Z> (lsh, varptrs, vals);
        }
      } else {
        if (exact[0]) {
          interpolate<Z,O,Z> (lsh, varptrs, vals);
        } else {
          interpolate<O,O,Z> (lsh, varptrs, vals);
        }
      }
    } else {
      if (exact[1]) {
        if (exact[0]) {
          interpolate<Z,Z,O> (lsh, varptrs, vals);
        } else {
          interpolate<O,Z,O> (lsh, varptrs, vals);
        }
      } else {
        if (exact[0]) {
          interpolate<Z,O,O> (lsh, varptrs, vals);
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
      ASSERT (0);
    }
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  // Set up an interpolation starting from global coordinates
  fasterp_setup_t::
  fasterp_setup_t (cGH const * restrict const cctkGH,
                   fasterp_glocs_t const & locations,
                   int const order_)
    : order (order_)
  {
    // Some global properties
    int const npoints = locations.size();
    
    
    
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
      ASSERT (not ierr);
      
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
    
    setup (cctkGH, local_locations);
  }
  
  
  
  // Set up an interpolation starting from local coordinates
  fasterp_setup_t::
  fasterp_setup_t (cGH const * restrict const cctkGH,
                   fasterp_llocs_t const & locations,
                   int const order_)
    : order (order_)
  {
    setup (cctkGH, locations);
  }
  
  
  
  // Helper for setting up an interpolation
  void
  fasterp_setup_t::
  setup (cGH const * restrict const cctkGH,
         fasterp_llocs_t const & locations)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose) CCTK_VInfo (CCTK_THORNSTRING,
                             "Setting up interpolation for %d grid points",
                             int(locations.size()));
    
    // Some global properties
    int const npoints = locations.size();
    int const nprocs = CCTK_nProcs (cctkGH);
    // int const myproc = CCTK_MyProc (cctkGH);
    
    mrc_t::determine_mrc_info();
    int const maxmrc = mrc_t::get_max_ind();
    
    MPI_Comm & comm_world = * (MPI_Comm *) GetMPICommWorld (cctkGH);
    
    
    
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
      ASSERT (not ierr);
      delta.AT(m) /= Carpet::maxspacereflevelfact;
    }
    
    // Calculate refinement levels, components, and integer grid point
    // indices
    if (verbose) CCTK_INFO ("Mapping points onto components");
    vector<fasterp_iloc_t> ilocs (npoints);
    vector<int> proc (npoints);
    fill_with_poison (proc);
    vector<int> nlocs (nprocs, 0);
    int n_nz_nlocs = 0;
    ASSERT (Carpet::is_level_mode());
    for (int n=0; n<npoints; ++n) {
      int const m = locations.maps.AT(n);
      rvect const pos (locations.coords[0].AT(n),
                       locations.coords[1].AT(n),
                       locations.coords[2].AT(n));
      
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
      ASSERT (rl>=0 and c>=0);
      
      ibbox const & ext =
        Carpet::vdd.AT(m)->boxes.AT(Carpet::mglevel).AT(rl).AT(c).exterior;
      rvect dpos = rpos - rvect(ipos);
      
      // Convert from Carpet indexing to grid point indexing
      ASSERT (all (ipos % ext.stride() == ivect(0)));
      ipos /= ext.stride();
      dpos /= rvect(ext.stride());
      ASSERT (all (abs(dpos) <= rvect(0.5)));
      
      ivect const ind = ipos - ext.lower() / ext.stride();
      ivect const lsh = ext.shape() / ext.stride();
      int const ind3d = index (lsh, ind);
      // TODO: ASSERT that there are enough ghost zones
      
      // Store result
      fasterp_iloc_t & iloc = ilocs.AT(n);
      iloc.mrc    = mrc_t(m, rl, c);
      iloc.ind3d  = ind3d;
      iloc.offset = dpos;
      
      // Find source processor
      int const p = Carpet::vhh.AT(m)->processor(rl,c);
      proc.AT(n) = p;
      if (nlocs.AT(p) == 0) ++n_nz_nlocs;
      ++ nlocs.AT(p);
      
      // Output
      if (veryverbose) {
        cout << "Point #" << n << " at " << pos << ": iloc " << iloc << endl;
      }
    }
    
    
    
    // Find mapping from processors to "processor indices": It may be
    // too expensive to store data for all processors, so we store
    // only data about those processors with which we actually
    // communicate.
    if (verbose) CCTK_INFO ("Determine set of communicating processors");
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
      ASSERT (pp == n_nz_nlocs);
      ASSERT (offset == npoints);
      recv_descr.npoints = npoints;
    }
    
    // Create a mapping "index" from location index, as specified by
    // the user, to the index order in which the data are received
    // from the other processors.
    if (verbose) {
      CCTK_INFO ("Determine inter-processor gather index permutation");
    }
    {
      vector<int> index (recv_descr.procs.size());
      fill_with_poison (index);
      for (int pp=0; pp<int(recv_descr.procs.size()); ++pp) {
        index.AT(pp) = recv_descr.procs.AT(pp).offset;
      }
      recv_descr.index.resize (npoints);
      for (int n=0; n<npoints; ++n) {
        int const p = proc.AT(n);
        int const pp = recv_descr.procinds.AT(p);
        ASSERT (pp >= 0);
        recv_descr.index.AT(n) = index.AT(pp);
        ++index.AT(pp);
      }
      for (int pp=0; pp<int(recv_descr.procs.size()); ++pp) {
        int const recv_npoints = index.AT(pp) - recv_descr.procs.AT(pp).offset;
        ASSERT (recv_npoints == recv_descr.procs.AT(pp).npoints);
      }
#ifndef NDEBUG
      vector<bool> received (npoints);
      fill (received, false);
#pragma omp parallel for
      for (int n=0; n<npoints; ++n) {
        assert (not received.AT(n));
        received.AT(n) = true;
      }
#pragma omp parallel for
      for (int n=0; n<npoints; ++n) {
        assert (received.AT(n));
      }
#endif
    }
    
    
    
    // Count the number of points which have to be sent to other
    // processors, and exchange this information with MPI
    if (verbose) {
      CCTK_INFO ("Count and exchange number of communicated grid points");
    }
    vector<int> recv_npoints (nprocs), recv_offsets (nprocs);
    fill (recv_npoints, 0);
    fill_with_poison (recv_offsets);
    {
      int offset = 0;
      for (int pp=0; pp<int(recv_descr.procs.size()); ++pp) {
        int const p = recv_descr.procs.AT(pp).p;
        recv_npoints.AT(p) = recv_descr.procs.AT(pp).npoints;
        recv_offsets.AT(p) = offset;
        offset += recv_npoints.AT(p);
      }
      ASSERT (offset == npoints);
    }
    vector<int> send_npoints (nprocs), send_offsets (nprocs);
    fill_with_poison (send_npoints);
    fill_with_poison (send_offsets);
    MPI_Alltoall (&recv_npoints.front(), 1, MPI_INT,
                  &send_npoints.front(), 1, MPI_INT,
                  comm_world);
    int npoints_send = 0;
    for (int p=0; p<nprocs; ++p) {
      send_offsets.AT(p) = npoints_send;
      npoints_send += send_npoints.AT(p);
    }
    
    // Scatter the location information into a send-array, and
    // exchange it with MPI
    vector<fasterp_iloc_t> scattered_ilocs(npoints);
#pragma omp parallel for
    for (int n=0; n<npoints; ++n) {
      scattered_ilocs.AT(recv_descr.index.AT(n)) = ilocs.AT(n);
    }
    vector<fasterp_iloc_t> gathered_ilocs(npoints_send);
    MPI_Alltoallv
      (&scattered_ilocs.front(), &recv_npoints.front(), &recv_offsets.front(),
       fasterp_iloc_t::mpi_datatype(),
       &gathered_ilocs.front(),  &send_npoints.front(), &send_offsets.front(),
       fasterp_iloc_t::mpi_datatype(),
       comm_world);
    
    
    
    // Fill in send descriptors
    send_descr.npoints = npoints_send;
    
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
          send_proc.offset  = send_offsets.AT(p);
          send_proc.npoints = send_npoints.AT(p);
          ++pp;
        }
      }
      ASSERT (pp == n_nz_recv_npoints);
    }
    
    
    
    // Calculate stencil coefficients
    if (verbose) CCTK_INFO ("Calculate stencil coefficients");
    
    for (size_t pp=0; pp<send_descr.procs.size(); ++pp) {
      send_proc_t & send_proc = send_descr.procs.AT(pp);
      
      vector<int> mrc2comp (maxmrc, -1);
      vector<int> comp2mrc (maxmrc);
      fill_with_poison (comp2mrc);
      int comps = 0;
      
      vector<int> npoints_comp (maxmrc);
      fill (npoints_comp, 0);
      
      // TODO: parallelise with OpenMP
      for (int n=0; n<send_proc.npoints; ++n) {
        fasterp_iloc_t const & iloc = gathered_ilocs.AT(send_proc.offset + n);
        int const mrc = iloc.mrc.get_ind();
        if (mrc2comp.AT(mrc) == -1) {
          mrc2comp.AT(mrc) = comps;
          comp2mrc.AT(comps) = mrc;
          ++ comps;
        }
        ++ npoints_comp.AT(mrc);
      }
      send_proc.comps.resize(comps);
      
      int offset = 0;
      for (int comp=0; comp<comps; ++comp) {
        send_comp_t & send_comp = send_proc.comps.AT(comp);
        int const mrc = comp2mrc.AT(comp);
        send_comp.mrc = mrc;
        send_comp.locs.reserve (npoints_comp.AT(mrc));
        
        mrc_t const themrc (mrc);
        int const m  = themrc.m;
        int const rl = themrc.rl;
        int const c  = themrc.c;
        ibbox const & ext =
          Carpet::vdd.AT(m)->boxes.AT(Carpet::mglevel).AT(rl).AT(c).exterior;
        send_comp.lsh = ext.shape() / ext.stride();
        
        send_comp.offset = offset;
        send_comp.npoints = npoints_comp.AT(mrc);
        offset += send_comp.npoints;
      }
      ASSERT (offset == send_proc.npoints);
      
      send_proc.index.resize (send_proc.npoints);
      fill_with_poison (send_proc.index);
#pragma omp parallel for
      for (int n=0; n<send_proc.npoints; ++n) {
        fasterp_iloc_t const & iloc = gathered_ilocs.AT(send_proc.offset + n);
        int const mrc = iloc.mrc.get_ind();
        int const comp = mrc2comp.AT(mrc);
        send_comp_t & send_comp = send_proc.comps.AT(comp);
        
        send_proc.index.AT(n) = send_comp.offset + send_comp.locs.size();
        
        fasterp_src_loc_t sloc;
        sloc.calc_stencil (iloc, send_comp.lsh, order);
#pragma omp critical
        send_comp.locs.push_back (sloc);
      }
      
      for (int comp=0; comp<comps; ++comp) {
        send_comp_t & send_comp = send_proc.comps.AT(comp);
        assert (int(send_comp.locs.size()) == send_comp.npoints);
      }
      
    } // for pp
    
    if (verbose) CCTK_INFO ("Done.");
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
    DECLARE_CCTK_PARAMETERS;
    
    size_t const nvars = varinds.size();
    if (verbose) CCTK_VInfo (CCTK_THORNSTRING,
                             "Interpolating %d variables", int(nvars));
    ASSERT (values.size() == nvars);
    for (size_t v=0; v<values.size(); ++v) {
      int const vi = varinds.AT(v);
      ASSERT (vi >= 0 and vi <= CCTK_NumVars());
      int const gi = CCTK_GroupIndexFromVarI (vi);
      ASSERT (gi >= 0);
      cGroup group;
      {
        int const ierr = CCTK_GroupData (gi, &group);
        ASSERT (not ierr);
      }
      ASSERT (group.grouptype == CCTK_GF);
      ASSERT (group.vartype == CCTK_VARIABLE_REAL);
      ASSERT (group.dim == dim);
#if 0
      // Cannot be called in level mode
      cGroupDynamicData dyndata;
      {
        int const ierr = CCTK_GroupDynamicData (cctkGH, gi, &dyndata);
        ASSERT (not ierr);
      }
      ASSERT (dyndata.activetimelevels >= 1);
#endif
      ASSERT (values.AT(v) != NULL);
    }
    
    MPI_Comm & comm_world = * (MPI_Comm *) GetMPICommWorld (cctkGH);
    int const mpi_tag = 0;
    
    // Post Irecvs
    if (verbose) CCTK_INFO ("Posting MPI_Irecvs");
    vector<CCTK_REAL> recv_points (recv_descr.npoints * nvars);
    fill_with_poison (recv_points);
    vector<MPI_Request> recv_reqs (recv_descr.procs.size());
    for (size_t pp=0; pp<recv_descr.procs.size(); ++pp) {
      recv_proc_t const & recv_proc = recv_descr.procs.AT(pp);
      
      MPI_Irecv (& recv_points.AT(recv_proc.offset * nvars),
                 recv_proc.npoints * nvars,
                 dist::datatype<CCTK_REAL>(), recv_proc.p, mpi_tag,
                 comm_world, & recv_reqs.AT(pp));
    }
    
    // Interpolate data and post Isends
    if (verbose) CCTK_INFO ("Interpolating and posting MPI_Isends");
    vector<CCTK_REAL> send_points (send_descr.npoints * nvars);
    fill_with_poison (send_points);
    vector<MPI_Request> send_reqs (send_descr.procs.size());
    for (size_t pp=0; pp<send_descr.procs.size(); ++pp) {
      send_proc_t const & send_proc = send_descr.procs.AT(pp);
      
      vector<CCTK_REAL>::iterator
        destit = send_points.begin() + send_proc.offset * nvars;
      vector<CCTK_REAL>::iterator const
        endit = destit + send_proc.npoints * nvars;
      
      for (size_t comp=0; comp<send_proc.comps.size(); ++comp) {
        send_comp_t const & send_comp = send_proc.comps.AT(comp);
        
        int const m  = send_comp.mrc.m;
        int const rl = send_comp.mrc.rl;
        int const c  = send_comp.mrc.c;
        int const tl = 0;
        
        // Get pointers to 3D data
        vector<CCTK_REAL const *> varptrs (nvars);
        for (size_t v=0; v<nvars; ++v) {
          int const gi = CCTK_GroupIndexFromVarI (varinds.AT(v));
          ASSERT (gi >= 0);
          int const vi = varinds.AT(v) - CCTK_FirstVarIndexI (gi);
          ASSERT (vi >= 0);
          varptrs.AT(v) =
            (CCTK_REAL const *)
            (* Carpet::arrdata.AT(gi).AT(m).data.AT(vi))
            (tl, rl, c, Carpet::mglevel)->storage();
          ASSERT (varptrs.AT(v));
        }
        
#pragma omp parallel for
        for (int n=0; n<int(send_comp.locs.size()); ++n) {
          vector<CCTK_REAL>::iterator const destit2 = destit + n * nvars;
          ASSERT (destit2 + nvars <= endit);
          send_comp.locs.AT(n).interpolate
            (send_comp.lsh, order, varptrs, & * destit2);
        }
        destit += nvars * send_comp.locs.size();
        
      } // for comps
      
      ASSERT (destit == endit);
      
      // Gather send points
      vector<CCTK_REAL> gathered_send_points (send_proc.npoints * nvars);
      fill_with_poison (gathered_send_points);
#pragma omp parallel for
      for (int n=0; n<send_proc.npoints; ++n) {
        int const nn = send_proc.offset + send_proc.index.AT(n);
        for (size_t v=0; v<nvars; ++v) {
          gathered_send_points.AT(n * nvars + v) =
            send_points.AT(nn * nvars + v);
        }
      }
#pragma omp parallel for
      for (int n=0; n<int(send_proc.npoints * nvars); ++n) {
        int const nn = send_proc.offset * nvars + n;
        send_points.AT(nn) = gathered_send_points.AT(n);
      }
      // TODO: Use a scatter index list instead of a gather index
      // list, and combine the scattering with the calculation
      
      MPI_Isend (& send_points.AT(send_proc.offset * nvars),
                 send_proc.npoints * nvars,
                 dist::datatype<CCTK_REAL>(), send_proc.p, mpi_tag,
                 comm_world, & send_reqs.AT(pp));
    } // for pp
    
    // Wait for Irecvs to complete
    if (verbose) CCTK_INFO ("Waiting for MPI_Irevcs to complete");
    MPI_Waitall (recv_reqs.size(), & recv_reqs.front(), MPI_STATUSES_IGNORE);
    
    // Gather data
    if (verbose) CCTK_INFO ("Gathering data");
#pragma omp parallel for
    for (int n=0; n<recv_descr.npoints; ++n) {
      size_t const nn = recv_descr.index.AT(n);
      for (size_t v=0; v<nvars; ++v) {
        values.AT(v)[n] = recv_points.AT(nn * nvars + v);
#ifndef NDEBUG
        recv_points.AT(nn * nvars + v) = poison;
#endif
      }
    }
    
    // Wait for Isends to complete
    if (verbose) CCTK_INFO ("Waiting for MPI_Isends to complete");
    MPI_Waitall (send_reqs.size(), & send_reqs.front(), MPI_STATUSES_IGNORE);
    
    if (verbose) CCTK_INFO ("Done.");
  }
  
  
  
} // namespace CarpetInterp2
