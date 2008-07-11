#ifndef FASTERP_HH
#define FASTERP_HH

#include <cstdlib>
#include <vector>

#include <cctk.h>

#include <defs.hh>
#include <dist.hh>
#include <typeprops.hh>
#include <vect.hh>



namespace CarpetInterp2 {
  
  using namespace std;
  
  
  
  int const dim = 3;
  
  // An interpolation point descriptor requires (3 * (max_order+1) +
  // 1) double precision values of memory
  int const max_order = 5;
  
  
  
  // Map C structures to MPI datatypes
  struct mpi_struct_descr_t {
    int          blocklength;
    MPI_Aint     displacement;
    MPI_Datatype type;
  };
  
  void
  create_mpi_datatype (size_t count,
                       mpi_struct_descr_t const descr[],
                       MPI_Datatype & newtype);
  
  
  
  // A global location, given by its global coordinates
  struct fasterp_glocs_t {
    vector<CCTK_REAL> coords[dim];
    fasterp_glocs_t (size_t const n)
    {
      for (int d=0; d<dim; ++d) {
        coords[d].resize(n);
      }
    }
    size_t size () const { return coords[0].size(); }
  };
  
  // A local location, given by map and local coordinates
  struct fasterp_llocs_t {
    vector<int> maps;
    vector<CCTK_REAL> coords[dim];
    fasterp_llocs_t (size_t const n)
    {
      maps.resize(n);
      for (int d=0; d<dim; ++d) {
        coords[d].resize(n);
      }
    }
    size_t size () const { return maps.size(); }
  };
  
  // An integer location, given by map, refinement level, and
  // component
  struct fasterp_iloc_t {
    int m, rl, c;
    
    // ivect idx;
    int ind3d;
    rvect offset;               // in terms of grid points
    
    static MPI_Datatype mpi_datatype ();
  };
  
  
  
  struct fasterp_dest_loc_t {
    // ivect idx;
    int ind3d;                  // destination grid point index
  };
  
  class fasterp_src_loc_t {
    CCTK_REAL coeffs[dim][max_order+1]; // interpolation coefficients
    bvect exact;
    
    // ivect idx;
    int ind3d;                  // source grid point offset
    
  public:
    void
    calc_stencil (fasterp_iloc_t const & iloc,
                  int order);
    void
    interpolate (ivect const & lsh,
                 int order,
                 vector<CCTK_REAL const *> const & varptrs,
                 CCTK_REAL * restrict vals)
      const;
    
  private:
    template <int O>
    void
    interpolate (ivect const & lsh,
                 vector<CCTK_REAL const *> const & varptrs,
                 CCTK_REAL * restrict vals)
      const;
    template <int O0, int O1, int O2>
    void
    interpolate (ivect const & lsh,
                 vector<CCTK_REAL const *> const & varptrs,
                 CCTK_REAL * restrict vals)
      const;
  };
  
  
  
  // A receive descriptor, describing what is received from other
  // processors
  struct recv_proc_t {
    int p;                      // sending processor
    int offset;
    int npoints;                // total number of received points
  };

  struct recv_descr_t {
    vector<recv_proc_t> procs;
    vector<int> procinds;
    int npoints;                // total number of received points
    
    vector<int> index;          // gather index list
  };
  
  // A send descriptor; describing what to send to other processors
  struct send_comp_t {
    // This structure does not exist for all components -- components
    // which are not accessed are not described, making this a sparse
    // data structure.  The field c contains the component number.
    vector<fasterp_src_loc_t> locs;
    int c;                      // source component
  };
  
  struct send_rl_t {
    vector<send_comp_t> comps;
    vector<int> compinds;
  };
  
  struct send_map_t {
    vector<send_rl_t> rls;
  };

  struct send_proc_t {
    // This structure does not exist for all processors -- processors
    // with which there is no communication are not described, making
    // this a sparse data structure.  The field p contains the
    // processor number.
    vector<send_map_t> maps;
    int p;                      // receiving processor
    int offset;
    int npoints;                // total number of sent points
  };
  
  struct send_descr_t {
    vector<send_proc_t> procs;
    // vector<int> procinds;
    int npoints;                // total number of sent points
  };
  
  
  
  class fasterp_setup_t {
    recv_descr_t recv_descr;
    send_descr_t send_descr;
    int order;
    
  public:
    fasterp_setup_t (cGH const * restrict cctkGH,
                     fasterp_glocs_t const & locations,
                     int order);
    
    ~ fasterp_setup_t ();
    
    void 
    interpolate (cGH const * restrict cctkGH,
                 vector<int> const & varinds,
                 vector<CCTK_REAL *> & values)
      const;
  };
  
  
  
} // namespace CarpetInterp2

#endif  // #define FASTERP_HH
