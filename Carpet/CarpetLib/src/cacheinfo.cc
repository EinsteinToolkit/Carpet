#include "cacheinfo.hh"

#include <cctk_Parameters.h>

#include <vectors.h>



template<int D>
vect<int,D>
pad_shape (bbox<int,D> const& extent)
{
  assert (all (extent.shape() >= 0));
  return pad_shape(extent.shape() / extent.stride());
}



template<int D>
vect<int,D>
pad_shape (vect<int,D> const& shape)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (all(shape>=0));
  
  static bool have_cacheinfo = false;
  static vector<cacheinfo_t> cacheinfo;
  if (not have_cacheinfo) {
    // Ignore L1 caches that are probably too small to be useful (e.g.
    // on Intel or AMD processors)
    // TODO: make this a parameter
    if (D1size >= 128*1024) {
      cacheinfo.push_back(cacheinfo_t(D1size, D1linesize, D1assoc));
    }
#if 0
    // TODO: this is too simplistic:
    // Add page size as a cache
    if (pagesize>0) {
      cacheinfo.push_back(cacheinfo_t(pagesize));
    }
#endif
    if (L2size>0) {
      cacheinfo.push_back(cacheinfo_t(L2size, L2linesize, L2assoc));
    }
    if (L3size>0) {
      cacheinfo.push_back(cacheinfo_t(L3size, L3linesize, L3assoc));
    }
    if (TLB_D1entries>0) {
      ptrdiff_t const TLB_D1size = TLB_D1entries * TLB_D1pagesize * TLB_D1assoc;
      cacheinfo.push_back(cacheinfo_t(TLB_D1size, TLB_D1pagesize, TLB_D1assoc));
    }
    if (TLB_L2entries>0) {
      ptrdiff_t const TLB_L2size = TLB_L2entries * TLB_L2pagesize * TLB_L2assoc;
      cacheinfo.push_back(cacheinfo_t(TLB_L2size, TLB_L2pagesize, TLB_L2assoc));
    }
    
    // TODO: sort caches by their sizes
    for (size_t n=0; n<cacheinfo.size(); ++n) {
      cacheinfo_t const& ci = cacheinfo.at(n);
      if (n>0) {
        // Ensure that the cache size is larger than the next lower
        // cache size
        assert (ci.size() > cacheinfo.at(n-1).size());
        // Ensure that the cache line size is evenly divided by the
        // next lower cache line size
        assert (ci.linesize() % cacheinfo.at(n-1).linesize() == 0);
        assert (ci.stride() > cacheinfo.at(n-1).stride());
      }
    } // for cacheinfo
    
    have_cacheinfo = true;
  } // if not have_cacheinfo
  
  vect<int,D> padded_shape;
  int accumulated_npoints = 1;
  for (int d=0; d<D; ++d) {
    int npoints = shape[d];
    
    if (d == 0) {
#if VECTORISE && VECTORISE_ALIGNED_ARRAYS
      // Pad array to a multiple of the vector size. Note that this is
      // a hard requirement, so that we can emit aligned load/store
      // operations.
      npoints = align_up (npoints, CCTK_REAL_VEC_SIZE);
#endif
      if (vector_size > 0) {
        npoints = align_up (npoints, vector_size);
      }
    }
    
    for (size_t n=0; n<cacheinfo.size(); ++n) {
      cacheinfo_t const& ci = cacheinfo.at(n);
      
      // Pad array in this direction to a multiple of this cache line
      // size
      assert (ci.linesize() % sizeof(CCTK_REAL) == 0);
      int const linesize = ci.linesize() / sizeof(CCTK_REAL);
      assert (is_power_of_2(linesize));
      if (npoints * accumulated_npoints >= linesize) {
        // The extent is at least one cache line long: round up to the
        // next full cache line
        npoints = align_up (npoints, linesize);
      } else {
#if 0
        // The extent is less than one cache line long: Ensure that
        // the array size divides the cache line size evenly by
        // rounding to the next power of 2
        // NOTE: This is disabled, since this would align everything
        // to powers of 2.
        npoints = next_power_of_2(npoints);
#endif
      }
      
      // Avoid multiples of the cache stride
      assert (ci.stride() % sizeof(CCTK_REAL) == 0);
      int const stride = ci.stride() / sizeof(CCTK_REAL);
      if (npoints * accumulated_npoints % stride == 0) {
        assert (linesize < stride);
        npoints += linesize;
      }
      
    } // for cacheinfo
    
    padded_shape[d] = npoints;
    accumulated_npoints *= npoints;
  }
  assert (prod (padded_shape) == accumulated_npoints);
  
  // self-check
  for (int d=0; d<D; ++d) {
    assert (padded_shape[d] >= shape[d]);
#if VECTORISE && VECTORISE_ALIGNED_ARRAYS
    if (d == 0) {
      assert (padded_shape[d] % CCTK_REAL_VEC_SIZE == 0);
    }
#endif
    if (vector_size > 0) {
      if (d == 0) {
        assert (padded_shape[d] % vector_size == 0);
      }
    }
    
    // TODO: add self-checks for the other requirements as well
  }
  
  if (verbose) {
    ostringstream buf;
    buf << "padding " << shape << " to " << padded_shape;
    CCTK_INFO (buf.str().c_str());
  }
  
  return padded_shape;
}



template vect<int,3> pad_shape (bbox<int,3> const& extent);
template vect<int,3> pad_shape (vect<int,3> const& shape);

template vect<int,4> pad_shape (bbox<int,4> const& extent);
template vect<int,4> pad_shape (vect<int,4> const& shape);
