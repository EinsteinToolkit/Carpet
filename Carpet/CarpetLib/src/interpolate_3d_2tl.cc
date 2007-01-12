#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "operator_prototypes.hh"
#include "typeprops.hh"

using namespace std;



namespace CarpetLib {


  
#define SRCIND3(i,j,k)                                  \
  index3 (srcioff + (i), srcjoff + (j), srckoff + (k),  \
          srciext, srcjext, srckext)
#define DSTIND3(i,j,k)                                  \
  index3 (dstioff + (i), dstjoff + (j), dstkoff + (k),  \
          dstiext, dstjext, dstkext)
  
  
  
  template <typename T>
  void
  interpolate_3d_2tl (T const * restrict const src1,
                      CCTK_REAL const t1,
                      T const * restrict const src2,
                      CCTK_REAL const t2,
                      ivect3 const & srcext,
                      T * restrict const dst,
                      CCTK_REAL const t,
                      ivect3 const & dstext,
                      ibbox3 const & srcbbox,
                      ibbox3 const & dstbbox,
                      ibbox3 const & regbbox)
  {
    typedef typename typeprops<T>::real RT;
    
    
    
#if 0
    // This is already guaranteed by bbox
    if (any (srcbbox.stride() == 0 or
             dstbbox.stride() == 0 or
             regbbox.stride() == 0))
    {
      CCTK_WARN (0, "Internal error: stride is zero");
    }
#endif
    
    if (any (srcbbox.stride() != regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    if (any (srcbbox.stride() != dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
#if 0
    // This needs to be allowed for cell centring
    if (any (srcbbox.lower() % srcbbox.stride() != 0 or
             dstbbox.lower() % dstbbox.stride() != 0 or
             regbbox.lower() % regbbox.stride() != 0))
    {
      CCTK_WARN (0, "Internal error: array origins are not integer multiples of the strides");
    }
#endif
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
#if 0
    // This is already guaranteed by bbox
    if (any ((srcbbox.upper() - srcbbox.lower()) % srcbbox.stride() != 0 or
             (dstbbox.upper() - dstbbox.lower()) % dstbbox.stride() != 0 or
             (regbbox.upper() - regbbox.lower()) % regbbox.stride() != 0))
    {
      CCTK_WARN (0, "Internal error: array extents are not integer multiples of the strides");
    }
#endif
    
    if (not regbbox.is_contained_in(srcbbox) or
        not regbbox.is_contained_in(dstbbox))
    {
      CCTK_WARN (0, "Internal error: region extent is not contained in array extent");
    }
    
    if (any (srcext != srcbbox.shape() / srcbbox.stride() or
             dstext != dstbbox.shape() / dstbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all ((regbbox.lower() - srcbbox.lower()) % srcbbox.stride() == 0));
    ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / srcbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
    ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();
    
    
    
    size_t const srciext = srcext[0];
    size_t const srcjext = srcext[1];
    size_t const srckext = srcext[2];
    
    size_t const dstiext = dstext[0];
    size_t const dstjext = dstext[1];
    size_t const dstkext = dstext[2];
    
    size_t const regiext = regext[0];
    size_t const regjext = regext[1];
    size_t const regkext = regext[2];
    
    size_t const srcioff = srcoff[0];
    size_t const srcjoff = srcoff[1];
    size_t const srckoff = srcoff[2];
    
    size_t const dstioff = dstoff[0];
    size_t const dstjoff = dstoff[1];
    size_t const dstkoff = dstoff[2];
    
    
    
    // Linear (first order) interpolation
    
    RT const eps = 1.0e-10;
    if (abs (t1 - t2) < eps) {
      CCTK_WARN (0, "Internal error: arrays have same time");
    }
    if (t < min (t1, t2) - eps or t > max (t1, t2) + eps) {
      CCTK_WARN (0, "Internal error: extrapolation in time");
    }
    
    RT const s1fac = (t - t2) / (t1 - t2);
    RT const s2fac = (t - t1) / (t2 - t1);
    
    
    
    // Loop over region
    for (size_t k=0; k<regkext; ++k) {
      for (size_t j=0; j<regjext; ++j) {
        for (size_t i=0; i<regiext; ++i) {
          
          dst [DSTIND3(i, j, k)] =
            + s1fac * src1 [SRCIND3(i, j, k)]
            + s2fac * src2 [SRCIND3(i, j, k)];
          
        }
      }
    }
    
  }
  
  
  
  template
  void
  interpolate_3d_2tl (CCTK_REAL const * restrict const src1,
                      CCTK_REAL const t1,
                      CCTK_REAL const * restrict const src2,
                      CCTK_REAL const t2,
                      ivect3 const & srcext,
                      CCTK_REAL * restrict const dst,
                      CCTK_REAL const t,
                      ivect3 const & dstext,
                      ibbox3 const & srcbbox,
                      ibbox3 const & dstbbox,
                      ibbox3 const & regbbox);
  
  template
  void
  interpolate_3d_2tl (CCTK_COMPLEX const * restrict const src1,
                      CCTK_REAL const t1,
                      CCTK_COMPLEX const * restrict const src2,
                      CCTK_REAL const t2,
                      ivect3 const & srcext,
                      CCTK_COMPLEX * restrict const dst,
                      CCTK_REAL const t,
                      ivect3 const & dstext,
                      ibbox3 const & srcbbox,
                      ibbox3 const & dstbbox,
                      ibbox3 const & regbbox);
  
  
  
} // namespace CarpetLib
