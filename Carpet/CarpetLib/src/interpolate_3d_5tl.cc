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
  interpolate_3d_5tl (T const * restrict const src1,
                      CCTK_REAL const t1,
                      T const * restrict const src2,
                      CCTK_REAL const t2,
                      T const * restrict const src3,
                      CCTK_REAL const t3,
                      T const * restrict const src4,
                      CCTK_REAL const t4,
                      T const * restrict const src5,
                      CCTK_REAL const t5,
                      ivect3 const & restrict srcext,
                      T * restrict const dst,
                      CCTK_REAL const t,
                      ivect3 const & restrict dstext,
                      ibbox3 const & restrict srcbbox,
                      ibbox3 const & restrict dstbbox,
                      ibbox3 const & restrict regbbox)
  {
    typedef typename typeprops<T>::real RT;
    
    
    
    if (any (srcbbox.stride() != regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    if (any (srcbbox.stride() != dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
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
    
    
    
    // Quadratic (second order) interpolation
    
    RT const eps = 1.0e-10;
    
    if (abs (t1 - t2) < eps or abs (t1 - t3) < eps or abs (t1 - t4) < eps or
        abs (t1 - t5) < eps or abs (t2 - t3) < eps or abs (t2 - t4) < eps or
        abs (t2 - t5) < eps or abs (t3 - t4) < eps or abs (t3 - t5) < eps or
        abs (t4 - t5) < eps)
    {
      CCTK_WARN (0, "Internal error: arrays have same time");
    }
    if (t < min (min (min (min (t1, t2), t3), t4), t5) - eps or
        t > max (max (max (max (t1, t2), t3), t4), t5) + eps)
    {
      CCTK_WARN (0, "Internal error: extrapolation in time");
    }
    
    RT const s1fac = (t - t2) * (t - t3) * (t - t4) * (t - t5) / ((t1 - t2) * (t1 - t3) * (t1 - t4) * (t1 - t5));
    RT const s2fac = (t - t1) * (t - t3) * (t - t4) * (t - t5) / ((t2 - t1) * (t2 - t3) * (t2 - t4) * (t2 - t5));
    RT const s3fac = (t - t1) * (t - t2) * (t - t4) * (t - t5) / ((t3 - t1) * (t3 - t2) * (t3 - t4) * (t3 - t5));
    RT const s4fac = (t - t1) * (t - t2) * (t - t3) * (t - t5) / ((t4 - t1) * (t4 - t2) * (t4 - t3) * (t4 - t5));
    RT const s5fac = (t - t1) * (t - t2) * (t - t3) * (t - t4) / ((t5 - t1) * (t5 - t2) * (t5 - t3) * (t5 - t4));
    
    
    
    // Loop over region
    for (size_t k=0; k<regkext; ++k) {
      for (size_t j=0; j<regjext; ++j) {
        for (size_t i=0; i<regiext; ++i) {
          
          dst [DSTIND3(i, j, k)] =
            + s1fac * src1 [SRCIND3(i, j, k)]
            + s2fac * src2 [SRCIND3(i, j, k)]
            + s3fac * src3 [SRCIND3(i, j, k)]
            + s4fac * src4 [SRCIND3(i, j, k)]
            + s5fac * src5 [SRCIND3(i, j, k)];
          
        }
      }
    }
    
  }
  
  
  
#define INSTANTIATE(T)                                  \
  template                                              \
  void                                                  \
  interpolate_3d_5tl (T const * restrict const src1,    \
                      CCTK_REAL const t1,               \
                      T const * restrict const src2,    \
                      CCTK_REAL const t2,               \
                      T const * restrict const src3,    \
                      CCTK_REAL const t3,               \
                      T const * restrict const src4,    \
                      CCTK_REAL const t4,               \
                      T const * restrict const src5,    \
                      CCTK_REAL const t5,               \
                      ivect3 const & restrict srcext,   \
                      T * restrict const dst,           \
                      CCTK_REAL const t,                \
                      ivect3 const & restrict dstext,   \
                      ibbox3 const & restrict srcbbox,  \
                      ibbox3 const & restrict dstbbox,  \
                      ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib