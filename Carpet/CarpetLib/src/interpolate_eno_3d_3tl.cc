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
  static inline
  T
  min3 (T const & x, T const & y, T const & z)
  {
    return min (x, min (y, z));
  }
  
  template <typename T>
  static inline
  T
  max3 (T const & x, T const & y, T const & z)
  {
    return max (x, max (y, z));
  }
  
  
  
  template <typename T>
  void
  interpolate_eno_3d_3tl (T const * restrict const src1,
                          CCTK_REAL const t1,
                          T const * restrict const src2,
                          CCTK_REAL const t2,
                          T const * restrict const src3,
                          CCTK_REAL const t3,
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
    
    if (abs (t1 - t2) < eps or abs (t1 - t3) < eps or abs (t2 - t3) < eps) {
      CCTK_WARN (0, "Internal error: arrays have same time");
    }
    if (t < min3 (t1, t2, t3) - eps or t > max3 (t1, t2, t3) + eps) {
      CCTK_WARN (0, "Internal error: extrapolation in time");
    }
    
    RT const s1fac3 = (t - t2) * (t - t3) / ((t1 - t2) * (t1 - t3));
    RT const s2fac3 = (t - t1) * (t - t3) / ((t2 - t1) * (t2 - t3));
    RT const s3fac3 = (t - t1) * (t - t2) / ((t3 - t1) * (t3 - t2));
    
    RT const s1fac2_12 = (t - t2) / (t1 - t2);
    RT const s2fac2_12 = (t - t1) / (t2 - t1);
    
    RT const s2fac2_23 = (t - t3) / (t2 - t3);
    RT const s3fac2_23 = (t - t2) / (t3 - t2);
    
    bool const use_12 =
      t >= min (t1, t2) - eps and t <= max (t1, t2) + eps;
    bool const use_23 =
      t >= min (t2, t3) - eps and t <= max (t2, t3) + eps;
    assert (use_12 or use_23);
    
    
    
    // Loop over region
    for (size_t k=0; k<regkext; ++k) {
      for (size_t j=0; j<regjext; ++j) {
        for (size_t i=0; i<regiext; ++i) {
          
          T const s1 = src1 [SRCIND3(i, j, k)];
          T const s2 = src2 [SRCIND3(i, j, k)];
          T const s3 = src3 [SRCIND3(i, j, k)];
          
          T d = s1fac3 * s1 + s2fac3 * s2 + s3fac3 * s3;
          
          // if ((d - max3 (s1, s2, s3)) * (d - min3 (s1, s2, s3)) < 0)
          if (d > max3 (s1, s2, s3) or d < min3 (s1, s2, s3)) {
            if (use_12) {
              d = s1fac2_12 * s1 + s2fac2_12 * s2;
            } else {
              d = s2fac2_23 * s2 + s3fac2_23 * s3;
            }
          }
          
          dst [DSTIND3(i, j, k)] = d;
          
        }
      }
    }
    
  }
  
  
  
#ifdef HAVE_CCTK_COMPLEX8
  template <>
  void
  interpolate_eno_3d_3tl (CCTK_COMPLEX8 const * restrict const src1,
                          CCTK_REAL const t1,
                          CCTK_COMPLEX8 const * restrict const src2,
                          CCTK_REAL const t2,
                          CCTK_COMPLEX8 const * restrict const src3,
                          CCTK_REAL const t3,
                          ivect3 const & restrict srcext,
                          CCTK_COMPLEX8 * restrict const dst,
                          CCTK_REAL const t,
                          ivect3 const & restrict dstext,
                          ibbox3 const & restrict srcbbox,
                          ibbox3 const & restrict dstbbox,
                          ibbox3 const & restrict regbbox)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "ENO for complex numbers is not supported");
  }
#endif
  
#ifdef HAVE_CCTK_COMPLEX16
  template <>
  void
  interpolate_eno_3d_3tl (CCTK_COMPLEX16 const * restrict const src1,
                          CCTK_REAL const t1,
                          CCTK_COMPLEX16 const * restrict const src2,
                          CCTK_REAL const t2,
                          CCTK_COMPLEX16 const * restrict const src3,
                          CCTK_REAL const t3,
                          ivect3 const & restrict srcext,
                          CCTK_COMPLEX16 * restrict const dst,
                          CCTK_REAL const t,
                          ivect3 const & restrict dstext,
                          ibbox3 const & restrict srcbbox,
                          ibbox3 const & restrict dstbbox,
                          ibbox3 const & restrict regbbox)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "ENO for complex numbers is not supported");
  }
#endif
  
#ifdef HAVE_CCTK_COMPLEX32
  template <>
  void
  interpolate_eno_3d_3tl (CCTK_COMPLEX32 const * restrict const src1,
                          CCTK_REAL const t1,
                          CCTK_COMPLEX32 const * restrict const src2,
                          CCTK_REAL const t2,
                          CCTK_COMPLEX32 const * restrict const src3,
                          CCTK_REAL const t3,
                          ivect3 const & restrict srcext,
                          CCTK_COMPLEX32 * restrict const dst,
                          CCTK_REAL const t,
                          ivect3 const & restrict dstext,
                          ibbox3 const & restrict srcbbox,
                          ibbox3 const & restrict dstbbox,
                          ibbox3 const & restrict regbbox)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "ENO for complex numbers is not supported");
  }
#endif
  
  
  
#define INSTANTIATE(T)                                          \
  template                                                      \
  void                                                          \
  interpolate_eno_3d_3tl (T const * restrict const src1,        \
                          CCTK_REAL const t1,                   \
                          T const * restrict const src2,        \
                          CCTK_REAL const t2,                   \
                          T const * restrict const src3,        \
                          CCTK_REAL const t3,                   \
                          ivect3 const & restrict srcext,       \
                          T * restrict const dst,               \
                          CCTK_REAL const t,                    \
                          ivect3 const & restrict dstext,       \
                          ibbox3 const & restrict srcbbox,      \
                          ibbox3 const & restrict dstbbox,      \
                          ibbox3 const & restrict regbbox);
#define CARPET_NO_COMPLEX
#include "instantiate"
#undef CARPET_NO_COMPLEX
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib
