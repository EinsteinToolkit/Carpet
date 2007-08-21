// See also Hern, "Numerical Relativity and Inhomogeneous
// Cosmologies", gr-qc/0004036, section 3.2, pp. 29 ff.; especially
// the last equation on page 37.



#include <algorithm>
#include <cassert>
#include <cmath>

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
  
  
  
  // Convert from the "standard" form of the grid function to the
  // "primitive" version, i.e., the antiderivative
  
  template <typename T>
  void
  prolongate_3d_cc_rf2_std2prim (T const * restrict const src,
                                 ivect3 const & restrict srcext,
                                 T * restrict const dst,
                                 ivect3 const & restrict dstext,
                                 ibbox3 const & restrict srcbbox,
                                 ibbox3 const & restrict dstbbox,
                                 ibbox3 const & restrict regbbox)
  {
    DECLARE_CCTK_PARAMETERS;
    
    typedef typename typeprops<T>::real RT;
    T (* const fromreal) (RT) = typeprops<T>::fromreal;
    
    
    
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
    assert (all (regbbox.stride() % 2 == 0));
    assert (all ((regbbox.lower() - srcbbox.lower() + regbbox.stride() / 2) %
                 regbbox.stride() == 0));
    ivect3 const srcoff =
      (regbbox.lower() - srcbbox.lower() + regbbox.stride() / 2) /
      regbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
    ivect3 const dstoff =
      (regbbox.lower() - dstbbox.lower()) / regbbox.stride();
    
    
    
    int const srciext = srcext[0];
    int const srcjext = srcext[1];
    int const srckext = srcext[2];
    
    int const dstiext = dstext[0];
    int const dstjext = dstext[1];
    int const dstkext = dstext[2];
    
    int const regiext = regext[0];
    int const regjext = regext[1];
    int const regkext = regext[2];
    
    int const srcioff = srcoff[0];
    int const srcjoff = srcoff[1];
    int const srckoff = srcoff[2];
    
    int const dstioff = dstoff[0];
    int const dstjoff = dstoff[1];
    int const dstkoff = dstoff[2];
    
    
    
    T const zero = fromreal (0);
    
    
    
#pragma omp parallel for
    for (int k=0; k<regkext; ++k) {
      for (int j=0; j<regjext; ++j) {
        for (int i=0; i<regiext; ++i) {
          if (i==0 or j==0 or k==0) {
            dst [DSTIND3(i, j, k)] = zero;
          } else {
            // // 1D
            // dst [DSTIND1(i)] =
            //   + dst [DSTIND1(i-1)]
            //   + src [SRCIND1(i-1)];
            // // 2D
            // dst [DSTIND2(i, j, k)] =
            //   + dst [DSTIND2(i-1, j)]
            //   + dst [DSTIND2(i, j-1)]
            //   - dst [DSTIND2(i-1, j-1)]
            //   + src [SRCIND2(i-1, j-1)];
            // 3D
            dst [DSTIND3(i, j, k)] =
              + dst [DSTIND3(i-1, j, k)]
              + dst [DSTIND3(i, j-1, k)]
              + dst [DSTIND3(i, j, k-1)]
              - dst [DSTIND3(i, j-1, k-1)]
              - dst [DSTIND3(i-1, j, k-1)]
              - dst [DSTIND3(i-1, j-1, k)]
              + dst [DSTIND3(i-1, j-1, k-1)]
              + src [SRCIND3(i-1, j-1, k-1)];
          }
        }
      }
    }
    
  }
  
  
  
#define INSTANTIATE(T)                                                  \
  template                                                              \
  void                                                                  \
  prolongate_3d_cc_rf2_std2prim (T const * restrict const src,          \
                                 ivect3 const & restrict srcext,        \
                                 T * restrict const dst,                \
                                 ivect3 const & restrict dstext,        \
                                 ibbox3 const & restrict srcbbox,       \
                                 ibbox3 const & restrict dstbbox,       \
                                 ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
  // Convert from the "primitive" form of the grid function to the
  // "standard" version
  
  template <typename T>
  void
  prolongate_3d_cc_rf2_prim2std (T const * restrict const src,
                                 ivect const & restrict srcext,
                                 T * restrict const dst,
                                 ivect const & restrict dstext,
                                 ibbox3 const & restrict srcbbox,
                                 ibbox3 const & restrict dstbbox,
                                 ibbox3 const & restrict regbbox)
  {
    DECLARE_CCTK_PARAMETERS;
    
    
    
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
    assert (all (regbbox.stride() % 2 == 0));
    assert (all ((regbbox.lower() - srcbbox.lower() - regbbox.stride() / 2) %
                 regbbox.stride() == 0));
    ivect3 const srcoff =
      (regbbox.lower() - srcbbox.lower() - regbbox.stride() / 2) /
      regbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
    ivect3 const dstoff =
      (regbbox.lower() - dstbbox.lower()) / regbbox.stride();
    
    
    
    int const srciext = srcext[0];
    int const srcjext = srcext[1];
    int const srckext = srcext[2];
    
    int const dstiext = dstext[0];
    int const dstjext = dstext[1];
    int const dstkext = dstext[2];
    
    int const regiext = regext[0];
    int const regjext = regext[1];
    int const regkext = regext[2];
    
    int const srcioff = srcoff[0];
    int const srcjoff = srcoff[1];
    int const srckoff = srcoff[2];
    
    int const dstioff = dstoff[0];
    int const dstjoff = dstoff[1];
    int const dstkoff = dstoff[2];
    
    
    
#pragma omp parallel for
    for (int k=0; k<regkext; ++k) {
      for (int j=0; j<regjext; ++j) {
        for (int i=0; i<regiext; ++i) {
          dst [DSTIND3(i, j, k)] = reffact2 *
            (- src [SRCIND3(i, j+1, k+1)]
             - src [SRCIND3(i+1, j, k+1)]
             - src [SRCIND3(i+1, j+1, k)]
             + src [SRCIND3(i+1, j, k)]
             + src [SRCIND3(i, j+1, k)]
             + src [SRCIND3(i, j, k+1)]
             - src [SRCIND3(i, j, k)]
             + src [SRCIND3(i+1, j+1, k+1)]);
        }
      }
    }
    
  }
  
  
  
#define INSTANTIATE(T)                                                  \
  template                                                              \
  void                                                                  \
  prolongate_3d_cc_rf2_prim2std (T const * restrict const src,          \
                                 ivect const & restrict srcext,         \
                                 T * restrict const dst,                \
                                 ivect const & restrict dstext,         \
                                 ibbox3 const & restrict srcbbox,       \
                                 ibbox3 const & restrict dstbbox,       \
                                 ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib
