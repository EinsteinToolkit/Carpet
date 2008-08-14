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
  
  
  
#define SRCIND3(i,j,k)                          \
  index3 (i, j, k,                              \
          srciext, srcjext, srckext)
#define DSTIND3(i,j,k)                          \
  index3 (i, j, k,                              \
          dstiext, dstjext, dstkext)
  
  
  
  // 1D interpolation coefficients
  static
  int const ncoeffs = 10;
  
  template <typename RT>
  static inline
  RT
  coeff (int const i)
  {
    static const RT coeffs[] = {
      -   35/RT(65536.0),
         385/RT(65536.0),
      -  495/RT(16384.0),
        1617/RT(16384.0),
      - 8085/RT(32768.0),
       24255/RT(32768.0),
        8085/RT(16384.0),
      - 1155/RT(16384.0),
         693/RT(65536.0),
      -   55/RT(65536.0)
    };
#ifdef CARPET_DEBUG
    assert (sizeof coeffs / sizeof *coeffs == ncoeffs);
    assert (i>=0 and i<ncoeffs);
#endif
    return coeffs[i];
  }
  
  
  
  // 0D "interpolation"
  template <typename T>
  static inline
  T
  interp0 (T const * restrict const p)
  {
    return * p;
  }
  
  // 1D interpolation
  template <typename T>
  static inline
  T
  interp1 (T const * restrict const p,
           size_t const d1)
  {
    typedef typename typeprops<T>::real RT;
    T res = typeprops<T>::fromreal (0);
    for (int i=0; i<ncoeffs; ++i) {
      res += coeff<RT>(i) * interp0<T> (p + i*d1);
    }
    return res;
  }
  
  // 2D interpolation
  template <typename T>
  static inline
  T
  interp2 (T const * restrict const p,
           size_t const d1,
           size_t const d2)
  {
    typedef typename typeprops<T>::real RT;
    T res = typeprops<T>::fromreal (0);
    for (int i=0; i<ncoeffs; ++i) {
      res += coeff<RT>(i) * interp1<T> (p + i*d2, d1);
    }
    return res;
  }
  
  // 3D interpolation
  template <typename T>
  static inline
  T
  interp3 (T const * restrict const p,
           size_t const d1,
           size_t const d2,
           size_t const d3)
  {
    typedef typename typeprops<T>::real RT;
    T res = typeprops<T>::fromreal (0);
    for (int i=0; i<ncoeffs; ++i) {
      res += coeff<RT>(i) * interp2<T> (p + i*d3, d1, d2);
    }
    return res;
  }
  
  
  
  template <typename T>
  void
  prolongate_3d_o9_rf2 (T const * restrict const src,
                        ivect3 const & restrict srcext,
                        T * restrict const dst,
                        ivect3 const & restrict dstext,
                        ibbox3 const & restrict srcbbox,
                        ibbox3 const & restrict dstbbox,
                        ibbox3 const & restrict regbbox)
  {
    typedef typename typeprops<T>::real RT;
    
    
    
    if (any (srcbbox.stride() <= regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    if (any (srcbbox.stride() != reffact2 * dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: source strides are not twice the destination strides");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all ((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
    ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / regbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
    ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();
    
    
    
    bvect3 const needoffsetlo = srcoff % reffact2 != 0 or regext > 1;
    bvect3 const needoffsethi = (srcoff + regext - 1) % reffact2 != 0 or regext > 1;
    ivect3 const offsetlo = either (needoffsetlo, 5, 0);
    ivect3 const offsethi = either (needoffsethi, 5, 0);
    
    
    
    if (not regbbox.expand(offsetlo, offsethi).is_contained_in(srcbbox) or
        not regbbox                           .is_contained_in(dstbbox))
    {
      CCTK_WARN (0, "Internal error: region extent is not contained in array extent");
    }
    
    if (any (srcext != srcbbox.shape() / srcbbox.stride() or
             dstext != dstbbox.shape() / dstbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes");
    }
    
    
    
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
    
    
    
    size_t const fi = srcioff % 2;
    size_t const fj = srcjoff % 2;
    size_t const fk = srckoff % 2;
    
    size_t const i0 = srcioff / 2;
    size_t const j0 = srcjoff / 2;
    size_t const k0 = srckoff / 2;
    
    
    
    size_t const srcdi = 1;
    size_t const srcdj = srcdi * srciext;
    size_t const srcdk = srcdj * srcjext;
    
    
    
    // Loop over fine region
    // Label scheme: l 8 fk fj fi
    
    size_t is, js, ks;
    size_t id, jd, kd;
    size_t i, j, k;
    
    // begin k loop
    k = 0;
    ks = k0;
    kd = dstkoff;
    if (fk == 0) goto l80;
    goto l81;
    
    // begin j loop
   l80:
    j = 0;
    js = j0;
    jd = dstjoff;
    if (fj == 0) goto l800;
    goto l801;
    
    // begin i loop
   l800:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8000;
    goto l8001;
    
    // kernel
   l8000:
    dst[DSTIND3(id,jd,kd)] = interp0<T> (& src[SRCIND3(is,js,ks)]);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8001;
    goto l900;
    
    // kernel
   l8001:
    dst[DSTIND3(id,jd,kd)] = interp1<T> (& src[SRCIND3(is-4,js,ks)], srcdi);
    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8000;
    goto l900;
    
    // end i loop
   l900:
    j = j+1;
    jd = jd+1;
    if (j < regjext) goto l801;
    goto l90;
    
    // begin i loop
   l801:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8010;
    goto l8011;
    
    // kernel
   l8010:
    dst[DSTIND3(id,jd,kd)] = interp1<T> (& src[SRCIND3(is,js-4,ks)], srcdj);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8011;
    goto l901;
    
    // kernel
   l8011:
    dst[DSTIND3(id,jd,kd)] =
      interp2<T> (& src[SRCIND3(is-4,js-4,ks)], srcdi, srcdj);
    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8010;
    goto l901;
    
    // end i loop
   l901:
    j = j+1;
    jd = jd+1;
    js = js+1;
    if (j < regjext) goto l800;
    goto l90;
    
    // end j loop
   l90:
    k = k+1;
    kd = kd+1;
    if (k < regkext) goto l81;
    goto l9;
    
    // begin j loop
   l81:
    j = 0;
    js = j0;
    jd = dstjoff;
    if (fj == 0) goto l810;
    goto l811;
    
    // begin i loop
   l810:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8100;
    goto l8101;
    
    // kernel
   l8100:
    dst[DSTIND3(id,jd,kd)] = interp1<T> (& src[SRCIND3(is,js,ks-4)], srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8101;
    goto l910;
    
    // kernel
   l8101:
    dst[DSTIND3(id,jd,kd)] =
      interp2<T> (& src[SRCIND3(is-4,js,ks-4)], srcdi, srcdj);
    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8100;
    goto l910;
    
    // end i loop
   l910:
    j = j+1;
    jd = jd+1;
    if (j < regjext) goto l811;
    goto l91;
    
    // begin i loop
   l811:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8110;
    goto l8111;
    
    // kernel
   l8110:
    dst[DSTIND3(id,jd,kd)] =
      interp2<T> (& src[SRCIND3(is,js-4,ks-4)], srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8111;
    goto l911;
    
    // kernel
   l8111:
    {
      dst[DSTIND3(id,jd,kd)] =
        interp3<T> (& src[SRCIND3(is-4,js-4,ks-4)], srcdi, srcdj, srcdk);
    }
    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8110;
    goto l911;
    
    // end i loop
   l911:
    j = j+1;
    jd = jd+1;
    js = js+1;
    if (j < regjext) goto l810;
    goto l91;
    
    // end j loop
   l91:
    k = k+1;
    kd = kd+1;
    ks = ks+1;
    if (k < regkext) goto l80;
    goto l9;
    
    // end k loop
   l9:;
    
  }
  
  
  
#define INSTANTIATE(T)                                          \
  template                                                      \
  void                                                          \
  prolongate_3d_o9_rf2 (T const * restrict const src,           \
                        ivect3 const & restrict srcext,         \
                        T * restrict const dst,                 \
                        ivect3 const & restrict dstext,         \
                        ibbox3 const & restrict srcbbox,        \
                        ibbox3 const & restrict dstbbox,        \
                        ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib
