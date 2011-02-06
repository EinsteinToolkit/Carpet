#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "gdata.hh"
#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;



namespace CarpetLib {
  

  
#define SRCIND3(i,j,k) index3 (i, j, k, srciext, srcjext, srckext)
#define DSTIND3(i,j,k) index3 (i, j, k, dstiext, dstjext, dstkext)
  
  
  
    
  namespace coeffs_3d_cc_rf2 {

  // 1D interpolation coefficients
  
  template<typename RT, int ORDER>
  struct coeffs1d {
    static const RT coeffs[];
    
    static ptrdiff_t const ncoeffs = ORDER+1;
    static ptrdiff_t const imin    = - ncoeffs/2 + 1;
    static ptrdiff_t const imax    = imin + ncoeffs;
    
    static inline
    RT
    get (int const di, int const i)
    {
      ptrdiff_t const j = di == 0 ? i - imin : imax-1 - (i-imin);
#ifdef CARPET_DEBUG
      assert (di == 0 or di == 1);
      assert (ncoeffs == sizeof coeffs / sizeof *coeffs);
      assert (j>=0 and j<ncoeffs);
#endif
      return coeffs[j];
    }
    
    // Test coefficients
    static
    void
    test ()
    {
      static bool tested = false;
      if (tested) return;
      tested = true;
      
      assert (ncoeffs == sizeof coeffs / sizeof *coeffs);
      
      // Test all orders and offsets
      bool error = false;
      for (int n=0; n<=ORDER; ++n) {
        for (int di=0; di<2; ++di) {
          RT res = RT(0.0);
          for (int i=imin; i<imax; ++i) {
            CCTK_REAL const dx = ORDER % 2;
            RT const x = RT(CCTK_REAL(i) - (di==0 ? 0.75 : 1.25 - dx));
            RT const y = ipow (x, n);
            res += get(di,i) * y;
          }
          RT const x0 = RT(0.0);
          RT const y0 = ipow (x0, n);
          if (not (good::abs (res - y0) < 1.0e-12)) {
            RT rt;
            ostringstream buf;
            buf << "Error in prolongate_3d_cc_rf2::coeffs_3d_cc_rf2\n"
                << "   RT=" << typestring(rt) << "\n"
                << "   ORDER=" << ORDER << "\n"
                << "   n=" << n << "\n"
                << "   di=" << di << "\n"
                << "   y0=" << y0 << ", res=" << res;
            CCTK_WARN (CCTK_WARN_ALERT, buf.str().c_str());
            error = true;
          }
        } // for di
      }   // for n
      if (error) {
        CCTK_WARN (CCTK_WARN_ABORT, "Aborting.");
      }
    }
  };
  
  
   
#define TYPECASE(N,RT)                          \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,0>::coeffs[] = {         \
    +1 / RT(1.0)                                \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,1>::coeffs[] = {         \
    +1 / RT(4.0),                               \
    +3 / RT(4.0)                                \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,2>::coeffs[] = {         \
    + 5 / RT(32.0),                             \
    +30 / RT(32.0),                             \
    - 3 / RT(32.0)                              \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,3>::coeffs[] = {         \
    -  5 / RT(128.0),                           \
    + 35 / RT(128.0),                           \
    +105 / RT(128.0),                           \
    -  7 / RT(128.0)                            \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,4>::coeffs[] = {         \
    -  45 / RT(2048.0),                         \
    + 420 / RT(2048.0),                         \
    +1890 / RT(2048.0),                         \
    - 252 / RT(2048.0)                          \
    +  35 / RT(2048.0)                          \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,5>::coeffs[] = {         \
    +  63 / RT(8192.0),                         \
    - 495 / RT(8192.0),                         \
    +2310 / RT(8192.0),                         \
    +6930 / RT(8192.0),                         \
    - 639 / RT(8192.0),                         \
    +  77 / RT(8192.0)                          \
  };

#define CARPET_NO_COMPLEX
#include "typecase.hh"
#undef CARPET_NO_COMPLEX
#undef TYPECASE
  
  } // namespace coeffs_3d_cc_rf2
  
  using namespace coeffs_3d_cc_rf2;
  
  
  
  // 0D "interpolation"
  template <typename T, int ORDER>
  static inline
  T
  interp0 (T const * restrict const p)
  {
    return * p;
  }
  
  // 1D interpolation
  template <typename T, int ORDER, int di>
  static inline
  T
  interp1 (T const * restrict const p,
           size_t const d1)
  {
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    T res = typeprops<T>::fromreal (0);
    for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
      res += coeffs::get(di,i) * interp0<T,ORDER> (p + i*d1);
    }
    return res;
  }
  
  // 2D interpolation
  template <typename T, int ORDER, int di, int dj>
  static inline
  T
  interp2 (T const * restrict const p,
           size_t const d1,
           size_t const d2)
  {
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    T res = typeprops<T>::fromreal (0);
    for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
      res += coeffs::get(dj,i) * interp1<T,ORDER,di> (p + i*d2, d1);
    }
    return res;
  }
  
  // 3D interpolation
  template <typename T, int ORDER, int di, int dj, int dk>
  static inline
  T
  interp3 (T const * restrict const p,
           size_t const d1,
           size_t const d2,
           size_t const d3)
  {
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    T res = typeprops<T>::fromreal (0);
    for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
      res += coeffs::get(dk,i) * interp2<T,ORDER,di,dj> (p + i*d3, d1, d2);
    }
    return res;
  }
  
  
  
  template <typename T, int ORDER>
  void
  prolongate_3d_cc_rf2 (T const * restrict const src,
                        ivect3 const & restrict srcext,
                        T * restrict const dst,
                        ivect3 const & restrict dstext,
                        ibbox3 const & restrict srcbbox,
                        ibbox3 const & restrict dstbbox,
                        ibbox3 const & restrict regbbox)
  {
    static_assert (ORDER>=0, "ORDER must be non-negative");
    
    typedef typename typeprops<T>::real RT;
    coeffs1d<RT,ORDER>::test ();
    
    
    
    if (any (srcbbox.stride() <= regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    if (any (srcbbox.stride() != reffact2 * dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: source strides are not twice the destination strides");
    }
    
    if (any (dstbbox.stride() % 2 != 0)) {
      CCTK_WARN (0, "Internal error: destination strides are not even");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all ((regbbox.lower() - srcbbox.lower() + regbbox.stride() / 2) % regbbox.stride() == 0));
    ivect3 const srcoff = (regbbox.lower() - srcbbox.lower() + regbbox.stride() / 2) / regbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
    ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();
    
    
    
    // The name "needoffset" does not make sense for cell centring
    bvect3 const needoffsetlo = srcoff % reffact2 != 0;
    bvect3 const needoffsethi = (srcoff + regext - 1) % reffact2 != 0;
    // This is probably wrong for odd orders
    ivect3 const offsetlo =
      ORDER%2!=0 ? ORDER/2 : either (needoffsetlo, ORDER/2+1, ORDER/2);
    ivect3 const offsethi =
      ORDER%2!=0 ? ORDER/2 : either (needoffsethi, ORDER/2+1, ORDER/2);
    
    
    
    if (not regbbox.expand(offsetlo, offsethi).is_contained_in(srcbbox) or
        not regbbox                           .is_contained_in(dstbbox))
    {
      cerr << "ORDER=" << ORDER << "\n"
           << "offsetlo=" << offsetlo << "\n"
           << "offsethi=" << offsethi << "\n"
           << "regbbox=" << regbbox << "\n"
           << "dstbbox=" << dstbbox << "\n"
           << "regbbox.expand=" << regbbox.expand(offsetlo, offsethi) << "\n"
           << "srcbbox=" << srcbbox << "\n";
      CCTK_WARN (0, "Internal error: region extent is not contained in array extent");
    }
    
    if (any (srcext != gdata::allocated_memory_shape(srcbbox.shape() / srcbbox.stride()) or
             dstext != gdata::allocated_memory_shape(dstbbox.shape() / dstbbox.stride())))
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
    
    
    
    // size_t const srcdi = SRCIND3(1,0,0) - SRCIND3(0,0,0);
    assert (SRCIND3(1,0,0) - SRCIND3(0,0,0) == 1);
    size_t const srcdi = 1;
    size_t const srcdj = SRCIND3(0,1,0) - SRCIND3(0,0,0);
    size_t const srcdk = SRCIND3(0,0,1) - SRCIND3(0,0,0);
    
    
    
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
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,0,0> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8001;
    goto l900;
    
    // kernel
  l8001:
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,1,0,0> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
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
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,1,0> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8011;
    goto l901;
    
    // kernel
  l8011:
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,1,1,0> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
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
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,0,1> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8101;
    goto l910;
    
    // kernel
  l8101:
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,1,0,1> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
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
      interp3<T,ORDER,0,1,1> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8111;
    goto l911;
    
    // kernel
  l8111:
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,1,1,1> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
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
  
  
  
#define TYPECASE(N,T)                                           \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_cc_rf2<T,0> (T const * restrict const src,      \
                             ivect3 const & restrict srcext,    \
                             T * restrict const dst,            \
                             ivect3 const & restrict dstext,    \
                             ibbox3 const & restrict srcbbox,   \
                             ibbox3 const & restrict dstbbox,   \
                             ibbox3 const & restrict regbbox);  \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_cc_rf2<T,1> (T const * restrict const src,      \
                             ivect3 const & restrict srcext,    \
                             T * restrict const dst,            \
                             ivect3 const & restrict dstext,    \
                             ibbox3 const & restrict srcbbox,   \
                             ibbox3 const & restrict dstbbox,   \
                             ibbox3 const & restrict regbbox);  \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_cc_rf2<T,2> (T const * restrict const src,      \
                             ivect3 const & restrict srcext,    \
                             T * restrict const dst,            \
                             ivect3 const & restrict dstext,    \
                             ibbox3 const & restrict srcbbox,   \
                             ibbox3 const & restrict dstbbox,   \
                             ibbox3 const & restrict regbbox);  \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_cc_rf2<T,3> (T const * restrict const src,      \
                             ivect3 const & restrict srcext,    \
                             T * restrict const dst,            \
                             ivect3 const & restrict dstext,    \
                             ibbox3 const & restrict srcbbox,   \
                             ibbox3 const & restrict dstbbox,   \
                             ibbox3 const & restrict regbbox);  \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_cc_rf2<T,4> (T const * restrict const src,      \
                             ivect3 const & restrict srcext,    \
                             T * restrict const dst,            \
                             ivect3 const & restrict dstext,    \
                             ibbox3 const & restrict srcbbox,   \
                             ibbox3 const & restrict dstbbox,   \
                             ibbox3 const & restrict regbbox);  \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_cc_rf2<T,5> (T const * restrict const src,      \
                             ivect3 const & restrict srcext,    \
                             T * restrict const dst,            \
                             ivect3 const & restrict dstext,    \
                             ibbox3 const & restrict srcbbox,   \
                             ibbox3 const & restrict dstbbox,   \
                             ibbox3 const & restrict regbbox);

#include "typecase.hh"
#undef TYPECASE
  
  
  
} // namespace CarpetLib
