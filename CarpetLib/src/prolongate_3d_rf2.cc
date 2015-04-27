#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "vectors.h"
#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;



namespace CarpetLib {
  
  
  
#define SRCIND3(i,j,k)                          \
  index3 (i, j, k,                              \
          srcipadext, srcjpadext, srckpadext,   \
          srciext, srcjext, srckext)
#define DSTIND3(i,j,k)                          \
  index3 (i, j, k,                              \
          dstipadext, dstjpadext, dstkpadext,   \
          dstiext, dstjext, dstkext)
#define SRCOFF3(i,j,k)                          \
  offset3 (i, j, k,                             \
           srcipadext, srcjpadext, srckpadext)
#define DSTOFF3(i,j,k)                          \
  offset3 (i, j, k,                             \
           dstipadext, dstjpadext, dstkpadext)
  
  
  
  namespace coeffs_3d_rf2 {

  // 1D interpolation coefficients
  
  template<typename RT, int ORDER>
  struct coeffs1d {
    static const RT coeffs[]
    CCTK_ATTRIBUTE_ALIGNED(CCTK_REAL_VEC_SIZE * CCTK_REAL_PRECISION);
    
    static ptrdiff_t const ncoeffs = ORDER+1;
    static ptrdiff_t const imin    = - ncoeffs/2 + 1;
    static ptrdiff_t const imax    = imin + ncoeffs;
    
    static inline
    RT const&
    get (ptrdiff_t const i)
    {
      static_assert (ncoeffs == sizeof coeffs / sizeof *coeffs,
                     "coefficient array has wrong size");
#ifdef CARPET_DEBUG
      assert (i>=imin and i<imax);
#endif
      ptrdiff_t const j = i - imin;
#ifdef CARPET_DEBUG
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
      
      static_assert (ncoeffs == sizeof coeffs / sizeof *coeffs,
                     "coefficient array has wrong size");
      
      // Do not test integer operators (they should be disabled
      // anyway)
      if (fabs (RT(0.5) - 0.5) > 1.0e-5) return;
      
      // Test all orders
      bool error = false;
      for (int n=0; n<=ORDER; ++n) {
        RT res = RT(0.0);
        for (ptrdiff_t i=imin; i<imax; ++i) {
          RT const x = RT(i) - RT(0.5);
          RT const y = ipow (x, n);
          res += get(i) * y;
        }
        RT const x0 = RT(0.0);
        RT const y0 = ipow (x0, n);
        // Allow losing 3 digits:
        CCTK_REAL const eps = RT(1.0e+3) * numeric_limits<RT>::epsilon();
        if (not (fabs (res - y0) < eps)) {
          RT rt;
          ostringstream buf;
          buf << "Error in prolongate_3d_rf2::coeffs_3d_rf2\n"
              << "   RT=" << typestring(rt) << "\n"
              << "   ORDER=" << ORDER << "\n"
              << "   n=" << n << "\n"
              << "   y0=" << y0 << ", res=" << res;
          CCTK_WARN (CCTK_WARN_ALERT, buf.str().c_str());
          error = true;
        }
      } // for n
      if (error) {
        CCTK_WARN (CCTK_WARN_ABORT, "Aborting.");
      }
    }
  };
  
  
  
#define TYPECASE(N,RT)                          \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,1>::coeffs[] = {         \
    +1 / RT(2.0),                               \
    +1 / RT(2.0)                                \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,3>::coeffs[] = {         \
    -1 / RT(16.0),                              \
    +9 / RT(16.0),                              \
    +9 / RT(16.0),                              \
    -1 / RT(16.0)                               \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,5>::coeffs[] = {         \
    +  3 / RT(256.0),                           \
    - 25 / RT(256.0),                           \
    +150 / RT(256.0),                           \
    +150 / RT(256.0),                           \
    - 25 / RT(256.0),                           \
    +  3 / RT(256.0)                            \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,7>::coeffs[] = {         \
      -   5 / RT(2048.0),                       \
      +  49 / RT(2048.0),                       \
      - 245 / RT(2048.0),                       \
      +1225 / RT(2048.0),                       \
      +1225 / RT(2048.0),                       \
      - 245 / RT(2048.0),                       \
      +  49 / RT(2048.0),                       \
      -   5 / RT(2048.0)                        \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,9>::coeffs[] = {         \
      +   35 / RT(65536.0),                     \
      -  405 / RT(65536.0),                     \
      + 2268 / RT(65536.0),                     \
      - 8820 / RT(65536.0),                     \
      +39690 / RT(65536.0),                     \
      +39690 / RT(65536.0),                     \
      - 8820 / RT(65536.0),                     \
      + 2268 / RT(65536.0),                     \
      -  405 / RT(65536.0),                     \
      +   35 / RT(65536.0)                      \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1d<RT,11>::coeffs[] = {        \
      -    63 / RT(524288.0),                   \
      +   847 / RT(524288.0),                   \
      -  5445 / RT(524288.0),                   \
      + 22869 / RT(524288.0),                   \
      - 76230 / RT(524288.0),                   \
      +320166 / RT(524288.0),                   \
      +320166 / RT(524288.0),                   \
      - 76230 / RT(524288.0),                   \
      + 22869 / RT(524288.0),                   \
      -  5445 / RT(524288.0),                   \
      +   847 / RT(524288.0),                   \
      -    63 / RT(524288.0)                    \
  };
  
#define CARPET_NO_COMPLEX
#define CARPET_NO_INT
#include "typecase.hh"
#undef TYPECASE
    
    
    
    extern "C"
    void CarpetLib_test_prolongate_3d_rf2 (CCTK_ARGUMENTS)
    {
      DECLARE_CCTK_ARGUMENTS;
      
#ifdef CCTK_REAL_PRECISION_4
#  define TYPECASE(N,RT)                        \
      coeffs1d<RT,1>::test();                   \
      coeffs1d<RT,3>::test();                   \
      coeffs1d<RT,5>::test();                   \
      coeffs1d<RT,7>::test();                   \
      coeffs1d<RT,9>::test();
#else
#  define TYPECASE(N,RT)                        \
      coeffs1d<RT,1>::test();                   \
      coeffs1d<RT,3>::test();                   \
      coeffs1d<RT,5>::test();                   \
      coeffs1d<RT,7>::test();                   \
      coeffs1d<RT,9>::test();                   \
      coeffs1d<RT,11>::test();
#endif
#define CARPET_NO_COMPLEX
#define CARPET_NO_INT
#include "typecase.hh"
#undef TYPECASE
    }
    
  } // namespace coeffs_3d_rf2
  
  using namespace coeffs_3d_rf2;
  
  
  
  // 0D "interpolation"
  template <typename T, int ORDER>
  static inline
  T const&
  interp0 (T const * restrict const p)
  {
    return * p;
  }
  
  // 1D interpolation
  template <typename T, int ORDER, int di>
  static inline
  T
  interp1 (T const * restrict const p,
           ptrdiff_t const d1)
  {
    typedef typeprops<T> typ;
    typedef typename typ::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    if (di == 0) {
      return interp0<T,ORDER> (p);
    } else {
      if (d1 == 1) {
#if 0
        T res = typ::fromreal (0);
        for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
          res += coeffs::get(i) * interp0<T,ORDER> (p + i);
        }
        return res;
#endif
        typedef vecprops<T> VP;
        typedef typename VP::vector_t VT;
        ptrdiff_t i = coeffs::imin;
        T res = typ::fromreal (0);
        if (coeffs::ncoeffs >= ptrdiff_t(VP::size())) {
          VT vres =
            VP::mul(VP::load(typ::fromreal(coeffs::get(i))),
                    VP::loadu(interp0<T,ORDER> (p + i)));
          i += VP::size();
#if defined(__INTEL_COMPILER)
          // Unroll the loop manually to help the Intel compiler
          // (This manual unrolling hurts with other compilers, e.g. PGI)
          assert (coeffs::ncoeffs / VP::size() <= 12);
          switch (coeffs::ncoeffs / VP::size()) {
            // Note that all case statements fall through
          case 12:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 11:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 10:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 9:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 8:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 7:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 6:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 5:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 4:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 3:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          case 2:
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
            i += VP::size();
          }
#else
          for (; i + VP::size() <= ptrdiff_t(coeffs::imax); i += VP::size()) {
            vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                            VP::loadu(interp0<T,ORDER> (p + i)),
                            vres);
          }
#endif
          for (int d=0; d<ptrdiff_t(VP::size()); ++d) {
            res += VP::elt(vres,d);
          }
        }
        assert (i == (ptrdiff_t(coeffs::imax) -
                      ptrdiff_t(coeffs::ncoeffs % VP::size())));
        for (i = coeffs::imax - coeffs::ncoeffs % VP::size();
             i < coeffs::imax;
             ++ i)
        {
          res += coeffs::get(i) * interp0<T,ORDER> (p + i*d1);
        }
        return res;
      } else {
        T res = typ::fromreal (0);
        for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
          res += coeffs::get(i) * interp0<T,ORDER> (p + i*d1);
        }
        return res;
      }
    }
  }
  
  // 2D interpolation
  template <typename T, int ORDER, int di, int dj>
  static inline
  T
  interp2 (T const * restrict const p,
           ptrdiff_t const d1,
           ptrdiff_t const d2)
  {
    typedef typeprops<T> typ;
    typedef typename typ::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    if (dj == 0) {
      return interp1<T,ORDER,di> (p, d1);
    } else {
      T res = typeprops<T>::fromreal (0);
      for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
        res += coeffs::get(i) * interp1<T,ORDER,di> (p + i*d2, d1);
      }
      return res;
    }
  }
  
  // 3D interpolation
  template <typename T, int ORDER, int di, int dj, int dk>
  static inline
  T
  interp3 (T const * restrict const p,
           ptrdiff_t const d1,
           ptrdiff_t const d2,
           ptrdiff_t const d3)
  {
    typedef typeprops<T> typ;
    typedef typename typ::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    if (dk == 0) {
      return interp2<T,ORDER,di,dj> (p, d1, d2);
    } else {
      T res = typeprops<T>::fromreal (0);
      for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
        res += coeffs::get(i) * interp2<T,ORDER,di,dj> (p + i*d3, d1, d2);
      }
      return res;
    }
  }
  
  
  
  // Check interpolation index ranges
  template <typename T, int ORDER>
  static inline
  void
  check_indices0 ()
  {
  }
  
  template <typename T, int ORDER, int di>
  static inline
  void
  check_indices1 (ptrdiff_t const is,
                  ptrdiff_t const srciext)
  {
    static_assert (di==0 or di==1, "di must be 0 or 1");
#ifdef CARPET_DEBUG
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    if (di == 0) {
      assert (is >= 0);
      assert (is < srciext);
    } else {
      assert (is + coeffs::imin >= 0);
      assert (is + coeffs::imax <= srciext);
    }
    check_indices0<T,ORDER> ();
#endif
  }
  
  template <typename T, int ORDER, int di, int dj>
  static inline
  void
  check_indices2 (ptrdiff_t const is,
                  ptrdiff_t const js,
                  ptrdiff_t const srciext,
                  ptrdiff_t const srcjext)
  {
    static_assert (dj==0 or dj==1, "dj must be 0 or 1");
#ifdef CARPET_DEBUG
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    if (dj == 0) {
      assert (js >= 0);
      assert (js < srcjext);
    } else {
      assert (js + coeffs::imin >= 0);
      assert (js + coeffs::imax <= srcjext);
    }
    check_indices1<T,ORDER,di> (is, srciext);
#endif
  }
  
  template <typename T, int ORDER, int di, int dj, int dk>
  static inline
  void
  check_indices3 (ptrdiff_t const is,
                  ptrdiff_t const js,
                  ptrdiff_t const ks,
                  ptrdiff_t const srciext,
                  ptrdiff_t const srcjext,
                  ptrdiff_t const srckext)
  {
    static_assert (dk==0 or dk==1, "dk must be 0 or 1");
#ifdef CARPET_DEBUG
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    if (dk == 0) {
      assert (ks >= 0);
      assert (ks < srckext);
    } else {
      assert (ks + coeffs::imin >= 0);
      assert (ks + coeffs::imax <= srckext);
    }
    check_indices2<T,ORDER,di,dj> (is,js, srciext,srcjext);
#endif
  }
  
  
  
  template <typename T, int ORDER>
  static void interp_driver_i (T const * restrict const src,
                               ivect3 const & restrict srcpadext,
                               ivect3 const & restrict srcext,
                               T * restrict const dst,
                               ivect3 const & restrict dstpadext,
                               ivect3 const & restrict dstext,
                               ivect3 const & restrict srcoff,
                               ivect3 const & restrict dstoff,
                               ivect3 const & restrict regext)
  {
    ptrdiff_t const srcipadext = srcpadext[0];
    ptrdiff_t const srcjpadext = srcpadext[1];
    ptrdiff_t const srckpadext = srcpadext[2];
    
    ptrdiff_t const dstipadext = dstpadext[0];
    ptrdiff_t const dstjpadext = dstpadext[1];
    ptrdiff_t const dstkpadext = dstpadext[2];
    
    ptrdiff_t const srciext = srcext[0];
    ptrdiff_t const srcjext = srcext[1];
    ptrdiff_t const srckext = srcext[2];
    
    ptrdiff_t const dstiext = dstext[0];
    ptrdiff_t const dstjext = dstext[1];
    ptrdiff_t const dstkext = dstext[2];
    
    ptrdiff_t const regiext = regext[0];
    ptrdiff_t const regjext = regext[1];
    ptrdiff_t const regkext = regext[2];
    
    ptrdiff_t const srcioff = srcoff[0];
    ptrdiff_t const dstioff = dstoff[0];
    
    ptrdiff_t const fi = srcioff % 2;
    ptrdiff_t const i0 = srcioff / 2;
    
    // ptrdiff_t const srcdi = SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
    ptrdiff_t const srcdi = 1;
    assert (srcdi == SRCOFF3(1,0,0) - SRCOFF3(0,0,0));
    
    // Loop over fine region
    for (ptrdiff_t k = 0; k < regkext; ++k) {
      for (ptrdiff_t j = 0; j < regjext; ++j) {
        
        ptrdiff_t i = 0;
        ptrdiff_t is = i0;
        ptrdiff_t id = dstioff;
        
        if (fi != 0 && i < regiext) {
          check_indices1<T,ORDER,1> (is, srciext);
          dst[DSTIND3(id,j,k)] =
            interp1<T,ORDER,1> (& src[SRCIND3(is,j,k)], srcdi);
          i = i+1;
          id = id+1;
          is = is+1;
        }
        
        while (i+1 < regiext) {
          check_indices1<T,ORDER,0> (is, srciext);
          dst[DSTIND3(id,j,k)] =
            interp1<T,ORDER,0> (& src[SRCIND3(is,j,k)], srcdi);
          i = i+1;
          id = id+1;
          
          check_indices1<T,ORDER,1> (is, srciext);
          dst[DSTIND3(id,j,k)] =
            interp1<T,ORDER,1> (& src[SRCIND3(is,j,k)], srcdi);
          i = i+1;
          id = id+1;
          is = is+1;
        }
        
        if (i < regiext) {
          check_indices1<T,ORDER,0> (is, srciext);
          dst[DSTIND3(id,j,k)] =
            interp1<T,ORDER,0> (& src[SRCIND3(is,j,k)], srcdi);
          i = i+1;
          id = id+1;
        }
        
        assert (i == regiext);
      }
    }
  }
  
  
  
  template <typename T, int ORDER>
  static void interp_driver_j (T const * restrict const src,
                               ivect3 const & restrict srcpadext,
                               ivect3 const & restrict srcext,
                               T * restrict const dst,
                               ivect3 const & restrict dstpadext,
                               ivect3 const & restrict dstext,
                               ivect3 const & restrict srcoff,
                               ivect3 const & restrict dstoff,
                               ivect3 const & restrict regext)
  {
    ptrdiff_t const srcipadext = srcpadext[0];
    ptrdiff_t const srcjpadext = srcpadext[1];
    ptrdiff_t const srckpadext = srcpadext[2];
    
    ptrdiff_t const dstipadext = dstpadext[0];
    ptrdiff_t const dstjpadext = dstpadext[1];
    ptrdiff_t const dstkpadext = dstpadext[2];
    
    ptrdiff_t const srciext = srcext[0];
    ptrdiff_t const srcjext = srcext[1];
    ptrdiff_t const srckext = srcext[2];
    
    ptrdiff_t const dstiext = dstext[0];
    ptrdiff_t const dstjext = dstext[1];
    ptrdiff_t const dstkext = dstext[2];
    
    ptrdiff_t const regiext = regext[0];
    ptrdiff_t const regjext = regext[1];
    ptrdiff_t const regkext = regext[2];
    
    ptrdiff_t const srcjoff = srcoff[1];
    ptrdiff_t const dstjoff = dstoff[1];
    
    ptrdiff_t const fj = srcjoff % 2;
    ptrdiff_t const j0 = srcjoff / 2;
    
    ptrdiff_t const srcdj = SRCOFF3(0,1,0) - SRCOFF3(0,0,0);
    
    // Loop over fine region
    for (ptrdiff_t k = 0; k < regkext; ++k) {
      
      ptrdiff_t j = 0;
      ptrdiff_t js = j0;
      ptrdiff_t jd = dstjoff;
      
      if (fj != 0 && j < regjext) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,1> (js, srcjext);
          dst[DSTIND3(i,jd,k)] =
            interp1<T,ORDER,1> (& src[SRCIND3(i,js,k)], srcdj);
        }
        j = j+1;
        jd = jd+1;
        js = js+1;
      }
      
      while (j+1 < regjext) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,0> (js, srcjext);
          dst[DSTIND3(i,jd,k)] =
            interp1<T,ORDER,0> (& src[SRCIND3(i,js,k)], srcdj);
        }
        j = j+1;
        jd = jd+1;
        
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,1> (js, srcjext);
          dst[DSTIND3(i,jd,k)] =
            interp1<T,ORDER,1> (& src[SRCIND3(i,js,k)], srcdj);
        }
        j = j+1;
        jd = jd+1;
        js = js+1;
      }
      
      if (j < regjext) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,0> (js, srcjext);
          dst[DSTIND3(i,jd,k)] =
            interp1<T,ORDER,0> (& src[SRCIND3(i,js,k)], srcdj);
        }
        j = j+1;
        jd = jd+1;
      }
      
      assert (j == regjext);
    }
  }
  
  
  
  template <typename T, int ORDER>
  static void interp_driver_k (T const * restrict const src,
                               ivect3 const & restrict srcpadext,
                               ivect3 const & restrict srcext,
                               T * restrict const dst,
                               ivect3 const & restrict dstpadext,
                               ivect3 const & restrict dstext,
                               ivect3 const & restrict srcoff,
                               ivect3 const & restrict dstoff,
                               ivect3 const & restrict regext)
  {
    ptrdiff_t const srcipadext = srcpadext[0];
    ptrdiff_t const srcjpadext = srcpadext[1];
    ptrdiff_t const srckpadext = srcpadext[2];
    
    ptrdiff_t const dstipadext = dstpadext[0];
    ptrdiff_t const dstjpadext = dstpadext[1];
    ptrdiff_t const dstkpadext = dstpadext[2];
    
    ptrdiff_t const srciext = srcext[0];
    ptrdiff_t const srcjext = srcext[1];
    ptrdiff_t const srckext = srcext[2];
    
    ptrdiff_t const dstiext = dstext[0];
    ptrdiff_t const dstjext = dstext[1];
    ptrdiff_t const dstkext = dstext[2];
    
    ptrdiff_t const regiext = regext[0];
    ptrdiff_t const regjext = regext[1];
    ptrdiff_t const regkext = regext[2];
    
    ptrdiff_t const srckoff = srcoff[2];
    ptrdiff_t const dstkoff = dstoff[2];
    
    ptrdiff_t const fk = srckoff % 2;
    ptrdiff_t const k0 = srckoff / 2;
    
    ptrdiff_t const srcdk = SRCOFF3(0,0,1) - SRCOFF3(0,0,0);
    
    // Loop over fine region
    ptrdiff_t k = 0;
    ptrdiff_t ks = k0;
    ptrdiff_t kd = dstkoff;
    
    if (fk != 0 && k < regkext) {
      for (ptrdiff_t j = 0; j < regjext; ++j) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,1> (ks, srckext);
          dst[DSTIND3(i,j,kd)] =
            interp1<T,ORDER,1> (& src[SRCIND3(i,j,ks)], srcdk);
        }
      }
      k = k+1;
      kd = kd+1;
      ks = ks+1;
    }
    
    while (k+1 < regkext) {
      for (ptrdiff_t j = 0; j < regjext; ++j) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,0> (ks, srckext);
          dst[DSTIND3(i,j,kd)] =
            interp1<T,ORDER,0> (& src[SRCIND3(i,j,ks)], srcdk);
        }
      }
      k = k+1;
      kd = kd+1;
      
      for (ptrdiff_t j = 0; j < regjext; ++j) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,1> (ks, srckext);
          dst[DSTIND3(i,j,kd)] =
            interp1<T,ORDER,1> (& src[SRCIND3(i,j,ks)], srcdk);
        }
      }
      k = k+1;
      kd = kd+1;
      ks = ks+1;
    }
    
    if (k < regkext) {
      for (ptrdiff_t j = 0; j < regjext; ++j) {
        for (ptrdiff_t i = 0; i < regiext; ++i) {
          check_indices1<T,ORDER,0> (ks, srckext);
          dst[DSTIND3(i,j,kd)] =
            interp1<T,ORDER,0> (& src[SRCIND3(i,j,ks)], srcdk);
        }
      }
      k = k+1;
      kd = kd+1;
    }
    
    assert (k == regkext);
  }
  
  
  
  template <typename T, int ORDER>
  static void interp_driver_3d (T const * restrict const src,
                                ivect3 const & restrict srcpadext,
                                ivect3 const & restrict srcext,
                                T * restrict const dst,
                                ivect3 const & restrict dstpadext,
                                ivect3 const & restrict dstext,
                                ivect3 const & restrict srcoff,
                                ivect3 const & restrict dstoff,
                                ivect3 const & restrict regext)
  {
    ptrdiff_t const srcipadext = srcpadext[0];
    ptrdiff_t const srcjpadext = srcpadext[1];
    ptrdiff_t const srckpadext = srcpadext[2];
    
    ptrdiff_t const dstipadext = dstpadext[0];
    ptrdiff_t const dstjpadext = dstpadext[1];
    ptrdiff_t const dstkpadext = dstpadext[2];
    
    ptrdiff_t const srciext = srcext[0];
    ptrdiff_t const srcjext = srcext[1];
    ptrdiff_t const srckext = srcext[2];
    
    ptrdiff_t const dstiext = dstext[0];
    ptrdiff_t const dstjext = dstext[1];
    ptrdiff_t const dstkext = dstext[2];
    
    ptrdiff_t const regiext = regext[0];
    ptrdiff_t const regjext = regext[1];
    ptrdiff_t const regkext = regext[2];
    
    ptrdiff_t const srcioff = srcoff[0];
    ptrdiff_t const srcjoff = srcoff[1];
    ptrdiff_t const srckoff = srcoff[2];
    
    ptrdiff_t const dstioff = dstoff[0];
    ptrdiff_t const dstjoff = dstoff[1];
    ptrdiff_t const dstkoff = dstoff[2];
    
    
    
    ptrdiff_t const fi = srcioff % 2;
    ptrdiff_t const fj = srcjoff % 2;
    ptrdiff_t const fk = srckoff % 2;
    
    ptrdiff_t const i0 = srcioff / 2;
    ptrdiff_t const j0 = srcjoff / 2;
    ptrdiff_t const k0 = srckoff / 2;
    
    
    
    // ptrdiff_t const srcdi = SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
    ptrdiff_t const srcdi = 1;
    assert (srcdi == SRCOFF3(1,0,0) - SRCOFF3(0,0,0));
    ptrdiff_t const srcdj = SRCOFF3(0,1,0) - SRCOFF3(0,0,0);
    ptrdiff_t const srcdk = SRCOFF3(0,0,1) - SRCOFF3(0,0,0);
    
    
    
    // Loop over fine region
    // Label scheme: l 8 fk fj fi
    
    ptrdiff_t is, js, ks;
    ptrdiff_t id, jd, kd;
    ptrdiff_t i, j, k;
    
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
    check_indices3<T,ORDER,0,0,0> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,0,0> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8001;
    goto l900;
    
    // kernel
  l8001:
    check_indices3<T,ORDER,1,0,0> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,1,0,0> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
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
    check_indices3<T,ORDER,0,1,0> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,1,0> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8011;
    goto l901;
    
    // kernel
  l8011:
    check_indices3<T,ORDER,1,1,0> (is,js,ks, srciext, srcjext, srckext);
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
    check_indices3<T,ORDER,0,0,1> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,0,1> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8101;
    goto l910;
    
    // kernel
  l8101:
    check_indices3<T,ORDER,1,0,1> (is,js,ks, srciext, srcjext, srckext);
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
    check_indices3<T,ORDER,0,1,1> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,1,1> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8111;
    goto l911;
    
    // kernel
  l8111:
    check_indices3<T,ORDER,1,1,1> (is,js,ks, srciext, srcjext, srckext);
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
  
  
  
  template <typename T, int ORDER>
  void
  prolongate_3d_rf2 (T const * restrict const src,
                     ivect3 const & restrict srcpadext,
                     ivect3 const & restrict srcext,
                     T * restrict const dst,
                     ivect3 const & restrict dstpadext,
                     ivect3 const & restrict dstext,
                     ibbox3 const & restrict srcbbox,
                     ibbox3 const & restrict dstbbox,
                     ibbox3 const & restrict,
                     ibbox3 const & restrict regbbox,
                     void * extraargs)
  {
    assert (not extraargs);
    
    static_assert (ORDER>=0 and ORDER % 2 == 1,
                   "ORDER must be non-negative and odd");
    
    typedef typename typeprops<T>::real RT;
    coeffs1d<RT,ORDER>::test();
    
    
    
    if (any (srcbbox.stride() <= regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_ERROR ("Internal error: strides disagree");
    }
    
    if (any (srcbbox.stride() != reffact2 * dstbbox.stride())) {
      CCTK_ERROR ("Internal error: source strides are not twice the destination strides");
    }
    
    if (any (srcbbox.lower() % srcbbox.stride() != 0)) {
      CCTK_ERROR ("Internal error: source bbox is not aligned with vertices");
    }
    if (any (dstbbox.lower() % dstbbox.stride() != 0)) {
      CCTK_ERROR ("Internal error: destination bbox is not aligned with vertices");
    }
    if (any (regbbox.lower() % regbbox.stride() != 0)) {
      CCTK_ERROR ("Internal error: prolongation region bbox is not aligned with vertices");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_ERROR ("Internal error: region extent is empty");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all ((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
    ivect3 const srcoff =
      (regbbox.lower() - srcbbox.lower()) / regbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
    ivect3 const dstoff =
      (regbbox.lower() - dstbbox.lower()) / regbbox.stride();
    
    
    
    bvect3 const needoffsetlo = srcoff % reffact2 != 0;
    bvect3 const needoffsethi = (srcoff + regext - 1) % reffact2 != 0;
    ivect3 const offsetlo =
      either (needoffsetlo, ORDER/2*reffact2+1,
              either (regext > 1, ORDER/2*reffact2, 0));
    ivect3 const offsethi =
      either (needoffsethi, ORDER/2*reffact2+1,
              either (regext > 1, ORDER/2*reffact2, 0));
    ibbox3 const src_regbbox_fine = regbbox.expand(offsetlo, offsethi);
    ibbox3 const src_regbbox = src_regbbox_fine.expanded_for (srcbbox);
    assert (src_regbbox.expanded_for (src_regbbox_fine) == src_regbbox_fine);
    ibbox3 const & dst_regbbox = regbbox;
    if (not src_regbbox.is_contained_in(srcbbox)) {
      CCTK_ERROR ("Internal error: source region extent is not contained in array extent");
    }
    if (not dst_regbbox.is_contained_in(dstbbox)) {
      CCTK_ERROR ("Internal error: destination region extent is not contained in array extent");
    }
    
    
    
#warning "TODO"
#if 1
    
    interp_driver_3d<T,ORDER> (src, srcpadext, srcext,
                               dst, dstpadext, dstext,
                               srcoff, dstoff, regext);
    
#else
    
    // Copy in the original source
    // TODO: Don't copy
    ivect3 srcdataext = src_regbbox.shape() / src_regbbox.stride();
    ivect3 srcdataoff = offsetlo;
    vector<T> srcdata (prod (srcdataext));
    {
      ptrdiff_t const srcipadext = srcpadext[0];
      ptrdiff_t const srcjpadext = srcpadext[1];
      ptrdiff_t const srckpadext = srcpadext[2];
      ptrdiff_t const srciext = srcext[0];
      ptrdiff_t const srcjext = srcext[1];
      ptrdiff_t const srckext = srcext[2];
      ptrdiff_t const srcioff = srcoff[0];
      ptrdiff_t const srcjoff = srcoff[1];
      ptrdiff_t const srckoff = srcoff[2];
      ptrdiff_t const i0 = srcioff / 2 - srcdataoff[0];
      ptrdiff_t const j0 = srcjoff / 2 - srcdataoff[1];
      ptrdiff_t const k0 = srckoff / 2 - srcdataoff[2];
      assert (i0 >= 0 && i0 + srcdataext[0] < srcext[0]);
      assert (j0 >= 0 && j0 + srcdataext[1] < srcext[1]);
      assert (k0 >= 0 && k0 + srcdataext[2] < srcext[2]);
      assert (all (srcdataext == src_regbbox.shape() / src_regbbox.stride()));
      assert (all (srcdataext <= srcext));
      for (int k=0; k<srcdataext[2]; ++k) {
        for (int j=0; j<srcdataext[1]; ++j) {
          for (int i=0; i<srcdataext[0]; ++i) {
            srcdata[i + srcdataext[0] * (j + srcdataext[1] * k)] =
              src[SRCIND3(i0+i, j0+j, k0+k)];
          }
        }
      }
    }
    
    bvect3 did_handle_direction = false;
    for (int ndirs = 0; ndirs < 3; ++ndirs) {
      // Handle longest direction last, since this is the most
      // efficient, as prolongation increases the number of points
      int dir = -1;
      for (int d = 0; d < 3; ++d) {
        if (!did_handle_direction[d] &&
            (dir < 0 || srcdataext[d] < srcdataext[dir]))
        {
          dir = d;
        }
      }
      
      // Allocate destination array
      ivect3 dstdataext = srcdataext;
      dstdataext[dir] = (dst_regbbox.shape() / dst_regbbox.stride())[dir];
      ivect3 dstdataoff = srcdataoff;
      dstdataoff[dir] = 0;
      vector<T> dstdata (prod (dstdataext));
      
      // Prolongate
      void (*interp_driver) (T const * restrict const src,
                             ivect3 const & restrict srcpadext,
                             ivect3 const & restrict srcext,
                             T * restrict const dst,
                             ivect3 const & restrict dstpadext,
                             ivect3 const & restrict dstext,
                             ivect3 const & restrict srcoff,
                             ivect3 const & restrict dstoff,
                             ivect3 const & restrict regext) = 0;
      switch (dir) {
      case 0: interp_driver = interp_driver_i<T,ORDER>; break;
      case 1: interp_driver = interp_driver_j<T,ORDER>; break;
      case 2: interp_driver = interp_driver_k<T,ORDER>; break;
      }
      interp_driver (&srcdata.front(), srcdataext, srcdataext,
                     &dstdata.front(), dstdataext, dstdataext,
                     srcdataoff, ivect3(0), dstdataext);
      
      // Switch arrays
      swap (srcdata, dstdata);  // should really be a move
      srcdataext = dstdataext;
      srcdataoff = dstdataoff;
      did_handle_direction[dir] = true;
    }
    
    // Copy out the result
    // TODO: Don't copy
    {
      ptrdiff_t const dstipadext = dstpadext[0];
      ptrdiff_t const dstjpadext = dstpadext[1];
      ptrdiff_t const dstkpadext = dstpadext[2];
      ptrdiff_t const dstiext = dstext[0];
      ptrdiff_t const dstjext = dstext[1];
      ptrdiff_t const dstkext = dstext[2];
      ptrdiff_t const dstioff = dstoff[0];
      ptrdiff_t const dstjoff = dstoff[1];
      ptrdiff_t const dstkoff = dstoff[2];
      assert (all (srcdataext == dst_regbbox.shape() / dst_regbbox.stride()));
      assert (all (srcdataext <= dstext));
      for (int k=0; k<srcdataext[2]; ++k) {
        for (int j=0; j<srcdataext[1]; ++j) {
          for (int i=0; i<srcdataext[0]; ++i) {
            dst[DSTIND3(dstioff+i, dstjoff+j, dstkoff+k)] =
              srcdata[i + srcdataext[0] * (j + srcdataext[1] * k)];
          }
        }
      }
    }
    
#endif
  }
  
  
  
#define TYPECASE(N,T)                                           \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_rf2<T,1> (T const * restrict const src,         \
                          ivect3 const & restrict srcpadext,    \
                          ivect3 const & restrict srcext,       \
                          T * restrict const dst,               \
                          ivect3 const & restrict dstpadext,    \
                          ivect3 const & restrict dstext,       \
                          ibbox3 const & restrict srcbbox,      \
                          ibbox3 const & restrict dstbbox,      \
                          ibbox3 const & restrict,              \
                          ibbox3 const & restrict regbbox,      \
                          void * extraargs);                    \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_rf2<T,3> (T const * restrict const src,         \
                          ivect3 const & restrict srcpadext,    \
                          ivect3 const & restrict srcext,       \
                          T * restrict const dst,               \
                          ivect3 const & restrict dstpadext,    \
                          ivect3 const & restrict dstext,       \
                          ibbox3 const & restrict srcbbox,      \
                          ibbox3 const & restrict dstbbox,      \
                          ibbox3 const & restrict,              \
                          ibbox3 const & restrict regbbox,      \
                          void * extraargs);                    \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_rf2<T,5> (T const * restrict const src,         \
                          ivect3 const & restrict srcpadext,    \
                          ivect3 const & restrict srcext,       \
                          T * restrict const dst,               \
                          ivect3 const & restrict dstpadext,    \
                          ivect3 const & restrict dstext,       \
                          ibbox3 const & restrict srcbbox,      \
                          ibbox3 const & restrict dstbbox,      \
                          ibbox3 const & restrict,              \
                          ibbox3 const & restrict regbbox,      \
                          void * extraargs);                    \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_rf2<T,7> (T const * restrict const src,         \
                          ivect3 const & restrict srcpadext,    \
                          ivect3 const & restrict srcext,       \
                          T * restrict const dst,               \
                          ivect3 const & restrict dstpadext,    \
                          ivect3 const & restrict dstext,       \
                          ibbox3 const & restrict srcbbox,      \
                          ibbox3 const & restrict dstbbox,      \
                          ibbox3 const & restrict,              \
                          ibbox3 const & restrict regbbox,      \
                          void * extraargs);                    \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_rf2<T,9> (T const * restrict const src,         \
                          ivect3 const & restrict srcpadext,    \
                          ivect3 const & restrict srcext,       \
                          T * restrict const dst,               \
                          ivect3 const & restrict dstpadext,    \
                          ivect3 const & restrict dstext,       \
                          ibbox3 const & restrict srcbbox,      \
                          ibbox3 const & restrict dstbbox,      \
                          ibbox3 const & restrict,              \
                          ibbox3 const & restrict regbbox,      \
                          void * extraargs);                    \
                                                                \
  template                                                      \
  void                                                          \
  prolongate_3d_rf2<T,11> (T const * restrict const src,        \
                           ivect3 const & restrict srcpadext,   \
                           ivect3 const & restrict srcext,      \
                           T * restrict const dst,              \
                           ivect3 const & restrict dstpadext,   \
                           ivect3 const & restrict dstext,      \
                           ibbox3 const & restrict srcbbox,     \
                           ibbox3 const & restrict dstbbox,     \
                           ibbox3 const & restrict,             \
                           ibbox3 const & restrict regbbox,     \
                           void * extraargs);

#define CARPET_NO_INT
#include "typecase.hh"
#undef TYPECASE
  
  
  
} // namespace CarpetLib
