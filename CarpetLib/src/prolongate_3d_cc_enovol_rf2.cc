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



//
// Grid point locations and their indices:
//
// global   0   4      12      20      28   |
// local    |   0       1       2       3   |
//          |   C       C       C       C   |
//          | f   f   f   f   f   f   f   f |
// local    | 0   1   2   3   4   5   6   7 |
// global   0 2   6  10  14  18  22  26  30 |
//
// Interpolation with even interpolation order:
//    offset zero (fine index = 2 * coarse index):
//       [1]
//    offset one (fine index = 2 * coarse index + 1):
//       [1]
// Interpolation with odd interpolation order:
//    offset zero:
//       [1 1 0]/2
//    offset one:
//       [0 1 1]/2
// The centres of these stencils are located at the coarse grid point
// corresponding to the fine grid point (the interpolation target)
// minus the offset. Example: fine grid 8 -> coarse grid 8, fine grid
// 12 -> also coarse grid 8.



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
           srciext, srcjext, srckext)
#define DSTOFF3(i,j,k)                          \
  offset3 (i, j, k,                             \
           dstiext, dstjext, dstkext)
  
  
    
  namespace coeffs_3d_cc_enovol_rf2 {
  
  // 1D interpolation coefficients
  
  template<typename RT, int ORDER, int OFFSET = 0>
  struct coeffs1dc {
    static const RT coeffs[];
  };
  
  template<typename RT, int ORDER>
  struct diffC1dc {
    static const RT coeffs[];
  };
  
  template<typename RT, int ORDER, int OFFSET = 0>
  struct coeffs1d {
    typedef coeffs1dc<RT,ORDER,OFFSET> coeffs_t;
    typedef diffC1dc<RT,ORDER> diff_t;
    
    static ptrdiff_t const ncoeffs = ORDER+1;
    static ptrdiff_t const imin    = OFFSET;
    static ptrdiff_t const imax    = imin + ncoeffs;
    // maximally left shifted stencil
    static ptrdiff_t const minimin    = -ORDER;
    // maximally right shifted stencil
    static ptrdiff_t const maximax    = ORDER+1;
    // maximal number of points needed
    static ptrdiff_t const maxncoeffs = maximax - minimin;
    
    
    static inline
    RT
    get (ptrdiff_t const i)
    {
      static_assert (ncoeffs == sizeof coeffs_t::coeffs / sizeof *coeffs_t::coeffs,
                     "coefficient array has wrong size");
#ifdef CARPET_DEBUG
      assert (i>=imin and i<imax);
#endif
      ptrdiff_t const j = i-imin;
#ifdef CARPET_DEBUG
      assert (j>=0 and j<ncoeffs);
#endif
      return coeffs_t::coeffs[j];
    }
    
    static inline
    RT
    diff (ptrdiff_t const i)
    {
      static_assert (ncoeffs == sizeof diff_t::coeffs / sizeof *diff_t::coeffs,
                     "coefficient array has wrong size");
#ifdef CARPET_DEBUG
      assert (i>=imin and i<imax);
#endif
      ptrdiff_t const j = i-imin;
#ifdef CARPET_DEBUG
      assert (j>=0 and j<ncoeffs);
#endif
      return diff_t::coeffs[j];
    }
    
    
    // Test coefficients
    static
    void
    test ()
    {
      static bool tested = false;
      if (tested) return;
      tested = true;
      
      static_assert (ncoeffs == sizeof coeffs_t::coeffs / sizeof *coeffs_t::coeffs,
                     "coefficient array has wrong size");
      
      // Test all orders
      bool error = false;
      for (int n=0; n<=ORDER; ++n) {
        RT res = RT(0.0);
        for (ptrdiff_t i=imin; i<imax; ++i) {
          RT const x0 = RT(CCTK_REAL(i));
          RT const x1 = RT(CCTK_REAL(i)+1.0);
          //cout << "x=" << 0.5*(x1+x0) << endl;
          RT const y = (ipow (x1, n+1) - ipow (x0, n+1)) / ((x1-x0)*(n+1));
          //cout << "y=" << y << endl;
          res += get(i) * y;
        }
        RT const x0 = RT(0.0);
        RT const x1 = RT(0.5);
        //cout << "x0=" << x0 << endl;
        RT const y0 = (ipow (x1, n+1) - ipow (x0, n+1)) / ((x1-x0)*(n+1));
        if (not (fabs (res - y0) < 1.0e-12)) {
          RT rt;
          ostringstream buf;
          buf << "Error in prolongate_3d_cc_rf2::coeffs_3d_cc_rf2\n"
              << "   RT=" << typestring(rt) << "\n"
              << "   ORDER=" << ORDER << "\n"
              << "   OFFSET=" << OFFSET << "\n"
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
  const RT coeffs1dc<RT,2,0>::coeffs[] = {      \
    +11 / RT(8.0),                              \
    - 4 / RT(8.0),                              \
    + 1 / RT(8.0)                               \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1dc<RT,2,-1>::coeffs[] = {     \
    + 1 / RT(8.0),                              \
    + 8 / RT(8.0),                              \
    - 1 / RT(8.0)                               \
  };                                            \
                                                \
  template<>                                    \
  const RT coeffs1dc<RT,2,-2>::coeffs[] = {     \
    - 1 / RT(8.0),                              \
    + 4 / RT(8.0),                              \
    + 5 / RT(8.0)                               \
  };                                            \
                                                \
                                                \
  template<>                                    \
  const RT diffC1dc<RT,2>::coeffs[] = {         \
    +1.0,                                       \
    -2.0,                                       \
    +1.0                                        \
  };                                            \
  
  
#define CARPET_NO_COMPLEX
#define CARPET_NO_INT
#include "typecase.hh"
#undef TYPECASE
  
  } // namespace coeffs_3d_cc_rf2_eno
  
  using namespace coeffs_3d_cc_enovol_rf2;
  
  /*
  template <typename T>
  static inline 
  T 
  minmod(const T a, const T b)
  {
     if (a * b < 0)
        return T(0);
     else if (fabs(a) < fabs(b))
        return a;
     else
        return b;
     return 0;
  }*/
  
  
  template <typename T, int ORDER, int di, int dj, int dk>
  static inline
  T
  interp3 (T const * restrict const p,
           size_t const d1,
           size_t const d2,
           size_t const d3)
  {
    if (d3 == 0) return *p; // end recursion when there is only a single point left

    static_assert (dk==0 or dk==1, "dk must be 0 or 1");
    typedef typename typeprops<T>::real RT;
    T res = typeprops<T>::fromreal (0);
    // get function values needed for the stencil
    T f[coeffs1d<RT,ORDER>::maxncoeffs];
    for (ptrdiff_t i=coeffs1d<RT,ORDER>::minimin; i<coeffs1d<RT,ORDER>::maximax; ++i) {
       f[i-coeffs1d<RT,ORDER>::minimin] = interp3<T,ORDER,0,di,dj> (p + i*d3, 0, d1, d2);
    }
    
    switch (ORDER)
    {
       case 2: {
	 const int shiftleft   = -2;
	 const int shiftcenter = -1;
	 const int shiftright  =  0;
	 typedef coeffs1d<RT,ORDER,shiftleft>   lcoeffs;
	 typedef coeffs1d<RT,ORDER,shiftcenter> ccoeffs;
	 typedef coeffs1d<RT,ORDER,shiftright>  rcoeffs;
	 T V[3] = { typeprops<T>::fromreal (0),
	            typeprops<T>::fromreal (0),
	            typeprops<T>::fromreal (0) };
	 // compute undivided differences for left-shifted stencil
	 for (ptrdiff_t i=lcoeffs::imin; i<lcoeffs::imax; ++i) {
	    V[0] += lcoeffs::diff(i) * f[i-lcoeffs::minimin];
	 }
	 // compute undivided differences for centered stencil
	 for (ptrdiff_t i=ccoeffs::imin; i<ccoeffs::imax; ++i) {
	    V[1] += ccoeffs::diff(i) * f[i-lcoeffs::minimin];
	 }
	 // compute undivided differences for right-shifted stencil
	 for (ptrdiff_t i=rcoeffs::imin; i<rcoeffs::imax; ++i) {
	    V[2] += rcoeffs::diff(i) * f[i-lcoeffs::minimin];
	 }
	 
	 // check that divided differences do not change sign: if so go back to first order!
	 if (V[0]*V[2] <= 0 || V[0]*V[1] <= 0 || V[1]*V[2] <= 0)
	 {
	   // switch back to zeroeth order!
           // TODO: come up with higher order fallback?
           res = f[-coeffs1d<RT,ORDER>::minimin];
	 } else {
           int min = 1;
           for (int i=0; i < 3; ++i)
              if (fabs(V[i]) < fabs(V[min])) min = i;
           
           switch (min) {
              case 0:
                 // use maximally left-shifted stencil since it is smoother
                 for (ptrdiff_t i=lcoeffs::imin; i<lcoeffs::imax; ++i) {
                    res += lcoeffs::get(i) * f[i-lcoeffs::minimin];
                 }
                 break;
              case 1:
                 // use centered stencil since it is smoother
                 for (ptrdiff_t i=ccoeffs::imin; i<ccoeffs::imax; ++i) {
                    res += ccoeffs::get(i) * f[i-lcoeffs::minimin];
                 }
                 break;
              case 2:
                 // use right-shifted stencil since it is smoother
                 for (ptrdiff_t i=rcoeffs::imin; i<rcoeffs::imax; ++i) {
                    res += rcoeffs::get(i) * f[i-lcoeffs::minimin];
                 }
                 break;
           }
         }
       }
       break;
       default: assert(0);
    }
    
    // check that result is reasonable!
    RT res1 = 2*f[-coeffs1d<RT,ORDER>::minimin] - res;
    if ((res - f[-1-coeffs1d<RT,ORDER>::minimin]) * (f[-coeffs1d<RT,ORDER>::minimin] - res) < 0 or
        (f[1-coeffs1d<RT,ORDER>::minimin] - res1) * (res1 - f[-coeffs1d<RT,ORDER>::minimin]) < 0)
    {
      // switch back to zeroeth order!
      // TODO: come up with higher order fallback?
      res = f[-coeffs1d<RT,ORDER>::minimin];
      res1 = 2*f[-coeffs1d<RT,ORDER>::minimin] - res;
    }
    
    if (dk) // this trivially enforces the conservative property
       res = res1;

    return res;
  }
  
  // Check interpolation index ranges
  template <typename T, int ORDER>
  static inline
  void
  check_indices0 ()
  {
  }
  
  template <typename T, int ORDER>
  static inline
  void
  check_indices1 (ptrdiff_t const is,
                  ptrdiff_t const srciext)
  {
#ifdef CARPET_DEBUG
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    assert (is + coeffs::minimin >= 0);
    assert (is + coeffs::maximax <= srciext);
    check_indices0<T,ORDER> ();
#endif
  }
  
  template <typename T, int ORDER>
  static inline
  void
  check_indices2 (ptrdiff_t const is,
                  ptrdiff_t const js,
                  ptrdiff_t const srciext,
                  ptrdiff_t const srcjext)
  {
#ifdef CARPET_DEBUG
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    assert (js + coeffs::minimin >= 0);
    assert (js + coeffs::maximax <= srcjext);
    check_indices1<T,ORDER> (is, srciext);
#endif
  }
  
  template <typename T, int ORDER>
  static inline
  void
  check_indices3 (ptrdiff_t const is,
                  ptrdiff_t const js,
                  ptrdiff_t const ks,
                  ptrdiff_t const srciext,
                  ptrdiff_t const srcjext,
                  ptrdiff_t const srckext)
  {
#ifdef CARPET_DEBUG
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs;
    assert (ks + coeffs::minimin >= 0);
    assert (ks + coeffs::maximax <= srckext);
    check_indices2<T,ORDER> (is,js, srciext,srcjext);
#endif
  }
  
  
  
  template <typename T, int ORDER, int di>
  static void
  testinterp ()
  {
    static bool tested = false;
    if (tested) return;
    tested = true;
    
    typedef typename typeprops<T>::real RT;
    typedef coeffs1d<RT,ORDER> coeffs_t;

    // get function values needed for the stencil
    T f[coeffs_t::maxncoeffs];
    for (ptrdiff_t i=coeffs_t::minimin; i<coeffs_t::maximax; ++i) {
        RT const x0 = RT(CCTK_REAL(i));
        RT const x1 = RT(CCTK_REAL(i)+1.0);
        //cout << "x=" << 0.5*(x1+x0) << endl;
        RT const y = (ipow (x1, ORDER) - ipow (x0, ORDER)) / ((x1-x0)*(ORDER));
        //cout << "y=" << y << endl;
       f[i-coeffs_t::minimin] = y;
    }
    T const res = interp3<T,ORDER,0,0,di> (f - coeffs_t::minimin, 0,0,1);
    
    RT const x0 = di * RT(0.5);
    RT const x1 = di * RT(0.5) + RT(0.5);
    //cout << "x0=" << x0 << endl;
    RT const y0 = (ipow (x1, ORDER) - ipow (x0, ORDER)) / ((x1-x0)*(ORDER));
    if (not (fabs (res - y0) < 1.0e-12)) {
      RT rt;
      ostringstream buf;
      buf << "Error in prolongate_3d_cc_rf2::coeffs_3d_cc_rf2\n"
          << "   RT=" << typestring(rt) << "\n"
          << "   ORDER=" << ORDER << "\n"
          << "   di=" << di << "\n"
          << "   y0=" << y0 << ", res=" << res;
      CCTK_WARN (CCTK_WARN_ALERT, buf.str().c_str());
      CCTK_WARN (CCTK_WARN_ABORT, "Aborting.");
    }
  }
  
  
  
  template <typename T, int ORDER>
  void
  prolongate_3d_cc_enovol_rf2 (T const * restrict const src,
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
    
    static_assert (ORDER>=0, "ORDER must be non-negative");
    
    typedef typename typeprops<T>::real RT;
    coeffs1d<RT,ORDER,0>::test ();
    coeffs1d<RT,ORDER,-1>::test ();
    coeffs1d<RT,ORDER,-2>::test ();

    //testinterp<RT,ORDER,0> ();
    //testinterp<RT,ORDER,1> ();
    
    
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
    
    
    
    // Determine the stencil radius (see diagram at the top of this
    // file)
    ivect3 offsetlo, offsethi;
    {
      assert (all ((regbbox.upper() - srcbbox.lower() + regbbox.stride() / 2) % regbbox.stride() == 0));
      for (int d=0; d<3; ++d) {
        offsetlo[d] = - coeffs1d<RT,ORDER>::minimin;
        offsethi[d] = coeffs1d<RT,ORDER>::maximax - 1;
      }
    }
    
    
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
    
    
    
    size_t const srcipadext = srcpadext[0];
    size_t const srcjpadext = srcpadext[1];
    size_t const srckpadext = srcpadext[2];
    
    size_t const dstipadext = dstpadext[0];
    size_t const dstjpadext = dstpadext[1];
    size_t const dstkpadext = dstpadext[2];
    
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
    
    // size_t const srcdi = SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
    size_t const srcdi = 1;
    assert (srcdi == SRCOFF3(1,0,0) - SRCOFF3(0,0,0));
    size_t const srcdj = SRCOFF3(0,1,0) - SRCOFF3(0,0,0);
    size_t const srcdk = SRCOFF3(0,0,1) - SRCOFF3(0,0,0);
    
    
    
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
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,0,0> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8001;
    goto l900;
    
    // kernel
  l8001:
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
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
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,1,0> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8011;
    goto l901;
    
    // kernel
  l8011:
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
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
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,0,1> (& src[SRCIND3(is,js,ks)], srcdi, srcdj, srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8101;
    goto l910;
    
    // kernel
  l8101:
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
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
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
    dst[DSTIND3(id,jd,kd)] =
      interp3<T,ORDER,0,1,1> (& src[SRCIND3(is,js,ks)], srcdi,srcdj,srcdk);
    i = i+1;
    id = id+1;
    if (i < regiext) goto l8111;
    goto l911;
    
    // kernel
  l8111:
    check_indices3<T,ORDER> (is,js,ks, srciext, srcjext, srckext);
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
  

  // Specialise for complex types
  
#define TYPECASE(N,T)                                                   \
                                                                        \
  template <>                                                           \
  void                                                                  \
  prolongate_3d_cc_enovol_rf2<T,2> (T const * restrict const src,       \
                                    ivect3 const & restrict srcpadext,  \
                                    ivect3 const & restrict srcext,     \
                                    T * restrict const dst,             \
                                    ivect3 const & restrict dstpadext,  \
                                    ivect3 const & restrict dstext,     \
                                    ibbox3 const & restrict srcbbox,    \
                                    ibbox3 const & restrict dstbbox,    \
                                    ibbox3 const & restrict,            \
                                    ibbox3 const & restrict regbbox,    \
                                    void * extraargs)                   \
  {                                                                     \
    CCTK_WARN (CCTK_WARN_ABORT,                                         \
               "ENO operators are not supported for CCTK_COMPLEX");     \
  }

#define CARPET_COMPLEX
#include "typecase.hh"
#undef TYPECASE                                                                 


  
#define TYPECASE(N,T)                                                   \
                                                                        \
  template                                                              \
  void                                                                  \
  prolongate_3d_cc_enovol_rf2<T,2> (T const * restrict const src,       \
                                    ivect3 const & restrict srcpadext,  \
                                    ivect3 const & restrict srcext,     \
                                    T * restrict const dst,             \
                                    ivect3 const & restrict dstpadext,  \
                                    ivect3 const & restrict dstext,     \
                                    ibbox3 const & restrict srcbbox,    \
                                    ibbox3 const & restrict dstbbox,    \
                                    ibbox3 const & restrict,            \
                                    ibbox3 const & restrict regbbox,    \
                                    void * extraargs);

#define CARPET_NO_INT
#define CARPET_NO_COMPLEX
#include "typecase.hh"
#undef TYPECASE

  
  
} // namespace CarpetLib
