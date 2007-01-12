#ifndef TYPEPROPS_HH
#define TYPEPROPS_HH

#include <cctk.h>



template <typename T>
struct typeprops {
  typedef T complex;
  typedef T real;
  static inline complex fromreal (real const x) { return x; }
};

#ifdef HAVE_CCTK_COMPLEX8
template <>
struct typeprops <CCTK_COMPLEX8> {
  typedef CCTK_COMPLEX8 complex;
  typedef CCTK_REAL4 real;
  static inline complex fromreal (real const x) { return CCTK_Cmplx8 (x, 0); }
};
#endif

#ifdef HAVE_CCTK_COMPLEX16
template <>
struct typeprops <CCTK_COMPLEX16> {
  typedef CCTK_COMPLEX16 complex;
  typedef CCTK_REAL8 real;
  static inline complex fromreal (real const x) { return CCTK_Cmplx16 (x, 0); }
};
#endif

#ifdef HAVE_CCTK_COMPLEX32
template <>
struct typeprops <CCTK_COMPLEX32> {
  typedef CCTK_COMPLEX32 complex;
  typedef CCTK_REAL16 real;
  static inline complex fromreal (real const x) { return CCTK_Cmplx32 (x, 0); }
};
#endif



#endif // #ifndef TYPEPROPS_HH
