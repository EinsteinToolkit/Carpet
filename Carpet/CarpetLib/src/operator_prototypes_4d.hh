#ifndef OPERATOR_PROTOTYPES_4D
#define OPERATOR_PROTOTYPES_4D

#include <cstdlib>

#include <cctk.h>

#include "defs.hh"
#include "bbox.hh"
#include "vect.hh"

#include "operator_prototypes.hh"



namespace CarpetLib {
  
  using namespace std;
  
  
  
  static inline
  size_t
  index4 (size_t const i, size_t const j, size_t const k, size_t const kl,
          size_t const exti, size_t const extj, size_t const extk, size_t const extl)
    CCTK_ATTRIBUTE_CONST;
  static inline
  size_t
  index4 (size_t const i, size_t const j, size_t const k, size_t const l,
          size_t const exti, size_t const extj, size_t const extk, size_t const extl)
  {
#ifdef CARPET_DEBUG
    assert (static_cast <ptrdiff_t> (i) >= 0 and i < exti);
    assert (static_cast <ptrdiff_t> (j) >= 0 and j < extj);
    assert (static_cast <ptrdiff_t> (k) >= 0 and k < extk);
    assert (static_cast <ptrdiff_t> (l) >= 0 and l < extl);
#endif
    
    return i + exti * (j + extj * (k + extk * l));
  }
  
  
  
  static int const dim4 = 4;
  
  typedef vect <bool, dim4> bvect4;
  typedef vect <int, dim4> ivect4;
  typedef bbox <int, dim4> ibbox4;
  
  
  
  template <typename T>
  void
  copy_4d (T const * restrict const src,
           ivect4 const & restrict srcext,
           T * restrict const dst,
           ivect4 const & restrict dstext,
           ibbox4 const & restrict srcbbox,
           ibbox4 const & restrict dstbbox,
           ibbox4 const & restrict regbbox);
  
  
  
  template <typename T>
  void
  prolongate_4d_o1_rf2 (T const * restrict const src,
                        ivect4 const & restrict srcext,
                        T * restrict const dst,
                        ivect4 const & restrict dstext,
                        ibbox4 const & restrict srcbbox,
                        ibbox4 const & restrict dstbbox,
                        ibbox4 const & restrict regbbox);
  
  
  
  template <typename T>
  void
  restrict_4d_rf2 (T const * restrict const src,
                   ivect4 const & restrict srcext,
                   T * restrict const dst,
                   ivect4 const & restrict dstext,
                   ibbox4 const & restrict srcbbox,
                   ibbox4 const & restrict dstbbox,
                   ibbox4 const & restrict regbbox);
  
  
  
} // namespace CarpetLib



#endif // #ifndef OPERATOR_PROTOTYPES_4D