#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>

#include "gdata.hh"
#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;



namespace CarpetLib {


  
#define SRCIND3(i,j,k)                                  \
  index3 (srcioff + (i), srcjoff + (j), srckoff + (k),  \
          srciext, srcjext, srckext)
#define DSTIND3(i,j,k)                                  \
  index3 (dstioff + (i), dstjoff + (j), dstkoff + (k),  \
          dstiext, dstjext, dstkext)
  
  
  
  // 0D "restriction"
  template <typename T>
  struct restrict0 {
    static inline T call (T const * restrict const p)
    {
      return * p;
    }
  };
  
  // 1D restriction
  template <typename T, int centi>
  struct restrict1 {
    static inline T call (T const * restrict const p,
                          size_t const d1);
  };
  template <typename T>
  struct restrict1<T,0> {
    static inline T call (T const * restrict const p,
                          size_t const d1)
    {
      T const res = restrict0<T>::call (p);
      return res;
    }
  };
  template <typename T>
  struct restrict1<T,1> {
    static inline T call (T const * restrict const p,
                          size_t const d1)
    {
      typedef typename typeprops<T>::real RT;
      RT const one  = 1;
      RT const half = one/2;
      T const res =
        half * restrict0<T>::call (p     ) +
        half * restrict0<T>::call (p + d1);
      return res;
    }
  };
  
  // 2D restriction
  template <typename T, int centi, int centj>
  struct restrict2 {
    static inline T call (T const * restrict const p,
                          size_t const d1, size_t const d2);
  };
  template <typename T, int centi>
  struct restrict2<T,centi,0> {
    static inline T call (T const * restrict const p,
                          size_t const d1, size_t const d2)
    {
      T const res = restrict1<T,centi>::call (p, d1);
      return res;
    }
  };
  template <typename T, int centi>
  struct restrict2<T,centi,1> {
    static inline T call (T const * restrict const p,
                          size_t const d1, size_t const d2)
    {
      typedef typename typeprops<T>::real RT;
      RT const one  = 1;
      RT const half = one/2;
      T const res =
        half * restrict1<T,centi>::call (p     , d1) +
        half * restrict1<T,centi>::call (p + d2, d1);
      return res;
    }
  };
  
  // 3D restriction
  template <typename T, int centi, int centj, int centk>
  struct restrict3 {
    static inline T call (T const * restrict const p,
                          size_t const d1, size_t const d2, size_t const d3);
  };
  template <typename T, int centi, int centj>
  struct restrict3<T,centi,centj,0> {
    static inline T call (T const * restrict const p,
                          size_t const d1, size_t const d2, size_t const d3)
    {
      T const res = restrict2<T,centi,centj>::call (p, d1, d2);
      return res;
    }
  };
  template <typename T, int centi, int centj>
  struct restrict3<T,centi,centj,1> {
    static inline T call (T const * restrict const p,
                          size_t const d1, size_t const d2, size_t const d3)
    {
      typedef typename typeprops<T>::real RT;
      RT const one  = 1;
      RT const half = one/2;
      T const res =
        half * restrict2<T,centi,centj>::call (p     , d1, d2) +
        half * restrict2<T,centi,centj>::call (p + d3, d1, d2);
      return res;
    }
  };
  
  
  
  template <typename T, int centi, int centj, int centk>
  void
  restrict_3d_vc_rf2 (T const * restrict const src,
                      ivect3 const & restrict srcext,
                      T * restrict const dst,
                      ivect3 const & restrict dstext,
                      ibbox3 const & restrict srcbbox,
                      ibbox3 const & restrict dstbbox,
                      ibbox3 const & restrict regbbox)
  {
    DECLARE_CCTK_PARAMETERS;
    
    typedef typename typeprops<T>::real RT;
    
    // false: vertex centered, true: cell centered
    ivect const icent (centi, centj, centk);
    assert (all (icent == 0 or icent == 1));
    bvect const cent (icent);
    
    
    
    if (any (srcbbox.stride() >= regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    if (any (reffact2 * srcbbox.stride() != dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: destination strides are not twice the source strides");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
    if (not regbbox.expanded_for(srcbbox).is_contained_in(srcbbox) or
        not regbbox.is_contained_in(dstbbox))
    {
      CCTK_WARN (0, "Internal error: region extent is not contained in array extent");
    }
    
    if (any (srcext != gdata::allocated_memory_shape(srcbbox.shape() / srcbbox.stride()) or
             dstext != gdata::allocated_memory_shape(dstbbox.shape() / dstbbox.stride())))
    {
      CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all (srcbbox.stride() % either (cent, 2, 1) == 0));
    if (not (all ((regbbox.lower() - srcbbox.lower() -
                   either (cent, srcbbox.stride() / 2, 0)) %
                  srcbbox.stride() == 0)))
    {
      cout << "restrict_3d_vc_rf2.cc\n";
      cout << "regbbox=" << regbbox << "\n";
      cout << "srcbbox=" << srcbbox << "\n";
      cout << "cent=" << cent << "\n";
    }
    assert (all ((regbbox.lower() - srcbbox.lower() -
                  either (cent, srcbbox.stride() / 2, 0)) %
                 srcbbox.stride() == 0));
    ivect3 const srcoff =
      (regbbox.lower() - srcbbox.lower() -
       either (cent, srcbbox.stride() / 2, 0)) /
      srcbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
    ivect3 const dstoff =
      (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();
    
    
    
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
    
    size_t const srcdi = SRCIND3(1,0,0) - SRCIND3(0,0,0);
    size_t const srcdj = SRCIND3(0,1,0) - SRCIND3(0,0,0);
    size_t const srcdk = SRCIND3(0,0,1) - SRCIND3(0,0,0);
    
    
    
    // Loop over coarse region
    for (int k=0; k<regkext; ++k) {
      for (int j=0; j<regjext; ++j) {
        for (int i=0; i<regiext; ++i) {
          
          dst [DSTIND3(i, j, k)] =
            restrict3<T,centi,centj,centk>::call
            (& src[SRCIND3(2*i, 2*j, 2*k)], srcdi, srcdj, srcdk);
          
        }
      }
    }
    
  }
  
  
  
#define TYPECASE(N,T)                                                   \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,0,0,0> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,0,0,1> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,0,1,0> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,0,1,1> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,1,0,0> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,1,0,1> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,1,1,0> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);        \
  template                                                              \
  void                                                                  \
  restrict_3d_vc_rf2<T,1,1,1> (T const * restrict const src,            \
                               ivect3 const & restrict srcext,          \
                               T * restrict const dst,                  \
                               ivect3 const & restrict dstext,          \
                               ibbox3 const & restrict srcbbox,         \
                               ibbox3 const & restrict dstbbox,         \
                               ibbox3 const & restrict regbbox);
#include "typecase.hh"
#undef TYPECASE
  
  
  
} // namespace CarpetLib
