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
  
  
  
  template <typename T>
  void
  restrict_3d_cc_rf2 (T const * restrict const src,
                      ivect3 const & restrict srcext,
                      T * restrict const dst,
                      ivect3 const & restrict dstext,
                      ibbox3 const & restrict srcbbox,
                      ibbox3 const & restrict dstbbox,
                      ibbox3 const & restrict regbbox)
  {
    DECLARE_CCTK_PARAMETERS;
    
    typedef typename typeprops<T>::real RT;
    
    
    
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
    
    if (any (srcext != srcbbox.shape() / srcbbox.stride() or
             dstext != dstbbox.shape() / dstbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all (srcbbox.stride() % 2 == 0));
    assert (all ((regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) %
                 srcbbox.stride() == 0));
    ivect3 const srcoff =
      (regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) /
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
    
    
    
    RT const one = 1;
    
    RT const f1 = one/2;
    RT const f2 = one/2;
    
    
    
    // Loop over coarse region
#pragma omp parallel for
    for (int k=0; k<regkext; ++k) {
      for (int j=0; j<regjext; ++j) {
        for (int i=0; i<regiext; ++i) {
          
          dst [DSTIND3(i, j, k)] =
            + f1*f1*f1 * src [SRCIND3(2*i  , 2*j  , 2*k  )]
            + f2*f1*f1 * src [SRCIND3(2*i+1, 2*j  , 2*k  )]
            + f1*f2*f1 * src [SRCIND3(2*i  , 2*j+1, 2*k  )]
            + f2*f2*f1 * src [SRCIND3(2*i+1, 2*j+1, 2*k  )]
            + f1*f1*f2 * src [SRCIND3(2*i  , 2*j  , 2*k+1)]
            + f2*f1*f2 * src [SRCIND3(2*i+1, 2*j  , 2*k+1)]
            + f1*f2*f2 * src [SRCIND3(2*i  , 2*j+1, 2*k+1)]
            + f2*f2*f2 * src [SRCIND3(2*i+1, 2*j+1, 2*k+1)];
          
        }
      }
    }
    
  }
  
  
  
#define INSTANTIATE(T)                                  \
  template                                              \
  void                                                  \
  restrict_3d_cc_rf2 (T const * restrict const src,     \
                      ivect3 const & restrict srcext,   \
                      T * restrict const dst,           \
                      ivect3 const & restrict dstext,   \
                      ibbox3 const & restrict srcbbox,  \
                      ibbox3 const & restrict dstbbox,  \
                      ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib
