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
  void
  restrict_3d_rf2 (T const * restrict const src,
                   ivect3 const & restrict srcext,
                   T * restrict const dst,
                   ivect3 const & restrict dstext,
                   ibbox3 const & restrict srcbbox,
                   ibbox3 const & restrict dstbbox,
                   ibbox3 const & restrict regbbox)
  {
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
    
    
    
    // Loop over coarse region
    for (size_t k=0; k<regkext; ++k) {
      for (size_t j=0; j<regjext; ++j) {
        for (size_t i=0; i<regiext; ++i) {
          
          dst [DSTIND3(i, j, k)] = src [SRCIND3(2*i, 2*j, 2*k)];
          
        }
      }
    }
    
  }
  
  
  
  template
  void
  restrict_3d_rf2 (CCTK_INT const * restrict const src,
                   ivect3 const & restrict srcext,
                   CCTK_INT * restrict const dst,
                   ivect3 const & restrict dstext,
                   ibbox3 const & restrict srcbbox,
                   ibbox3 const & restrict dstbbox,
                   ibbox3 const & restrict regbbox);
  
  template
  void
  restrict_3d_rf2 (CCTK_REAL const * restrict const src,
                   ivect3 const & restrict srcext,
                   CCTK_REAL * restrict const dst,
                   ivect3 const & restrict dstext,
                   ibbox3 const & restrict srcbbox,
                   ibbox3 const & restrict dstbbox,
                   ibbox3 const & restrict regbbox);
  
  template
  void
  restrict_3d_rf2 (CCTK_COMPLEX const * restrict const src,
                   ivect3 const & restrict srcext,
                   CCTK_COMPLEX * restrict const dst,
                   ivect3 const & restrict dstext,
                   ibbox3 const & restrict srcbbox,
                   ibbox3 const & restrict dstbbox,
                   ibbox3 const & restrict regbbox);
  
  
  
} // namespace CarpetLib
