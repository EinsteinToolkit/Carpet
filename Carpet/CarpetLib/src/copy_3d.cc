#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

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
  copy_3d (T const * restrict const src,
           ivect3 const & restrict srcext,
           T * restrict const dst,
           ivect3 const & restrict dstext,
           ibbox3 const & restrict srcbbox,
           ibbox3 const & restrict dstbbox,
           ibbox3 const & restrict regbbox)
  {
    if (any (srcbbox.stride() != regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      cout << "copy_3d.cc:" << endl
           << "srcbbox=" << srcbbox << endl
           << "dstbbox=" << dstbbox << endl
           << "regbbox=" << regbbox << endl;
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
    assert (all ((regbbox.lower() - srcbbox.lower()) % srcbbox.stride() == 0));
    ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / srcbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
    ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();
    
    
    
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
    
    
    
    // Loop over region
#pragma omp parallel for
    for (int k=0; k<regkext; ++k) {
      for (int j=0; j<regjext; ++j) {
        for (int i=0; i<regiext; ++i) {
          
          dst [DSTIND3(i, j, k)] = src [SRCIND3(i, j, k)];
          
        }
      }
    }
    
  }
  
  
  
#define INSTANTIATE(T)                          \
  template                                      \
  void                                          \
  copy_3d (T const * restrict const src,        \
           ivect3 const & restrict srcext,      \
           T * restrict const dst,              \
           ivect3 const & restrict dstext,      \
           ibbox3 const & restrict srcbbox,     \
           ibbox3 const & restrict dstbbox,     \
           ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib
