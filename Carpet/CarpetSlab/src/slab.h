/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.h,v 1.1 2002/10/24 10:53:48 schnetter Exp $ */

#ifndef CARPETSLAB_H
#define CARPETSLAB_H

#include "cctk.h"

#ifdef __cplusplus
namespace CarpetSlab {
  extern "C" {
#endif
    
    int Hyperslab_GetHyperslab (cGH* const GH,
				const int target_proc,
				const int vindex,
				const int vtimelvl,
				const int hdim,
				const int global_startpoint [/*vdim*/],
				const int directions [/*vdim*/],
				const int lengths [/*hdim*/],
				const int downsample [/*hdim*/],
				void** const hdata,
				int hsize [/*hdim*/]);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetSlab */
#endif

#endif /* !defined(CARPETSLAB_H) */
