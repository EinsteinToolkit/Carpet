// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.hh,v 1.2 2003/11/05 16:18:39 schnetter Exp $

#ifndef CARPETSLAB_HH
#define CARPETSLAB_HH

#include "cctk.h"

#include "slab.h"

namespace CarpetSlab {
  
  // Non-standard interface -- don't use
  void FillSlab (const cGH* const cgh,
		 const int dest_proc,
		 const int n,
		 const int tl,
		 const int hdim,
		 const int origin[/*vdim*/],
		 const int dirs[/*hdim*/],
		 const int stride[/*hdim*/],
		 const int length[/*hdim*/],
                 void* const hdata);
  
  // Non-standard interface -- don't use
  void* GetSlab (const cGH* const cgh,
		 const int dest_proc,
		 const int n,
		 const int tl,
		 const int hdim,
		 const int origin[/*vdim*/],
		 const int dirs[/*hdim*/],
		 const int stride[/*hdim*/],
		 const int length[/*hdim*/]);
  
} // namespace CarpetSlab

#endif // !defined(CARPETSLAB_HH)
