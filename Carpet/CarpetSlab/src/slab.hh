// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.hh,v 1.1 2002/10/24 10:53:48 schnetter Exp $

#ifndef CARPETSLAB_HH
#define CARPETSLAB_HH

#include "cctk.h"

#include "slab.h"

namespace CarpetSlab {
  
  void* GetSlab (cGH* const cgh,
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
