// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/Attic/carpetslab.hh,v 1.4 2002/09/01 14:52:30 schnetter Exp $

#ifndef CARPETSLAB_HH
#define CARPETSLAB_HH

#include "cctk.h"

#include "carpetslab.h"

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
