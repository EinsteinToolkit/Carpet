// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/Attic/carpetslab.hh,v 1.3 2001/03/10 20:55:09 eschnett Exp $

#ifndef CARPETSLAB_HH
#define CARPETSLAB_HH

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
