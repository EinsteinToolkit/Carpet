// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/Attic/carpetslab.h,v 1.1 2001/03/07 13:01:11 eschnett Exp $

#include "cctk.h"

int Hyperslab_GetLocalHyperslab (cGH* GH,
				 int vindex,
				 int vtimelvl,
				 int hdim,
				 const int global_startpoint [/*vdim*/],
				 const int directions [/*vdim*/],
				 const int lengths [/*hdim*/],
				 const int downsample [/*hdim*/],
				 void** hdata,
				 int hsize [/*hdim*/],
				 int ghsize [/*hdim*/],
				 int hoffset [/*hdim*/]);

int Hyperslab_GetHyperslab (cGH* GH,
			    int target_proc,
			    int vindex,
			    int vtimelvl,
			    int hdim,
			    const int global_startpoint [/*vdim*/],
			    const int directions [/*vdim*/],
			    const int lengths [/*hdim*/],
			    const int downsample [/*hdim*/],
			    void** hdata,
			    int hsize [/*hdim*/]);
