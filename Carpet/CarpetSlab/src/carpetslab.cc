// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/Attic/carpetslab.cc,v 1.10 2001/07/04 12:29:52 schnetter Exp $

#include <alloca.h>
#include <assert.h>
#include <stdlib.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/bboxset.hh"
#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/gdata.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "carpetslab.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/Attic/carpetslab.cc,v 1.10 2001/07/04 12:29:52 schnetter Exp $";



namespace CarpetSlab {
  
  using namespace Carpet;
  
  
  
  void* GetSlab (cGH* const cgh,
		 const int dest_proc,
		 const int n,
		 const int ti,
		 const int hdim,
		 const int origin[/*vdim*/],
		 const int dirs[/*hdim*/],
		 const int stride[/*hdim*/],
		 const int length[/*hdim*/])
  {
    // Check global state
    assert (reflevel>=0);
    assert (mglevel>=0);
    
    // Save global state
    int saved_component = component;
    component = -1;
    
    // Check Cactus grid hierarchy
    assert (cgh);
    
    // Check destination processor
    assert (dest_proc>=-1 && dest_proc<CCTK_nProcs(cgh));
    
    // Check variable index
    assert (n>=0 && n<CCTK_NumVars());
    
    // Get info about variable
    const int group = CCTK_GroupIndexFromVarI(n);
    assert (group>=0);
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int var = n - n0;
    assert (var>=0);
    
    // Get info about group
    cGroup gp;
    CCTK_GroupData (group, &gp);
    assert (gp.dim<=dim);
    assert (CCTK_QueryGroupStorageI(cgh, group));
    const int typesize = CCTK_VarTypeSize(gp.vartype);
    assert (typesize>0);
    
    // Check dimension
    assert (hdim>=0 && hdim<=gp.dim);
    
    // Check timelevel
    const int num_tl = gp.numtimelevels;
    assert (ti>=0 && ti<num_tl);
    const int tl = -ti;
    
    // Check origin
//     for (int d=0; d<dim; ++d) {
//       assert (origin[d]>=0 && origin[d]<=sizes[d]);
//     }
    
    // Check directions
    for (int dd=0; dd<hdim; ++dd) {
      assert (dirs[dd]>=1 && dirs[dd]<=dim);
    }
    
    // Check stride
    for (int dd=0; dd<hdim; ++dd) {
      assert (stride[dd]>0);
    }
    
    // Check length
    for (int dd=0; dd<hdim; ++dd) {
      assert (length[dd]>=0);
    }
    
    // Check extent
//     for (int dd=0; dd<hdim; ++dd) {
//       assert (origin[dirs[dd]-1] + length[dd] <= sizes[dirs[dd]]);
//     }
    
    // Get insider information about variable
    const gh<dim>* myhh;
    const dh<dim>* mydd;
    const generic_gf<dim>* myff;
    switch (gp.grouptype) {
    case CCTK_SCALAR:
      abort();
    case CCTK_ARRAY:
      assert (group < (int)arrdata.size());
      myhh = arrdata[group].hh;
      mydd = arrdata[group].dd;
      assert (var < (int)arrdata[group].data.size());
      myff = arrdata[group].data[var];
      break;
    case CCTK_GF:
      myhh = hh;
      mydd = Carpet::dd;	// dd is int dd above for old C++ compilers
      assert (group < (int)gfdata.size());
      assert (var < (int)gfdata[group].data.size());
      myff = gfdata[group].data[var];
      break;
    default:
      abort();
    }
    assert (myhh);
    assert (mydd);
    assert (myff);
    
    // Detemine collecting processor
    const int collect_proc = dest_proc<0 ? 0 : dest_proc;
    
    // Determine own rank
    const int rank = CCTK_MyProc(cgh);
    
    // Calculate global size
    int totalsize = 1;
    for (int dd=0; dd<hdim; ++dd) {
      totalsize *= length[dd];
    }
    
    // Allocate memory
    void* hdata = 0;
    if (dest_proc==-1 || rank==dest_proc) {
      hdata = malloc(totalsize * typesize);
      assert (hdata);
    }
    
    if (hh->components(reflevel) > 0) {
      
      // Only temporarily
      component = 0;
      
      // Get sample data
      const generic_data<dim>* mydata;
      mydata = (*myff)(tl, reflevel, component, mglevel);
      
      // Stride of data in memory
      const vect<int,dim> str = mydata->extent().stride();
      
      // Stride of collected data
      vect<int,dim> hstr = str;
      for (int dd=0; dd<hdim; ++dd) {
	hstr[dirs[dd]-1] *= stride[dd];
      }
      
      // Lower bound of collected data
      vect<int,dim> hlb(0);
      for (int d=0; d<gp.dim; ++d) {
	hlb[d] = origin[d] * str[d];
      }
      
      // Upper bound of collected data
      vect<int,dim> hub = hlb;
      for (int dd=0; dd<hdim; ++dd) {
	hub[dirs[dd]-1] += (length[dd]-1) * hstr[dirs[dd]-1];
      }
      
      // Calculate extent to collect
      const bbox<int,dim> hextent (hlb, hub, hstr);
      assert (hextent.num_points() == totalsize);
      
      // Create collector data object
      void* myhdata = rank==collect_proc ? hdata : 0;
      generic_data<dim>* const alldata = mydata->make_typed();
      alldata->allocate (hextent, collect_proc, myhdata);
      
      // Done with the temporary stuff
      mydata = 0;
      component = -1;
      
      // Loop over all components, copying data from them
      assert (component == -1);
      for (component=0; component<hh->components(reflevel); ++component) {
	
	// Get data object
	mydata = (*myff)(tl, reflevel, component, mglevel);
	
	// Calculate overlapping extents
	const bboxset<int,dim> myextents
	  = ((mydd->boxes[reflevel][component][mglevel].sync_not
	      | mydd->boxes[reflevel][component][mglevel].interior)
	     & hextent);
	
	// Loop over overlapping extents
	for (bboxset<int,dim>::const_iterator ext_iter = myextents.begin();
	     ext_iter != myextents.end();
	     ++ext_iter) {
	  
	  // Copy data
	  alldata->copy_from (mydata, *ext_iter);
	  
	}
	
      } // Loop over components
      component = -1;
      
      // Copy result to all processors
      if (dest_proc == -1) {
	for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
	  if (proc != collect_proc) {
	    
	    void* myhdata = rank==proc ? hdata : 0;
	    generic_data<dim>* const tmpdata = mydata->make_typed();
	    tmpdata->allocate (alldata->extent(), proc, myhdata);
	    tmpdata->copy_from (alldata, alldata->extent());
	    delete tmpdata;
	    
	  }
	}
      } // Copy result
      
      delete alldata;
      
    } // if components>0
    
    // Restore global state
    component = saved_component;
    
    // Success
    return hdata;
  }
  
  
  
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
			      int hsize [/*hdim*/])
  {
    const int gpdim = CCTK_GroupDimFromVarI(vindex);
    assert (gpdim>=1 && gpdim<=dim);
    
    // Check some arguments
    assert (hdim>=0 && hdim<=dim);
    
    // Check output arguments
    assert (hdata);
    assert (hsize);
    
    // Calculate more convenient representation of the direction
    int* const dirs = (int*)alloca(hdim * sizeof(int));
    // The following if statement is written according to the
    // definition of "dir".
    if (hdim==1) {
      // 1-dimensional hyperslab
      int mydir = 0;
      for (int d=0; d<dim; ++d) {
	if (directions[d]!=0) {
	  mydir = d+1;
	  break;
	}
      }
      assert (mydir>0);
      for (int d=0; d<dim; ++d) {
	if (d == mydir-1) {
	  assert (directions[d]!=0);
	} else {
	  assert (directions[d]==0);
	}
      }
      dirs[0] = mydir;
    } else if (hdim==dim) {
      // dim-dimensional hyperslab
      for (int dd=0; dd<hdim; ++dd) {
	dirs[dd] = dd+1;
      }
    } else if (hdim==2) {
      // 2-dimensional hyperslab with dim==3
      assert (dim==3);
      int mydir = 0;
      for (int d=0; d<dim; ++d) {
	if (directions[d]==0) {
	  mydir = d+1;
	  break;
	}
      }
      assert (mydir>0);
      for (int d=0; d<dim; ++d) {
	if (d == mydir-1) {
	  assert (directions[d]==0);
	} else {
	  assert (directions[d]!=0);
	}
      }
      int dd=0;
      for (int d=0; d<dim; ++d) {
	if (d != mydir-1) {
	  dirs[dd] = d+1;
	  ++dd;
	}
      }
      assert (dd==hdim);
    } else {
      abort();
    }
    
#if 0
    // Invert directions
    int invdir[dim];
    for (int d=0; d<dim; ++d) {
      invdir[d] = 0;
    }
    for (int dd=0; dd<hdim; ++dd) {
      invdir[dirs[dd]-1] = dd+1;
    }
#endif
    
    // Calculate lengths
    for (int dd=0; dd<hdim; ++dd) {
      if (lengths[dd]<0) {
	const int totlen = *CCTK_ArrayGroupSizeI(GH, dirs[dd]-1, vindex);
	assert (totlen>=0);
	// Partial argument check
	assert (global_startpoint[dirs[dd]-1]>=0);
	assert (global_startpoint[dirs[dd]-1]<=totlen);
	assert (downsample[dd]>0);
	hsize[dd] = (totlen - global_startpoint[dirs[dd]-1]) / downsample[dd];
      } else {
	hsize[dd] = lengths[dd];
      }
      assert (hsize[dd]>=0);
    }
    
    // Get the slab
    *hdata = GetSlab (GH,
		      target_proc,
		      vindex,
		      vtimelvl,
		      hdim,
		      global_startpoint,
		      dirs,
		      downsample,
		      hsize);
    
    // Return with success
    return 1;
  }
  
  
  
} // namespace CarpetSlab
