// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.cc,v 1.10 2003/06/18 18:28:08 schnetter Exp $

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "cctk.h"

#include "bbox.hh"
#include "bboxset.hh"
#include "dh.hh"
#include "gdata.hh"
#include "gh.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "slab.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.cc,v 1.10 2003/06/18 18:28:08 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetSlab_slab_cc);
}



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
    
    if (gp.grouptype==CCTK_GF && reflevel==-1) {
      CCTK_WARN (0, "It is not possible to use hyperslabbing for a grid function in global mode");
    }
    const int rl = gp.grouptype==CCTK_GF ? reflevel : 0;
    
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
    const ggf<dim>* myff;
    assert (group < (int)arrdata.size());
    myhh = arrdata[group].hh;
    assert (myhh);
    mydd = arrdata[group].dd;
    assert (mydd);
    assert (var < (int)arrdata[group].data.size());
    myff = arrdata[group].data[var];
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
      memset (hdata, 0, totalsize * typesize);
    }
    
    // Get sample data
    const gdata<dim>* mydata;
    mydata = (*myff)(tl, rl, 0, 0);
    
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
    gdata<dim>* const alldata = mydata->make_typed();
    alldata->allocate (hextent, collect_proc, myhdata);
    
    // Done with the temporary stuff
    mydata = 0;
    
    // Loop over all components, copying data from them
    BEGIN_LOCAL_COMPONENT_LOOP (cgh, gp.grouptype) {
      
      const int c = gp.grouptype==CCTK_GF ? component : 0;
      
      // Get data object
      mydata = (*myff)(tl, rl, c, mglevel);
      
      // Calculate overlapping extents
      const bboxset<int,dim> myextents
        = ((mydd->boxes[rl][c][mglevel].sync_not
            | mydd->boxes[rl][c][mglevel].interior)
           & hextent);
      
      // Loop over overlapping extents
      for (bboxset<int,dim>::const_iterator ext_iter = myextents.begin();
           ext_iter != myextents.end();
           ++ext_iter) {
        
        // Copy data
        alldata->copy_from (mydata, *ext_iter);
        
      }
      
    } END_LOCAL_COMPONENT_LOOP;
    
    // Copy result to all processors
    if (dest_proc == -1) {
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          
          void* myhdata = rank==proc ? hdata : 0;
          gdata<dim>* const tmpdata = mydata->make_typed();
          tmpdata->allocate (alldata->extent(), proc, myhdata);
          tmpdata->copy_from (alldata, alldata->extent());
          delete tmpdata;
          
        }
      }
    } // Copy result
    
    delete alldata;
    
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
    const int vdim = CCTK_GroupDimFromVarI(vindex);
    assert (vdim>=1 && vdim<=dim);
    
    // Check some arguments
    assert (hdim>=0 && hdim<=dim);
    
    // Check output arguments
    assert (hdata);
    assert (hsize);
    
    // Calculate more convenient representation of the direction
    int dirs[dim];              // should really be dirs[hdim]
    // The following if statement is written according to the
    // definition of "dir".
    if (hdim==1) {
      // 1-dimensional hyperslab
      int mydir = 0;
      for (int d=0; d<vdim; ++d) {
	if (directions[d]!=0) {
	  mydir = d+1;
	  break;
	}
      }
      assert (mydir>0);
      for (int d=0; d<vdim; ++d) {
	if (d == mydir-1) {
	  assert (directions[d]!=0);
	} else {
	  assert (directions[d]==0);
	}
      }
      dirs[0] = mydir;
    } else if (hdim==vdim) {
      // vdim-dimensional hyperslab
      for (int d=0; d<vdim; ++d) {
	dirs[d] = d+1;
      }
    } else if (hdim==2) {
      // 2-dimensional hyperslab with vdim==3
      assert (vdim==3);
      int mydir = 0;
      for (int d=0; d<vdim; ++d) {
	if (directions[d]==0) {
	  mydir = d+1;
	  break;
	}
      }
      assert (mydir>0);
      for (int d=0; d<vdim; ++d) {
	if (d == mydir-1) {
	  assert (directions[d]==0);
	} else {
	  assert (directions[d]!=0);
	}
      }
      int dd=0;
      for (int d=0; d<vdim; ++d) {
	if (d != mydir-1) {
	  dirs[dd] = d+1;
	  ++dd;
	}
      }
      assert (dd==hdim);
    } else {
      assert (0);
    }
    // Fill remaining length
    for (int d=vdim; d<dim; ++d) {
      dirs[d] = d+1;
    }
    
    // Calculate lengths
    for (int dd=0; dd<hdim; ++dd) {
      if (lengths[dd]<0) {
	int gsh[dim];
	int ierr = CCTK_GroupgshVI(GH, dim, gsh, vindex);
	assert (!ierr);
	const int totlen = gsh[dirs[dd]-1];
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
