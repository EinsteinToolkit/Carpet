// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.cc,v 1.17 2004/04/16 11:52:18 schnetter Exp $

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "cctk.h"

#include "util_Table.h"

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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/slab.cc,v 1.17 2004/04/16 11:52:18 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetSlab_slab_cc);
}



namespace CarpetSlab {
  
  using namespace Carpet;
  
  
  
  // Mapping object
  // (just store the mapping)
  struct mapping {
    int vindex;
    int hdim;
    int * origin;               // [vdim]
    int * dirs;                 // [hdim]
    int * stride;               // [hdim]
    int * length;               // [hdim]
  };
  
  
  
  int StoreMapping (mapping * const mp)
  {
    int const table = Util_TableCreate (UTIL_TABLE_FLAGS_DEFAULT);
    assert (table>=0);
    int const ierr = Util_TableSetPointer (table, mp, "mapping");
    assert (ierr>=0);
    return table;
  }
  
  mapping * RetrieveMapping (int const table)
  {
    CCTK_POINTER mp;
    int const ierr = Util_TableGetPointer (table, &mp, "mapping");
    assert (ierr>=0);
    return (mapping *)mp;
  }
  
  void DeleteMapping (int const table)
  {
    int const ierr = Util_TableDestroy (table);
    assert (ierr>=0);
  }
  
  
  
  void FillSlab (const cGH* const cgh,
		 const int dest_proc,
		 const int n,
		 const int ti,
		 const int hdim,
		 const int origin[/*vdim*/],
		 const int dirs[/*hdim*/],
		 const int stride[/*hdim*/],
		 const int length[/*hdim*/],
                 void* const hdata)
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
      CCTK_WARN (0, "It is not possible to use hyperslabbing for a grid function in global mode (use singlemap mode instead)");
    }
    const int rl = gp.grouptype==CCTK_GF ? reflevel : 0;
    
    if (gp.grouptype==CCTK_GF && Carpet::map==-1) {
      CCTK_WARN (0, "It is not possible to use hyperslabbing for a grid function in level mode (use singlemap mode instead)");
    }
    const int m = gp.grouptype==CCTK_GF ? Carpet::map : 0;
    
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
    myhh = arrdata.at(group).at(m).hh;
    assert (myhh);
    mydd = arrdata.at(group).at(m).dd;
    assert (mydd);
    assert (var < (int)arrdata.at(group).at(m).data.size());
    myff = arrdata.at(group).at(m).data.at(var);
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
    assert (hdata);
    if (dest_proc==-1 || rank==dest_proc) {
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
    assert (hextent.size() == totalsize);
    
    // Create collector data object
    void* myhdata = rank==collect_proc ? hdata : 0;
    gdata<dim>* const alldata = mydata->make_typed(-1);
    alldata->allocate (hextent, collect_proc, myhdata);
    
    // Done with the temporary stuff
    mydata = 0;
    
    for (comm_state<dim> state; !state.done(); state.step()) {
      
      // Loop over all components, copying data from them
      BEGIN_LOCAL_COMPONENT_LOOP (cgh, gp.grouptype) {
        
        // Get data object
        mydata = (*myff)(tl, rl, component, mglevel);
        
        // Calculate overlapping extents
        const bboxset<int,dim> myextents
          = ((mydd->boxes.at(rl).at(component).at(mglevel).sync_not
              | mydd->boxes.at(rl).at(component).at(mglevel).interior)
             & hextent);
        
        // Loop over overlapping extents
        for (bboxset<int,dim>::const_iterator ext_iter = myextents.begin();
             ext_iter != myextents.end();
             ++ext_iter) {
          
          // Copy data
          alldata->copy_from (state, mydata, *ext_iter);
          
        }
        
      } END_LOCAL_COMPONENT_LOOP;
      
    } // for step
    
    // Copy result to all processors
    if (dest_proc == -1) {
      vector<gdata<dim>*> tmpdata(CCTK_nProcs(cgh));
      vector<comm_state<dim> > state;
      
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          void* myhdata = rank==proc ? hdata : 0;
          tmpdata.at(proc) = mydata->make_typed(-1);
          tmpdata.at(proc)->allocate (alldata->extent(), proc, myhdata);
          tmpdata.at(proc)->copy_from (state.at(proc), alldata, alldata->extent());
        }
      }
      
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          tmpdata.at(proc)->copy_from (state.at(proc), alldata, alldata->extent());
        }
      }
      
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          tmpdata.at(proc)->copy_from (state.at(proc), alldata, alldata->extent());
          delete tmpdata.at(proc);
        }
      }
      
    } // Copy result
    
    delete alldata;
  }
  
  
  
  void* GetSlab (const cGH* const cgh,
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
      CCTK_WARN (0, "It is not possible to use hyperslabbing for a grid function in global mode (use singlemap mode instead)");
    }
    const int rl = gp.grouptype==CCTK_GF ? reflevel : 0;
    
    if (gp.grouptype==CCTK_GF && Carpet::map==-1) {
      CCTK_WARN (0, "It is not possible to use hyperslabbing for a grid function in level mode (use singlemap mode instead)");
    }
    const int m = gp.grouptype==CCTK_GF ? Carpet::map : 0;
    
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
    myhh = arrdata.at(group).at(m).hh;
    assert (myhh);
    mydd = arrdata.at(group).at(m).dd;
    assert (mydd);
    assert (var < (int)arrdata.at(group).at(m).data.size());
    myff = arrdata.at(group).at(m).data.at(var);
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
      assert (0);
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
    assert (hextent.size() == totalsize);
    
    // Create collector data object
    void* myhdata = rank==collect_proc ? hdata : 0;
    gdata<dim>* const alldata = mydata->make_typed(-1);
    alldata->allocate (hextent, collect_proc, myhdata);
    
    // Done with the temporary stuff
    mydata = 0;
    
    for (comm_state<dim> state; !state.done(); state.step()) {
      
      // Loop over all components, copying data from them
      BEGIN_LOCAL_COMPONENT_LOOP (cgh, gp.grouptype) {
        
        // Get data object
        mydata = (*myff)(tl, rl, component, mglevel);
        
        // Calculate overlapping extents
        const bboxset<int,dim> myextents
          = ((mydd->boxes.at(rl).at(component).at(mglevel).sync_not
              | mydd->boxes.at(rl).at(component).at(mglevel).interior)
             & hextent);
        
        // Loop over overlapping extents
        for (bboxset<int,dim>::const_iterator ext_iter = myextents.begin();
             ext_iter != myextents.end();
             ++ext_iter) {
          
          // Copy data
          alldata->copy_from (state, mydata, *ext_iter);
          
        }
        
      } END_LOCAL_COMPONENT_LOOP;
      
    } // for step
    
    // Copy result to all processors
    if (dest_proc == -1) {
      vector<gdata<dim>*> tmpdata(CCTK_nProcs(cgh));
      vector<comm_state<dim> > state;
      
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          void* myhdata = rank==proc ? hdata : 0;
          tmpdata.at(proc) = mydata->make_typed(-1);
          tmpdata.at(proc)->allocate (alldata->extent(), proc, myhdata);
          tmpdata.at(proc)->copy_from (state.at(proc), alldata, alldata->extent());
        }
      }
      
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          tmpdata.at(proc)->copy_from (state.at(proc), alldata, alldata->extent());
        }
      }
      
      for (int proc=0; proc<CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          tmpdata.at(proc)->copy_from (state.at(proc), alldata, alldata->extent());
          delete tmpdata.at(proc);
        }
      }
      
    } // Copy result
    
    delete alldata;
    
    // Success
    return hdata;
  }
  
  
  
  CCTK_INT CarpetSlab_Get (CCTK_POINTER_TO_CONST const cctkGH_,
                           CCTK_INT const mapping_handle,
                           CCTK_INT const proc,
                           CCTK_INT const vindex,
                           CCTK_INT const timelevel,
                           CCTK_INT const hdatatype,
                           void * const hdata)
  {
    cGH const * const cctkGH = (cGH const *) cctkGH_;
    
    // Check arguments
    assert (cctkGH);
    assert (mapping_handle>=0);
    assert (proc==-1 || proc>=0 && proc<CCTK_nProcs(cctkGH));
    assert (vindex>=0 && vindex<CCTK_NumVars());
    assert (timelevel>=0);
    assert (hdatatype>=0);
    assert (hdata);
    
    // Get mapping
    const mapping * const mp = RetrieveMapping (mapping_handle);
    assert (mp);
    
    // Calculate total size
    size_t size = 1;
    for (size_t d=0; d<(size_t)mp->hdim; ++d) {
      size *= mp->length[d];
    }
    
    // Get type size
    size_t const sz = CCTK_VarTypeSize (hdatatype);
    assert (sz>0);
    
    // Forward call
    FillSlab (cctkGH, proc, vindex, timelevel,
              mp->hdim, mp->origin, mp->dirs, mp->stride, mp->length, hdata);
    
    return 0;
  }
  
  
  
  CCTK_INT CarpetSlab_GetList (CCTK_POINTER_TO_CONST const cctkGH_,
                               CCTK_INT const mapping_handle,
                               CCTK_INT const num_arrays,
                               CCTK_INT const * const procs,
                               CCTK_INT const * const vindices,
                               CCTK_INT const * const timelevels,
                               CCTK_INT const * const hdatatypes,
                               void * const * const hdata,
                               CCTK_INT * const retvals)
  {
    cGH const * const cctkGH = (cGH const *) cctkGH_;
    
    // Check arguments
    assert (cctkGH);
    assert (mapping_handle>=0);
    assert (num_arrays>=0);
    assert (procs);
    assert (vindices);
    assert (timelevels);
    assert (hdatatypes);
    assert (hdata);
    assert (retvals);
    
    // Remember whether there were errors
    bool everyting_okay = true;
    
    // Loop over all slabs
    for (int n=0; n<num_arrays; ++n) {
      // Forward call
      retvals[n] = CarpetSlab_Get (cctkGH, mapping_handle, procs[n],
                                   vindices[n], timelevels[n], hdatatypes[n],
                                   hdata[n]);
      everyting_okay = everyting_okay && retvals[n];
    }
    
    return everyting_okay ? 0 : -1;
  }
  
  
  
  CCTK_INT CarpetSlab_LocalMappingByIndex (CCTK_POINTER_TO_CONST const cctkGH_,
                                           CCTK_INT const vindex,
                                           CCTK_INT const hdim,
                                           CCTK_INT const * const direction,
                                           CCTK_INT const * const origin,
                                           CCTK_INT const * const extent,
                                           CCTK_INT const * const downsample_,
                                           CCTK_INT const table_handle,
                                           CCTK_INT (* const conversion_fn) (CCTK_INT const nelems,
                                                                             CCTK_INT const src_stride,
                                                                             CCTK_INT const dst_stride,
                                                                             CCTK_INT const src_type,
                                                                             CCTK_INT const dst_type,
                                                                             void const * const from,
                                                                             void * const to),
                                           CCTK_INT * const hsize_local,
                                           CCTK_INT * const hsize_global,
                                           CCTK_INT * const hoffset_global)
  {
    CCTK_WARN (0, "not implemented");
    return 0;
  }
  
  
  
  CCTK_INT CarpetSlab_GlobalMappingByIndex (CCTK_POINTER_TO_CONST const cctkGH_,
                                            CCTK_INT const vindex,
                                            CCTK_INT const hdim,
                                            CCTK_INT const * const direction,
                                            CCTK_INT const * const origin,
                                            CCTK_INT const * const extent,
                                            CCTK_INT const * const downsample_,
                                            CCTK_INT const table_handle,
                                            CCTK_INT (* const conversion_fn) (CCTK_INT const nelems,
                                                                              CCTK_INT const src_stride,
                                                                              CCTK_INT const dst_stride,
                                                                              CCTK_INT const src_type,
                                                                              CCTK_INT const dst_type,
                                                                              void const * const from,
                                                                              void * const to),
                                            CCTK_INT * const hsize)
  {
    cGH const * const cctkGH = (cGH const *) cctkGH_;
    
    // Check arguments
    assert (cctkGH);
    assert (vindex>=0 && vindex<CCTK_NumVars());
    assert (hdim>=0 && hdim<=dim);
    assert (direction);
    assert (origin);
    assert (extent);
    assert (table_handle>=0);
    assert (hsize);
    
    // Get more information
    int const vdim = CCTK_GroupDimFromVarI (vindex);
    assert (vdim>=0 && vdim<dim);
    assert (hdim<=vdim);
    
    // Not implemented
    assert (! conversion_fn);
    
    // Allocate memory
    mapping * mp = new mapping;
    mp->origin = new int[vdim];
    mp->dirs   = new int[hdim];
    mp->stride = new int[hdim];
    mp->length = new int[hdim];
    
    // Calculate more convenient representation of the direction
    int dirs[dim];              // should really be dirs[hdim]
    // The following if statement is written according to the
    // definition of "dir".
    if (hdim==1) {
      // 1-dimensional hyperslab
      int mydir = 0;
      for (int d=0; d<vdim; ++d) {
	if (direction[d]!=0) {
	  mydir = d+1;
	  break;
	}
      }
      assert (mydir>0);
      for (int d=0; d<vdim; ++d) {
	if (d == mydir-1) {
	  assert (direction[d]!=0);
	} else {
	  assert (direction[d]==0);
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
	if (direction[d]==0) {
	  mydir = d+1;
	  break;
	}
      }
      assert (mydir>0);
      for (int d=0; d<vdim; ++d) {
	if (d == mydir-1) {
	  assert (direction[d]==0);
	} else {
	  assert (direction[d]!=0);
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
    vector<CCTK_INT> downsample(hdim);
    for (int dd=0; dd<hdim; ++dd) {
      if (extent[dd]<0) {
	int gsh[dim];
	int ierr = CCTK_GroupgshVI(cctkGH, dim, gsh, vindex);
	assert (!ierr);
	const int totlen = gsh[dirs[dd]-1];
	assert (totlen>=0);
	// Partial argument check
	assert (origin[dirs[dd]-1]>=0);
	assert (origin[dirs[dd]-1]<=totlen);
        downsample[dd] = downsample_ ? downsample_[dd] : 1;
	assert (downsample[dd]>0);
	hsize[dd] = (totlen - origin[dirs[dd]-1]) / downsample[dd];
      } else {
	hsize[dd] = extent[dd];
      }
      assert (hsize[dd]>=0);
    }
    
    // Store information
    mp->vindex = vindex;
    mp->hdim = hdim;
    for (size_t d=0; d<(size_t)hdim; ++d) {
      mp->origin[d] = origin[d];
      mp->dirs[d]   = dirs[d];
      mp->stride[d] = downsample[d];
      mp->length[d] = hsize[d];
    }
    
    return 0;
  }
  
  
  
  CCTK_INT CarpetSlab_FreeMapping (CCTK_INT const mapping_handle)
  {
    // Check arguments
    assert (mapping_handle>=0);
    
    // Get mapping
    mapping * mp = RetrieveMapping (mapping_handle);
    assert (mp);
    
    // Delete storage
    DeleteMapping (mapping_handle);
    
    // Free memory
    delete [] mp->origin;
    delete [] mp->dirs;
    delete [] mp->stride;
    delete [] mp->length;
    delete mp;
    
    return 0;
  }
  
  
  
  int Hyperslab_GetHyperslab (const cGH* const GH,
			      const int target_proc,
			      const int vindex,
			      const int vtimelvl,
			      const int hdim,
			      const int global_startpoint [/*vdim*/],
			      const int directions [/*vdim*/],
			      const int lengths [/*hdim*/],
			      const int downsample_ [/*hdim*/],
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
    vector<int> downsample(hdim);
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
        downsample[dd] = downsample_ ? downsample_[dd] : 1;
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
		      &downsample[0],
		      hsize);
    
    // Return with success
    return 1;
  }
  
  
  
} // namespace CarpetSlab
