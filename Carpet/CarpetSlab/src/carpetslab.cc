// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetSlab/src/Attic/carpetslab.cc,v 1.3 2001/03/07 13:01:11 eschnett Exp $

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/gdata.hh"
#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/dist.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "carpetslab.hh"



namespace CarpetSlab {
  
  using namespace Carpet;
  
  
  
  int Hyperslab_GetLocalHyperslab (cGH* cgh,
				   const int n,
				   const int tl,
				   const int hdim,
				   const int origin [/*vdim*/],
				   const int dir [/*vdim*/],
				   const int len [/*hdim*/],
				   const int downsample [/*hdim*/],
				   void** const hdata,
				   int hsize [/*hdim*/],
				   int ghsize [/*hdim*/],
				   int hoffset [/*hdim*/])
  {
    // check current status
    assert (mglevel>=0);
    assert (reflevel>=0);
    assert (component>=0);	// local := this component
    
    // check arguments
    assert (n>=0 && n<CCTK_NumVars());
    assert (tl>=0);
    assert (hdim>0 && hdim<=dim);
    // the following assertion is too strict (better allow arbitrary values)
    {
      const int group = CCTK_GroupIndexFromVarI(n);
      assert (group>=0);
      switch (CCTK_GroupTypeFromVarI(n)) {
      case CCTK_SCALAR:
	abort();
      case CCTK_ARRAY: {
	assert (group<(int)arrdata.size());
// 	assert (arrdata[group].hh->is_local(reflevel, component));
// 	const bbox<int,dim> ext =
// 	  arrdata[group].dd->boxes[reflevel][component][mglevel].exterior;
// 	for (int d=0; d<dim; ++d) {
// 	  assert (origin[d] >= ext.lower()[d] / ext.stride()[d]
// 		  && origin[d] <= ext.upper()[d] / ext.stride()[d]);
// 	}
	break;
      }
      case CCTK_GF: {
	assert (group<(int)gfdata.size());
// 	assert (hh->is_local(reflevel, component));
// 	const bbox<int,dim> ext =
// 	  dd->boxes[reflevel][component][mglevel].exterior;
// 	for (int d=0; d<dim; ++d) {
// 	  assert (origin[d] >= ext.lower()[d] / ext.stride()[d]
// 		  && origin[d] <= ext.upper()[d] / ext.stride()[d]);
// 	}
	break;
      }
      default:
	abort();
      }
    }
    // the following assertion is too strict (better allow arbitrary values)
    for (int d=0; d<dim; ++d) assert (dir[d]==0 || dir[d]==1);
    // the following assertion is too strict (better allow arbitrary values)
    for (int d=0; d<hdim; ++d) assert (downsample[d]>0);
    assert (hdata);
    assert (hsize);
    assert (ghsize);
    assert (hoffset);
    
    // get variable info
    const int group = CCTK_GroupIndexFromVarI(n);
    assert (group>=0);
    const int var = n - CCTK_FirstVarIndexI(group);
    assert (var>=0);
    
    // get grid hierarchy data hierarchy and grid function
    gh<dim>* myhh;
    dh<dim>* mydd;
    generic_gf<dim>* myff;
    switch (CCTK_GroupTypeFromVarI(n)) {
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
      mydd = dd;
      assert (group < (int)gfdata.size());
      assert (var < (int)gfdata[group].data.size());
      myff = gfdata[group].data[var];
      break;
    default:
      abort();
    }
    
    // get data
    const generic_data<dim>* mydata
      = (*myff)(tl, reflevel, component, mglevel);
    
    // get local and global bounding boxes
    assert (reflevel < (int)mydd->boxes.size());
    assert (component < (int)mydd->boxes[reflevel].size());
    assert (mglevel < (int)mydd->boxes[reflevel][component].size());
    const bbox<int,dim> locbox
      = mydd->boxes[reflevel][component][mglevel].exterior;
    const bbox<int,dim> globox
      = mydd->bases[reflevel][mglevel].exterior;
    
    // calculate more convenient representation of the direction
    vect<int,dim> stride[hdim];
    // the following if statement is written according to the
    // definition of "dir".
    if (hdim==1) {
      for (int d=0; d<dim; ++d) stride[0][d] = dir[d] * downsample[0];
    } else if (hdim==dim) {
      for (int dd=0; dd<hdim; ++dd) {
	for (int d=0; d<dim; ++d) stride[dd][d] = d==dd ? downsample[dd] : 0;
      }
    } else if (hdim==2) {
      assert (dim==3);
      if (dir[0]==0) {
	assert (dir[1]!=0 && dir[2]!=0);
	stride[0] = vect<int,dim>::dir(1);
	stride[1] = vect<int,dim>::dir(2);
      } else if (dir[1]==0) {
	assert (dir[0]!=0 && dir[2]!=0);
	stride[0] = vect<int,dim>::dir(0);
	stride[1] = vect<int,dim>::dir(2);
      } else if (dir[2]==0) {
	assert (dir[0]!=0 && dir[1]!=0);
	stride[0] = vect<int,dim>::dir(0);
	stride[1] = vect<int,dim>::dir(1);
      } else {
	abort();
      }
      for (int dd=0; dd<hdim; ++dd) stride[dd] *= downsample[dd];
    } else {
      abort();
    }
    for (int dd=0; dd<hdim; ++dd) stride[dd] *= locbox.stride();
    
    // local lower bound
    vect<int,dim> lbound;
    for (int d=0; d<dim; ++d) lbound[d] = origin[d];
    lbound *= locbox.stride();
    
    // local upper bound
    vect<int,dim> ubound = lbound;
    for (int dd=0; dd<hdim; ++dd) {
      if (len[dd]<0) {
	assert (any(stride[dd]>0));
	while (all(ubound + stride[dd] <= locbox.upper())) {
	  ubound += stride[dd];
	}
      } else {
	ubound += stride[dd] * (len[dd]-1);
      }
    }
    
    lbound = max(lbound, locbox.lower());
    ubound = min(ubound, locbox.upper());
    
    // local size
    int total_hsize = 1;
    for (int dd=0; dd<hdim; ++dd) {
      hsize[dd] = 0;
      assert (any(stride[dd]>0));
      while (all(lbound + stride[dd] * hsize[dd] <= ubound)) ++hsize[dd];
      total_hsize *= hsize[dd];
    }
    assert (total_hsize>=0);
    
    // sanity check
    if (total_hsize>0) {
      vect<int,dim> tmp = lbound;
      for (int dd=0; dd<hdim; ++dd) {
	tmp += stride[dd] * (hsize[dd]-1);
      }
      assert (all(tmp == ubound));
    }
    
    // global lower bound
    vect<int,dim> glbound;
    for (int d=0; d<dim; ++d) glbound[d] = origin[d];
    glbound *= globox.stride();
    
    // global upper bound
    vect<int,dim> gubound = glbound;
    for (int dd=0; dd<hdim; ++dd) {
      if (len[dd]<0) {
	assert (any(stride[dd]>0));
	while (all(gubound + stride[dd] <= globox.upper())) {
	  gubound += stride[dd];
	}
      } else {
	gubound += stride[dd] * (len[dd]-1);
      }
    }
    
    glbound = max(glbound, globox.lower());
    gubound = min(gubound, globox.upper());
    
    // global size
    int total_ghsize = 1;
    for (int dd=0; dd<hdim; ++dd) {
      ghsize[dd] = 0;
      assert (any(stride[dd]>0));
      while (all(glbound + stride[dd] * ghsize[dd] <= gubound)) {
	++ghsize[dd];
      }
      total_ghsize *= ghsize[dd];
    }
    assert (total_ghsize>=0);
    
    // sanity check
    if (total_ghsize>0) {
      vect<int,dim> tmp = glbound;
      for (int dd=0; dd<hdim; ++dd) {
	tmp += stride[dd] * (ghsize[dd]-1);
      }
      assert (all(tmp == gubound));
    }
    
    // local to global offset
    for (int dd=0; dd<hdim; ++dd) {
      hoffset[dd] = 0;
      assert (any(stride[dd]>0));
      while (all(glbound + stride[dd] * (hoffset[dd]+1) <= lbound)) {
	++hoffset[dd];
      }
    }
    
    // sanity check
    {
      vect<int,dim> tmp = glbound;
      for (int dd=0; dd<hdim; ++dd) {
	tmp += stride[dd] * hoffset[dd];
      }
      assert (all(tmp == lbound));
    }
    
    // bail out if this component is on another processor
    if (! myhh->is_local(reflevel, component)) {
      *hdata = 0;
      for (int dd=0; dd<hdim; ++dd) {
	hsize[dd] = 0;
	hoffset[dd] = 0;
      }
      return -1;
    }
    
    // allocate the memory
    *hdata = malloc(total_hsize * CCTK_VarTypeSize(CCTK_VarTypeI(n)));
    assert (*hdata);
    
    if (total_hsize>0) {
      
      // copy the data to user memory
      char* const       dest = (char*)*hdata;
      const char* const src  = (const char*)mydata->storage();
      const int         sz   = CCTK_VarTypeSize(CCTK_VarTypeI(n));
      
      int dest_index[hdim];
      for (int dd=0; dd<hdim; ++dd) dest_index[dd] = 0;
      for (;;) {
	
	vect<int,dim> src_index = lbound;
	for (int dd=0; dd<hdim; ++dd) src_index += stride[dd] * dest_index[dd];
	
	int di = 0;
	for (int dd=hdim-1; dd>=0; --dd) di = di * hsize[dd] + dest_index[dd];
	
	const int si = mydata->offset(src_index);
	
	memcpy(dest + sz*di, src + sz*si, sz);
	
	for (int dd=0; dd<hdim; ++dd) {
	  ++dest_index[dd];
	  if (dest_index[dd]<hsize[dd]) break;
	  dest_index[dd]=0;
	  if (dd==hdim-1) goto done;
	}
      }
    done: ;
      
    } // if total_hsize>0
    
    return 0;
  }
  
  
  
  int Hyperslab_GetHyperslab (cGH* cgh,
			      const int target_proc,
			      const int n,
			      const int tl,
			      const int hdim,
			      const int origin [/*vdim*/],
			      const int dir [/*vdim*/],
			      const int len [/*hdim*/],
			      const int downsample [/*hdim*/],
			      void** const hdata,
			      int hsize [/*hdim*/])
  {
    // check current status
    assert (mglevel>=0);
    assert (reflevel>=0);
    
    // in order to work, this requires that all processors have the
    // same number of components
    const int saved_component = component;
    component = -1;
    
    // check arguments
    assert (n>=0 && n<CCTK_NumVars());
    assert (tl>=0);
    assert (hdim>0 && hdim<=dim);
    // the following assertion is too strict (better allow arbitrary values)
    {
      // TODO: make sure that origin is within the extent of this
      // refinement / multigrid level
      // (but no such extent is stored in dh)
      // (it is now; use it!)
      const int group = CCTK_GroupIndexFromVarI(n);
      assert (group>=0);
      switch (CCTK_GroupTypeFromVarI(n)) {
      case CCTK_SCALAR:
	abort();
      case CCTK_ARRAY:
	assert (group<(int)arrdata.size());
	break;
      case CCTK_GF:
	assert (group<(int)gfdata.size());
	break;
      default:
	abort();
      }
    }
    // the following assertion is too strict (better allow arbitrary values)
    for (int d=0; d<dim; ++d) assert (dir[d]==0 || dir[d]==1);
    // the following assertion is too strict (better allow arbitrary values)
    for (int d=0; d<hdim; ++d) assert (downsample[d]>0);
    assert (hdata);
    assert (hsize);
    
    int collect_proc = target_proc;
    if (collect_proc<0) collect_proc = 0;
    
    assert (hh->components(reflevel)>0);
    *hdata = 0;
    for (int dd=0; dd<hdim; ++dd) hsize[dd] = 0;
    int totalhsize = 0;
    
    const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
    
    MPI_Datatype type;
    switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)				\
    case N: {					\
      assert (sz == sizeof(T));			\
      T dummy;					\
      type = dist::datatype(dummy);		\
      break;					\
    }
#include "Carpet/Carpet/src/typecase"
#undef TYPECASE
    default:
      abort();
    }
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    
    // loop over all components
    for (component=0; component<hh->components(reflevel); ++component) {
      
      void* myhdata;
      int myhsize[hdim], ghsize[hdim], hoffset[hdim];
      
      const int retval = Hyperslab_GetLocalHyperslab
	(cgh, n, tl, hdim, origin, dir, len, downsample,
	 &myhdata, myhsize, ghsize, hoffset);
      
      if (hh->is_local(reflevel,component)) {
	assert (retval == 0);
      } else {
	assert (retval == -1);
      }
      
      int mytotalsize = 1;
      for (int dd=0; dd<hdim; ++dd) mytotalsize *= myhsize[dd];
      
      if (component==0) {
	if (target_proc<0 || rank == target_proc) {
	  
	  for (int dd=0; dd<hdim; ++dd) hsize[dd] = ghsize[dd];
	  
	  totalhsize = 1;
	  for (int dd=0; dd<hdim; ++dd) totalhsize *= hsize[dd];
	  
	  if (rank == collect_proc) {
	    *hdata = malloc(totalhsize * sz);
	    assert (*hdata);
	  } else {
	    *hdata = 0;
	  }
	
	}
      }
      
      if (!myhdata && rank == collect_proc) {
	MPI_Status status;
	assert (mytotalsize==0);
	const int src = hh->proc(reflevel, component);
	MPI_Recv (&mytotalsize, 1, MPI_INT, src, 2001,
		  dist::comm, &status);
	myhdata = malloc(mytotalsize * sz);
	assert (myhdata);
	MPI_Recv (myhdata, mytotalsize, type, src, 2001,
		  dist::comm, &status);
      } else if (myhdata && rank != collect_proc) {
	MPI_Send (&mytotalsize, 1, MPI_INT, collect_proc, 2001, dist::comm);
	MPI_Send (myhdata, mytotalsize, type, collect_proc, 2001, dist::comm);
	free (myhdata);
	myhdata = 0;
	mytotalsize = 0;
      }
      
      if (myhdata>0) {
	assert (rank == collect_proc);
	
	if (mytotalsize>0) {
	  int dest_index[hdim], src_index[hdim];
	  for (int dd=0; dd<hdim; ++dd) dest_index[dd] = hoffset[dd];
	  for (int dd=0; dd<hdim; ++dd) src_index[dd] = 0;
	  for (;;) {
	    int di=0;
	    for (int dd=hdim-1; dd>=0; --dd) {
	      di = di * hsize[dd] + dest_index[dd];
	    }
	    int si=0;
	    for (int dd=hdim-1; dd>=0; --dd) {
	      si = si * myhsize[dd] + src_index[dd];
	    }
	    
	    memcpy ((char*)*hdata + sz*di, (char*)myhdata + sz*si, sz);
	    
	    for (int dd=0; dd<hdim; ++dd) {
	      ++dest_index[dd];
	      ++src_index[dd];
	      if (src_index[dd] < myhsize[dd]) break;
	      dest_index[dd] = hoffset[dd];
	      src_index[dd] = 0;
	      if (dd==hdim-1) goto done;
	    }
	  }
	done: ;
	  
	} // if mytotalsize>0
	
	free (myhdata);
	myhdata = 0;
      } else {
	assert (rank != collect_proc);
      }
      
    }
    
    if (target_proc<0) {
      MPI_Bcast (hdata, totalhsize*sz, MPI_BYTE,
		 collect_proc, CarpetMPICommunicator());
    }
    
    component = saved_component;
    
    return 0;
  }
  
} // namespace CarpetSlab
