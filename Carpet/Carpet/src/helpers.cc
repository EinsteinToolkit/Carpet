#include <assert.h>
#include <stdlib.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/dist.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/helpers.cc,v 1.1 2001/07/04 12:29:48 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  int Barrier (cGH* cgh)
  {
    MPI_Barrier (dist::comm);
    return 0;
  }
  
  
  
  int Exit (cGH* cgh, int retval)
  {
    CCTK_Barrier (cgh);
    dist::finalize();
    exit (retval);
  }
  
  int Abort (cGH* cgh, int retval)
  {
    MPI_Abort (dist::comm, retval);
    abort ();
  }
  
  
  
  int MyProc (cGH* cgh)
  {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    return rank;
  }
  
  int nProcs (cGH* cgh)
  {
    int size;
    MPI_Comm_size (dist::comm, &size);
    return size;
  }
  
  
  
  const int* ArrayGroupSizeB (cGH* cgh, int dir, int group,
			      const char* groupname)
  {
    static const int zero = 0, one = 1;
    
    if (component == -1) {
      // global routine
      return &zero;
    }
    
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    const int gpdim = CCTK_GroupDimI(group);
    
    assert (dir>=0 && dir<gpdim);
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      
      const int var = CCTK_FirstVarIndexI(group);
      assert (var>=0 && var<CCTK_NumVars());
      
      switch (CCTK_GroupTypeFromVarI(var)) {
	
      case CCTK_SCALAR:
	assert (group<(int)scdata.size());
	return &one;
	
      case CCTK_ARRAY:
	assert (group<(int)arrdata.size());
	return &arrdata[group].size[dir];
	
      case CCTK_GF:
	assert (group<(int)gfdata.size());
	return &gfsize[dir];
	
      default:
	abort();
      }
      
    } else {
      
      // no storage
      return &zero;
      
    }
  }
  
  
  
  MPI_Comm CarpetMPICommunicator ()
  {
    return dist::comm;
  }
  
  
  
  int mintl (const checktimes where, const int num_tl)
  {
    assert (num_tl>0);
    switch (where) {
    case currenttime:
      return 0;
    case currenttimebutnotifonly:
      // don't include current time if there is only one time level
      return num_tl>1 ? 0: 1;
    case allbutlasttime:
      // do include current time if there is only one time level
      return 0;
    case allbutcurrenttime:
      return 1;
    case alltimes:
      return 0;
    default:
      abort();
    }
  }
  
  int maxtl (const checktimes where, const int num_tl)
  {
    assert (num_tl>0);
    switch (where) {
    case currenttime:
      return 0;
    case currenttimebutnotifonly:
      return 0;
    case allbutlasttime:
      return num_tl-2;
    case allbutcurrenttime:
      return num_tl-1;
    case alltimes:
      return num_tl-1;
    default:
      abort();
    }
  }
  
  
  
  void Checkpoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (verbose) {
      va_list args;
      char msg[1000];
      va_start (args, fmt);
      vsnprintf (msg, sizeof(msg), fmt, args);
      va_end (args);
      CCTK_INFO (msg);
    }
    if (barriers) {
      MPI_Barrier (dist::comm);
    }
  }
  
  
  
  void UnsupportedVarType (const int vindex)
  {
    assert (vindex>=0 && vindex<CCTK_NumVars());
    CCTK_VWarn
      (0, __LINE__, __FILE__, CCTK_THORNSTRING,
       "Carpet does not support the type of the variable \"%s\".\n"
       "Either enable support for this type, "
       "or change the type of this variable.", CCTK_FullName(vindex));
  }
  
  
  
  void set_reflevel (cGH* cgh, const int rl)
  {
    // Check
    assert (rl>=0 && rl<hh->reflevels());
    assert (component == -1);
    
    // Change
    reflevel = rl;
    const bbox<int,dim>& base = hh->baseextent;
    reflevelfact = (int)floor(pow((double)hh->reffact, reflevel)+0.5);
    cgh->cctk_delta_time = base_delta_time / reflevelfact;
    for (int d=0; d<dim; ++d) {
      cgh->cctk_gsh[d]
	= ((base.shape() / base.stride()
	    + dd->lghosts + dd->ughosts)[d] - 1) * reflevelfact + 1;
      cgh->cctk_levfac[d] = reflevelfact;
    }
  }
  
  
  
  void set_mglevel (cGH* cgh, const int ml)
  {
    assert (ml==0);
    assert (component==-1);
    mglevel = ml;
    cgh->cctk_convlevel = mglevel;
  }
  
  
  
  void set_component (cGH* cgh, const int c)
  {
    assert (c==-1 || (c>=0 && c<hh->components(reflevel)));
    component = c;
      
    if (component == -1) {
      // Global mode -- no component is active
      
      // Set Cactus parameters to pseudo values
      for (int d=0; d<dim; ++d) {
	cgh->cctk_lsh[d]      = 0xdeadbeef;
	cgh->cctk_bbox[2*d  ] = 0xdeadbeef;
	cgh->cctk_bbox[2*d+1] = 0xdeadbeef;
	cgh->cctk_lbnd[d]     = 0xdeadbeef;
	cgh->cctk_ubnd[d]     = 0xdeadbeef;
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = 0xdeadbeef;
	}
      }
      
      // Set local grid function and array sizes
      for (int d=0; d<dim; ++d) {
	gfsize[d] = 0xdeadbeef;
      }
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	const int n0 = CCTK_FirstVarIndexI(group);
	switch (CCTK_GroupTypeFromVarI(n0)) {
	case CCTK_SCALAR:
	  break;
	case CCTK_ARRAY:
	  for (int d=0; d<dim; ++d) {
	    arrdata[group].size[d] = 0xdeadbeef;
	  }
	  break;
	case CCTK_GF:
	  break;
	default:
	  abort();
	}
      }
      
      // Set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n); ++tl) {
	  
	  const int group = CCTK_GroupIndexFromVarI(n);
	  assert (group>=0);
	  const int var   = n - CCTK_FirstVarIndexI(group);
	  assert (var>=0);
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // Group has storage
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR:
	      // Scalar variables can be accessed
	      assert (group<(int)scdata.size());
	      assert (var<(int)scdata[group].data.size());
	      assert (reflevel<(int)scdata[group].data[var].size());
	      assert (tl<(int)scdata[group].data[var][reflevel].size());
	      cgh->data[n][tl] = scdata[group].data[var][reflevel][tl];
	      break;
	      
	    case CCTK_ARRAY:
	    case CCTK_GF:
	      // Arrays and grid functions cannot be accessed
	      cgh->data[n][tl] = 0;
	      break;
	      
	    default:
	      abort();
	    }
	    
	  } else {
	    // Group has no storage
	    
	    cgh->data[n][tl] = 0;
	    
	  }
	  
	} // for tl
      }	// for n
      
    } else {
      // Local mode -- a component is active
      
      // Set Cactus parameters
      for (int d=0; d<dim; ++d) {
	const bbox<int,dim>& ext
	  = dd->boxes[reflevel][component][mglevel].exterior;
	cgh->cctk_lsh[d] = (ext.shape() / ext.stride())[d];
	cgh->cctk_lbnd[d] = (ext.lower() / ext.stride())[d];
	cgh->cctk_ubnd[d] = (ext.upper() / ext.stride())[d];
	assert (cgh->cctk_lsh[d]>=0 && cgh->cctk_lsh[d]<=cgh->cctk_gsh[d]);
	assert (cgh->cctk_lbnd[d]>=0 && cgh->cctk_ubnd[d]<cgh->cctk_gsh[d]);
	assert (cgh->cctk_lbnd[d]<=cgh->cctk_ubnd[d]+1);
	// No outer boundaries on the finer grids
	cgh->cctk_bbox[2*d  ]
	  = reflevel==0 && cgh->cctk_lbnd[d] == 0;
	cgh->cctk_bbox[2*d+1]
	  = reflevel==0 && cgh->cctk_ubnd[d] == cgh->cctk_gsh[d]-1;
#if 0
	// Do allow outer boundaries on the finer grids (but this is
	// generally inconsistent -- c. f. periodicity)
	const bbox<int,dim>& base = hh->baseextent;
	cgh->cctk_bbox[2*d  ] = (ext.lower() < base.lower())[d];
	cgh->cctk_bbox[2*d+1] = (ext.upper() > base.upper())[d];
#endif
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  // TODO: support staggering
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
	}
      }
      
      // Set local grid function and array sizes
      const bbox<int,dim> ext =
	dd->boxes[reflevel][component][mglevel].exterior;
      for (int d=0; d<dim; ++d) {
	gfsize[d] = ext.shape()[d] / ext.stride()[d];
      }
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	const int n0 = CCTK_FirstVarIndexI(group);
	switch (CCTK_GroupTypeFromVarI(n0)) {
	case CCTK_SCALAR:
	  break;
	case CCTK_ARRAY: {
	  assert (reflevel < (int)arrdata[group].dd->boxes.size());
	  assert (component < (int)arrdata[group].dd->boxes[reflevel].size());
	  assert (mglevel
		  < (int)arrdata[group].dd->boxes[reflevel][component].size());
	  const bbox<int,dim> ext =
	    arrdata[group].dd->boxes[reflevel][component][mglevel].exterior;
	  for (int d=0; d<dim; ++d) {
	    arrdata[group].size[d] = ext.shape()[d] / ext.stride()[d];
	  }
	}
	case CCTK_GF:
	  break;
	default:
	  abort();
	}
      }
      
      // Set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	
	const int group = CCTK_GroupIndexFromVarI(n);
	assert (group>=0);
	const int var   = n - CCTK_FirstVarIndexI(group);
	assert (var>=0);
	const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	assert (num_tl>0);
	
	for (int tl=0; tl<num_tl; ++tl) {
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // Group has storage
	    
	    switch (CCTK_GroupTypeFromVarI(n)) {
	      
	    case CCTK_SCALAR:
	      assert (group<(int)scdata.size());
	      assert (var<(int)scdata[group].data.size());
	      assert (reflevel<(int)scdata[group].data[var].size());
	      assert (tl<(int)scdata[group].data[var][reflevel].size());
	      cgh->data[n][tl] = scdata[group].data[var][reflevel][tl];
	      break;
	      
	    case CCTK_ARRAY:
	      assert (group<(int)arrdata.size());
	      assert (var<(int)arrdata[group].data.size());
	      cgh->data[n][tl]
		= ((*arrdata[group].data[var])
		   (-tl, reflevel, component, mglevel)->storage());
	      break;
	      
	    case CCTK_GF:
	      assert (group<(int)gfdata.size());
	      assert (var<(int)gfdata[group].data.size());
	      cgh->data[n][tl]
		= ((*gfdata[group].data[var])
		   (-tl, reflevel, component, mglevel)->storage());
	      break;
	      
	    default:
	      abort();
	    }
	    if (CCTK_GroupTypeFromVarI(n)==CCTK_SCALAR
		|| hh->is_local(reflevel,component)) {
	      assert (cgh->data[n][tl]);
	    } else {
	      assert (! cgh->data[n][tl]);
	    }
	    
	  } else {
	    // Group has no storage
	    
	    cgh->data[n][tl] = 0;
	    
	  } // if ! has storage
	  
	} // for tl
      } // for n
      
    }
  }
  
} // namespace Carpet
