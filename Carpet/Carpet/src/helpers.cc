#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/dist.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/helpers.cc,v 1.15 2002/01/09 13:56:26 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  int Barrier (const cGH* cgh)
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
  
  
  
  int MyProc (const cGH* cgh)
  {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    return rank;
  }
  
  int nProcs (const cGH* cgh)
  {
    int size;
    MPI_Comm_size (dist::comm, &size);
    return size;
  }
  
  
  
  MPI_Comm CarpetMPIComm ()
  {
    return dist::comm;
  }
  
  
  
  MPI_Datatype CarpetMPIDatatype (const int vartype)
  {
    switch (vartype) {
#define TYPECASE(N,T)				\
    case N: {					\
      T dummy;					\
      return dist::datatype(dummy);		\
    }
#include "typecase"
#undef TYPECASE
    default:
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Carpet does not support the variable type %d.", vartype);
    }
    // notreached
    return MPI_CHAR;
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
  
  
  
  void Waypoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (verbose || veryverbose) {
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
  
  void Checkpoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (veryverbose) {
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
    reflevelfact = ipow(hh->reffact, reflevel);
    cgh->cctk_delta_time = base_delta_time / reflevelfact;
    for (int d=0; d<dim; ++d) {
      cgh->cctk_gsh[d]
	= ((base.shape() / base.stride()
	    + dd->lghosts + dd->ughosts)[d] - 1) * reflevelfact + 1;
      cgh->cctk_levfac[d] = reflevelfact;
    }
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      const bbox<int,dim>& base = arrdata[group].hh->baseextent;
      for (int d=0; d<dim; ++d) {
	((int*)arrdata[group].info.gsh)[d]
	  = (((base.shape() / base.stride()
	       + arrdata[group].dd->lghosts + arrdata[group].dd->ughosts)[d]
	      - 1)
	     * reflevelfact + 1);
      }
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
    assert (reflevel>=0 && reflevel<hh->reflevels());
    assert (c==-1 || (c>=0 && c<hh->components(reflevel)));
    component = c;
    assert (component==-1
	    || (mglevel>=0 && mglevel<hh->mglevels(reflevel,component)));
    
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
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	for (int d=0; d<dim; ++d) {
	  ((int*)arrdata[group].info.lsh)[d]      = 0xdeadbeef;
	  ((int*)arrdata[group].info.bbox)[2*d  ] = 0xdeadbeef;
	  ((int*)arrdata[group].info.bbox)[2*d+1] = 0xdeadbeef;
	  ((int*)arrdata[group].info.lbnd)[d]     = 0xdeadbeef;
	  ((int*)arrdata[group].info.ubnd)[d]     = 0xdeadbeef;
	}
      }
      
      // Set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	for (int tl=0; tl<CCTK_NumTimeLevelsFromVarI(n); ++tl) {
	  
	  const int group = CCTK_GroupIndexFromVarI(n);
	  assert (group>=0);
	  const int var   = n - CCTK_FirstVarIndexI(group);
	  assert (var>=0);
	  
	  if (CCTK_GroupTypeI(group) != CCTK_SCALAR) {
	    // Arrays and grid functions cannot be accessed
	    
	    cgh->data[n][tl] = 0;
	    
	  } else {
	    // Scalars can be accessed
	    
	    if (CCTK_QueryGroupStorageI(cgh, group)) {
	      // Group has storage
	      
	      assert (group<(int)arrdata.size());
	      assert (var<(int)arrdata[group].data.size());
	      assert (arrdata[group].data[var]);
	      const int c = CCTK_MyProc(cgh);
	      assert (hh->is_local(reflevel,c));
	      cgh->data[n][tl]
		= ((*arrdata[group].data[var])
		   (-tl, reflevel, c, mglevel)->storage());
	      assert (cgh->data[n][tl]);
	      
	    } else {
	      // Group has no storage
	      
	      cgh->data[n][tl] = 0;
	      
	    } // if ! has storage
	    
	  } // if group type is SCALAR
	  
	} // for tl
      }	// for n
      
    } else {
      // Local mode -- a component is active
      
      // Set Cactus parameters
      {
	assert (reflevel < (int)dd->boxes.size());
	assert (component < (int)dd->boxes[reflevel].size());
	assert (mglevel < (int)dd->boxes[reflevel][component].size());
	const bbox<int,dim>& ext
	  = dd->boxes[reflevel][component][mglevel].exterior;
	for (int d=0; d<dim; ++d) {
	  cgh->cctk_lsh[d] = (ext.shape() / ext.stride())[d];
	  cgh->cctk_lbnd[d] = (ext.lower() / ext.stride())[d];
	  cgh->cctk_ubnd[d] = (ext.upper() / ext.stride())[d];
	  assert (cgh->cctk_lsh[d]>=0 && cgh->cctk_lsh[d]<=cgh->cctk_gsh[d]);
	  assert (cgh->cctk_lbnd[d]>=0 && cgh->cctk_ubnd[d]<cgh->cctk_gsh[d]);
	  assert (cgh->cctk_lbnd[d]<=cgh->cctk_ubnd[d]+1);
	  // Do allow outer boundaries on the finer grids
	  cgh->cctk_bbox[2*d  ] = cgh->cctk_lbnd[d] == 0;
	  cgh->cctk_bbox[2*d+1] = cgh->cctk_ubnd[d] == cgh->cctk_gsh[d]-1;
	  for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	    // TODO: support staggering
	    cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
	  }
	}
      }
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	assert (reflevel < (int)arrdata[group].dd->boxes.size());
	assert (component < (int)arrdata[group].dd->boxes[reflevel].size());
	assert (mglevel < (int)arrdata[group].dd->boxes[reflevel][component].size());
	const bbox<int,dim>& ext
	  = arrdata[group].dd->boxes[reflevel][component][mglevel].exterior;
	for (int d=0; d<dim; ++d) {
	  ((int*)arrdata[group].info.lsh)[d]
	    = (ext.shape() / ext.stride())[d];
	  ((int*)arrdata[group].info.lbnd)[d]
	    = (ext.lower() / ext.stride())[d];
	  ((int*)arrdata[group].info.ubnd)[d]
	    = (ext.upper() / ext.stride())[d];
	  assert (arrdata[group].info.lsh[d]>=0
		  && arrdata[group].info.lsh[d]<=arrdata[group].info.gsh[d]);
	  assert (arrdata[group].info.lbnd[d]>=0
		  && arrdata[group].info.ubnd[d]<arrdata[group].info.gsh[d]);
	  assert (arrdata[group].info.lbnd[d]<=arrdata[group].info.ubnd[d]+1);
	  // Do allow outer boundaries on the finer grids
	  ((int*)arrdata[group].info.bbox)[2*d  ]
	    = arrdata[group].info.lbnd[d] == 0;
	  ((int*)arrdata[group].info.bbox)[2*d+1]
	    = arrdata[group].info.ubnd[d] == arrdata[group].info.gsh[d]-1;
	}
      }
      
      // Set Cactus pointers to data
      for (int n=0; n<CCTK_NumVars(); ++n) {
	
	const int group = CCTK_GroupIndexFromVarI(n);
	assert (group>=0);
	const int var = n - CCTK_FirstVarIndexI(group);
	assert (var>=0);
	const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	assert (num_tl>0);
	
	for (int tl=0; tl<num_tl; ++tl) {
	  
	  if (CCTK_QueryGroupStorageI(cgh, group)) {
	    // Group has storage
	    
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    assert (arrdata[group].data[var]);
	    cgh->data[n][tl]
	      = ((*arrdata[group].data[var])
		 (-tl, reflevel, component, mglevel)->storage());
	    if (hh->is_local(reflevel,component)) {
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
