#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/dist.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/helpers.cc,v 1.27 2002/10/24 10:39:39 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_helpers_cc);
}



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
    // if there is no cgh yet, assume nothing has been initialised
    // yet, and don't use dist::comm
    int rank;
    MPI_Comm_rank (cgh ? dist::comm : MPI_COMM_WORLD, &rank);
    return rank;
  }
  
  int nProcs (const cGH* cgh)
  {
    // if there is no cgh yet, assume nothing has been initialised
    // yet, and don't use dist::comm
    int size;
    MPI_Comm_size (cgh ? dist::comm : MPI_COMM_WORLD, &size);
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
      return num_tl>1 ? 0 : 1;
    case allbutlasttime:
      // do include current time if there is only one time level
      return 0;
    case allbutcurrenttime:
      return 1;
    case alltimes:
      return 0;
    default:
      assert (0);
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
      assert (0);
    }
  }
  
  
  
  void Waypoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (verbose || veryverbose) {
      va_list args;
      char msg[1000];
      va_start (args, fmt);
      vsnprintf (msg, sizeof msg, fmt, args);
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
      vsnprintf (msg, sizeof msg, fmt, args);
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
    assert (rl>=0 && rl<maxreflevels && rl<hh->reflevels());
    assert (mglevel == -1);
    assert (component == -1);
    
    // Change
    reflevel = rl;
    reflevelfact = ipow(reffact, reflevel);
    vect<int,dim>::ref(cgh->cctk_levfac) = reflevelfact;
  }
  
  
  
  void set_mglevel (cGH* cgh, const int ml)
  {
    // Check
    assert (ml==-1 || (ml>=0 && ml<mglevels));
    assert (component == -1);
    
    // Change
    mglevel = ml;
    mglevelfact = ipow(mgfact, mglevel);
    cgh->cctk_convlevel = mglevel;
    
    // Set gsh
    if (mglevel == -1) {
      
      cgh->cctk_delta_time = 0xdeadbeef;
      vect<int,dim>::ref(cgh->cctk_gsh) = 0xdeadbeef;
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	vect<int,dim>::ref((int*)arrdata[group].info.gsh) = 0xdeadbeef;
      }
      
    } else {
      
      const bbox<int,dim>& base = hh->baseextent;
      cgh->cctk_delta_time = base_delta_time / reflevelfact * mglevelfact;
      assert (all(base.shape() % base.stride() == 0));
      assert (all((base.shape() / base.stride()) % mglevelfact == 0));
      vect<int,dim>::ref(cgh->cctk_gsh)
	= (base.shape() / base.stride() - 1) / mglevelfact * reflevelfact + 1;
      for (int group=0; group<CCTK_NumGroups(); ++group) {
	const bbox<int,dim>& base = arrdata[group].hh->baseextent;
	vect<int,dim>::ref((int*)arrdata[group].info.gsh)
	  = ((base.shape() / base.stride() - 1)
	     / mglevelfact * reflevelfact + 1);
      }
      
    } // if mglevel != -1
  }
  
  
  
  void set_component (cGH* cgh, const int c)
  {
    assert (reflevel>=0 && reflevel<hh->reflevels());
    assert (c==-1 || (c>=0 && c<hh->components(reflevel)));
    component = c;
    assert (component==-1
	    || (mglevel>=0 && mglevel<hh->mglevels(reflevel,component)));
    
    
    
    // Set Cactus parameters
    if (component == -1) {
      // Global mode -- no component is active
      
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
      
    } else {
      // Local mode
      
      assert (reflevel>=0 && reflevel < (int)dd->boxes.size());
      assert (component>=0 && component < (int)dd->boxes[reflevel].size());
      assert (mglevel>=0
	      && mglevel < (int)dd->boxes[reflevel][component].size());
      const bbox<int,dim>& ext
	= dd->boxes[reflevel][component][mglevel].exterior;
      for (int d=0; d<dim; ++d) {
	cgh->cctk_lsh[d] = (ext.shape() / ext.stride())[d];
	cgh->cctk_lbnd[d] = (ext.lower() / ext.stride())[d];
	cgh->cctk_ubnd[d] = (ext.upper() / ext.stride())[d];
	assert (cgh->cctk_lsh[d]>=0 && cgh->cctk_lsh[d]<=cgh->cctk_gsh[d]);
// 	assert (cgh->cctk_lbnd[d]>=0 && cgh->cctk_ubnd[d]<cgh->cctk_gsh[d]);
	assert (cgh->cctk_ubnd[d]-cgh->cctk_lbnd[d]+1 == cgh->cctk_lsh[d]);
	assert (cgh->cctk_lbnd[d]<=cgh->cctk_ubnd[d]+1);
	cgh->cctk_bbox[2*d  ] = hh->outer_boundaries[reflevel][component][d][0];
	cgh->cctk_bbox[2*d+1] = hh->outer_boundaries[reflevel][component][d][1];
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  // TODO: support staggering
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
	}
      }
      
    } // if local mode
    
    
    
    // Set Cactus parameters for all groups
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      if (mglevel == -1
	  || (CCTK_GroupTypeI(group) == CCTK_GF
	      && component == -1)) {
	// Global mode for a grid function: not active
	
	for (int d=0; d<dim; ++d) {
	  ((int*)arrdata[group].info.lsh)[d]      = 0xdeadbeef;
	  ((int*)arrdata[group].info.bbox)[2*d  ] = 0xdeadbeef;
	  ((int*)arrdata[group].info.bbox)[2*d+1] = 0xdeadbeef;
	  ((int*)arrdata[group].info.lbnd)[d]     = 0xdeadbeef;
	  ((int*)arrdata[group].info.ubnd)[d]     = 0xdeadbeef;
	}
	
      } else {
	// Local mode, or array or scalar: active
	
	int rl, c;
	if (CCTK_GroupTypeI(group) == CCTK_GF) {
	  rl = reflevel;
	  c = component;
	} else {
	  rl = 0;
	  c = CCTK_MyProc(cgh);
	}
	assert (rl>=0 && rl < (int)arrdata[group].dd->boxes.size());
	assert (c>=0 && c < (int)arrdata[group].dd->boxes[rl].size());
	assert (mglevel>=0 && mglevel < (int)dd->boxes[rl][c].size());
	const bbox<int,dim>& bext = arrdata[group].hh->baseextent;
	const bbox<int,dim>& iext = arrdata[group].hh->extents[rl][c][mglevel];
	const bbox<int,dim>& ext
	  = arrdata[group].dd->boxes[rl][c][mglevel].exterior;
	for (int d=0; d<dim; ++d) {
	  ((int*)arrdata[group].info.lsh)[d]
	    = (ext.shape() / ext.stride())[d];
	  ((int*)arrdata[group].info.lbnd)[d]
	    = (ext.lower() / ext.stride())[d];
	  ((int*)arrdata[group].info.ubnd)[d]
	    = (ext.upper() / ext.stride())[d];
	  assert (arrdata[group].info.lsh[d]>=0
		  && arrdata[group].info.lsh[d]<=arrdata[group].info.gsh[d]);
// 	  assert (arrdata[group].info.lbnd[d]>=0
// 		  && arrdata[group].info.ubnd[d]<arrdata[group].info.gsh[d]);
	  assert (arrdata[group].info.ubnd[d]-arrdata[group].info.lbnd[d]+1
		  == arrdata[group].info.lsh[d]);
	  assert (arrdata[group].info.lbnd[d]<=arrdata[group].info.ubnd[d]+1);
	  // Do not allow outer boundaries on the finer grids
	  if (rl==0) {
	    assert (iext.lower()[d] >= bext.lower()[d]);
	    assert (iext.upper()[d] <= bext.upper()[d]);
	    ((int*)arrdata[group].info.bbox)[2*d  ]
	      = iext.lower()[d] == bext.lower()[d];
	    ((int*)arrdata[group].info.bbox)[2*d+1]
	      = iext.upper()[d] == bext.upper()[d];
	  } else {
	    ((int*)arrdata[group].info.bbox)[2*d  ] = 0;
	    ((int*)arrdata[group].info.bbox)[2*d+1] = 0;
	  }
	} // for d
	
      } // if local mode
      
    } // for group
    
    
    
    // Set Cactus pointers to data
    for (int n=0; n<CCTK_NumVars(); ++n) {
      
      const int group = CCTK_GroupIndexFromVarI(n);
      assert (group>=0);
      const int var = n - CCTK_FirstVarIndexI(group);
      assert (var>=0);
      const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
      assert (num_tl>0);
      
      for (int tl=0; tl<num_tl; ++tl) {
	
	if (mglevel != -1
	    && CCTK_QueryGroupStorageI(cgh, group)) {
	  // Group has storage
	  
	  if (CCTK_GroupTypeI(group) == CCTK_GF) {
	    // It is a grid function
	    
	    if (component == -1) {
	      // Global mode
	      
	      // Grid functions cannot be accessed in global mode
	      cgh->data[n][tl] = 0;
	      
	    } else {
	      // local mode
	      
	      assert (reflevel>=0 && reflevel < (int)dd->boxes.size());
	      assert (component>=0
		      && component < (int)dd->boxes[reflevel].size());
	      assert (mglevel>=0
		      && mglevel < (int)dd->boxes[reflevel][component].size());
	      
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
	      
	    } // if global mode
	    
	  } else {
	    // Scalars and arrays can always be accessed
	    
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    assert (arrdata[group].data[var]);
	    const int rl = 0;
	    const int c = CCTK_MyProc(cgh);
	    assert (rl>=0 && rl<(int)arrdata[group].dd->boxes.size());
	    assert (c>=0 && c<(int)arrdata[group].dd->boxes[rl].size());
	    assert (mglevel>=0 && mglevel < (int)dd->boxes[rl][c].size());
	    assert (hh->is_local(reflevel,c));
	    cgh->data[n][tl]
	      = ((*arrdata[group].data[var]) (-tl, rl, c, mglevel)->storage());
	    assert (cgh->data[n][tl]);
	    
	  } // if scalar or array
	  
	} else {
	  // Group has no storage
	  
	  cgh->data[n][tl] = 0;
	  
	} // if has no storage
	
      } // for tl
    } // for n
    
  }
  
} // namespace Carpet
