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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/helpers.cc,v 1.31 2003/05/12 13:48:05 schnetter Exp $";
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
    return -999;
  }
  
  int Abort (cGH* cgh, int retval)
  {
    MPI_Abort (dist::comm, retval);
    abort ();
    return -999;
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
    return -999;
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
    return -999;
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
    assert (rl==-1 || (rl>=0 && rl<maxreflevels && rl<maxreflevels));
    assert (mglevel == -1);
    assert (component == -1);
    
    // Change
    reflevel = rl;
    if (reflevel == -1) {
      // global mode
      reflevelfact = 0xdeadbeef;
    } else {
      // level mode or local mode
      reflevelfact = ipow(reffact, reflevel);
    }
    vect<int,dim>::ref(cgh->cctk_levfac) = reflevelfact;
  }
  
  
  
  void set_mglevel (cGH* cgh, const int ml)
  {
    // Check
    assert (ml==-1 || (ml>=0 && ml<mglevels));
    assert (reflevel>=0 && reflevel<hh->reflevels());
    assert (component == -1);
    
    // Change
    mglevel = ml;
    mglevelfact = ipow(mgfact, mglevel);
    cgh->cctk_convlevel = mglevel;
    
    // Set gsh
    if (mglevel == -1) {
      
      mglevelfact = 0xdeadbeef;
      cgh->cctk_convlevel = 0xdeadbeef;
      
      cgh->cctk_delta_time = 0xdeadbeef;
      for (int d=0; d<dim; ++d) {
        cgh->cctk_origin_space[d] = 0xdeadbeef;
      }
      
      vect<int,dim>::ref(cgh->cctk_gsh) = 0xdeadbeef;
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          vect<int,dim>::ref((int*)arrdata[group].info.gsh) = 0xdeadbeef;
        }
      }
      
    } else {
      
      mglevelfact = ipow(mgfact, mglevel);
      cgh->cctk_convlevel = mglevel;
      
      // TODO: set cctk_time here as well
      cgh->cctk_delta_time = base_delta_time / reflevelfact * mglevelfact;
      
      {
        const bbox<int,dim>& baseext = dd->bases[reflevel][mglevel].exterior;
        for (int d=0; d<dim; ++d) {
          cgh->cctk_origin_space[d] = base_origin_space[d] + cgh->cctk_delta_space[d] / cgh->cctk_levfac[d] * baseext.lower()[d] / baseext.stride()[d];
        }
      }
      
      const bbox<int,dim>& baseext = dd->bases[reflevel][mglevel].exterior;
      vect<int,dim>::ref(cgh->cctk_gsh)	= baseext.shape() / baseext.stride();
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          const bbox<int,dim>& baseext = arrdata[group].dd->bases[reflevel][mglevel].exterior;
          for (int d=0; d<dim; ++d) {
            ((int*)arrdata[group].info.gsh)[d] = (baseext.shape() / baseext.stride())[d];
          }
        }
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
      // Level mode -- no component is active
      
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
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          
          for (int d=0; d<arrdata[group].info.dim; ++d) {
            ((int*)arrdata[group].info.lsh)[d]      = 0xdeadbeef;
            ((int*)arrdata[group].info.bbox)[2*d  ] = 0xdeadbeef;
            ((int*)arrdata[group].info.bbox)[2*d+1] = 0xdeadbeef;
            ((int*)arrdata[group].info.lbnd)[d]     = 0xdeadbeef;
            ((int*)arrdata[group].info.ubnd)[d]     = 0xdeadbeef;
          }
          
          const int firstvar = CCTK_FirstVarIndexI (group);
          const int numvars = CCTK_NumVarsInGroupI (group);
          const int num_tl = CCTK_NumTimeLevelsFromVarI (firstvar);
          
          assert (group<(int)arrdata.size());
          for (int var=0; var<numvars; ++var) {
            assert (var<(int)arrdata[group].data.size());
            for (int tl=0; tl<num_tl; ++tl) {
              cgh->data[firstvar+var][tl] = 0;
            }
          }
          
        } // if grouptype
      } // for var
      
    } else {
      // Local mode
      
      assert (reflevel>=0 && reflevel < (int)dd->boxes.size());
      assert (component>=0 && component < (int)dd->boxes[reflevel].size());
      assert (mglevel>=0 && mglevel < (int)dd->boxes[reflevel][component].size());
      
      const bbox<int,dim>& baseext = dd->bases[reflevel][mglevel].exterior;
      const vect<vect<bool,2>,dim>& obnds = hh->outer_boundaries[reflevel][component];
      const bbox<int,dim>& ext = dd->boxes[reflevel][component][mglevel].exterior;
      
      for (int d=0; d<dim; ++d) {
        
	cgh->cctk_lsh[d] = (ext.shape() / ext.stride())[d];
	cgh->cctk_lbnd[d] = ((ext.lower() - baseext.lower()) / ext.stride())[d];
	cgh->cctk_ubnd[d] = ((ext.upper() - baseext.lower()) / ext.stride())[d];
	cgh->cctk_bbox[2*d  ] = obnds[d][0];
	cgh->cctk_bbox[2*d+1] = obnds[d][1];
	for (int stg=0; stg<CCTK_NSTAGGER; ++stg) {
	  // TODO: support staggering
	  cgh->cctk_lssh[CCTK_LSSH_IDX(stg,d)] = cgh->cctk_lsh[d];
	}
        
	assert (cgh->cctk_lsh[d]>=0 && cgh->cctk_lsh[d]<=cgh->cctk_gsh[d]);
 	assert (cgh->cctk_lbnd[d]>=0 && cgh->cctk_ubnd[d]<cgh->cctk_gsh[d]);
	assert (cgh->cctk_ubnd[d]-cgh->cctk_lbnd[d]+1 == cgh->cctk_lsh[d]);
	assert (cgh->cctk_lbnd[d]<=cgh->cctk_ubnd[d]+1);
        
      }
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_GroupTypeI(group) == CCTK_GF) {
          
          assert (reflevel>=0 && reflevel < (int)arrdata[group].dd->boxes.size());
          assert (component>=0 && component < (int)arrdata[group].dd->boxes[reflevel].size());
          assert (mglevel>=0 && mglevel < (int)arrdata[group].dd->boxes[reflevel][component].size());
          
          const bbox<int,dim>& baseext = arrdata[group].dd->bases[reflevel][mglevel].exterior;
          const vect<vect<bool,2>,dim>& obnds = arrdata[group].hh->outer_boundaries[reflevel][component];
          const bbox<int,dim>& ext = arrdata[group].dd->boxes[reflevel][component][mglevel].exterior;
          
          for (int d=0; d<arrdata[group].info.dim; ++d) {
            
            ((int*)arrdata[group].info.lsh)[d]  = (ext.shape() / ext.stride())[d];
            ((int*)arrdata[group].info.lbnd)[d] = ((ext.lower() - baseext.lower()) / ext.stride())[d];
            ((int*)arrdata[group].info.ubnd)[d] = ((ext.upper() - baseext.lower()) / ext.stride())[d];
            ((int*)arrdata[group].info.bbox)[2*d  ] = obnds[d][0];
            ((int*)arrdata[group].info.bbox)[2*d+1] = obnds[d][1];
            
            assert (arrdata[group].info.lsh[d]>=0 && arrdata[group].info.lsh[d]<=arrdata[group].info.gsh[d]);
            assert (arrdata[group].info.lbnd[d]>=0 && arrdata[group].info.ubnd[d]<arrdata[group].info.gsh[d]);
            assert (arrdata[group].info.ubnd[d]-arrdata[group].info.lbnd[d]+1 == arrdata[group].info.lsh[d]);
            assert (arrdata[group].info.lbnd[d]<=arrdata[group].info.ubnd[d]+1);
          }
          
          const int firstvar = CCTK_FirstVarIndexI (group);
          const int numvars = CCTK_NumVarsInGroupI (group);
          const int num_tl = CCTK_NumTimeLevelsFromVarI (firstvar);
          
          assert (hh->is_local(reflevel,component));
          
          assert (group<(int)arrdata.size());
          for (int var=0; var<numvars; ++var) {
            assert (var<(int)arrdata[group].data.size());
            for (int tl=0; tl<num_tl; ++tl) {
              cgh->data[firstvar+var][tl] = ((*arrdata[group].data[var]) (-tl, reflevel, component, mglevel)->storage());
            }
          }
          
        } // if grouptype
      } // for group
      
    } // if local mode
    
  }
  
  
  
  // This is a temporary measure to call a schedule group from a
  // global mode function.
  int CallScheduleGroup (cGH * const cgh, const char * const group)
  {
    assert (component == -1);
    CCTK_ScheduleTraverse (group, cgh, CallFunction);
    return 0;
  }
  
  
  
  // This is a temporary measure to call a local mode function from a
  // global mode function.  A more elegant way would be to reuse the
  // CallFunction stuff, or function aliasing.  Is there a way for the
  // user to get at the cFunctionData structure?
  int CallLocalFunction (cGH * const cgh,
                         void (* const function) (cGH * const cgh))
  {
    BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
      function (cgh);
    } END_LOCAL_COMPONENT_LOOP(cgh);
    return 0;
  }
  
} // namespace Carpet
