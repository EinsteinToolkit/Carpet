#include <assert.h>
#include <stdlib.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/gf.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Storage.cc,v 1.1 2001/07/04 12:29:47 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  int EnableGroupStorage (cGH* cgh, const char* groupname)
  {
    Checkpoint ("%*sEnableGroupStorage %s", 2*reflevel, "", groupname);
    
    // TODO: Enabling storage for one refinement level has to enable
    // it for all other refinement levels as well.  Disabling must
    // wait until all refinement levels have been disabled.
    
    // TODO: Invent a mode "reflevel==-1" that is global, i. e. has
    // effect for all refinement levels.  This mode is used during
    // INITIAL, and en-/disabling storage in it is also global.
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was enabled previously
      return 1;
    }
    
    // There is a difference between the Cactus time levels and the
    // Carpet time levels.  If there are n time levels, then the
    // Cactus time levels are numbered 0 ... n-1, with the current
    // time level being n-1.  In Carpet, the time levels are numbered
    // -(n-1) ... 0, where the current time level is always 0.
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tmin = 1 - num_tl;
    const int tmax = 0;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      assert (scdata[group].data.size()==0
	      || scdata[group].data[0].size()==0
	      || scdata[group].data[0][0].size()==0
	      || scdata[group].data[0][0][0] == 0);
      for (int var=0; var<(int)scdata[group].data.size(); ++var) {
	for (int rl=0; rl<(int)scdata[group].data[var].size(); ++rl) {
	  for (int tl=0; tl<(int)scdata[group].data[var][rl].size(); ++tl) {
	    const int n = n0 + var;
	    switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	    case N:					\
	      scdata[group].data[var][rl][tl] = new T;	\
	      break;
#include "typecase"
#undef TYPECASE
	    default:
	      UnsupportedVarType(n);
	    } // switch
	  } // for
	} // for
      }	// for
      break;
      
    case CCTK_ARRAY:
      assert (arrdata[group].data.size()==0
	      || arrdata[group].data[0] == 0);
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	const int n = n0 + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)							\
	case N:								\
	  arrdata[group].data[var] = new gf<T,dim>			\
	    (CCTK_VarName(n), *arrdata[group].tt, *arrdata[group].dd,	\
	     tmin, tmax);						\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  UnsupportedVarType(n);
	} // switch
      }	// for
      break;
      
    case CCTK_GF:
      assert (gfdata[group].data.size()==0
	      || gfdata[group].data[0] == 0);
      for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	const int n = n0 + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)						\
	case N:							\
	  gfdata[group].data[var] = new gf<T,dim>		\
	    (CCTK_VarName(n), *tt, *(dh<dim>*)dd, tmin, tmax);	\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  UnsupportedVarType(n);
	} // switch
      }	// for
      break;
      
    default:
      abort();
    }
    
    // Reinitialise Cactus variables
    set_component (cgh, component);
    PoisonGroup (cgh, group, alltimes);
    
    // storage was not enabled previously
    return 0;
  }
  
  
  
  int DisableGroupStorage (cGH* cgh, const char* groupname)
  {
    Checkpoint ("%*sDisableGroupStorage %s", 2*reflevel, "", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was disabled previously
      return 0;
    }
    
    // XXX
    CCTK_WARN (1, "Cannot disable storage -- storage management is not yet consistent for FMR");
    return 1;
    
    const int n0 = CCTK_FirstVarIndexI(group);
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR:
      if (scdata[group].data.size()==0
	  || scdata[group].data[0].size()==0
	  || scdata[group].data[0][0].size()==0
	  || scdata[group].data[0][0][0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)scdata[group].data.size(); ++var) {
	const int n = n0 + var;
	for (int rl=0; rl<(int)scdata[group].data[var].size(); ++rl) {
	  for (int tl=0; tl<(int)scdata[group].data[var][rl].size(); ++tl) {
	    switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)						\
	    case N:						\
	      delete (T*)scdata[group].data[var][rl][tl];	\
	      break;
#include "typecase"
#undef TYPECASE
	    default:
	      UnsupportedVarType(n);
	    }
	    scdata[group].data[var][rl][tl] = 0;
	  }
	}
      }
      break;
      
    case CCTK_ARRAY:
      if (arrdata[group].data.size()==0 || arrdata[group].data[0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	case N:						\
	  delete (gf<T,dim>*)arrdata[group].data[var];	\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  UnsupportedVarType(n);
	} // switch
	arrdata[group].data[var] = 0;
      }	// for
      break;
      
    case CCTK_GF:
      if (gfdata[group].data.size()==0
	  || gfdata[group].data[0] == 0) {
	// group already has no storage
	break;
      }
      for (int var=0; var<(int)gfdata[group].data.size(); ++var) {
	const int n = CCTK_FirstVarIndexI(group) + var;
	switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)					\
	case N:						\
	  delete (gf<T,dim>*)gfdata[group].data[var];	\
	  break;
#include "typecase"
#undef TYPECASE
	default:
	  UnsupportedVarType(n);
	} // switch
	gfdata[group].data[var] = 0;
      }	// for
      break;
      
    default:
      abort();
    }
    
    // Reinitialise Cactus variables
    set_component (cgh, component);
    
    // storage was not disabled previously
    return 1;
  }
  
  
  
  int QueryGroupStorageB (cGH* cgh, int group, const char* groupname)
  {
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    const int n = CCTK_FirstVarIndexI(group);
    assert (n>=0 && n<CCTK_NumVars());
    const int var = 0;
    
    switch (CCTK_GroupTypeFromVarI(n)) {
      
    case CCTK_SCALAR: {
      assert (group<(int)scdata.size());
      assert (var<(int)scdata[group].data.size());
      assert (reflevel<(int)scdata[group].data[var].size());
      const int tl=0;
      assert (tl<(int)scdata[group].data[var][reflevel].size());
      return scdata[group].data[var][reflevel][tl] != 0;
    }
      
    case CCTK_ARRAY: {
      assert (group<(int)arrdata.size());
      assert (var<(int)arrdata[group].data.size());
      return arrdata[group].data[var] != 0;
    }
      
    case CCTK_GF: {
      assert (group<(int)gfdata.size());
      assert (var<(int)gfdata[group].data.size());
      return gfdata[group].data[var] != 0;
    }
      
    default:
      abort();
    }
  }
  
} // namespace Carpet
