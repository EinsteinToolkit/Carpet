#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/gf.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Storage.cc,v 1.15 2003/05/02 14:20:26 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Storage_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  int EnableGroupStorage (cGH* cgh, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("%*sEnableGroupStorage %s", 2*reflevel, "", groupname);
    
    // TODO: Enabling storage for one refinement level has to enable
    // it for all other refinement levels as well.  Disabling must
    // wait until all refinement levels have been disabled.
    
    // TODO: Invent a mode "reflevel==-1" that is global, i.e. has
    // effect for all refinement levels.  This mode is used during
    // INITIAL, and en-/disabling storage in it is also global.
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was enabled previously
      return 1;
    }
    
    // Check whether this group has transfer operators
    if (! arrdata[group].do_transfer) {
      const int grouptype = CCTK_GroupTypeI(group);
      if (grouptype == CCTK_GF) {
        const int var = CCTK_FirstVarIndexI(group);
        const int vartype = CCTK_VarTypeI(var);
        const char * vartypename = CCTK_VarTypeName(vartype);
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "(Allocating storage for Cactus group \"%s\".)  Note: This group (which has the variable type %s) will be neither prolongated nor restricted.",
                    groupname, vartypename);
      }
    }
    
    // There is a difference between the Cactus time levels and the
    // Carpet time levels.  If there are n time levels, then the
    // Cactus time levels are numbered 0 ... n-1, with the current
    // time level being 0.  In Carpet, the time levels are numbered
    // -(n-1) ... 0, where the current time level is also 0.
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tmin = 1 - num_tl;
    const int tmax = 0;
    const int my_prolongation_order_time
      = num_tl==1 ? 0 : prolongation_order_time;
    
    if (max_refinement_levels > 1) {
      if (num_tl <= my_prolongation_order_time) {
	CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "There are not enough time levels for the desired temporal prolongation order in the grid function group \"%s\".  With Carpet::prolongation_order_time=%d, you need at least %d time levels.",
		    CCTK_GroupName(group),
		    prolongation_order_time, prolongation_order_time+1);
      }
    }
    
    assert (arrdata[group].data.size()==0
	    || arrdata[group].data[0] == 0);
    for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
      const int n = n0 + var;
      switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)							\
      case N:								\
	arrdata[group].data[var] = new gf<T,dim>			\
	  (CCTK_VarName(n), *arrdata[group].tt, *arrdata[group].dd,	\
	   tmin, tmax, my_prolongation_order_time);			\
	break;
#include "typecase"
#undef TYPECASE
      default:
	UnsupportedVarType(n);
      } // switch
    } // for
    
    // Reinitialise Cactus variables
    if (component!=-1) set_component (cgh, component);
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
    
    assert (arrdata[group].data.size());
    assert (arrdata[group].data[0]);
    for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
      const int n = n0 + var;
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
    } // for
    
    // Reinitialise Cactus variables
    set_component (cgh, component);
    
    // storage was not disabled previously
    return 1;
  }
  
  
  
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname)
  {
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    const int n = CCTK_FirstVarIndexI(group);
    assert (n>=0 && n<CCTK_NumVars());
    const int var = 0;
    
    assert (group<(int)arrdata.size());
    assert (var<(int)arrdata[group].data.size());
    return arrdata[group].data[var] != 0;
  }
  
  
  
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname)
  {
    static const int zero = 0;
    
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (component == -1) {
      // global routine
      return &zero;
    }
    
    const int gpdim = arrdata[group].info.dim;
    
    assert (dir>=0 && dir<gpdim);
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      
      const int var = CCTK_FirstVarIndexI(group);
      assert (var>=0 && var<CCTK_NumVars());
      
      assert (group<(int)arrdata.size());
      return &arrdata[group].info.lsh[dir];
      
    } else {
      
      // no storage
      return &zero;
      
    }
  }
  
  
  
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data)
  {
    assert (group>=0 && group<CCTK_NumGroups());
    *data = arrdata[group].info;
    return 0;
  }
  
} // namespace Carpet
