#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Comm.cc,v 1.12 2002/10/12 13:02:52 schnetter Exp $";

CCTK_FILEVERSION(Carpet_Comm_cc)



namespace Carpet {
  
  using namespace std;
  
  
  
  int SyncGroup (cGH* cgh, const char* groupname)
  {
    if (hh->local_components(reflevel) > 1) assert (component == -1);
    
    Checkpoint ("%*sSyncGroup %s", 2*reflevel, "", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot synchronise group \"%s\" because it has no storage",
		  groupname);
      return -1;
    }
    
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tl = 0;
    
    assert (group<(int)arrdata.size());
    for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
	if (reflevel>0) {
	  for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	    // use the current time here (which may be modified by the
	    // user)
 	    const CCTK_REAL time = cgh->cctk_time / base_delta_time;
	    if (false) {
	      const CCTK_REAL time1 = tt->time (tl, reflevel, mglevel);
	      const CCTK_REAL time2 = cgh->cctk_time / base_delta_time;
	      assert (fabs(time1 - time2) < 1e-10);
	    }
	    
	    arrdata[group].data[var]->ref_bnd_prolongate
	      (tl, reflevel, c, mglevel, time);
	  }
	}
      }
      for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	arrdata[group].data[var]->sync (tl, reflevel, c, mglevel);
      }
    }
    
    return 0;
  }
  
  
  
  int EnableGroupComm (cGH* cgh, const char* groupname)
  {
    // Communication is always enabled
    return 0;
  }
  
  int DisableGroupComm (cGH* cgh, const char* groupname)
  {
    // Communication is always enabled
    return -1;
  }
  
} // namespace Carpet
