#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Checksum.cc,v 1.8 2002/10/24 10:39:38 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Checksum_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void CalculateChecksums (cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("%*sCalculateChecksums", 2*reflevel, "");
    
    checksums.resize(CCTK_NumVars());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	
	const int n = CCTK_FirstVarIndexI(group) + var;
	const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	assert (num_tl>0);
	const int min_tl = mintl(where, num_tl);
	const int max_tl = maxtl(where, num_tl);
	
	checksums[n].resize(maxreflevels);
	checksums[n][reflevel].resize(num_tl);
	for (int tl=min_tl; tl<=max_tl; ++tl) {
	  checksums[n][reflevel][tl].resize(hh->components(reflevel));
	  BEGIN_COMPONENT_LOOP(cgh) {
	    checksums[n][reflevel][tl][component].valid = false;
	  } END_COMPONENT_LOOP(cgh);
	}
      }
    }
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
	  assert (sz>0);
	  const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	  assert (num_tl>0);
	  const int min_tl = mintl(where, num_tl);
	  const int max_tl = maxtl(where, num_tl);
	  
	  for (int tl=min_tl; tl<=max_tl; ++tl) {
	    BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
	      const int gpdim = arrdata[group].info.dim;
	      int np = 1;
	      for (int d=0; d<gpdim; ++d) {
		np *= *CCTK_ArrayGroupSizeI(cgh, d, group);
	      }
	      const void* data = cgh->data[n][tl];
	      int chk = 0;
	      for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
		chk += ((const int*)data)[i];
	      }
	      checksums[n][reflevel][tl][component].sum = chk;
	      checksums[n][reflevel][tl][component].valid = true;
	    } END_LOCAL_COMPONENT_LOOP(cgh);
	  } // for tl
	} // for var
      }	// if has storage
    } // for group
  }
  
  
  
  void CheckChecksums (cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("%*sCheckChecksums", 2*reflevel, "");
    
    assert ((int)checksums.size()==CCTK_NumVars());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n));
	  assert (sz>0);
	  const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	  assert (num_tl>0);
	  const int min_tl = mintl(where, num_tl);
	  const int max_tl = maxtl(where, num_tl);
	  
	  assert ((int)checksums[n].size()==maxreflevels);
	  assert ((int)checksums[n][reflevel].size()==num_tl);
	  
	  for (int tl=min_tl; tl<=max_tl; ++tl) {
	    
	    bool unexpected_change = false;
	    
	    assert ((int)checksums[n][reflevel][tl].size()
		    == hh->components(reflevel));
	    BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
	      if (checksums[n][reflevel][tl][component].valid) {
		const int gpdim = arrdata[group].info.dim;
		int np = 1;
		for (int d=0; d<gpdim; ++d) {
		  np *= *CCTK_ArrayGroupSizeI(cgh, d, group);
		}
		const void* data = cgh->data[n][tl];
		int chk = 0;
		for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
		  chk += ((const int*)data)[i];
		}
		unexpected_change
		  = (unexpected_change
		     || chk != checksums[n][reflevel][tl][component].sum);
	      }
	    } END_LOCAL_COMPONENT_LOOP(cgh);
	    
	    if (unexpected_change) {
	      char* fullname = CCTK_FullName(n);
	      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
			  "Timelevel %d of the variable \"%s\" has changed unexpectedly.",
			  tl, fullname);
	      free (fullname);
	    }
	    
	  } // for tl
	  
	} // for var
      }	// if has storage
    } // for group
  }
  
} // namespace Carpet
