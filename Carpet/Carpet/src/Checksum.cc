#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Checksum.cc,v 1.13 2004/01/25 14:57:27 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Checksum_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void CalculateChecksums (const cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("CalculateChecksums");
    
    checksums.resize(maxreflevels);
    checksums[reflevel].resize(mglevels);
    checksums[reflevel][mglevel].resize(CCTK_NumGroups());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
        const int grouptype = CCTK_GroupTypeI(group);
        if (reflevel == 0 || grouptype == CCTK_GF) {
          checksums[reflevel][mglevel][group].resize(arrdata[group].size());
          BEGIN_MAP_LOOP(cgh, grouptype) {
            checksums[reflevel][mglevel][group][map].resize(arrdata[group][map].hh->components(reflevel));
            BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype) {
              const int nvars = CCTK_NumVarsInGroupI(group);
              checksums[reflevel][mglevel][group][map][component].resize(nvars);
              if (nvars > 0) {
                
                const int n0 = CCTK_FirstVarIndexI(group);
                const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
                assert (sz>0);
                
                ivect size(1);
                const int gpdim = groupdata[group].info.dim;
                for (int d=0; d<gpdim; ++d) {
                  size[d] = groupdata[group].info.lsh[d];
                }
                const int np = prod(size);
                
                const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
                assert (num_tl>0);
                const int min_tl = mintl(where, num_tl);
                const int max_tl = maxtl(where, num_tl);
                
                for (int var=0; var<nvars; ++var) {
                  checksums[reflevel][mglevel][group][map][component][var].resize(num_tl);
                  for (int tl=min_tl; tl<=max_tl; ++tl) {
                    
                    const int n = n0 + var;
                    const void* data = cgh->data[n][tl];
                    unsigned int chk = 0;
                    for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
                      chk += ((const unsigned int*)data)[i];
                    }
                    
                    checksums[reflevel][mglevel][group][map][component][var][tl].sum = chk;
                    checksums[reflevel][mglevel][group][map][component][var][tl].valid = true;
                    
                  } // for tl
                } // for var
              } // if group has vars
            } END_LOCAL_COMPONENT_LOOP;
          } END_MAP_LOOP;
        } // if grouptype fits
      } // if storage
    } // for group
    
  }
  
  
  
  void CheckChecksums (const cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("CheckChecksums");
    
    assert ((int)checksums.size()==maxreflevels);
    assert ((int)checksums[reflevel].size()==mglevels);
    assert ((int)checksums[reflevel][mglevel].size()==CCTK_NumGroups());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
        const int grouptype = CCTK_GroupTypeI(group);
        if (reflevel == 0 || grouptype == CCTK_GF) {
          assert (checksums[reflevel][mglevel][group].size()==arrdata[group].size());
          BEGIN_MAP_LOOP(cgh, grouptype) {
            assert ((int)checksums[reflevel][mglevel][group][map].size()==arrdata[group][map].hh->components(reflevel));
            BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype) {
              const int nvars = CCTK_NumVarsInGroupI(group);
              assert ((int)checksums[reflevel][mglevel][group][map][component].size()==nvars);
              if (nvars > 0) {
                
                const int n0 = CCTK_FirstVarIndexI(group);
                const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
                assert (sz>0);
                
                ivect size(1);
                const int gpdim = groupdata[group].info.dim;
                for (int d=0; d<gpdim; ++d) {
                  size[d] = groupdata[group].info.lsh[d];
                }
                const int np = prod(size);
                
                const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
                assert (num_tl>0);
                const int min_tl = mintl(where, num_tl);
                const int max_tl = maxtl(where, num_tl);
                
                for (int var=0; var<nvars; ++var) {
                  assert ((int)checksums[reflevel][mglevel][group][map][component][var].size()==num_tl);
                  for (int tl=min_tl; tl<=max_tl; ++tl) {
                    if (checksums[reflevel][mglevel][group][map][component][var][tl].valid) {
                      
                      const int n = n0 + var;
                      const void* data = cgh->data[n][tl];
                      unsigned int chk = 0;
                      for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
                        chk += ((const unsigned int*)data)[i];
                      }
                      const bool unexpected_change =
                        chk != checksums[reflevel][mglevel][group][map][component][var][tl].sum;
                      
                      if (unexpected_change) {
                        char* fullname = CCTK_FullName(n);
                        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                                    "Timelevel %d, component %d, refinement level %d of the variable \"%s\" has changed unexpectedly",
                                    tl, component, reflevel, fullname);
                        free (fullname);
                      }
                      
                    } // if valid
                  } // for tl
                } // for var
              } // if group has vars
            } END_LOCAL_COMPONENT_LOOP;
          } END_MAP_LOOP;
        } // if grouptype fits
      } // if storage
    } // for group
    
  }
  
} // namespace Carpet
