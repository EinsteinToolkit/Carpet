#include <cassert>
#include <cstdlib>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // Checksum information
  struct ckdesc {
    bool valid;
    unsigned int sum;
  };
  
  // Helper class
  struct ckdesc4 {
    vector<vector<vector<vector<ckdesc> > > > a; // [m][c][var][tl]
  };
  
  // Checksum information
  vector<vector<vector<ckdesc4> > > checksums; // [rl][ml][group]
  
  
  
  // The parameter where specifies which time levels should be
  // poisoned.  what specifies what kind of grid variables should be
  // poisoned.
  void CalculateChecksums (const cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("CalculateChecksums");
    
    checksums.resize(maxreflevels);
    checksums.at(reflevel).resize(mglevels);
    checksums.at(reflevel).at(mglevel).resize(CCTK_NumGroups());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
        const int grouptype = CCTK_GroupTypeI(group);
        if (reflevel == 0 or grouptype == CCTK_GF) {
          checksums.at(reflevel).at(mglevel).at(group).a.resize(arrdata.at(group).size());
          BEGIN_MAP_LOOP(cgh, grouptype) {
            checksums.at(reflevel).at(mglevel).at(group).a.at(map).resize(arrdata.at(group).at(map).hh->components(reflevel));
            BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype) {
              const int nvars = CCTK_NumVarsInGroupI(group);
              checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).resize(nvars);
              if (nvars > 0) {
                
                const int n0 = CCTK_FirstVarIndexI(group);
                const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
                assert (sz>0);
                
                ivect size(1);
                const int gpdim = groupdata.at(group).info.dim;
                for (int d=0; d<gpdim; ++d) {
                  size[d] = groupdata.at(group).info.lsh[d];
                }
                const int np = prod(size);
                
                const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
                assert (num_tl>0);
                const int min_tl = min_timelevel(where, num_tl);
                const int max_tl = max_timelevel(where, num_tl);
                
                for (int var=0; var<nvars; ++var) {
                  checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).at(var).resize(num_tl);
                  for (int tl=min_tl; tl<=max_tl; ++tl) {
                    
                    const int n = n0 + var;
                    const void* data = cgh->data[n][tl];
                    unsigned int chk = 0;
                    for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
                      chk += ((const unsigned int*)data)[i];
                    }
                    
                    checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).at(var).at(tl).sum = chk;
                    checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).at(var).at(tl).valid = true;
                    
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
    assert ((int)checksums.at(reflevel).size()==mglevels);
    assert ((int)checksums.at(reflevel).at(mglevel).size()==CCTK_NumGroups());
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
        const int grouptype = CCTK_GroupTypeI(group);
        if (reflevel == 0 or grouptype == CCTK_GF) {
          assert (checksums.at(reflevel).at(mglevel).at(group).a.size()==arrdata.at(group).size());
          BEGIN_MAP_LOOP(cgh, grouptype) {
            assert ((int)checksums.at(reflevel).at(mglevel).at(group).a.at(map).size()==arrdata.at(group).at(map).hh->components(reflevel));
            BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype) {
              const int nvars = CCTK_NumVarsInGroupI(group);
              assert ((int)checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).size()==nvars);
              if (nvars > 0) {
                
                const int n0 = CCTK_FirstVarIndexI(group);
                const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
                assert (sz>0);
                
                ivect size(1);
                const int gpdim = groupdata.at(group).info.dim;
                for (int d=0; d<gpdim; ++d) {
                  size[d] = groupdata.at(group).info.lsh[d];
                }
                const int np = prod(size);
                
                const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
                assert (num_tl>0);
                const int min_tl = min_timelevel(where, num_tl);
                const int max_tl = max_timelevel(where, num_tl);
                
                for (int var=0; var<nvars; ++var) {
                  assert ((int)checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).at(var).size()==num_tl);
                  for (int tl=min_tl; tl<=max_tl; ++tl) {
                    if (checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).at(var).at(tl).valid) {
                      
                      const int n = n0 + var;
                      const void* data = cgh->data[n][tl];
                      unsigned int chk = 0;
                      for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
                        chk += ((const unsigned int*)data)[i];
                      }
                      const bool unexpected_change =
                        chk != checksums.at(reflevel).at(mglevel).at(group).a.at(map).at(component).at(var).at(tl).sum;
                      
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
