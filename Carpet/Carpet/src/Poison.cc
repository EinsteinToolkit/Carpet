#include <cassert>
#include <cstdlib>
#include <cstring>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // The parameter where specifies which time levels should be
  // poisoned.  what specifies what kind of grid variables should be
  // poisoned.
  void Poison (const cGH* cgh, const checktimes where, const int what)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! poison_new_timelevels) return;
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
        int const grouptype = CCTK_GroupTypeI (group);
        if (what == 0 or
            (what == CCTK_GF and grouptype == CCTK_GF) or
            (what == CCTK_ARRAY and (grouptype == CCTK_ARRAY or
                                     grouptype == CCTK_SCALAR)))
        {
          PoisonGroup (cgh, group, where);
        }
      } // if has storage
    } // for group
  }
  
  
  
  void PoisonGroup (const cGH* cgh, const int group, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (group>=0 and group<CCTK_NumGroups());
    
    if (! poison_new_timelevels) return;
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      char * const groupname = CCTK_GroupName(group);
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot poison group \"%s\" because it has no storage",
		  groupname);
      free (groupname);
      return;
    }
    
    const int nvar = CCTK_NumVarsInGroupI(group);
    if (nvar == 0) return;
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
    assert (sz>0);
    
    const int num_tl = CCTK_ActiveTimeLevelsVI(cgh, n0);
    assert (num_tl>0);
    const int min_tl = min_timelevel(where, num_tl);
    const int max_tl = max_timelevel(where, num_tl);
    
    if (min_tl <= max_tl) {
      
      {
        char * const groupname = CCTK_GroupName(group);
        Checkpoint ("PoisonGroup \"%s\"", groupname);
        free (groupname);
      }
      
      const int grouptype = CCTK_GroupTypeI(group);
      
      BEGIN_MAP_LOOP(cgh, grouptype) {
        BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype) {
          
          ivect size(1);
          const int gpdim = groupdata.at(group).info.dim;
          for (int d=0; d<gpdim; ++d) {
            size[d] = groupdata.at(group).info.lsh[d];
          }
          const int np = prod(size);
          
          for (int var=0; var<nvar; ++var) {
            const int n = n0 + var;
            for (int tl=min_tl; tl<=max_tl; ++tl) {
              memset (cgh->data[n][tl], poison_value, np*sz);
            } // for tl
          } // for var
          
        } END_LOCAL_COMPONENT_LOOP;
      } END_MAP_LOOP;
      
    } // if tl
  }
  
  
  
  void PoisonCheck (const cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! check_for_poison) return;
    
    Checkpoint ("PoisonCheck");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      const int nvar = CCTK_NumVarsInGroupI(group);
      if (nvar > 0 && CCTK_QueryGroupStorageI(cgh, group)) {
        
        const int grouptype = CCTK_GroupTypeI(group);
        const int n0 = CCTK_FirstVarIndexI(group);
        assert (n0>=0);
        const int tp = CCTK_VarTypeI(n0);
        const int gpdim = groupdata.at(group).info.dim;
        
        const int num_tl = CCTK_ActiveTimeLevelsVI(cgh, n0);
        assert (num_tl>0);
        const int min_tl = min_timelevel(where, num_tl);
        const int max_tl = max_timelevel(where, num_tl);
        
        BEGIN_MAP_LOOP(cgh, grouptype) {
          BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype) {
            
            ivect size(1);
            for (int d=0; d<gpdim; ++d) {
              size[d] = groupdata.at(group).info.lsh[d];
            }
            const int np = prod(size);
            
            for (int var=0; var<nvar; ++var) {
              const int n = n0 + var;
              
              for (int tl=min_tl; tl<=max_tl; ++tl) {
                
                const void* const data = cgh->data[n][tl];
                int numpoison=0;
                for (int k=0; k<size[2]; ++k) {
                  for (int j=0; j<size[1]; ++j) {
                    for (int i=0; i<size[0]; ++i) {
                      const int idx = i + size[0] * (j + size[1] * k);
                      bool poisoned=false;
                      switch (tp) {
#define TYPECASE(N,T)                                                   \
                      case N: {                                         \
                        T worm;                                         \
                        memset (&worm, poison_value, sizeof worm);      \
                        const T & val = ((const T*)data)[idx];          \
                        poisoned = memcmp (&worm, &val, sizeof worm) == 0; \
                        break;                                          \
                      }
#include "typecase"
#undef TYPECASE
                      default:
                        UnsupportedVarType(n);
                      }
                      if (poisoned) {
                        ++numpoison;
                        if (max_poison_locations==-1
                            or numpoison<=max_poison_locations) {
                          char* fullname = CCTK_FullName(n);
                          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                                      "At iteration %d: timelevel %d, component %d, map %d, refinement level %d of the variable \"%s\" contains poison at [%d,%d,%d]",
                                      cgh->cctk_iteration,
                                      tl, component, map, reflevel,
                                      fullname, i,j,k);
                          free (fullname);
                        }
                      } // if poisoned
                    } // for i
                  } // for j
                } // for k
                if (max_poison_locations!=-1 and numpoison>max_poison_locations) {
                  char* fullname = CCTK_FullName(n);
                  CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "At iteration %d: timelevel %d, component %d, map %d, refinement level %d of the variable \"%s\" contains poison at %d of %d locations; not all locations were printed",
                              cgh->cctk_iteration,
                              tl, component, map, reflevel,
                              fullname, numpoison, np);
                  free (fullname);
                } else if (numpoison>0) {
                  CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                              "Found poison at %d of %d locations",
                              numpoison, np);
                }
                
              } // for tl
            } // for var
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
        
      } // if has storage
    } // for group
  }
  
} // namespace Carpet
