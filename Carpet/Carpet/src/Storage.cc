#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <defs.hh>
#include <dh.hh>
#include <gf.hh>
#include <operators.hh>
#include <typeprops.hh>

#include <carpet.hh>



namespace Carpet {
  
  using namespace std;
  
  
  
  static int
  GroupStorageCrease (const cGH* cctkGH, int n_groups, const int* groups,
                      const int* tls, int* status,
                      const bool inc);
  
  static void
  GroupStorageCheck (cGH const * const cctkGH,
                     int const group,
                     int const ml, int const rl);
  
  
  
  int
  GroupStorageCrease (const cGH* cctkGH, int n_groups, const int* groups,
                      const int* tls, int* status,
                      const bool inc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH);
    assert (n_groups >= 0);
    assert (groups);
    assert (tls);
    for (int n=0; n<n_groups; ++n) {
      if (groups[n] < 0 or groups[n] >= CCTK_NumGroups()) {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group index %d is illegal", groups[n]);
        return -1;
      }
      assert (groups[n] >= 0 and groups[n] < CCTK_NumGroups());
      // TODO: tls[n] can also be -1; in that case, all time
      // levels should be activated / deactivated
      assert (tls[n] >= 0);
    }
    
    bool const can_do = is_meta_mode() or is_global_mode() or is_level_mode();
    bool const all_ml = is_meta_mode();
    int const min_ml = all_ml ? 0        : mglevel;
    int const max_ml = all_ml ? mglevels : mglevel+1;
    
    int min_num_timelevels = INT_MAX;
    
    for (int n=0; n<n_groups; ++n) {
      int const group = groups[n];
      
      if(storage_verbose)
      {
        char * const groupname = CCTK_GroupName (group);
        assert (groupname);
        Checkpoint ("  %s: %screase to %d",
                    groupname, inc ? "in" : "de", tls[n]);
        free (groupname);
      }
      
      cGroup gp;
      check (not CCTK_GroupData (group, &gp));
      
      bool const all_rl = is_meta_mode() or is_global_mode();
      bool const is_array = gp.grouptype != CCTK_GF;
      int const min_rl = is_array ? 0 : all_rl ? 0         : reflevel;
      int const max_rl = is_array ? 1 : all_rl ? reflevels : reflevel+1;
      
      int const firstvarindex = CCTK_FirstVarIndexI (group);
      assert (gp.numvars == 0
              or (firstvarindex >= 0 and firstvarindex < CCTK_NumVars()));
      
      // Check an assumption
      if (not gp.vectorgroup) assert (gp.vectorlength == 1);
      
      // Allocate the time levels
      for (int ml=min_ml; ml<max_ml; ++ml) {
        for (int rl=min_rl; rl<max_rl; ++rl) {
          
          // Record previous number of allocated time levels
          if (status) {
            // Note: This remembers only the last level
            status[n] = groupdata.AT(group).activetimelevels.AT(ml).AT(rl);
          }
          
          // Only do something if the number of time levels actually
          // needs to be changed -- do nothing otherwise
          
          const bool do_increase
            =     inc and tls[n] > groupdata.AT(group).activetimelevels.AT(ml).AT(rl);
          const bool do_decrease
            = not inc and tls[n] < groupdata.AT(group).activetimelevels.AT(ml).AT(rl);
          if (do_increase or do_decrease) {
            
            if (not can_do) {
              char * const groupname = CCTK_GroupName (group);
              char const * const modestring
                = (is_meta_mode() ? "meta"
                   : is_global_mode() ? "global"
                   : is_level_mode() ? "level"
                   : is_singlemap_mode() ? "singlemap"
                   : is_local_mode() ? "local"
                   : NULL);
              CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "Cannot change storage for group \"%s\" in %s mode",
                          groupname, modestring);
              free (groupname);
            }
            assert (can_do);
            
            // No storage change in local mode
            // TODO: This this seems superfluous, given the test above
            if (gp.grouptype == CCTK_GF) {
              assert ((map == -1 or maps == 1)
                      and (component == -1
                           or vhh.AT(0)->local_components(reflevel) == 1));
            }
            
            // Set the new number of active time levels
            groupdata.AT(group).activetimelevels.AT(ml).AT(rl) = tls[n];
            
            for (int m=0; m<(int)arrdata.AT(group).size(); ++m) {
              for (int var=0; var<gp.numvars; ++var) {
#ifdef CCTK_HAVE_CONTIGUOUS_GROUPS
                bool const contiguous = gp.contiguous;
#else
                bool const contiguous = false;
#endif
                const int vectorindex
                  = (contiguous
                     ? var
                     : gp.vectorgroup ? var % gp.vectorlength : 0);
                const int vectorlength
                  = (contiguous
                     ? gp.numvars
                     : gp.vectorgroup ? gp.vectorlength : 1);
                assert (vectorindex>=0 and vectorindex<gp.numvars);
                assert (vectorlength>0 and vectorlength<=gp.numvars);
                ggf* const vectorleader
                  = (vectorindex>0
                     ? arrdata.AT(group).AT(m).data.AT(var - vectorindex)
                     : NULL);
                const int varindex = firstvarindex + var;
                // TODO: allocate these in SetupGH, and after recomposing
                if (not arrdata.AT(group).AT(m).data.AT(var)) {
                  switch (specific_cactus_type(gp.vartype)) {
#define TYPECASE(N,T)                                                   \
                    case N:                                             \
                      arrdata.AT(group).AT(m).data.AT(var) = new gf<T>  \
                      (varindex,                                        \
                       groupdata.AT(group).transport_operator,          \
                       *arrdata.AT(group).AT(m).tt,                     \
                       *arrdata.AT(group).AT(m).dd,                     \
                       prolongation_order_time,                         \
                       vectorlength, vectorindex, (gf<T>*)vectorleader); \
                    break;
#include "typecase.hh"
#undef TYPECASE
                  default:
                    UnsupportedVarType (varindex);
                  } // switch gp.vartype
                } // if not allocated
                
                arrdata.AT(group).AT(m).data.AT(var)->set_timelevels
                  (ml, rl, tls[n]);
                
                // Set the data pointers for grid arrays
                if (gp.grouptype != CCTK_GF) {
                  assert (rl==0 and m==0);
                  for (int tl=0; tl<gp.numtimelevels; ++tl) {
                    if (tl < groupdata.AT(group).activetimelevels.AT(ml).AT(rl)) {
                      ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
                      void *const ptr = ff->data_pointer(tl, 0, 0, 0)->storage();
                      cctkGH->data[varindex][tl] = ptr;
                    } else {
                      cctkGH->data[varindex][tl] = NULL;
                    }
                  }
                } // if grouptype != GF
              
              } // for var
            } // for m
            
          } // if really change the number of active time levels
          
          // Complain if there are not enough active time levels
          GroupStorageCheck (cctkGH, group, ml, rl);
          
          // Record (minimum of) current number of time levels
          min_num_timelevels =
            min(min_num_timelevels,
                groupdata.AT(group).activetimelevels.AT(ml).AT(rl));
          
        } // for rl
      } // for ml
      
    } // for n
    if (min_num_timelevels == INT_MAX) {
      min_num_timelevels = 0;
    }
    
    return do_allow_past_timelevels ? 
      min_num_timelevels : min(1,min_num_timelevels);
  }
  
  
  
  int
  GroupStorageIncrease (const cGH* cctkGH, int n_groups, const int* groups,
                        const int* tls, int* status)
  {
    DECLARE_CCTK_PARAMETERS

    if(storage_verbose) {
      Checkpoint ("GroupStorageIncrease");
    }
    return
      GroupStorageCrease (cctkGH, n_groups, groups, tls, status, true);
  }
  
  
  
  int
  GroupStorageDecrease (const cGH* cctkGH, int n_groups, const int* groups,
                        const int* tls, int* status)
  {
    DECLARE_CCTK_PARAMETERS

    if(storage_verbose) {
      Checkpoint ("GroupStorageDecrease");
    }
    return
      GroupStorageCrease (cctkGH, n_groups, groups, tls, status, false);
  }
  
  
  
  int
  EnableGroupStorage (const cGH* cctkGH, const char* groupname)
  {
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 and group<CCTK_NumGroups());
    const int tls = CCTK_MaxTimeLevelsGI(group);
    int status;
    GroupStorageIncrease (cctkGH, 1, &group, &tls, &status);
    // Return whether storage was allocated previously
    return status;
  }
  
  
  
  int
  DisableGroupStorage (const cGH* cctkGH, const char* groupname)
  {
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 and group<CCTK_NumGroups());
    const int tls = 0;
    int status;
    GroupStorageDecrease (cctkGH, 1, &group, &tls, &status);
    // Return whether storage was allocated previously
    return status;
  }
  
  
  
  int
  QueryGroupStorageB (const cGH* cctkGH, int group, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    if (group<0 or group>=CCTK_NumGroups()) {
      CCTK_WARN (1, "QueryGroupStorage: illegal group specified");
      return -1;
    }
    int const grouptype = CCTK_GroupTypeI (group);
    if (is_meta_mode() or is_global_mode()) {
      if (grouptype == CCTK_GF) return -2;
    }
    if (groupdata.size() == 0) return -3;
    int const rl = grouptype == CCTK_GF ? reflevel : 0;
    // Return whether storage is allocated
    CCTK_INT const deadbeef = get_deadbeef();
    assert (groupdata.AT(group).activetimelevels.AT(mglevel).AT(rl) != deadbeef);
    return groupdata.AT(group).activetimelevels.AT(mglevel).AT(rl) > 0;
  }
  
  
  
  const int*
  ArrayGroupSizeB (const cGH* cctkGH, int dir, int group, const char* groupname)
  {
    static const int zero = 0;
    static const int error = 0;
    
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 and group<CCTK_NumGroups());
    
    if (mglevel == -1) {
      return &error;            // meta mode
    }
    
    const int gptype = CCTK_GroupTypeI (group);
    if (gptype == CCTK_GF and map == -1) {
      return &error;            // global or level mode for a GF
    }
    
    const int gpdim = groupdata.AT(group).info.dim;
    assert (dir>=0 and dir<gpdim);
    
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {
      
      return &groupdata.AT(group).info.lsh[dir];
      
    } else {
      
      // no storage
      return &zero;
      
    }
  }
  
  
  
  int GroupDynamicData (const cGH* cctkGH, int group, cGroupDynamicData* data)
  {
    // Return values:
    //  0 for success
    // -1 if given pointer to data structure is NULL
    // -3 if given GH pointer is invalid
    // (-77 if group has zero variables)
    // -78 if group does not exist
    if (not cctkGH) return -3;
    if (not (group>=0 and group<CCTK_NumGroups())) return -78;
    if (not data) return -1;
    assert (cctkGH);
    assert (group>=0 and group<CCTK_NumGroups());
    assert (data);
    *data = groupdata.AT(group).info;
    return 0;
  }
  
  
  
  void
  GroupStorageCheck (cGH const * const cctkGH,
                     int const group,
                     int const ml, int const rl)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (not do_warn_about_storage) return;
    if (max_refinement_levels == 1) return;
    
    cGroup gp;
    check (not CCTK_GroupData (group, & gp));
    
    if (gp.grouptype == CCTK_GF) {
      if (groupdata.AT(group).transport_operator != op_none and
          groupdata.AT(group).transport_operator != op_sync and
          groupdata.AT(group).transport_operator != op_restrict and
          groupdata.AT(group).transport_operator != op_copy)
      {
        if (groupdata.AT(group).activetimelevels.AT(ml).AT(rl) != 0 and
            (groupdata.AT(group).activetimelevels.AT(ml).AT(rl) <
             prolongation_order_time+1))
        {
          static vector<bool> didwarn;
          int const numgroups = CCTK_NumGroups();
          if ((int)didwarn.size() < numgroups) {
            didwarn.resize (numgroups, false);
          }
          if (not didwarn.AT(group)) {
            // Warn only once per group
            didwarn.AT(group) = true;
            char * const groupname = CCTK_GroupName (group);
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "There are not enough time levels for the desired temporal prolongation order in the grid function group \"%s\".  With Carpet::prolongation_order_time=%d, you need at least %d time levels.",
                        groupname,
                        (int)prolongation_order_time,
                        (int)(prolongation_order_time+1));
            free (groupname);
          }
        }
      }
    }
  }
  
  void
  GroupsStorageCheck (cGH const * const cctkGH)
  {
    for (int group = 0; group < CCTK_NumGroups(); ++ group) {
      for (int ml = 0; ml < mglevels; ++ ml) {
        for (int rl = 0; rl < reflevels; ++ rl) {
          GroupStorageCheck (cctkGH, group, ml, rl);
        }
      }
    }
  }

} // namespace Carpet
