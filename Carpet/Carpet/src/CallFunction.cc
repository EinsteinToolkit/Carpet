#include <algorithm>
#include <cassert>
#include <cstdlib>

#include "cctk.h"
#include "cctki_GHExtensions.h"

#include "gh.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  static void SyncGroupsInScheduleBlock( cFunctionData* attribute, cGH* cgh );
  /// Traverse one function on all components of one refinement level
  /// of one multigrid level.
  int CallFunction (void* function, ///< the function to call
                    cFunctionData* attribute, ///< attributes of the function
                    void* data) ///< ???
  {
//     Checkpoint ("Starting CallFunction...");
    
    cGH* cgh = (cGH*)data;
    
    assert (!! attribute->meta
            + !! attribute->global
            + !! attribute->level
            + !! attribute->singlemap
            + !! attribute->local <= 1);
    
    assert (!! attribute->loop_global
            + !! attribute->loop_level
            + !! attribute->loop_singlemap
            + !! attribute->loop_local <= 1);
    
    if (attribute->meta or is_meta_mode()) {
      // Convtest operation
      
      if (do_meta_mode) {
        if (attribute->loop_local) {
          BEGIN_META_MODE(cgh) {
            BEGIN_MGLEVEL_LOOP(cgh) {
              BEGIN_REFLEVEL_LOOP(cgh) {
                BEGIN_MAP_LOOP(cgh, CCTK_GF) {
                  BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
                    Checkpoint ("Meta time local mode call at %s to %s::%s",
                                attribute->where, attribute->thorn, attribute->routine);
                    const int res = CCTK_CallFunction (function, attribute, data);
                    assert (res==0);
                  } END_LOCAL_COMPONENT_LOOP;
                } END_MAP_LOOP;
                SyncGroupsInScheduleBlock( attribute, cgh );
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_singlemap) {
          BEGIN_META_MODE(cgh) {
            BEGIN_MGLEVEL_LOOP(cgh) {
              BEGIN_REFLEVEL_LOOP(cgh) {
                BEGIN_MAP_LOOP(cgh, CCTK_GF) {
                  Checkpoint ("Meta time singlemap mode call at %s to %s::%s",
                              attribute->where, attribute->thorn, attribute->routine);
                  const int res = CCTK_CallFunction (function, attribute, data);
                  assert (res==0);
                } END_MAP_LOOP;
                SyncGroupsInScheduleBlock( attribute, cgh );
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_level) {
          BEGIN_META_MODE(cgh) {
            BEGIN_MGLEVEL_LOOP(cgh) {
              BEGIN_REFLEVEL_LOOP(cgh) {
                Checkpoint ("Meta time level mode call at %s to %s::%s",
                            attribute->where, attribute->thorn, attribute->routine);
                const int res = CCTK_CallFunction (function, attribute, data);
                assert (res==0);
                SyncGroupsInScheduleBlock( attribute, cgh );
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_global) {
          BEGIN_META_MODE(cgh) {
            BEGIN_MGLEVEL_LOOP(cgh) {
              Checkpoint ("Meta time global mode call at %s to %s::%s",
                          attribute->where, attribute->thorn, attribute->routine);
              const int res = CCTK_CallFunction (function, attribute, data);
              assert (res==0);
              SyncGroupsInScheduleBlock( attribute, cgh );
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else {
          BEGIN_META_MODE(cgh) {
            Checkpoint ("Meta mode call at %s to %s::%s",
                        attribute->where, attribute->thorn, attribute->routine);
            const int res = CCTK_CallFunction (function, attribute, data);
            assert (res==0);
            SyncGroupsInScheduleBlock( attribute, cgh );
          } END_META_MODE;
        }
      }
      
    } else if (attribute->global or is_global_mode()) {
      // Global operation: call once
      
      assert (! attribute->loop_meta);
      
      if (do_global_mode) {
        if (attribute->loop_local) {
          BEGIN_GLOBAL_MODE(cgh) {
            BEGIN_REFLEVEL_LOOP(cgh) {
              BEGIN_MAP_LOOP(cgh, CCTK_GF) {
                BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
                  Checkpoint ("Global time local mode call at %s to %s::%s",
                              attribute->where, attribute->thorn, attribute->routine);
                  const int res = CCTK_CallFunction (function, attribute, data);
                  assert (res==0);
                } END_LOCAL_COMPONENT_LOOP;
              } END_MAP_LOOP;
              SyncGroupsInScheduleBlock( attribute, cgh );
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else if (attribute->loop_singlemap) {
          BEGIN_GLOBAL_MODE(cgh) {
            BEGIN_REFLEVEL_LOOP(cgh) {
              BEGIN_MAP_LOOP(cgh, CCTK_GF) {
                Checkpoint ("Global time singlemap mode call at %s to %s::%s",
                            attribute->where, attribute->thorn, attribute->routine);
                const int res = CCTK_CallFunction (function, attribute, data);
                assert (res==0);
              } END_MAP_LOOP;
              SyncGroupsInScheduleBlock( attribute, cgh );
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else if (attribute->loop_level) {
          BEGIN_GLOBAL_MODE(cgh) {
            BEGIN_REFLEVEL_LOOP(cgh) {
              Checkpoint ("Global time level mode call at %s to %s::%s",
                          attribute->where, attribute->thorn, attribute->routine);
              const int res = CCTK_CallFunction (function, attribute, data);
              assert (res==0);
              SyncGroupsInScheduleBlock( attribute, cgh );
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else {
          BEGIN_GLOBAL_MODE(cgh) {
            Checkpoint ("Global mode call at %s to %s::%s",
                        attribute->where, attribute->thorn, attribute->routine);
            const int res = CCTK_CallFunction (function, attribute, data);
            assert (res==0);
            SyncGroupsInScheduleBlock( attribute, cgh );
          } END_GLOBAL_MODE;
        }
      }
      
    } else if (attribute->level) {
      // Level operation: call once per refinement level
      
      assert (! attribute->loop_meta);
      assert (! attribute->loop_global);
      
      if (attribute->loop_local) {
        BEGIN_MAP_LOOP(cgh, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
            Checkpoint ("Level time local mode call at %s to %s::%s",
                        attribute->where, attribute->thorn, attribute->routine);
            const int res = CCTK_CallFunction (function, attribute, data);
            assert (res==0);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } else if (attribute->loop_singlemap) {
        BEGIN_MAP_LOOP(cgh, CCTK_GF) {
          Checkpoint ("Level time singlemap mode call at %s to %s::%s",
                      attribute->where, attribute->thorn, attribute->routine);
          const int res = CCTK_CallFunction (function, attribute, data);
          assert (res==0);
        } END_MAP_LOOP;
      } else {
        Checkpoint ("Level mode call at %s to %s::%s",
                    attribute->where, attribute->thorn, attribute->routine);
        const int res = CCTK_CallFunction (function, attribute, data);
        assert (res==0);
      }
      SyncGroupsInScheduleBlock( attribute, cgh );
      
    } else if (attribute->singlemap) {
      // Single map operation: call once per refinement level and map
      
      assert (! attribute->loop_meta);
      assert (! attribute->loop_global);
      assert (! attribute->loop_level);
      
      if (attribute->loop_local) {
        BEGIN_MAP_LOOP(cgh, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
            Checkpoint ("Singlemap time local mode call at %s to %s::%s",
                        attribute->where, attribute->thorn, attribute->routine);
            const int res = CCTK_CallFunction (function, attribute, data);
            assert (res==0);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } else {
        BEGIN_MAP_LOOP(cgh, CCTK_GF) {
          Checkpoint ("Singlemap mode call at %s to %s::%s",
                      attribute->where, attribute->thorn, attribute->routine);
          const int res = CCTK_CallFunction (function, attribute, data);
          assert (res==0);
        } END_MAP_LOOP;
      }
      SyncGroupsInScheduleBlock( attribute, cgh );
      
    } else {
      // Local operation: call once per component
      
      assert (! attribute->loop_meta);
      assert (! attribute->loop_global);
      assert (! attribute->loop_level);
      assert (! attribute->loop_singlemap);
      
      BEGIN_MAP_LOOP(cgh, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
          Checkpoint ("Local mode call at %s to %s::%s",
                      attribute->where, attribute->thorn, attribute->routine);
          const int res = CCTK_CallFunction (function, attribute, data);
          assert (res==0);
        } END_LOCAL_COMPONENT_LOOP;
      }	END_MAP_LOOP;
      SyncGroupsInScheduleBlock( attribute, cgh );
      
    }
    
//     Checkpoint ("done with CallFunction.");
    
    // The return value indicates whether the grid functions have been
    // synchronised.
    // 0: let the flesh do the synchronisation
    // 1: we did the synchronisation
    return 1;
  }

  struct typed_group {
    int vartype;
    vector<int> members;
  };

  void SyncGroupsInScheduleBlock( cFunctionData* attribute, cGH* cgh )
  {
    // check if there is anything to do
    if (attribute->n_SyncGroups <= 0) return;

    vector<group_set> groups;

    // sort all grid variables into sets of the same vartype
    for (int g = 0; g < attribute->n_SyncGroups; g++) {
      // skip empty groups
      const int group = attribute->SyncGroups[g];
      if (CCTK_NumVarsInGroupI (group) <= 0) continue;

      group_set newset;
      const int firstvar = CCTK_FirstVarIndexI (group);
      newset.vartype = CCTK_VarTypeI (firstvar);
      assert (newset.vartype >= 0);
      int c;
      for (c = 0; c < groups.size(); c++) {
        if (newset.vartype == groups[c].vartype) {
          break;
        }
      }
      if (c == groups.size()) {
        groups.push_back (newset);
      }
      groups[c].members.push_back (group);
    }

    for (int c = 0; c < groups.size(); c++) {
      SyncGroups (cgh, groups[c]);
    }
  }
  
} // namespace Carpet
