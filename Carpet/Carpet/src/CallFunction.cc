#include <algorithm>
#include <cassert>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctki_GHExtensions.h>

#include <gh.hh>

#include "carpet.hh"
#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  static void
  SyncGroupsInScheduleBlock (cFunctionData * attribute, cGH * cctkGH);
  
  /// Traverse one function on all components of one refinement level
  /// of one multigrid level.
  int CallFunction (void * function, ///< the function to call
                    cFunctionData * attribute, ///< attributes of the function
                    void * data) ///< private data for CCTK_CallFunction
  {
    DECLARE_CCTK_PARAMETERS;
    
    static Timer total_timer (timerSet(), "CallFunction");
    static Timer user_timer  (timerSet(), "CallFunction::thorns");
    static Timer sync_timer  (timerSet(), "CallFunction::syncs");
    
    total_timer.start();
    
    cGH * cctkGH = static_cast<cGH *> (data);
    
    assert (not not attribute->meta +
            not not attribute->global +
            not not attribute->level +
            not not attribute->singlemap +
            not not attribute->local
            <= 1);
    
    assert (not not attribute->loop_global +
            not not attribute->loop_level +
            not not attribute->loop_singlemap +
            not not attribute->loop_local
            <= 1);
    
    if (attribute->meta or is_meta_mode()) {
      // Convtest operation
      
      if (do_meta_mode) {
        if (attribute->loop_local) {
          BEGIN_META_MODE(cctkGH) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                  BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                    Checkpoint ("Meta time local mode call at %s to %s::%s",
                                attribute->where,
                                attribute->thorn, attribute->routine);
                    user_timer.start();
                    const int res = CCTK_CallFunction
                      (function, attribute, data);
                    user_timer.stop();
                    assert (res==0);
                  } END_LOCAL_COMPONENT_LOOP;
                } END_MAP_LOOP;
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_singlemap) {
          BEGIN_META_MODE(cctkGH) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                  Checkpoint ("Meta time singlemap mode call at %s to %s::%s",
                              attribute->where,
                              attribute->thorn, attribute->routine);
                  user_timer.start();
                  const int res = CCTK_CallFunction (function, attribute, data);
                  user_timer.stop();
                  assert (res==0);
                } END_MAP_LOOP;
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_level) {
          BEGIN_META_MODE(cctkGH) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                Checkpoint ("Meta time level mode call at %s to %s::%s",
                            attribute->where,
                            attribute->thorn, attribute->routine);
                user_timer.start();
                const int res = CCTK_CallFunction (function, attribute, data);
                user_timer.stop();
                assert (res==0);
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_global) {
          BEGIN_META_MODE(cctkGH) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              Checkpoint ("Meta time global mode call at %s to %s::%s",
                          attribute->where,
                          attribute->thorn, attribute->routine);
              user_timer.start();
              const int res = CCTK_CallFunction (function, attribute, data);
              user_timer.stop();
              assert (res==0);
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else {
          BEGIN_META_MODE(cctkGH) {
            Checkpoint ("Meta mode call at %s to %s::%s",
                        attribute->where,
                        attribute->thorn, attribute->routine);
            user_timer.start();
            const int res = CCTK_CallFunction (function, attribute, data);
            user_timer.stop();
            assert (res==0);
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        }
      }
      
    } else if (attribute->global or is_global_mode()) {
      // Global operation: call once
      
      assert (not attribute->loop_meta);
      
      if (do_global_mode) {
        if (attribute->loop_local) {
          BEGIN_GLOBAL_MODE(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                  Checkpoint ("Global time local mode call at %s to %s::%s",
                              attribute->where,
                              attribute->thorn, attribute->routine);
                  user_timer.start();
                  const int res = CCTK_CallFunction (function, attribute, data);
                  user_timer.stop();
                  assert (res==0);
                } END_LOCAL_COMPONENT_LOOP;
              } END_MAP_LOOP;
              sync_timer.start();
              SyncGroupsInScheduleBlock (attribute, cctkGH) ;
              sync_timer.stop();
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else if (attribute->loop_singlemap) {
          BEGIN_GLOBAL_MODE(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                Checkpoint ("Global time singlemap mode call at %s to %s::%s",
                            attribute->where,
                            attribute->thorn, attribute->routine);
                user_timer.start();
                const int res = CCTK_CallFunction (function, attribute, data);
                user_timer.stop();
                assert (res==0);
              } END_MAP_LOOP;
              sync_timer.start();
              SyncGroupsInScheduleBlock (attribute, cctkGH) ;
              sync_timer.stop();
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else if (attribute->loop_level) {
          BEGIN_GLOBAL_MODE(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              Checkpoint ("Global time level mode call at %s to %s::%s",
                          attribute->where,
                          attribute->thorn, attribute->routine);
              user_timer.start();
              const int res = CCTK_CallFunction (function, attribute, data);
              user_timer.stop();
              assert (res==0);
              sync_timer.start();
              SyncGroupsInScheduleBlock (attribute, cctkGH) ;
              sync_timer.stop();
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else {
          BEGIN_GLOBAL_MODE(cctkGH) {
            Checkpoint ("Global mode call at %s to %s::%s",
                        attribute->where,
                        attribute->thorn, attribute->routine);
            user_timer.start();
            const int res = CCTK_CallFunction (function, attribute, data);
            user_timer.stop();
            assert (res==0);
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              sync_timer.start();
              SyncGroupsInScheduleBlock (attribute, cctkGH) ;
              sync_timer.stop();
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        }
      }
      
    } else if (attribute->level) {
      // Level operation: call once per refinement level
      
      assert (not attribute->loop_meta);
      assert (not attribute->loop_global);
      
      if (attribute->loop_local) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            Checkpoint ("Level time local mode call at %s to %s::%s",
                        attribute->where,
                        attribute->thorn, attribute->routine);
            user_timer.start();
            const int res = CCTK_CallFunction (function, attribute, data);
            user_timer.stop();
            assert (res==0);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } else if (attribute->loop_singlemap) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          Checkpoint ("Level time singlemap mode call at %s to %s::%s",
                      attribute->where,
                      attribute->thorn, attribute->routine);
          user_timer.start();
          const int res = CCTK_CallFunction (function, attribute, data);
          user_timer.stop();
          assert (res==0);
        } END_MAP_LOOP;
      } else {
        Checkpoint ("Level mode call at %s to %s::%s",
                    attribute->where,
                    attribute->thorn, attribute->routine);
        user_timer.start();
        const int res = CCTK_CallFunction (function, attribute, data);
        user_timer.stop();
        assert (res==0);
      }
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    } else if (attribute->singlemap) {
      // Single map operation: call once per refinement level and map
      
      assert (not attribute->loop_meta);
      assert (not attribute->loop_global);
      assert (not attribute->loop_level);
      
      if (attribute->loop_local) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            Checkpoint ("Singlemap time local mode call at %s to %s::%s",
                        attribute->where,
                        attribute->thorn, attribute->routine);
            user_timer.start();
            const int res = CCTK_CallFunction (function, attribute, data);
            user_timer.stop();
            assert (res==0);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } else {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          Checkpoint ("Singlemap mode call at %s to %s::%s",
                      attribute->where,
                      attribute->thorn, attribute->routine);
          user_timer.start();
          const int res = CCTK_CallFunction (function, attribute, data);
          user_timer.stop();
          assert (res==0);
        } END_MAP_LOOP;
      }
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    } else {
      // Local operation: call once per component
      
      assert (not attribute->loop_meta);
      assert (not attribute->loop_global);
      assert (not attribute->loop_level);
      assert (not attribute->loop_singlemap);
      
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          Checkpoint ("Local mode call at %s to %s::%s",
                      attribute->where,
                      attribute->thorn, attribute->routine);
          user_timer.start();
          const int res = CCTK_CallFunction (function, attribute, data);
          user_timer.stop();
          assert (res==0);
        } END_LOCAL_COMPONENT_LOOP;
      }	END_MAP_LOOP;
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    }
    
    if (schedule_barriers) CCTK_Barrier (cctkGH);
    
    total_timer.stop();
    
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

  void SyncGroupsInScheduleBlock (cFunctionData* attribute, cGH* cctkGH)
  {
    // check if there is anything to do
    if (attribute->n_SyncGroups <= 0) return;

    // remove all empty groups from the list
    vector<int> groups;
    for (int g = 0; g < attribute->n_SyncGroups; g++) {
      const int group = attribute->SyncGroups[g];
      if (CCTK_NumVarsInGroupI (group) > 0) {
        groups.push_back (group);
      }
    }

    SyncProlongateGroups (cctkGH, groups);
  }
  
} // namespace Carpet
