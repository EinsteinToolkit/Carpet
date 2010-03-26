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
  CallScheduledFunction (char const * restrict time_and_mode,
                         void * function,
                         cFunctionData * attribute,
                         void * data,
                         Timer & user_timer);
  
  static void
  SyncGroupsInScheduleBlock (cFunctionData * attribute, cGH * cctkGH);
  
  /// Traverse one function on all components of one refinement level
  /// of one multigrid level.
  int CallFunction (void * function, ///< the function to call
                    cFunctionData * attribute, ///< attributes of the function
                    void * data) ///< private data for CCTK_CallFunction
  {
    DECLARE_CCTK_PARAMETERS;
    
    static Timer total_timer ("CallFunction");
    static Timer user_timer  ("CallFunction::thorns");
    static Timer sync_timer  ("CallFunction::syncs");
    
    total_timer.start();
    
    cGH * cctkGH = static_cast<cGH *> (data);
    
    assert (int (not not attribute->meta) +
            int (not not attribute->meta_early) +
            int (not not attribute->meta_late) +
            int (not not attribute->global) +
            int (not not attribute->global_early) +
            int (not not attribute->global_late) +
            int (not not attribute->level) +
            int (not not attribute->singlemap) +
            int (not not attribute->local)
            <= 1);
    
    assert (not not attribute->loop_global +
            not not attribute->loop_level +
            not not attribute->loop_singlemap +
            not not attribute->loop_local
            <= 1);
    
    if (attribute->meta or attribute->meta_early or attribute->meta_late or
        is_meta_mode())
    {
      // Convtest operation
      
      if ((attribute->meta and do_meta_mode) or
          (attribute->meta_early and do_early_meta_mode) or
          (attribute->meta_late and do_late_meta_mode) or
          is_meta_mode())
      {
        if (attribute->loop_local) {
          BEGIN_META_MODE(cctkGH) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                  BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                    CallScheduledFunction
                      ("Meta time local mode",
                       function, attribute, data, user_timer);
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
                  CallScheduledFunction
                    ("Meta time singlemap mode",
                     function, attribute, data, user_timer);
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
                CallScheduledFunction
                  ("Meta time level mode",
                   function, attribute, data, user_timer);
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else if (attribute->loop_global) {
          BEGIN_META_MODE(cctkGH) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              CallScheduledFunction
                ("Meta time global mode",
                 function, attribute, data, user_timer);
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                sync_timer.start();
                SyncGroupsInScheduleBlock (attribute, cctkGH) ;
                sync_timer.stop();
              } END_REFLEVEL_LOOP;
            } END_MGLEVEL_LOOP;
          } END_META_MODE;
        } else {
          BEGIN_META_MODE(cctkGH) {
            CallScheduledFunction
              ("Meta mode",
               function, attribute, data, user_timer);
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
      
    } else if (attribute->global or attribute->global_early or
               attribute->global_late or is_global_mode())
    {
      // Global operation: call once
      
      if ((attribute->global and do_global_mode) or
          (attribute->global_early and do_early_global_mode) or
          (attribute->global_late and do_late_global_mode) or
          is_global_mode())
      {
        if (attribute->loop_local) {
          BEGIN_GLOBAL_MODE(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                  CallScheduledFunction
                    ("Global time local mode",
                     function, attribute, data, user_timer);
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
                CallScheduledFunction
                  ("Global time singlemap mode",
                   function, attribute, data, user_timer);
              } END_MAP_LOOP;
              sync_timer.start();
              SyncGroupsInScheduleBlock (attribute, cctkGH) ;
              sync_timer.stop();
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else if (attribute->loop_level) {
          BEGIN_GLOBAL_MODE(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              CallScheduledFunction
                ("Global time level mode",
                 function, attribute, data, user_timer);
              sync_timer.start();
              SyncGroupsInScheduleBlock (attribute, cctkGH) ;
              sync_timer.stop();
            } END_REFLEVEL_LOOP;
          } END_GLOBAL_MODE;
        } else {
          BEGIN_GLOBAL_MODE(cctkGH) {
            CallScheduledFunction
              ("Global mode",
               function, attribute, data, user_timer);
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
      
      if (attribute->loop_local) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            CallScheduledFunction
              ("Level time local mode",
               function, attribute, data, user_timer);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } else if (attribute->loop_singlemap) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          CallScheduledFunction
            ("Level time singlemap mode",
             function, attribute, data, user_timer);
        } END_MAP_LOOP;
      } else {
        CallScheduledFunction
          ("Level mode",
           function, attribute, data, user_timer);
      }
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    } else if (attribute->singlemap) {
      // Single map operation: call once per refinement level and map
      
      if (attribute->loop_local) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            CallScheduledFunction
              ("Singlemap time local mode",
               function, attribute, data, user_timer);
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } else {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
          CallScheduledFunction
            ("Singlemap mode",
             function, attribute, data, user_timer);
        } END_MAP_LOOP;
      }
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    } else {
      // Local operation: call once per component
      
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          CallScheduledFunction
            ("Local mode",
             function, attribute, data, user_timer);
        } END_LOCAL_COMPONENT_LOOP;
      }	END_MAP_LOOP;
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    }
    
    if (schedule_barriers) {
      static unsigned int magic = 0xe8932329UL; // a random starting value
      unsigned int recv = magic;
      Checkpoint ("About to Bcast %u", magic);
      MPI_Bcast (& recv, 1, MPI_UNSIGNED, 0, dist::comm());
      Checkpoint ("Finished Bcast");
      if (recv != magic) {
        CCTK_WARN (CCTK_WARN_ABORT,
                   "Inconsistent communication schedule: not all processes return to CallFunction at the same time");
      }
      ++ magic;
      Checkpoint ("About to Barrier");
      MPI_Barrier (dist::comm());
      Checkpoint ("Finished Barrier");
    }
    
    total_timer.stop();
    
    // The return value indicates whether the grid functions have been
    // synchronised.
    // 0: let the flesh do the synchronisation
    // 1: we did the synchronisation
    return 1;
  }
  
  
  
  void
  CallScheduledFunction (char const * restrict const time_and_mode,
                         void * const function,
                         cFunctionData * const attribute,
                         void * const data,
                         Timer & user_timer)
  {
    cGH const * const cctkGH = static_cast <cGH const *> (data);
    Checkpoint ("%s call at %s to %s::%s",
                time_and_mode,
                attribute->where,
                attribute->thorn, attribute->routine);
    CallBeforeRoutines (cctkGH, function, attribute, data);
    user_timer.start();
    int const res = CCTK_CallFunction (function, attribute, data);
    user_timer.stop();
    CallAfterRoutines (cctkGH, function, attribute, data);
    assert (res==0);
  }
  
  
  
  void SyncGroupsInScheduleBlock (cFunctionData* attribute, cGH* cctkGH)
  {
    // check if there is anything to do
    if (attribute->n_SyncGroups <= 0) return;

    // remove all empty groups from the list
    vector<int> groups;
    groups.reserve (attribute->n_SyncGroups);
    for (int g = 0; g < attribute->n_SyncGroups; g++) {
      const int group = attribute->SyncGroups[g];
      if (CCTK_NumVarsInGroupI (group) > 0) {
        groups.push_back (group);
      }
    }

    SyncProlongateGroups (cctkGH, groups);
  }
  
} // namespace Carpet
