#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctki_GHExtensions.h>

#include <gh.hh>

#include <carpet.hh>
#include <Timers.hh>

#include "adler32.hh"



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
    static Timer user_timer  ("thorns");
    static Timer sync_timer  ("syncs");
    
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
                BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
                  BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                    CallScheduledFunction
                      ("Meta time local mode",
                       function, attribute, data, user_timer);
                  } END_LOCAL_COMPONENT_LOOP;
                } END_LOCAL_MAP_LOOP;
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
              BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
                BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                  CallScheduledFunction
                    ("Global time local mode",
                     function, attribute, data, user_timer);
                } END_LOCAL_COMPONENT_LOOP;
              } END_LOCAL_MAP_LOOP;
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
        BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            CallScheduledFunction
              ("Level time local mode",
               function, attribute, data, user_timer);
          } END_LOCAL_COMPONENT_LOOP;
        } END_LOCAL_MAP_LOOP;
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
        BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            CallScheduledFunction
              ("Singlemap time local mode",
               function, attribute, data, user_timer);
          } END_LOCAL_COMPONENT_LOOP;
        } END_LOCAL_MAP_LOOP;
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
      
      BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          CallScheduledFunction
            ("Local mode",
             function, attribute, data, user_timer);
        } END_LOCAL_COMPONENT_LOOP;
      }	END_LOCAL_MAP_LOOP;
      sync_timer.start();
      SyncGroupsInScheduleBlock (attribute, cctkGH) ;
      sync_timer.stop();
      
    }
    
    if (schedule_barriers) {
      // Create an ID that is almost unique for this scheduled
      // function call
      stringstream buf;
      buf << attribute->meta
          << attribute->meta_early
          << attribute->meta_late
          << attribute->global
          << attribute->global_early
          << attribute->global_late
          << attribute->level
          << attribute->singlemap
          << attribute->local << "\n";
      buf << attribute->where << "\n";
      buf << attribute->thorn << "\n";
      buf << attribute->routine << "\n";
      string const str = buf.str();
      int const id = adler32(str.c_str(), str.length());
      Carpet::NamedBarrier (NULL, id, "Carpet::CallFunction");
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
    int const skip = CallBeforeRoutines (cctkGH, function, attribute, data);
    if (not skip) {
      Timer timer(attribute->routine);
      
      // Save the time step size
      CCTK_REAL const saved_cctk_delta_time = cctkGH->cctk_delta_time;
      
      user_timer.start();
      timer.start();
      if (CCTK_IsFunctionAliased("Accelerator_PreCallFunction")) {
        Timer pre_timer("PreCall");
        pre_timer.start();
        Accelerator_PreCallFunction(cctkGH, attribute);
        pre_timer.stop();
      }
      int const res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      if (CCTK_IsFunctionAliased("Accelerator_PostCallFunction")) {
        Timer post_timer("PostCall");
        post_timer.start();
        Accelerator_PostCallFunction(cctkGH, attribute);
        post_timer.stop();
      }
      timer.stop();
      user_timer.stop();
      
      // Manage the time step size. If the time step size changes
      // during initialisation, assume it is thorn Time, and update
      // the time hierarchy. If it changes during evolution, assume it
      // is MoL, and do nothing.
      if (cctkGH->cctk_iteration == 0 and
          cctkGH->cctk_delta_time != saved_cctk_delta_time)
      {
        // The user changed cctk_delta_time during initialisation --
        // update our internals and the time hierarchy
        bool const is_global =
          attribute->meta         or
          attribute->meta_early   or
          attribute->meta_late    or
          attribute->global       or
          attribute->global_early or
          attribute->global_late;
        delta_time =
          cctkGH->cctk_delta_time / mglevelfact *
          (is_global ? 1.0 : timereflevelfact);
        for (int ml=0; ml<mglevels; ++ml) {
          for (int rl=0; rl<reflevels; ++rl) {
            // Update the time delta
            CCTK_REAL const dt =
              delta_time / timereffacts.AT(rl) * ipow(mgfact, ml);
            tt->set_delta(ml,rl,dt);
            CCTK_REAL const t0 = tt->get_time(ml,rl,0);
            // Update the times of the past timelevels
            for (int tl=1; tl<timelevels; ++tl) {
              CCTK_REAL const t = t0 - tl * dt;
              tt->set_time(ml,rl,tl,t);
            }
          }
        }
      }
      
    }
    CallAfterRoutines (cctkGH, function, attribute, data);
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
