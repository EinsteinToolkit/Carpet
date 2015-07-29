#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>
#include <iterator>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Requirements.hh>

#include <Timer.hh>

#include <ggf.hh>
#include <gh.hh>

#include <carpet.hh>

//#undef KNARFDEBL
//#define KNARFDEBL 2

namespace Carpet {

  using namespace std;

  static string intvectString(const vector<int> vect) {
    std::stringstream result;
    std::copy(vect.begin(), vect.end(), std::ostream_iterator<int>(result, " "));
    return result.str();
  }

  static void ProlongateGroupBoundaries (const cGH* cctkGH,
                                         const vector<int>& groups);


  // Carpet's overload function for CCTK_SyncGroupsByDirI()
  // which synchronises a set of groups given by their indices.
  // Synchronisation of individual directions is not (yet?) implemented.
  //
  // returns the number of groups successfully synchronised
  //        -1 if a group in the set doesn't have storage assigned
  int SyncGroupsByDirI (const cGH* cctkGH,
                        int num_groups,
                        const int *groups,
                        const int *directions)
  {
    int group, retval = 0;
    vector<int> groups_set;
    groups_set.reserve(num_groups);

    // individual directions aren't supported (yet?)
    if (directions != NULL) {
      CCTK_WARN (0, "Carpet doesn't support synchronisation of individual "
                    "directions");
    }

    for (group = 0; group < num_groups; group++) {
      if (CCTK_NumVarsInGroupI (groups[group]) > 0) {
        groups_set.push_back (groups[group]);
      }
    }

    if (groups_set.size() > 0) {
      retval = SyncProlongateGroups (cctkGH, groups_set);
      if (retval == 0) {
        retval = groups_set.size();
      }
    }

    return retval;
  }


  // synchronises ghostzones and prolongates boundaries of a set of groups
  //
  // returns 0 for success and -1 if the set contains a group with no storage
  int SyncProlongateGroups_ (const cGH* cctkGH, const vector<int>& groups,
                            cFunctionData const* function_data)
  {
    int retval = 0;
    DECLARE_CCTK_PARAMETERS;

    assert (groups.size() > 0);

    // check consistency of all groups:
    // create a new set with empty and no-storage groups removed
    vector<int> goodgroups;
    goodgroups.reserve (groups.size());
    for (size_t group = 0; group < groups.size(); group++) {
      const int g = groups.AT(group);
      const int grouptype = CCTK_GroupTypeI (g);
      char* const groupname = CCTK_GroupName (g);
      Checkpoint ("SyncGroup \"%s\" iteration=%d time=%g",
                  groupname,
                  cctkGH->cctk_iteration, (double) cctkGH->cctk_time);

      if (grouptype == CCTK_GF) {
        if (reflevel == -1) {
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Cannot synchronise in global mode "
                      "(Tried to synchronise group \"%s\")",
                      groupname);
        }
        if (map != -1 and component == -1) {
          if (maps == 1) {
            CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Synchronising group \"%s\" in singlemap mode",
                        groupname);
          } else {
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Cannot synchronise in singlemap mode "
                        "(Tried to synchronise group \"%s\")",
                        groupname);
          }
        }
        if (component != -1) {
          if (maps == 1 and vhh.AT(map)->local_components(reflevel) == 1) {
            CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Synchronising group \"%s\" in local mode",
                        groupname);
          } else {
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Cannot synchronise in local mode "
                        "(Tried to synchronise group \"%s\")",
                        groupname);
          }
        }
      }

      if (not CCTK_QueryGroupStorageI (cctkGH, g)) {
        CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot synchronise group \"%s\" because it has no storage",
                    groupname);
        retval = -1;
      }
      else if (CCTK_NumVarsInGroupI (g) > 0) {
        goodgroups.push_back(g);
      }

      free (groupname);
    } // for g

    if (goodgroups.size() > 0) {
      
#ifdef REQUIREMENTS_HH
      Requirements::Sync(function_data, goodgroups,
                         cctkGH->cctk_iteration, reflevel, timelevel);
#endif
      
      // prolongate boundaries
      bool const local_do_prolongate = do_prolongate and not do_taper;

      // PTODO: Is this necessary?
      if(use_psamr and carpet_cctk_iteration>0)
      {
        if(reflevel>0 and (sync_during_time_integration or local_do_prolongate))
        {
          static Timers::Timer timer("Sync");
          timer.start();
          reflevel--;
          KNARFDEB(2, "KNARF preprolongate sync lvl %d on proc %d, groups: %s\n",
                      reflevel, CCTK_MyProc(NULL), intvectString(goodgroups).c_str());
          SyncGroups(cctkGH, goodgroups);
          timer.stop();
          reflevel++;
        }
      }

      if (local_do_prolongate) {
        static Timers::Timer timer ("Prolongate");
        timer.start();
        ProlongateGroupBoundaries (cctkGH, goodgroups);
        timer.stop();
        KNARFDEB(2, "KNARF do prolongate lvl %d on proc %d\n", reflevel, CCTK_MyProc(NULL));
      }
      else
        KNARFDEB(2, "KNARF skip prolongate lvl %d on proc %d\n", reflevel, CCTK_MyProc(NULL));

      // This was found to be necessary on Hopper, otherwise memory
      // seems to be overwritten while prolongating/syncronizing. It
      // looks liks this might be an MPI implementation issue, but
      // this is not clear. A barrier at this point seems to be a
      // sufficient workaround, and is now used on Hopper. For more
      // information about this ask Frank Loeffler
      // <knarf@cct.lsu.edu>.
#ifdef CARPET_MPI_BARRIER_PROLONGATE_SYNC
      Carpet::NamedBarrier(cctkGH,
                           8472211063, "CARPET_MPI_BARRIER_PROLONGATE_SYNC");
#endif

      if (sync_barriers)
      {
        static Timers::Timer barrier_timer ("ProlongateSyncBarrier");
        barrier_timer.start();
        CCTK_Barrier(cctkGH);
        barrier_timer.stop();
      }

      // synchronise ghostzones
      if (sync_during_time_integration or local_do_prolongate) {
        static Timers::Timer timer ("Sync");
        timer.start();
        KNARFDEB(2, "KNARF postprolongate sync lvl %d on proc %d, groups: %s\n",
                    reflevel, CCTK_MyProc(NULL), intvectString(goodgroups).c_str());
        SyncGroups (cctkGH, goodgroups);
        timer.stop();
      }
      
    }
    
    return retval;
  }

  // KNARF TODO
  int SyncProlongateGroups (const cGH* cctkGH, const vector<int>& groups,
                            cFunctionData const* function_data)
  {
    DECLARE_CCTK_PARAMETERS;
    // PTODO: Why check iteration?
    if (use_psamr and carpet_cctk_iteration != 0)
    {
      int save_reflevel = reflevel;
      reflevel = reflevels-1;
      // Prolongate the finest grid
      int retval = SyncProlongateGroups_(cctkGH, groups, function_data);
      reflevel = save_reflevel;

      if (carpet_next_coarse_lvl < 0)
        return retval;
      if (carpet_next_coarse_lvl > 0)
      {
        reflevel = carpet_next_coarse_lvl;
        // Prolongate one of the coarse grids
        retval = SyncProlongateGroups_(cctkGH, groups, function_data);
        reflevel = save_reflevel;
        // PTODO: We want to synchronize 'reflevel' as well - or do we?
        SyncGroups(cctkGH, groups);
      }
      else
      {
        reflevel = 0;
        SyncGroups(cctkGH, groups);
      }
      reflevel = save_reflevel;
      return retval;
    }
    return SyncProlongateGroups_(cctkGH, groups, function_data);
  }

  // Prolongate the boundaries of all CCTK_GF groups in the given set
  static void ProlongateGroupBoundaries (const cGH* cctkGH,
                                         const vector<int>& groups)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // use the current time here (which may be modified by the user)
    CCTK_REAL time = cctkGH->cctk_time;
    if (use_psamr)
      time = cctk_initial_time + carpet_level_iteration[reflevel] * delta_time / maxtimereflevelfact;

    if (reflevel == 0)
      return;
    
    Checkpoint ("ProlongateGroups");

    assert (groups.size() > 0);
    
    // TODO: Transmit only the regions that are actually needed
    if (CCTK_IsFunctionAliased("Accelerator_RequireValidData")) {
      // TODO: Require only as many time levels as needed
      // TODO: Handle time levels correctly
      vector<CCTK_INT> vis, rls, tls;
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int grouptype = CCTK_GroupTypeI (g);
        if (grouptype != CCTK_GF) {
          continue;
        }
        const int num_tl = prolongation_order_time + 1;
        const int var0 = CCTK_FirstVarIndexI(g);
        const int varn = CCTK_NumVarsInGroupI(g);
        for (int vi=var0; vi<var0+varn; ++vi) {
          for (int tl=0; tl<num_tl; ++tl) {
            vis.push_back(vi);
            rls.push_back(reflevel-1);
            tls.push_back(tl);
          }
        }
      }
      const CCTK_INT on_device = 0;
      Accelerator_RequireValidData
        (cctkGH,
         &vis.front(), &rls.front(), &tls.front(), vis.size(), on_device);
    }

    static vector<Timers::Timer*> timers;
    if (timers.empty()) {
      timers.push_back(new Timers::Timer("comm_state[0].create"));
      for (astate state = static_cast<astate>(0);
           state != state_done;
           state = static_cast<astate>(static_cast<int>(state)+1))
      {
        ostringstream name1;
        name1 << "comm_state[" << timers.size() << "]"
              << "." << tostring(state) << ".user";
        timers.push_back(new Timers::Timer(name1.str()));
        ostringstream name2;
        name2 << "comm_state[" << timers.size() << "]"
              << "." << tostring(state) << ".step";
        timers.push_back(new Timers::Timer(name2.str()));
      }
    }
    
    vector<Timers::Timer*>::iterator ti = timers.begin();
    (*ti)->start();
    for (comm_state state; not state.done(); state.step()) {
      (*ti)->stop(); ++ti; (*ti)->start();
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int grouptype = CCTK_GroupTypeI (g);
        if (grouptype != CCTK_GF) {
          continue;
        }
        assert (reflevel>=0 and reflevel<reflevels);
        const int active_tl =
          groupdata.AT(g).activetimelevels.AT(mglevel).AT(reflevel);
        assert (active_tl>=0);
        const int tl = active_tl > 1 ? timelevel : 0;
        
        for (int m = 0; m < (int)arrdata.AT(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
            ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
            gv->ref_bnd_prolongate_all
              (state, tl, reflevel, mglevel, time);
          }
        }
      }
      (*ti)->stop(); ++ti; (*ti)->start();
    }
    (*ti)->stop();
    ++ti; assert(ti == timers.end());
    
    // TODO: Transmit only the regions that were actually written
    // (i.e. the buffer zones)
    if (CCTK_IsFunctionAliased("Accelerator_NotifyDataModified")) {
      // TODO: Handle time levels correctly
      vector<CCTK_INT> vis, rls, tls;
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int grouptype = CCTK_GroupTypeI (g);
        if (grouptype != CCTK_GF) {
          continue;
        }
        const int var0 = CCTK_FirstVarIndexI(g);
        const int varn = CCTK_NumVarsInGroupI(g);
        for (int vi=var0; vi<var0+varn; ++vi) {
          vis.push_back(vi);
          rls.push_back(reflevel-1);
          tls.push_back(timelevel);
        }
      }
      const CCTK_INT on_device = 0;
      Accelerator_NotifyDataModified
        (cctkGH,
         &vis.front(), &rls.front(), &tls.front(), vis.size(), on_device);
    }
  }


  // synchronises a set of groups
  void SyncGroups (const cGH* cctkGH, const vector<int>& groups)
  {
    DECLARE_CCTK_PARAMETERS;

    Checkpoint ("SyncGroups");

    assert (groups.size() > 0);

    KNARFDEB(1, "KNARF SYNC, rl %d, proc %d, groups: %s\n", reflevel, CCTK_MyProc(NULL),
                intvectString(groups).c_str());
    
    if (CCTK_IsFunctionAliased("Accelerator_PreSync")) {
      vector<CCTK_INT> groups_(groups.size());
      for (size_t i=0; i<groups.size(); ++i) groups_[i] = groups[i];
      Accelerator_PreSync(cctkGH, &groups_.front(), groups_.size());
    }

    static vector<Timers::Timer*> timers;
    if (timers.empty()) {
      timers.push_back(new Timers::Timer("comm_state[0].create"));
      for (astate state = static_cast<astate>(0);
           state != state_done;
           state = static_cast<astate>(static_cast<int>(state)+1))
      {
        ostringstream name1;
        name1 << "comm_state[" << timers.size() << "]"
              << "." << tostring(state) << ".user";
        timers.push_back(new Timers::Timer(name1.str()));
        ostringstream name2;
        name2 << "comm_state[" << timers.size() << "]"
              << "." << tostring(state) << ".step";
        timers.push_back(new Timers::Timer(name2.str()));
      }
    }
    
    vector<Timers::Timer*>::iterator ti = timers.begin();
    (*ti)->start();
    for (comm_state state; not state.done(); state.step()) {
      (*ti)->stop(); ++ti; (*ti)->start();
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int grouptype = CCTK_GroupTypeI (g);
        const int ml = grouptype == CCTK_GF ? mglevel : 0;
        const int rl = grouptype == CCTK_GF ? reflevel : 0;
        const int active_tl = groupdata.AT(g).activetimelevels.AT(ml).AT(rl);
        assert (active_tl>=0);
        const int tl = active_tl > 1 ? timelevel : 0;
        for (int m = 0; m < (int)arrdata.AT(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
            arrdesc& array = arrdata.AT(g).AT(m);
            array.data.AT(v)->sync_all (state, tl, rl, ml);
          }
        }
      }
      (*ti)->stop(); ++ti; (*ti)->start();
    }
    (*ti)->stop();
    ++ti; assert(ti == timers.end());
    
    if (CCTK_IsFunctionAliased("Accelerator_PostSync")) {
      vector<CCTK_INT> groups_(groups.size());
      for (size_t i=0; i<groups.size(); ++i) groups_[i] = groups[i];
      Accelerator_PostSync(cctkGH, &groups_.front(), groups_.size());
    }
  }


  int EnableGroupComm (const cGH* cctkGH, const char* groupname)
  {
    // Communication is always enabled
    return 0;
  }

  int DisableGroupComm (const cGH* cctkGH, const char* groupname)
  {
    // Communication is always enabled
    return -1;
  }

} // namespace Carpet
