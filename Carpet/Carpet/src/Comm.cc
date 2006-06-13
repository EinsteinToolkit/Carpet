#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"



namespace Carpet {

  using namespace std;


  static void ProlongateGroupBoundaries (const cGH* cctkGH,
                                         CCTK_REAL initial_time,
                                         const vector<int>& groups);


  // Carpet's overload function for CCTK_SyncGroup()
  // which synchronises a single group
  //
  // returns 0 for success
  //        -1 if the group doesn't have storage assigned
  //        -2 if the given groupname is invalid
  int SyncGroup (const cGH* cctkGH, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    int retval;

    const int group = CCTK_GroupIndex(groupname);
    if (group >= 0) {
      assert (group < (int)arrdata.size());

      const vector<int> groups(1, group);
      retval = SyncProlongateGroups (cctkGH, groups);
    } else {
      retval = -2;
    }

    return retval;
  }


  // synchronises ghostzones and prolongates boundaries of a set of groups
  //
  // returns 0 for success and -1 if the set contains a group with no storage
  int SyncProlongateGroups (const cGH* cctkGH, const vector<int>& groups)
  {
    int retval = 0;
    DECLARE_CCTK_PARAMETERS;

    assert (groups.size() > 0);

    // check consistency of all groups:
    // create a new set with empty and no-storage groups removed
    vector<int> goodgroups;
    for (size_t g = 0; g < groups.size(); g++) {
      const int group = groups[g];
      const int grouptype = CCTK_GroupTypeI (group);
      char* groupname = CCTK_GroupName (group);
      Checkpoint ("SyncGroup \"%s\" time=%g",
                  groupname, (double) cctkGH->cctk_time);

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
          if (maps == 1 and vhh.at(map)->local_components(reflevel) == 1) {
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

      if (not CCTK_QueryGroupStorageI (cctkGH, group)) {
        CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot synchronise group \"%s\" because it has no storage",
                    groupname);
        retval = -1;
      }
      else if (CCTK_NumVarsInGroupI (group) > 0) {
        goodgroups.push_back(group);
      }

      free (groupname);
    }

    if (goodgroups.size() > 0) {
      // prolongate boundaries
      if (reflevel > 0) {
        if (do_prolongate) {
          bool local_do_prolongate;
          if (use_tapered_grids) {
            // TODO: Check iteration number instead
            CCTK_REAL mytime;
            CCTK_REAL parenttime;
            if (map == -1) {
              mytime = vtt.at(0)->time (0, reflevel, mglevel);
              parenttime = vtt.at(0)->time (0, reflevel - 1, mglevel);
            } else {
              mytime = vtt.at(map)->time (0, reflevel, mglevel);
              parenttime = vtt.at(map)->time (0, reflevel - 1, mglevel);
            }
            bool const in_sync =
              abs (mytime - parenttime) < 1.0e-10 * abs (delta_time);
            local_do_prolongate = in_sync;
          } else {
            local_do_prolongate = true;
          }
          if (local_do_prolongate) {
            ProlongateGroupBoundaries (cctkGH, cctk_initial_time, goodgroups);
          }
        }
      }

      // synchronise ghostzones
      SyncGroups (cctkGH, goodgroups);
    }

    return retval;
  }

  // Prolongate the boundaries of all CCTK_GF groups in the given set
  static void ProlongateGroupBoundaries (const cGH* cctkGH,
                                         CCTK_REAL initial_time,
                                         const vector<int>& groups)
  {
    DECLARE_CCTK_PARAMETERS;
    const int tl = 0;

    // use the current time here (which may be modified by the user)
    const CCTK_REAL time
      = (cctkGH->cctk_time - initial_time) / delta_time;

    for (comm_state state; not state.done(); state.step()) {
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups[group];
        const int grouptype = CCTK_GroupTypeI (g);
        if (grouptype != CCTK_GF) {
          continue;
        }
        assert (reflevel>=0 and reflevel<reflevels);

        for (int m = 0; m < (int)arrdata.at(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.at(g).at(m).data.size(); ++v) {
            ggf *const gv = arrdata.at(g).at(m).data.at(v);
            for (int c = 0; c < vhh.at(m)->components(reflevel); ++c) {
              gv->ref_bnd_prolongate (state, tl, reflevel, c, mglevel, time);
            }
          }
        }
      }
    }
  }


  // synchronises a set of groups
  void SyncGroups (const cGH* cctkGH, const vector<int>& groups)
  {
    DECLARE_CCTK_PARAMETERS;
    const int tl = 0;

    assert (groups.size() > 0);

    for (comm_state state; not state.done(); state.step()) {
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups[group];
        const int grouptype = CCTK_GroupTypeI (g);
        const int ml = grouptype == CCTK_GF ? mglevel : 0;
        const int rl = grouptype == CCTK_GF ? reflevel : 0;
        for (int m = 0; m < (int)arrdata.at(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.at(g).at(m).data.size(); ++v) {
            arrdesc& array = arrdata.at(g).at(m);
            for (int c = 0; c < array.hh->components(rl); ++c) {
              array.data.at(v)->sync (state, tl, rl, c, ml);
            }
          }
        }
      }
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
