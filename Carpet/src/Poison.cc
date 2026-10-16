#include <cassert>
#include <cstdlib>
#include <cstring>
#include <limits>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <Timer.hh>

#include <defs.hh>
#include <typeprops.hh>
#include "PreSync.hh"

#include <carpet.hh>

#include <loopcontrol.h>

namespace Carpet {

using namespace std;

// The parameter where specifies which time levels should be
// poisoned.  what specifies what kind of grid variables should be
// poisoned.
void Poison(cGH const *const cctkGH, checktimes const where, int const what) {
  DECLARE_CCTK_PARAMETERS;

  assert(what == 0 or what == CCTK_GF or what == CCTK_ARRAY);

  if (not poison_new_timelevels)
    return;

  Timers::Timer timer("Poison");
  timer.start();

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {
      int const grouptype = CCTK_GroupTypeI(group);
      if (what == 0 or (what == CCTK_GF and grouptype == CCTK_GF) or
          (what == CCTK_ARRAY and
           (grouptype == CCTK_ARRAY or grouptype == CCTK_SCALAR))) {
        PoisonGroup(cctkGH, group, where);
      }
    } // if has storage
  }   // for group
  timer.stop();
}

void PoisonGroup(cGH const *const cctkGH, int const group,
                 checktimes const where) {
  DECLARE_CCTK_PARAMETERS;

  assert(group >= 0 and group < CCTK_NumGroups());

  if (not poison_new_timelevels)
    return;

  if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Cannot poison group \"%s\" because it has no storage",
               groupname);
    free(groupname);
    return;
  }

  int const nvar = CCTK_NumVarsInGroupI(group);
  if (nvar == 0)
    return;
  int const n0 = CCTK_FirstVarIndexI(group);
  assert(n0 >= 0);
  int const sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
  assert(sz > 0);

  int const table = CCTK_GroupTagsTableI(group);
  assert(table >= 0);
  bool persistent;
  char buf[100];
  int const ilen = Util_TableGetString(table, sizeof buf, buf, "Persistent");
  if (ilen > 0) {
    if (CCTK_EQUALS(buf, "yes")) {
      persistent = true;
    } else if (CCTK_EQUALS(buf, "no")) {
      persistent = false;
    } else {
      assert(0);
    }
  } else if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // default
    persistent = true;
  } else {
    assert(0);
  }

  int const num_tl = CCTK_ActiveTimeLevelsVI(cctkGH, n0);
  assert(num_tl > 0);
  int const min_tl = min_timelevel(where, num_tl, persistent);
  int const max_tl = max_timelevel(where, num_tl, persistent);

  if (min_tl <= max_tl) {

    {
      char *const groupname = CCTK_GroupName(group);
      Checkpoint("PoisonGroup \"%s\"", groupname);
      free(groupname);
    }

    CCTK_INT const poison_value = get_poison_value();

    int const grouptype = CCTK_GroupTypeI(group);

    BEGIN_LOCAL_MAP_LOOP(cctkGH, grouptype) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, grouptype) {

        ivect size(1);
        int const gpdim = groupdata.AT(group).info.dim;
        for (int d = 0; d < gpdim; ++d) {
          size[d] = groupdata.AT(group).info.ash[d];
        }
        int const np = prod(size);

        for (int var = 0; var < nvar; ++var) {
          int const n = n0 + var;
          for (int tl = min_tl; tl <= max_tl; ++tl) {
            memset(cctkGH->data[n][tl], poison_value, size_t(np) * sz);
          } // for tl
        }   // for var
      }
      END_LOCAL_COMPONENT_LOOP;
    }
    END_LOCAL_MAP_LOOP;

    // mark as invalid
    int const rl = grouptype == CCTK_GF ? reflevel : 0;
    for (int var = 0; var < nvar; ++var) {
      int const map0 = 0;
      ggf *const ff = arrdata.AT(group).AT(map0).data.AT(var);
      assert(ff);
      for (int tl = min_tl; tl <= max_tl; ++tl) {
        ff->set_valid(mglevel, rl, tl, CCTK_VALID_NOWHERE);
      }
    }

  } // if tl
}

void PoisonCheck(cGH const *const cctkGH, checktimes const where) {
  DECLARE_CCTK_PARAMETERS;

  if (not check_for_poison)
    return;

  Checkpoint("PoisonCheck");
  Timers::Timer timer("PoisonCheck");
  timer.start();

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    int const nvar = CCTK_NumVarsInGroupI(group);
    if (nvar > 0 and CCTK_QueryGroupStorageI(cctkGH, group)) {

      int const grouptype = CCTK_GroupTypeI(group);
      int const n0 = CCTK_FirstVarIndexI(group);
      assert(n0 >= 0);
      int const tp = CCTK_VarTypeI(n0);
      int const gpdim = groupdata.AT(group).info.dim;

      int const table = CCTK_GroupTagsTableI(group);
      assert(table >= 0);
      bool persistent;
      char buf[100];
      int const ilen =
          Util_TableGetString(table, sizeof buf, buf, "Persistent");
      if (ilen > 0) {
        if (CCTK_EQUALS(buf, "yes")) {
          persistent = true;
        } else if (CCTK_EQUALS(buf, "no")) {
          persistent = false;
        } else {
          assert(0);
        }
      } else if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        // default
        persistent = true;
      } else {
        assert(0);
      }

      CCTK_INT const poison_value = get_poison_value();

      int const num_tl = CCTK_ActiveTimeLevelsVI(cctkGH, n0);
      assert(num_tl > 0);
      int const min_tl = min_timelevel(where, num_tl, persistent);
      int const max_tl = max_timelevel(where, num_tl, persistent);

      BEGIN_LOCAL_MAP_LOOP(cctkGH, grouptype) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, grouptype) {

          ivect size(1), asize(1);
          for (int d = 0; d < gpdim; ++d) {
            size[d] = groupdata.AT(group).info.lsh[d];
            asize[d] = groupdata.AT(group).info.ash[d];
          }
          int const np = prod(size);

          for (int var = 0; var < nvar; ++var) {
            int const n = n0 + var;

            for (int tl = min_tl; tl <= max_tl; ++tl) {

              const void *const data = cctkGH->data[n][tl];
              int numpoison = 0;
              for (int k = 0; k < size[2]; ++k) {
                for (int j = 0; j < size[1]; ++j) {
                  for (int i = 0; i < size[0]; ++i) {
                    int const idx = i + asize[0] * (j + asize[1] * k);
                    bool poisoned = false;
                    switch (specific_cactus_type(tp)) {
// cast to (void *) below to avoid warning about CCTK_COMPLEX being a
// non-trivial type. C++11 guarantees that std::complex<double> has the same
// memory layout as C99's complex type, which is (likely) a trivial type in
// memory.
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T worm;                                                                    \
    memset(static_cast<void *>(&worm), poison_value, sizeof worm);             \
    const T &val = ((const T *)data)[idx];                                     \
    poisoned = memcmp(&worm, &val, sizeof worm) == 0;                          \
    break;                                                                     \
  }
#include "typecase.hh"
#undef TYPECASE
                    default:
                      UnsupportedVarType(n);
                    }
                    if (poisoned) {
                      ++numpoison;
                      if (max_poison_locations == -1 or
                          numpoison <= max_poison_locations) {
                        char *fullname = CCTK_FullName(n);
                        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                                   "At iteration %d: timelevel %d, component "
                                   "%d, map %d, refinement level %d of the "
                                   "variable \"%s\" contains poison at "
                                   "[%d,%d,%d]",
                                   cctkGH->cctk_iteration, tl, component, map,
                                   reflevel, fullname, i, j, k);
                        free(fullname);
                      }
                    } // if poisoned
                  }   // for i
                }     // for j
              }       // for k
              if (max_poison_locations != -1 and
                  numpoison > max_poison_locations) {
                char *fullname = CCTK_FullName(n);
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "At iteration %d: timelevel %d, component %d, map "
                           "%d, refinement level %d of the variable \"%s\" "
                           "contains poison at %d of %d locations; not all "
                           "locations were printed",
                           cctkGH->cctk_iteration, tl, component, map, reflevel,
                           fullname, numpoison, np);
                free(fullname);
              } else if (numpoison > 0) {
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "Found poison at %d of %d locations", numpoison, np);
              }

            } // for tl
          }   // for var
        }
        END_LOCAL_COMPONENT_LOOP;
      }
      END_LOCAL_MAP_LOOP;

    } // if has storage
  }   // for group
  timer.stop();
}


namespace {

struct taper_shell_state {
  int shell;
  bool refilled;
};

// One counter per multigrid level and refinement level. A parent
// boundary prolongation clears the counter. Each integrator substep
// that does not refill advances it by one ghost-width shell.
taper_shell_state &taper_state(int const mg, int const rl) {
  static vector<vector<taper_shell_state>> states;
  if (int(states.size()) <= mg)
    states.resize(size_t(mg) + 1);
  if (int(states.at(mg).size()) <= rl) {
    taper_shell_state const fresh = {0, false};
    states.at(mg).resize(size_t(rl) + 1, fresh);
  }
  return states.at(mg).at(rl);
}

// Shell 0 is the ghost layer that touches a buffer. Later shells step
// inward by one ghost width. Those are the zones a tapered step must
// not read again until a parent prolongation refills them.
ibset taper_star_shell(int const shell) {
  assert(shell >= 0);
  dh const &dd = *vdd.AT(map);
  dh::light_dboxes const &light =
      dd.light_boxes.AT(mglevel).AT(reflevel).AT(component);
  dh::local_dboxes const &local =
      dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
  i2vect const &width = dd.ghost_widths.AT(reflevel);
  ibset const ghosts = ibset(light.exterior) - ibset(light.interior);
  if (shell == 0)
    return ghosts & local.buffers.expand(width);
  ibset grown = ghosts;
  for (int step = 0; step < shell - 1; ++step)
    grown = grown.expand(width);
  return local.buffers & (grown.expand(width) - grown);
}

bool group_uses_prolongation(int const group) {
  char *const name = CCTK_GroupName(group);
  int const table = CCTK_GroupTagsTableI(group);
  if (table < 0) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Group \"%s\" has no tag table", name ? name : "");
    free(name);
    std::abort();
  }
  char buf[200];
  int const ilen =
      Util_TableGetString(table, int(sizeof buf), buf, "Prolongation");
  if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    free(name);
    return true;
  }
  if (ilen < 0 or ilen >= int(sizeof buf)) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot read the Prolongation tag of group \"%s\"",
                name ? name : "");
    free(name);
    std::abort();
  }
  bool const prolongate = not CCTK_EQUALS(buf, "None");
  free(name);
  return prolongate;
}

void poison_star_shell(cGH const *const cctkGH, ibset const &region) {
  if (region.empty())
    return;
  CCTK_REAL const nval = numeric_limits<CCTK_REAL>::quiet_NaN();
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_GroupTypeI(group) != CCTK_GF)
      continue;
    if (not CCTK_QueryGroupStorageI(cctkGH, group))
      continue;
    if (not group_uses_prolongation(group))
      continue;
    int const nvar = CCTK_NumVarsInGroupI(group);
    int const n0 = CCTK_FirstVarIndexI(group);
    if (nvar <= 0 or n0 < 0)
      continue;
    for (int var = 0; var < nvar; ++var) {
      int const n = n0 + var;
      if (CCTK_VarTypeI(n) != CCTK_VARIABLE_REAL)
        continue;
      // Same time level ProlongateGroupBoundaries writes.
      int const active_tl = CCTK_ActiveTimeLevelsVI(cctkGH, n);
      if (active_tl <= 0)
        continue;
      int const tl = active_tl > 1 ? timelevel : 0;
      if (tl < 0 or tl >= active_tl)
        continue;
      void *const raw = cctkGH->data[n][tl];
      if (raw == 0) {
        char const *const name = CCTK_VarName(n);
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Variable \"%s\" has no data pointer at time level %d",
                    name ? name : "", tl);
        std::abort();
      }
      CCTK_REAL *const data = static_cast<CCTK_REAL *>(raw);
      LOOP_OVER_BSET(cctkGH, region, box, imin, imax) {
        CCTK_LOOP3(poison_taper_shell, i, j, k, imin[0], imin[1], imin[2],
                   imax[0], imax[1], imax[2], cctkGH->cctk_ash[0],
                   cctkGH->cctk_ash[1], cctkGH->cctk_ash[2]) {
          int const idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
          data[idx] = nval;
        }
        CCTK_ENDLOOP3(poison_taper_shell);
      }
      END_LOOP_OVER_BSET;
    }
  }
}

} // namespace

void NoteTaperHaloRefilled() {
  if (mglevel < 0 or reflevel < 0)
    return;
  taper_state(mglevel, reflevel).refilled = true;
}

extern "C" void CarpetPoisonTaper(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (not poison_taper)
    return;
  if (reflevel <= 0)
    return;
  if (cctk_dim != 3) {
    CCTK_ERROR("poison_taper supports three dimensions");
  }

  taper_shell_state &state = taper_state(mglevel, reflevel);
  if (state.refilled) {
    state.shell = 0;
    state.refilled = false;
  }
  // The last substep of a non-tapered step prolongates after this
  // routine. That refill replaces every shell, including the one this
  // substep just consumed.
  if (do_prolongate and not do_taper) {
    state.shell = 0;
    return;
  }

  int const shell = state.shell;
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      ibset const region = taper_star_shell(shell);
      poison_star_shell(cctkGH, region);
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_LOCAL_MAP_LOOP;
  state.shell = shell + 1;
}

} // namespace Carpet
