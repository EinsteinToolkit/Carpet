#include "all_state.hh"
#include "all_clauses.hh"
#include "clauses.hh"
#include "util.hh"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cstdlib>

#include <cstdlib>

using namespace std;

namespace Requirements {

// Ignore requirements in these variables; these variables are
// always considered valid
std::vector<bool> ignored_variables;

static void add_ignored_variable(int const id, const char *const opstring,
                                 void *const callback_arg) {
  std::vector<bool> &ivs = *static_cast<std::vector<bool> *>(callback_arg);
  ivs.AT(id) = true;
}

void all_state_t::setup(int const maps) {
  DECLARE_CCTK_PARAMETERS;
  assert(vars.empty());
  vars.resize(CCTK_NumVars());
  for (variables_t::iterator ivar = vars.begin(); ivar != vars.end(); ++ivar) {
    reflevels_t &rls = *ivar;
    int const vi = &*ivar - &*vars.begin();
    assert(rls.empty());
    // Allocate one refinement level initially
    int const nrls = 1;
    rls.resize(nrls);
    for (reflevels_t::iterator irl = rls.begin(); irl != rls.end(); ++irl) {
      maps_t &ms = *irl;
      assert(ms.empty());
      int const group_type = CCTK_GroupTypeFromVarI(vi);
      int const nms = group_type == CCTK_GF ? maps : 1;
      if (verbose) {
        char *const fullname = CCTK_FullName(vi);
        int const rl = &*irl - &*rls.begin();
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Setting up %d maps for variable %s(rl=%d)", nms, fullname,
                   rl);
        free(fullname);
      }
      ms.resize(nms);
      for (maps_t::iterator im = ms.begin(); im != ms.end(); ++im) {
        timelevels_t &tls = *im;
        assert(tls.empty());
        // Not allocating any time levels here
      }
    }
  }
  assert(ignored_variables.empty());
  ignored_variables.resize(CCTK_NumVars());
  const int iret =
      CCTK_TraverseString(ignore_these_variables, add_ignored_variable,
                          (void *)&ignored_variables, CCTK_GROUP_OR_VAR);
  assert(iret >= 0);
}

// Update internal data structures when Carpet changes the number of
// active timelevels for a group
void all_state_t::change_storage(vector<int> const &groups,
                                 vector<int> const &timelevels,
                                 int const reflevel) {
  DECLARE_CCTK_PARAMETERS;
  assert(groups.size() == timelevels.size());
  for (vector<int>::const_iterator igi = groups.begin(),
                                   itl = timelevels.begin();
       igi != groups.end(); ++igi, ++itl) {
    int const gi = *igi;
    int const tl = *itl;
    bool const is_array = CCTK_GroupTypeI(gi) != CCTK_GF;
    int const v0 = CCTK_FirstVarIndexI(gi);
    int const nv = CCTK_NumVarsInGroupI(gi);
    for (int vi = v0; vi < v0 + nv; ++vi) {
      reflevels_t &rls = vars.AT(vi);
      int const reflevels = int(rls.size());
      bool const all_rl = reflevel == -1;
      int const min_rl = is_array ? 0 : all_rl ? 0 : reflevel;
      int const max_rl = is_array ? 1 : all_rl ? reflevels : reflevel + 1;
      assert(min_rl >= 0 and max_rl <= reflevels);
      for (int rl = min_rl; rl < max_rl; ++rl) {
        maps_t &ms = rls.AT(rl);
        for (maps_t::iterator im = ms.begin(); im != ms.end(); ++im) {
          timelevels_t &tls = *im;
          int const ntls = int(tls.size());
          if (tl < ntls) {
            // Free some storage
            if (verbose) {
              char *const fullname = CCTK_FullName(vi);
              int const m = &*im - &*ms.begin();
              CCTK_VInfo(CCTK_THORNSTRING, "Decreasing storage to %d time "
                                           "levels for variable %s(rl=%d,m=%d)",
                         tl, fullname, rl, m);
              free(fullname);
            }
            tls.resize(tl);
          } else if (tl > ntls) {
            // Allocate new storage
            if (verbose) {
              char *const fullname = CCTK_FullName(vi);
              int const m = &*im - &*ms.begin();
              CCTK_VInfo(CCTK_THORNSTRING, "Increasing storage to %d time "
                                           "levels for variable %s(rl=%d,m=%d)",
                         tl, fullname, rl, m);
              free(fullname);
            }
            // The default constructor for gridpoint_t sets all
            // data to "invalid"
            tls.resize(tl);
          }
        }
      }
    }
  }
}

// Update internal data structures when Carpet regrids
void all_state_t::regrid(int const reflevels) {
  DECLARE_CCTK_PARAMETERS;
  assert(old_vars.empty());
  old_vars.resize(vars.size());

  int const ng = CCTK_NumGroups();
  for (int gi = 0; gi < ng; ++gi) {
    int const group_type = CCTK_GroupTypeI(gi);
    switch (group_type) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      // Grid arrays remain unchanged
      break;
    case CCTK_GF: {
      // Only grid functions are regridded
      int const v0 = CCTK_FirstVarIndexI(gi);
      int const nv = CCTK_NumVarsInGroupI(gi);
      for (int vi = v0; vi < v0 + nv; ++vi) {
        reflevels_t &rls = vars.AT(vi);
        reflevels_t &old_rls = old_vars.AT(vi);
        assert(old_rls.empty());
        swap(rls, old_rls);
        // Delete (unused) old refinement levels
        int const old_reflevels = int(old_rls.size());
        for (int rl = reflevels; rl < old_reflevels; ++rl) {
          maps_t &old_ms = old_rls.AT(rl);
          if (verbose) {
            char *const fullname = CCTK_FullName(vi);
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Deleting unused refinement level %d of variable %s", rl,
                       fullname);
            free(fullname);
          }
          old_ms.clear();
        }
        // Allocate new refinement levels
        rls.resize(reflevels);
        maps_t const &old_ms = old_rls.AT(0);
        int const old_maps = int(old_ms.size());
        int const maps = old_maps;
        for (int rl = old_reflevels; rl < reflevels; ++rl) {
          maps_t &ms = rls.AT(rl);
          if (verbose) {
            char *const fullname = CCTK_FullName(vi);
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Allocating new refinement level %d for variable %s", rl,
                       fullname);
            free(fullname);
          }
          ms.resize(maps);
          for (maps_t::iterator im = ms.begin(); im != ms.end(); ++im) {
            int const crl = 0;
            int const m = &*im - &*ms.begin();
            timelevels_t &tls = *im;
            assert(tls.empty());
            int const ntls = int(old_rls.AT(crl).AT(m).size());
            // Allocate undefined timelevels
            tls.resize(ntls);
          }
        }
      }
      break;
    }
    default:
      assert(0);
    }
  }
}

// Update internal data structures when Carpet recomposes
void all_state_t::recompose(int const iteration, int const reflevel,
                            valid::valid_t const where) {
  DECLARE_CCTK_PARAMETERS;
  int const ng = CCTK_NumGroups();
  location_t loc("recompose");
  loc.it = iteration;
  loc.rl = reflevel;
  for (int gi = 0; gi < ng; ++gi) {
    int const group_type = CCTK_GroupTypeI(gi);
    switch (group_type) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      // Grid arrays remain unchanged
      break;
    case CCTK_GF: {
      // Only grid functions are regridded
      int const v0 = CCTK_FirstVarIndexI(gi);
      int const nv = CCTK_NumVarsInGroupI(gi);
      for (int vi = v0; vi < v0 + nv; ++vi) {
        loc.vi = vi;
        reflevels_t &rls = vars.AT(vi);
        maps_t &ms = rls.AT(reflevel);
        reflevels_t &old_rls = old_vars.AT(vi);
        int const old_reflevels = int(old_rls.size());
        if (reflevel < old_reflevels) {
          // This refinement level is regridded
          maps_t &old_ms = old_rls.AT(reflevel);
          assert(not old_ms.empty());
          assert(ms.empty());
          swap(ms, old_ms);
          for (maps_t::iterator im = ms.begin(), old_im = old_ms.begin();
               im != ms.end(); ++im, ++old_im) {
            loc.m = &*im - &*ms.begin();
            timelevels_t &tls = *im;
            if (verbose) {
              char *const fullname = CCTK_FullName(vi);
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Recomposing variable %s(rl=%d,m=%d)", fullname,
                         reflevel, loc.m);
              free(fullname);
            }
            for (timelevels_t::iterator itl = tls.begin(); itl != tls.end();
                 ++itl) {
              loc.tl = &*itl - &*tls.begin();
              gridpoint_t &gp = *itl;
              switch (where) {
              case valid::nowhere:
                gp.set_interior(false, loc);
              // fall through
              case valid::interior:
                // Recomposing sets only the interior
                gp.set_boundary(false, loc);
                gp.set_ghostzones(false, loc);
                gp.set_boundary_ghostzones(false, loc);
              // fall through
              case valid::everywhere:
                // do nothing
                break;
              default:
                assert(0);
              }
            }
          }
          assert(old_ms.empty());
        } else {
          // This refinement level is new
          assert(where == valid::nowhere);
        }
      }
      break;
    }
    default:
      assert(0);
    }
  }
}

void all_state_t::regrid_free() {
  // Ensure all old maps have been recomposed
  for (variables_t::const_iterator ivar = old_vars.begin();
       ivar != old_vars.end(); ++ivar) {
    reflevels_t const &old_rls = *ivar;
    for (reflevels_t::const_iterator irl = old_rls.begin();
         irl != old_rls.end(); ++irl) {
      maps_t const &old_ms = *irl;
      assert(old_ms.empty());
    }
  }
  old_vars.clear();
}

// Update internal data structures when Carpet cycles timelevels
void all_state_t::cycle(int const reflevel) {
  int const ng = CCTK_NumGroups();
  for (int gi = 0; gi < ng; ++gi) {
    int const group_type = CCTK_GroupTypeI(gi);
    bool do_cycle;
    switch (group_type) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      // Grid arrays are cycled in global mode
      do_cycle = reflevel == -1;
      break;
    case CCTK_GF:
      // Grid functions are cycled in level mode
      do_cycle = reflevel >= 0;
      break;
    default:
      assert(0);
    }
    if (do_cycle) {
      // Translate global mode to refinement level 0
      int const rl = reflevel >= 0 ? reflevel : 0;
      int const v0 = CCTK_FirstVarIndexI(gi);
      int const nv = CCTK_NumVarsInGroupI(gi);
      for (int vi = v0; vi < v0 + nv; ++vi) {
        reflevels_t &rls = vars.AT(vi);
        maps_t &ms = rls.AT(rl);
        for (maps_t::iterator im = ms.begin(); im != ms.end(); ++im) {
          timelevels_t &tls = *im;
          int const ntl = int(tls.size());
          if (ntl >= 1) {
            // Only cycle variables with sufficient storage
            for (int tl = ntl - 1; tl > 0; --tl) {
              tls.AT(tl) = tls.AT(tl - 1);
            }
#if 0
              // The new time level is uninitialised
              // TODO: keep it valid to save time, since MoL will
              // copy it anyway?
              tls.AT(0) = gridpoint_t();
#endif
          }
        }
      }
    }
  }
}

// Check that the grid is in the required state before a given
// function is executed
void all_state_t::before_routine(cFunctionData const *const function_data,
                                 all_clauses_t &all_clauses,
                                 int const iteration, int const reflevel,
                                 int const map, int const timelevel,
                                 int const timelevel_offset) const {
  location_t loc("routine", function_data);
  loc.it = iteration;
  // Loop over all clauses
  clauses_t const &clauses = all_clauses.get_clauses(function_data);
  for (vector<clause_t>::const_iterator iclause = clauses.reads.begin();
       iclause != clauses.reads.end(); ++iclause) {
    clause_t const &clause = *iclause;
    for (vector<int>::const_iterator ivar = clause.vars.begin();
         ivar != clause.vars.end(); ++ivar) {
      int const vi = *ivar;
      loc.vi = vi;

      if (ignored_variables.AT(vi))
        continue;

      // Loop over all (refinement levels, maps, time levels)
      reflevels_t const &rls = vars.AT(vi);
      int const reflevels = int(rls.size());
      int min_rl, max_rl;
      if (clause.all_reflevels or reflevel == -1) {
        min_rl = 0;
        max_rl = reflevels;
      } else {
        min_rl = reflevel;
        max_rl = min_rl + 1;
      }
      for (int rl = min_rl; rl < max_rl; ++rl) {
        loc.rl = rl;

        maps_t const &ms = rls.AT(rl);
        int const maps = int(ms.size());
        int min_m, max_m;
        if (clause.all_maps or map == -1) {
          min_m = 0;
          max_m = maps;
        } else {
          min_m = map;
          max_m = min_m + 1;
        }
        for (int m = min_m; m < max_m; ++m) {
          loc.m = m;

          timelevels_t const &tls = ms.AT(m);
          int const ntls = int(tls.size());
          assert(ntls >= clause.min_num_timelevels());
          assert(timelevel != -1);
          for (int tl = 0; tl < ntls; ++tl) {
            loc.tl = tl;
            if (tl >= timelevel_offset and
                clause.active_on_timelevel(tl - timelevel_offset)) {
              gridpoint_t const &gp = tls.AT(tl);
              gp.check_state(clause, loc);
            }
          }
        }
      }
    }
  }
}

// Update internal data structures after a function has been
// executed to reflect the fact that some variables are now valid
void all_state_t::after_routine(cFunctionData const *const function_data,
                                all_clauses_t &all_clauses, int const iteration,
                                int const reflevel, int const map,
                                int const timelevel,
                                int const timelevel_offset) {
  location_t loc("routine", function_data);
  loc.it = iteration;
  // Loop over all clauses
  clauses_t const &clauses = all_clauses.get_clauses(function_data);
  for (vector<clause_t>::const_iterator iclause = clauses.writes.begin();
       iclause != clauses.writes.end(); ++iclause) {
    clause_t const &clause = *iclause;
    for (vector<int>::const_iterator ivar = clause.vars.begin();
         ivar != clause.vars.end(); ++ivar) {
      int const vi = *ivar;
      loc.vi = vi;

      // Loop over all (refinement levels, maps, time levels)
      reflevels_t &rls = vars.AT(vi);
      int const reflevels = int(rls.size());
      int min_rl, max_rl;
      if (clause.all_reflevels or reflevel == -1) {
        min_rl = 0;
        max_rl = reflevels;
      } else {
        min_rl = reflevel;
        max_rl = min_rl + 1;
      }
      for (int rl = min_rl; rl < max_rl; ++rl) {
        loc.rl = rl;

        maps_t &ms = rls.AT(rl);
        int const maps = int(ms.size());
        int min_m, max_m;
        if (clause.all_maps or map == -1) {
          min_m = 0;
          max_m = maps;
        } else {
          min_m = map;
          max_m = min_m + 1;
        }
        for (int m = min_m; m < max_m; ++m) {
          loc.m = m;

          timelevels_t &tls = ms.AT(m);
          int const ntls = int(tls.size());
          assert(ntls >= clause.min_num_timelevels());
          assert(timelevel != -1);
          for (int tl = 0; tl < ntls; ++tl) {
            loc.tl = tl;
            if (tl >= timelevel_offset and
                clause.active_on_timelevel(tl - timelevel_offset)) {
              gridpoint_t &gp = tls.AT(tl);
              // TODO: If this variable is both read and written
              // (i.e. if this is a projection), then only the
              // written region remains valid
              gp.update_state(clause, loc);
            }
          }
        }
      }
    }
  }
}

// Update internal data structures when Carpet syncs
void all_state_t::sync(cFunctionData const *const function_data,
                       vector<int> const &groups, int const iteration,
                       int const reflevel, int const timelevel) {
  location_t loc("synchronisation", function_data);
  loc.it = iteration;
  loc.rl = reflevel;
  loc.tl = timelevel;
  // Loop over all variables
  for (vector<int>::const_iterator igi = groups.begin(); igi != groups.end();
       ++igi) {
    int const gi = *igi;
    bool do_sync;
    int const group_type = CCTK_GroupTypeI(gi);
    switch (group_type) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      // Grid arrays are synced in global mode
      do_sync = reflevel == -1;
      break;
    case CCTK_GF:
      // Grid functions are synced in level mode
      do_sync = reflevel >= 0;
      break;
    default:
      assert(0);
    }
    if (do_sync) {
      // Translate global mode to refinement level 0
      int const rl = reflevel >= 0 ? reflevel : 0;
      int const v0 = CCTK_FirstVarIndexI(gi);
      int const nv = CCTK_NumVarsInGroupI(gi);
      for (int vi = v0; vi < v0 + nv; ++vi) {
        if (ignored_variables.AT(vi))
          continue;
        loc.vi = vi;

        reflevels_t &rls = vars.AT(vi);
        maps_t &ms = rls.AT(rl);
        int const maps = int(ms.size());
        for (int m = 0; m < maps; ++m) {
          loc.m = m;
          timelevels_t &tls = ms.AT(m);
          int const tl = timelevel;
          gridpoint_t &gp = tls.AT(tl);

          // Synchronising requires a valid interior
          if (not gp.interior()) {
            gp.report_error(loc, "interior");
          }

          // Synchronising (i.e. prolongating) requires valid data
          // on all time levels of the same map of the next
          // coarser refinement level
          if (rl > 0) {
            location_t cloc(loc);
            cloc.info = "prolongation";
            int const crl = rl - 1;
            cloc.rl = crl;
            maps_t const &cms = rls.AT(crl);
            timelevels_t const &ctls = cms.AT(m);
            // TODO: use prolongation_order_time instead?
            int const ctimelevels = int(ctls.size());
            for (int ctl = 0; ctl < ctimelevels; ++ctl) {
              cloc.tl = ctl;
              gridpoint_t const &cgp = ctls.AT(ctl);
              if (not(cgp.interior() and cgp.boundary() and cgp.ghostzones() and
                      cgp.boundary_ghostzones())) {
                cgp.report_error(cloc, "everywhere");
              }
            }
          }

          // Synchronising sets all ghost zones, and sets boundary
          // ghost zones if boundary zones are set
          if (gp.boundary()) {
            if (gp.ghostzones() and gp.boundary_ghostzones()) {
              gp.report_warning(loc, "ghostzones+boundary_ghostzones");
            }
          } else {
            if (gp.ghostzones()) {
              gp.report_warning(loc, "ghostzones");
            }
          }
          gp.set_ghostzones(true, loc);
          gp.set_boundary_ghostzones(gp.boundary(), loc);
        }
      }
    }
  }
}

// Update internal data structures when Carpet restricts
void all_state_t::restrict1(vector<int> const &groups, int const iteration,
                            int const reflevel) {
  location_t loc("restriction");
  loc.it = iteration;
  loc.rl = reflevel;
  // Loop over all variables
  for (vector<int>::const_iterator igi = groups.begin(); igi != groups.end();
       ++igi) {
    int const gi = *igi;
    bool do_restrict;
    int const group_type = CCTK_GroupTypeI(gi);
    switch (group_type) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      // Grid arrays are synced in global mode
      do_restrict = reflevel == -1;
      break;
    case CCTK_GF:
      // Grid functions are synced in level mode
      do_restrict = reflevel >= 0;
      break;
    default:
      assert(0);
    }
    if (do_restrict) {
      // Translate global mode to refinement level 0
      int const rl = reflevel >= 0 ? reflevel : 0;
      int const v0 = CCTK_FirstVarIndexI(gi);
      int const nv = CCTK_NumVarsInGroupI(gi);
      for (int vi = v0; vi < v0 + nv; ++vi) {
        loc.vi = vi;
        if (ignored_variables.AT(vi))
          continue;

        reflevels_t &rls = vars.AT(vi);
        int const reflevels = int(rls.size());
        maps_t &ms = rls.AT(rl);
        int const maps = int(ms.size());
        for (int m = 0; m < maps; ++m) {
          loc.m = m;
          timelevels_t &tls = ms.AT(m);
          int const tl = 0;
          loc.tl = tl;
          gridpoint_t &gp = tls.AT(tl);

          // Restricting requires a valid interior (otherwise we
          // cannot be sure that all of the interior is valid
          // afterwards)
          if (not gp.interior()) {
            gp.report_error(loc, "interior");
          }

          // Restricting requires valid data on the current time
          // level of the same map of the next finer refinement
          // level
          if (rl < reflevels - 1) {
            int const frl = rl + 1;
            maps_t const &fms = rls.AT(frl);
            timelevels_t const &ftls = fms.AT(m);
            int const ftl = 0;
            location_t floc(loc);
            floc.rl = frl;
            floc.tl = ftl;
            gridpoint_t const &fgp = ftls.AT(ftl);
            if (not(fgp.interior() and fgp.boundary() and fgp.ghostzones() and
                    fgp.boundary_ghostzones())) {
              fgp.report_error(floc, "everywhere");
            }
          }

          // Restricting fills (part of) the interior, but leaves
          // ghost zones and boundary zones undefined
          gp.set_boundary(false, loc);
          gp.set_ghostzones(false, loc);
          gp.set_boundary_ghostzones(false, loc);
        }
      }
    }
  }
}

void all_state_t::output(ostream &os) const {
  os << "all_state:" << std::endl;
  os << "vars:" << std::endl;
  os << vars << std::endl;
  os << "old_vars:" << std::endl;
  os << old_vars << std::endl;
}

template ostream &output(ostream &os,
                         const vector<all_state_t::timelevels_t> &v);
template ostream &output(ostream &os, const vector<all_state_t::maps_t> &v);
template ostream &output(ostream &os,
                         const vector<all_state_t::reflevels_t> &v);
template ostream &output(ostream &os,
                         const vector<all_state_t::variables_t> &v);

void all_state_t::invalidate(vector<int> const &vars1, int const reflevel,
                             int const map, int const timelevel) {
  // Loop over all variables
  for (vector<int>::const_iterator ivi = vars1.begin(); ivi != vars1.end();
       ++ivi) {
    int const vi = *ivi;
    reflevels_t &rls = vars.AT(vi);
    maps_t &ms = rls.AT(reflevel);
    timelevels_t &tls = ms.AT(map);
    // This time level is uninitialised
    tls.AT(timelevel) = gridpoint_t();
  }
}
}
