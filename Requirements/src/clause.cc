#include "clause.hh"
#include "util.hh"

#include <cctk.h>
#include <cctki_Schedule.h>

#include <cstdlib>
#include <cstring>
#include <vector>

using namespace std;

namespace Requirements {

// template ostream& output (ostream& os, const vector<clause_t>& v);

void clause_t::interpret_options(cFunctionData const *const function_data) {
  if (function_data->meta or function_data->meta_early or
      function_data->meta_late or function_data->global or
      function_data->global_early or function_data->global_late) {
    assert(not all_reflevels);
    all_reflevels = true;
  }
  if (function_data->level) {
    assert(not all_maps);
    all_maps = true;
  }
  // Ignore singlemap and local options
  // Ignore loop_* options
}

void clause_t::parse_clause(char const *const clause1) {
  char *const clause = strdup(clause1);
  char *p = clause;

  // Remove trailing "(...)" modifier, if any
  p = strchr(p, '(');
  if (p)
    *p = '\0';
  int const gi = CCTK_GroupIndex(clause);
  if (gi >= 0) {
    // A group
    int const v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    int const nv = CCTK_NumVarsInGroupI(gi);
    assert(nv >= 0);
    for (int vi = v0; vi < v0 + nv; ++vi) {
      vars.push_back(vi);
    }
  } else {
    // Not a group - should be a variable
    int const vi = CCTK_VarIndex(clause);
    if (vi < 0) {
      CCTK_VWarn(
          CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
          "could not obtain variable/group index for '%s' in clause '%s': %d",
          clause, clause1, vi);
    }
    assert(vi >= 0);
    vars.push_back(vi);
  }

  // Parse modifiers
  // TODO: Use CarpetLib parser for this
  // TODO: add user friendly error messages
  // TODO: teach the flesh about commas within the READS/WRITES block
  if (p) {
    ++p;
    for (;;) {
      size_t const len = strcspn(p, ";)");
      char const c = p[len];
      assert(c);
      p[len] = '\0';
      if (CCTK_EQUALS(p, "everywhere")) {
        assert(not everywhere and not interior and not boundary and
               not boundary_ghostzones);
        everywhere = true;
      } else if (CCTK_EQUALS(p, "interior")) {
        assert(not everywhere and not interior);
        interior = true;
      } else if (CCTK_EQUALS(p, "boundary")) {
        assert(not everywhere and not boundary);
        boundary = true;
      } else if (CCTK_EQUALS(p, "boundary_ghostzones")) {
        assert(not everywhere and not boundary_ghostzones);
        boundary_ghostzones = true;
      } else if (CCTK_EQUALS(p, "timelevel0")) {
        assert(not timelevel0 and not all_timelevels);
        timelevel0 = true;
      } else if (CCTK_EQUALS(p, "timelevel1")) {
        assert(not timelevel1 and not all_timelevels);
        timelevel1 = true;
      } else if (CCTK_EQUALS(p, "timelevel2")) {
        assert(not timelevel2 and not all_timelevels);
        timelevel2 = true;
      } else if (CCTK_EQUALS(p, "all_timelevels")) {
        // TODO: look at schedule group instead
        assert(not timelevel0 and not timelevel1 and not timelevel2 and
               not all_timelevels);
        all_timelevels = true;
      } else {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unknown modifier '%s' in clause '%s'", p, clause1);
      }
      if (c == ')')
        break;
      assert(c == ';');
      p += len + 1;
    }
  }

  free(clause);
}

int clause_t::min_num_timelevels() const {
  if (timelevel2)
    return 3;
  if (timelevel1)
    return 2;
  return 1;
}

bool clause_t::active_on_timelevel(int const tl) const {
  assert(tl >= 0);
  if (all_timelevels)
    return true;
  if (timelevel0 and tl == 0)
    return true;
  if (timelevel1 and tl == 1)
    return true;
  if (timelevel2 and tl == 2)
    return true;
  bool const no_timelevel_clause =
      not timelevel0 and not timelevel1 and not timelevel2;
  if (no_timelevel_clause and tl == 0)
    return true;
  return false;
}

void clause_t::output(ostream &os) const {
  char *const groupname = CCTK_GroupNameFromVarI(vars.AT(0));
  os << groupname;
  free(groupname);
  os << "{";
  for (vector<int>::const_iterator ivi = vars.begin(); ivi != vars.end();
       ++ivi) {
    if (ivi != vars.begin())
      os << ",";
    char *const fullname = CCTK_FullName(*ivi);
    os << fullname;
    free(fullname);
  }
  os << "}(";
  if (everywhere)
    os << "everywhere;";
  if (interior)
    os << "interior;";
  if (boundary)
    os << "boundary;";
  if (boundary_ghostzones)
    os << "boundary_ghostzones;";
  if (timelevel0)
    os << "timelevel0;";
  if (timelevel1)
    os << "timelevel1;";
  if (timelevel2)
    os << "timelevel2;";
  if (all_timelevels)
    os << "all_timelevels;";
  if (all_maps)
    os << "all_maps;";
  if (all_reflevels)
    os << "all_reflevels;";
  os << ")";
}
};
