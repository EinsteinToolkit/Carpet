#include <cmath>
#include <cstdlib>
#include <vector>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <ggf.hh>

#include <carpet.hh>
#include "functions.hh"
#include "hermite.hh"
#include "hermite_fill.hh"

namespace Carpet {

using namespace std;

bool hermite_prolongate = false;

// C linkage, outside the anonymous namespace, so these names are MoL's
// globals. A weak declaration stays unresolved when MoL is not linked.
extern "C" {
extern CCTK_INT *EvolvedVariableIndex __attribute__((weak));
extern CCTK_INT *RHSVariableIndex __attribute__((weak));
extern CCTK_INT MoLNumEvolvedVariables __attribute__((weak));
extern CCTK_INT MoLNumEvolvedVariablesSlow __attribute__((weak));
extern CCTK_INT MoLNumEvolvedArrayVariables __attribute__((weak));
extern CCTK_INT MoL_MSRK_Active(void) __attribute__((weak));
}

namespace {

struct Component {
  int lc;
  gdata *saved;
  gdata *scratch;
};

struct Slot {
  ggf *evolved;
  int ml;
  int parent_rl;
  int src_tl;
  vector<Component> components;
};

bool prolongates(operator_type const op) {
  return op != op_none and op != op_sync and op != op_restrict;
}

void require_mol_tables() {
  if (&EvolvedVariableIndex == nullptr or EvolvedVariableIndex == nullptr or
      &RHSVariableIndex == nullptr or RHSVariableIndex == nullptr or
      &MoLNumEvolvedVariables == nullptr or
      &MoLNumEvolvedVariablesSlow == nullptr or
      &MoLNumEvolvedArrayVariables == nullptr) {
    CCTK_ERROR("Hermite dense output needs MoL's evolved-variable tables.");
  }
  if (MoLNumEvolvedArrayVariables > 0) {
    CCTK_ERROR("Hermite dense output does not support evolved grid arrays.");
  }
  if (MoLNumEvolvedVariablesSlow > 0) {
    CCTK_ERROR("Hermite dense output does not support slow MoL variables.");
  }
  if (MoLNumEvolvedVariables < 0) {
    CCTK_ERROR("Hermite dense output saw a negative evolved-variable count.");
  }
}

int rhs_for(int const evolved) {
  for (int i = 0; i < MoLNumEvolvedVariables; ++i) {
    if (EvolvedVariableIndex[i] == evolved)
      return RHSVariableIndex[i];
  }
  return -1;
}

int field_offset(gdata const *const field, int const i, int const j,
                 int const k) {
  return i * field->stride()[0] + j * field->stride()[1] +
         k * field->stride()[2];
}

void require_same_layout(gdata const *const field, gdata const *const proto,
                         char const *const what) {
  if (field->elementsize() != proto->elementsize() or
      any(field->shape() != proto->shape()) or
      any(field->stride() != proto->stride())) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output found a %s grid whose layout does not "
                "match the evolved variable.",
                what);
  }
}

void fill_component(gdata *const dst, gdata const *const y_start,
                    gdata const *const y_end, gdata const *const dy_start,
                    gdata const *const dy_end, CCTK_REAL const theta,
                    CCTK_REAL const dt) {
  if (dst->elementsize() != static_cast<int>(sizeof(CCTK_REAL))) {
    CCTK_ERROR("Hermite dense output only supports real grid functions.");
  }
  require_same_layout(y_start, dst, "y(t_n)");
  require_same_layout(y_end, dst, "y(t_n+dt)");
  require_same_layout(dy_start, dst, "y'(t_n)");
  require_same_layout(dy_end, dst, "y'(t_n+dt)");
  hermite::Weights const w = hermite::weights(theta);
  CCTK_REAL const *const ys =
      static_cast<CCTK_REAL const *>(y_start->storage());
  CCTK_REAL const *const ye = static_cast<CCTK_REAL const *>(y_end->storage());
  CCTK_REAL const *const dys =
      static_cast<CCTK_REAL const *>(dy_start->storage());
  CCTK_REAL const *const dye =
      static_cast<CCTK_REAL const *>(dy_end->storage());
  CCTK_REAL *const out = static_cast<CCTK_REAL *>(dst->storage());
  ivect const shape = dst->shape();
  for (int k = 0; k < shape[2]; ++k) {
    for (int j = 0; j < shape[1]; ++j) {
      for (int i = 0; i < shape[0]; ++i) {
        int const off = field_offset(dst, i, j, k);
        out[off] = w.value_start * ys[off] + w.value_end * ye[off] +
                   dt * (w.slope_start * dys[off] + w.slope_end * dye[off]);
      }
    }
  }
}

ggf *variable_ggf(int const vi, int const map) {
  int const group = CCTK_GroupIndexFromVarI(vi);
  if (group < 0) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output cannot find the group of variable %d.",
                vi);
  }
  int const first = CCTK_FirstVarIndexI(group);
  int const var = vi - first;
  if (map < 0 or map >= (int)arrdata.AT(group).size() or var < 0 or
      var >= (int)arrdata.AT(group).AT(map).data.size()) {
    char *const name = CCTK_FullName(vi);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output has no grid for %s on map %d.",
                name ? name : "<unknown>", map);
    free(name);
  }
  return arrdata.AT(group).AT(map).data.AT(var);
}

bool multistep_slopes() {
  return &MoL_MSRK_Active != nullptr and MoL_MSRK_Active() != 0;
}

void prepare_variable(vector<Slot> &slots, ggf *const evolved, int const rhs,
                      int const map, CCTK_REAL const theta, CCTK_REAL const dt,
                      bool const multistep) {
  for (size_t s = 0; s < slots.size(); ++s) {
    if (slots.AT(s).evolved == evolved)
      return;
  }
  if (evolved->vectorlength != 1) {
    CCTK_ERROR("Hermite dense output does not handle vectorized grid functions.");
  }
  if (not prolongates(evolved->transport_operator)) {
    char *const name = CCTK_FullName(evolved->varindex);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output cannot fill %s because its prolongation "
                "operator does not prolongate boundaries.",
                name ? name : "<unknown>");
    free(name);
  }
  int const parent = reflevel - 1;
  int const ntl = evolved->timelevels(mglevel, parent);
  if (ntl < 3) {
    char *const name = CCTK_FullName(evolved->varindex);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output needs three time levels of %s on the "
                "parent to keep y(t_n) and y(t_n+dt) while it builds the "
                "prolongation field. Found %d.",
                name ? name : "<unknown>", ntl);
    free(name);
  }
  ggf *const rhs_gv = variable_ggf(rhs, map);
  int const rhs_levels = rhs_gv->timelevels(mglevel, parent);
  int const rhs_needed = multistep ? 3 : 2;
  if (rhs_levels < rhs_needed) {
    char *const name = CCTK_FullName(rhs);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output needs %d time levels of %s. Found %d. "
                "Set MoL::hermite_rhs_storage = yes. Multistep RK4 also "
                "needs MoL::msrk_rhs_storage = yes.",
                rhs_needed, name ? name : "<unknown>", rhs_levels);
    free(name);
  }
  hermite::SlopeLevels const slopes = hermite::slope_levels(multistep);

  Slot slot;
  slot.evolved = evolved;
  slot.ml = mglevel;
  slot.parent_rl = parent;
  slot.src_tl = ntl - 1;
  int const nlc = evolved->h.local_components(parent);
  slot.components.reserve(nlc);
  for (int lc = 0; lc < nlc; ++lc) {
    gdata *const y_end = evolved->data_pointer(0, parent, lc, mglevel);
    gdata *const y_start = evolved->data_pointer(1, parent, lc, mglevel);
    gdata *const dy_end =
        rhs_gv->data_pointer(slopes.end, parent, lc, mglevel);
    gdata *const dy_start =
        rhs_gv->data_pointer(slopes.start, parent, lc, mglevel);
    if (y_end == NULL or y_start == NULL or dy_end == NULL or
        dy_start == NULL or not y_end->has_storage() or
        not y_start->has_storage() or not dy_end->has_storage() or
        not dy_start->has_storage()) {
      CCTK_ERROR("Hermite dense output found a parent grid with no storage.");
    }
    gdata *const scratch = evolved->new_typed_data();
    scratch->allocate(y_end->extent(), y_end->padded_shape(),
                      y_end->padding_offset(), y_end->proc());
    fill_component(scratch, y_start, y_end, dy_start, dy_end, theta, dt);
    Component component;
    component.lc = lc;
    component.scratch = scratch;
    component.saved = evolved->exchange_data(slot.src_tl, parent, lc, mglevel,
                                             scratch);
    slot.components.push_back(component);
  }
  slots.push_back(slot);
}

vector<int> evolved_groups() {
  require_mol_tables();
  vector<int> groups;
  for (int i = 0; i < MoLNumEvolvedVariables; ++i) {
    int const vi = EvolvedVariableIndex[i];
    int const group = CCTK_GroupIndexFromVarI(vi);
    if (group < 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Hermite dense output cannot find the group of evolved "
                  "variable %d.",
                  vi);
    }
    bool seen = false;
    for (size_t g = 0; g < groups.size(); ++g) {
      if (groups.AT(g) == group)
        seen = true;
    }
    if (not seen)
      groups.push_back(group);
  }
  return groups;
}

} // namespace

struct HermiteFill::State {
  vector<Slot> slots;
};

HermiteFill::HermiteFill(cGH const *const cctkGH, vector<int> const &groups,
                         CCTK_REAL const time)
    : state_(new State) {
  if (not hermite_prolongate)
    return;
  if (reflevel <= 0)
    return;
  DECLARE_CCTK_PARAMETERS;
  if (not hermite_dense_output)
    return;

  require_mol_tables();
  bool any_evolved = false;
  for (size_t g = 0; g < groups.size(); ++g) {
    int const group = groups.AT(g);
    if (CCTK_GroupTypeI(group) != CCTK_GF)
      continue;
    int const first = CCTK_FirstVarIndexI(group);
    int const nvars = CCTK_NumVarsInGroupI(group);
    for (int v = 0; v < nvars; ++v) {
      if (rhs_for(first + v) >= 0)
        any_evolved = true;
    }
  }
  if (not any_evolved)
    return;

  if (tt == NULL)
    CCTK_ERROR("Hermite dense output has no Carpet time hierarchy.");
  int const parent = reflevel - 1;
  CCTK_REAL const t_end = tt->get_time(mglevel, parent, 0);
  CCTK_REAL const t_start = tt->get_time(mglevel, parent, 1);
  CCTK_REAL const dt = t_end - t_start;
  if (not(dt > 0.0)) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output needs a positive parent step. "
                "Level %d has t_n = %.17g and t_n+dt = %.17g.",
                parent, static_cast<double>(t_start),
                static_cast<double>(t_end));
  }
  bool const multistep = multistep_slopes();
  CCTK_REAL const tolerance =
      1.0e-10 * (1.0 + fabs(t_end) + fabs(t_start));
  bool const outside =
      time < t_start - tolerance or time > t_end + tolerance;
  if (outside and not multistep) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite dense output was asked for time %.17g outside the "
                "parent interval [%.17g, %.17g] on level %d.",
                static_cast<double>(time), static_cast<double>(t_start),
                static_cast<double>(t_end), parent);
  }
  if (not std::isfinite(static_cast<double>(time))) {
    CCTK_ERROR("Hermite dense output was asked for a non-finite time.");
  }
  CCTK_REAL const theta = (time - t_start) / dt;
  static bool announced_outside = false;
  if (outside and not announced_outside) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Hermite dense output is evaluating the cubic outside the "
               "parent step for a multistep RK4 stage.");
    announced_outside = true;
  }

  static bool announced = false;
  if (not announced) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Hermite dense output is filling fine ghost zones.");
    announced = true;
  }

  for (size_t g = 0; g < groups.size(); ++g) {
    int const group = groups.AT(g);
    if (CCTK_GroupTypeI(group) != CCTK_GF)
      continue;
    int const first = CCTK_FirstVarIndexI(group);
    int const nvars = CCTK_NumVarsInGroupI(group);
    for (int map = 0; map < (int)arrdata.AT(group).size(); ++map) {
      if ((int)arrdata.AT(group).AT(map).data.size() != nvars) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Hermite dense output found %d grids for group %d, which "
                    "has %d variables.",
                    (int)arrdata.AT(group).AT(map).data.size(), group, nvars);
      }
      for (int v = 0; v < nvars; ++v) {
        ggf *const evolved = arrdata.AT(group).AT(map).data.AT(v);
        int const rhs = rhs_for(first + v);
        if (rhs < 0)
          continue;
        prepare_variable(state_->slots, evolved, rhs, map, theta, dt,
                         multistep);
      }
    }
  }
}

HermiteFill::~HermiteFill() {
  for (size_t s = 0; s < state_->slots.size(); ++s) {
    Slot &slot = state_->slots.AT(s);
    for (size_t c = 0; c < slot.components.size(); ++c) {
      Component &component = slot.components.AT(c);
      gdata *const back = slot.evolved->exchange_data(
          slot.src_tl, slot.parent_rl, component.lc, slot.ml, component.saved);
      if (back != component.scratch) {
        CCTK_ERROR("Hermite dense output lost its scratch grid.");
      }
      delete component.scratch;
      component.scratch = NULL;
    }
  }
  delete state_;
  state_ = NULL;
}

bool HermiteFill::source_timelevel(ggf const *const evolved,
                                   int &src_tl) const {
  for (size_t s = 0; s < state_->slots.size(); ++s) {
    if (state_->slots.AT(s).evolved == evolved) {
      src_tl = state_->slots.AT(s).src_tl;
      return true;
    }
  }
  return false;
}

// The first right-hand side of a CSAMR step is evaluated at the step
// start. CallEvol has already stored the step end in cctk_time, and
// SyncGroupsByDirI has no time argument, so the fill temporarily uses
// the start. The level dt is delta_time / timereffacts[reflevel]. In
// level mode that equals cctk_delta_time / cctk_timefac: global mode
// stores the coarse dt, and level mode stores the level's time factor.
void Carpet_HermitePreStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (not hermite_dense_output or not use_csamr or not hermite_prolongate)
    return;
  if (reflevel <= 0 or do_taper)
    return;

  if (timereffacts.AT(reflevel) <= 0 or cctkGH->cctk_timefac <= 0) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite pre-step on level %d has time factor %d and "
                "cctk_timefac %.17g.",
                reflevel, timereffacts.AT(reflevel),
                static_cast<double>(cctkGH->cctk_timefac));
  }
  CCTK_REAL const level_dt =
      delta_time / static_cast<CCTK_REAL>(timereffacts.AT(reflevel));
  CCTK_REAL const gh_dt = cctkGH->cctk_delta_time /
                          static_cast<CCTK_REAL>(cctkGH->cctk_timefac);
  CCTK_REAL const tolerance = 1.0e-10 * (1.0 + fabs(level_dt));
  if (not(level_dt > 0.0) or fabs(gh_dt - level_dt) > tolerance) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Hermite pre-step on level %d needs one positive level dt. "
                "delta_time / time factor is %.17g and cctk_delta_time / "
                "cctk_timefac is %.17g.",
                reflevel, static_cast<double>(level_dt),
                static_cast<double>(gh_dt));
  }

  vector<int> const groups = evolved_groups();
  if (groups.empty())
    return;
  CCTK_REAL const step_end = cctkGH->cctk_time;
  cctkGH->cctk_time = step_end - level_dt;
  int const status = SyncGroupsByDirI(cctkGH, static_cast<int>(groups.size()),
                                      &groups.front(), NULL);
  cctkGH->cctk_time = step_end;
  if (status < 0) {
    CCTK_ERROR("Hermite pre-step could not synchronise an evolved group.");
  }
}

} // namespace Carpet
