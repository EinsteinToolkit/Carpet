#include <cctk.h>

#include <cstddef>
#include <cstring>
#include <limits>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <carpet.hh>
#include "timestat.hh"

namespace Carpet {

namespace {

struct ReplicatedSlot {
  int vi;
  int tl;
  std::vector<char> bytes;
};

std::vector<ReplicatedSlot> replicated_slots;

#ifdef _OPENMP
const std::vector<ReplicatedSlot> *active_view = nullptr;
bool overlap_worker = false;
bool overlap_active = false;
bool group_info_owned = false;
#pragma omp threadprivate(active_view, overlap_worker, overlap_active,        \
                          group_info_owned)
#else
const std::vector<ReplicatedSlot> *active_view = nullptr;
bool overlap_worker = false;
bool overlap_active = false;
bool group_info_owned = false;
#endif

#ifdef CCTK_MPI
MPI_Comm fine_comm = MPI_COMM_NULL;
MPI_Comm coarse_comm = MPI_COMM_NULL;
#endif

// Not threadprivate. The master fills it before the pair. The coarse
// thread assigns it into its own threadprivate groupdata.
std::vector<groupdesc> groupdata_snapshot;

// Boundary records the coarse cGH pointer. Freeing it each step makes
// the next step look like a new GH.
cGH *cached_coarse_gh = nullptr;

int *duplicate_ints(int const *const src, int const count, char const *const what) {
  if (count < 1)
    CCTK_VERROR("CSAMR group info %s has length %d.", what, count);
  if (src == nullptr)
    CCTK_VERROR("CSAMR group info is missing %s.", what);
  void *const raw = std::malloc(sizeof(int) * static_cast<std::size_t>(count));
  if (raw == nullptr)
    CCTK_VERROR("CSAMR could not copy group info %s.", what);
  std::memcpy(raw, src, sizeof(int) * static_cast<std::size_t>(count));
  return static_cast<int *>(raw);
}

void release_ints(int const *&ptr) {
  std::free(const_cast<int *>(ptr));
  ptr = nullptr;
}

std::size_t replicated_bytes(gdata const *const gd, int const vi) {
  if (gd == nullptr or not gd->has_storage() or gd->storage() == nullptr)
    CCTK_VERROR("CSAMR replicated snapshot found empty storage for variable %d.",
                vi);
  int const elem = CCTK_VarTypeSize(CCTK_VarTypeI(vi));
  if (elem <= 0)
    CCTK_VERROR("CSAMR replicated snapshot does not know the size of variable %d.",
                vi);
  if (gd->elementsize() != elem) {
    CCTK_VERROR("CSAMR replicated snapshot element size %d of variable %d does "
                "not match the Cactus size %d.",
                gd->elementsize(), vi, elem);
  }
  size_type const elements = gd->size();
  if (elements < 0)
    CCTK_VERROR("CSAMR replicated snapshot of variable %d has a negative size.",
                vi);
  unsigned long long const count = static_cast<unsigned long long>(elements);
  unsigned long long const width = static_cast<unsigned long long>(elem);
  if (width != 0 and
      count > std::numeric_limits<std::size_t>::max() / width) {
    CCTK_VERROR("CSAMR replicated snapshot of variable %d overflows.", vi);
  }
  return static_cast<std::size_t>(count * width);
}

gdata *replicated_gdata(int const vi, int const tl) {
  int const group = CCTK_GroupIndexFromVarI(vi);
  int const var = vi - CCTK_FirstVarIndexI(group);
  if (group < 0 or var < 0)
    CCTK_VERROR("CSAMR replicated variable %d has no group variable.", vi);
  ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
  if (ff == nullptr or ff->h.local_components(0) < 1)
    CCTK_VERROR("CSAMR replicated variable %d has no local component.", vi);
  return ff->data_pointer(tl, 0, 0, 0);
}

int *clone_ints(int const *const src, int const count, char const *const what) {
  if (count < 1)
    CCTK_VERROR("CSAMR cGH clone %s has length %d.", what, count);
  if (src == nullptr)
    CCTK_VERROR("CSAMR cGH clone is missing %s.", what);
  void *const raw = std::malloc(sizeof(int) * static_cast<std::size_t>(count));
  if (raw == nullptr)
    CCTK_VERROR("CSAMR cGH clone could not allocate %s.", what);
  std::memcpy(raw, src, sizeof(int) * static_cast<std::size_t>(count));
  return static_cast<int *>(raw);
}

CCTK_REAL *clone_reals(CCTK_REAL const *const src, int const count,
                       char const *const what) {
  if (count < 1)
    CCTK_VERROR("CSAMR cGH clone %s has length %d.", what, count);
  if (src == nullptr)
    CCTK_VERROR("CSAMR cGH clone is missing %s.", what);
  void *const raw =
      std::malloc(sizeof(CCTK_REAL) * static_cast<std::size_t>(count));
  if (raw == nullptr)
    CCTK_VERROR("CSAMR cGH clone could not allocate %s.", what);
  std::memcpy(raw, src, sizeof(CCTK_REAL) * static_cast<std::size_t>(count));
  return static_cast<CCTK_REAL *>(raw);
}

} // namespace

bool csamr_overlap_worker() { return overlap_worker; }

bool csamr_overlap_active() { return overlap_active; }

void csamr_set_overlap_worker(bool const on) { overlap_worker = on; }

void csamr_set_overlap_active(bool const on) { overlap_active = on; }

void csamr_use_replicated_snapshot(bool const on) {
  active_view = on ? &replicated_slots : nullptr;
}

void csamr_retarget_replicated(cGH *const cctkGH) {
  if (active_view == nullptr)
    return;
  if (cctkGH == nullptr)
    CCTK_ERROR("CSAMR replicated retarget was called without a GH.");
  for (ReplicatedSlot const &slot : *active_view) {
    gdata const *const gd = replicated_gdata(slot.vi, slot.tl);
    if (replicated_bytes(gd, slot.vi) != slot.bytes.size()) {
      CCTK_VERROR("CSAMR replicated snapshot of variable %d time level %d "
                  "changed size.",
                  slot.vi, slot.tl);
    }
    if (slot.bytes.empty())
      CCTK_VERROR("CSAMR replicated snapshot of variable %d is empty.", slot.vi);
    cctkGH->data[slot.vi][slot.tl] =
        const_cast<char *>(slot.bytes.data());
  }
}

void csamr_snapshot_replicated() {
  if (not replicated_slots.empty())
    CCTK_ERROR("CSAMR replicated snapshot was taken before the previous one "
               "was released.");
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    int const grouptype = CCTK_GroupTypeI(group);
    if (grouptype != CCTK_SCALAR and grouptype != CCTK_ARRAY)
      continue;
    if (groupdata.AT(group).activetimelevels.empty() or
        groupdata.AT(group).activetimelevels.AT(0).empty()) {
      CCTK_VERROR("CSAMR replicated snapshot has no time-level count for "
                  "group %d.",
                  group);
    }
    int const ntl = groupdata.AT(group).activetimelevels.AT(0).AT(0);
    if (ntl < 1)
      continue;
    int const nvars = CCTK_NumVarsInGroupI(group);
    int const first = CCTK_FirstVarIndexI(group);
    for (int var = 0; var < nvars; ++var) {
      int const vi = first + var;
      if (arrdata.AT(group).AT(0).data.AT(var) == nullptr or
          arrdata.AT(group).AT(0).data.AT(var)->h.local_components(0) < 1) {
        CCTK_VERROR("CSAMR replicated snapshot has no local component for "
                    "variable %d.",
                    vi);
      }
      for (int tl = 0; tl < ntl; ++tl) {
        gdata const *const gd = replicated_gdata(vi, tl);
        std::size_t const nbytes = replicated_bytes(gd, vi);
        if (nbytes == 0) {
          CCTK_VERROR("CSAMR replicated snapshot of variable %d time level %d "
                      "is empty.",
                      vi, tl);
        }
        ReplicatedSlot slot;
        slot.vi = vi;
        slot.tl = tl;
        slot.bytes.resize(nbytes);
        std::memcpy(slot.bytes.data(), gd->storage(), nbytes);
        replicated_slots.push_back(slot);
      }
    }
  }
}

void csamr_release_replicated() { replicated_slots.clear(); }

void csamr_share_groupdata() {
  if (groupdata.empty())
    CCTK_ERROR("CSAMR has no group info to share.");
  groupdata_snapshot = groupdata;
}

void csamr_release_shared_groupdata() { groupdata_snapshot.clear(); }

void csamr_privatize_group_info() {
  if (group_info_owned)
    CCTK_ERROR("CSAMR group info was privatized twice on one thread.");
  if (groupdata_snapshot.empty())
    CCTK_ERROR("CSAMR group info was not shared before the threaded pair.");
  // copyin left this vector empty on the coarse thread.
  groupdata = groupdata_snapshot;
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    cGroupDynamicData &info = groupdata.AT(group).info;
    info.gsh = duplicate_ints(info.gsh, dim, "gsh");
    info.lsh = duplicate_ints(info.lsh, dim, "lsh");
    info.ash = duplicate_ints(info.ash, dim, "ash");
    info.lbnd = duplicate_ints(info.lbnd, dim, "lbnd");
    info.ubnd = duplicate_ints(info.ubnd, dim, "ubnd");
    info.tile_min = duplicate_ints(info.tile_min, dim, "tile_min");
    info.tile_max = duplicate_ints(info.tile_max, dim, "tile_max");
    info.bbox = duplicate_ints(info.bbox, 2 * dim, "bbox");
    info.nghostzones = duplicate_ints(info.nghostzones, dim, "nghostzones");
  }
  group_info_owned = true;
}

void csamr_release_group_info() {
  if (not group_info_owned)
    return;
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    cGroupDynamicData &info = groupdata.AT(group).info;
    release_ints(info.gsh);
    release_ints(info.lsh);
    release_ints(info.ash);
    release_ints(info.lbnd);
    release_ints(info.ubnd);
    release_ints(info.tile_min);
    release_ints(info.tile_max);
    release_ints(info.bbox);
    release_ints(info.nghostzones);
  }
  group_info_owned = false;
}

namespace {

void *clone_var_pointers(cGH const *const src, int const vi) {
  int const ntl = CCTK_DeclaredTimeLevelsVI(vi);
  if (ntl < 1)
    CCTK_VERROR("CSAMR cGH clone variable %d declares %d time levels.", vi, ntl);
  if (src->data == nullptr or src->data[vi] == nullptr)
    CCTK_VERROR("CSAMR cGH clone is missing data for variable %d.", vi);
  void *const raw =
      std::malloc(sizeof(void *) * static_cast<std::size_t>(ntl));
  if (raw == nullptr)
    CCTK_VERROR("CSAMR cGH clone could not allocate data for variable %d.", vi);
  std::memcpy(raw, src->data[vi], sizeof(void *) * static_cast<std::size_t>(ntl));
  return raw;
}

} // namespace

cGH *csamr_clone_gh(cGH const *const src) {
  if (src == nullptr)
    CCTK_ERROR("CSAMR cGH clone was called without a GH.");
  if (cached_coarse_gh != nullptr)
    return cached_coarse_gh;
  void *const raw = std::malloc(sizeof(cGH));
  if (raw == nullptr)
    CCTK_ERROR("CSAMR could not allocate a cGH clone.");
  cGH *const gh = static_cast<cGH *>(raw);
  std::memcpy(gh, src, sizeof(cGH));
  gh->current_scheduled_function = nullptr;
  // extensions and identity stay shared. Everything the two updates
  // write gets its own allocation.
  int const n = src->cctk_dim > 0 ? src->cctk_dim : 1;
  gh->cctk_gsh = clone_ints(src->cctk_gsh, n, "gsh");
  gh->cctk_lsh = clone_ints(src->cctk_lsh, n, "lsh");
  gh->cctk_lbnd = clone_ints(src->cctk_lbnd, n, "lbnd");
  gh->cctk_ubnd = clone_ints(src->cctk_ubnd, n, "ubnd");
  gh->cctk_tile_min = clone_ints(src->cctk_tile_min, n, "tile_min");
  gh->cctk_tile_max = clone_ints(src->cctk_tile_max, n, "tile_max");
  gh->cctk_ash = clone_ints(src->cctk_ash, n, "ash");
  gh->cctk_to = clone_ints(src->cctk_to, n, "to");
  gh->cctk_from = clone_ints(src->cctk_from, n, "from");
  gh->cctk_levfac = clone_ints(src->cctk_levfac, n, "levfac");
  gh->cctk_levoff = clone_ints(src->cctk_levoff, n, "levoff");
  gh->cctk_levoffdenom = clone_ints(src->cctk_levoffdenom, n, "levoffdenom");
  gh->cctk_bbox = clone_ints(src->cctk_bbox, 2 * n, "bbox");
  gh->cctk_nghostzones = clone_ints(src->cctk_nghostzones, 2 * n, "nghostzones");
  gh->cctk_delta_space = clone_reals(src->cctk_delta_space, n, "delta_space");
  gh->cctk_origin_space = clone_reals(src->cctk_origin_space, n, "origin_space");

  int const nvars = CCTK_NumVars();
  int const nvar_slots = nvars > 0 ? nvars : 1;
  void *const data_raw =
      std::malloc(sizeof(void **) * static_cast<std::size_t>(nvar_slots));
  if (data_raw == nullptr)
    CCTK_ERROR("CSAMR cGH clone could not allocate the data table.");
  gh->data = static_cast<void ***>(data_raw);
  for (int vi = 0; vi < nvars; ++vi)
    gh->data[vi] = static_cast<void **>(clone_var_pointers(src, vi));

  int const ngroups = CCTK_NumGroups();
  int const ngroup_slots = ngroups > 0 ? ngroups : 1;
  void *const groups_raw = std::malloc(sizeof(cGHGroupData) *
                                      static_cast<std::size_t>(ngroup_slots));
  if (groups_raw == nullptr)
    CCTK_ERROR("CSAMR cGH clone could not allocate group data.");
  if (src->GroupData == nullptr)
    CCTK_ERROR("CSAMR cGH clone is missing group data.");
  std::memcpy(groups_raw, src->GroupData,
              sizeof(cGHGroupData) * static_cast<std::size_t>(ngroup_slots));
  gh->GroupData = static_cast<cGHGroupData *>(groups_raw);
  cached_coarse_gh = gh;
  return gh;
}

void csamr_free_gh(cGH *const gh) {
  if (gh == nullptr || gh == cached_coarse_gh)
    return;
  std::free(gh->cctk_gsh);
  std::free(gh->cctk_lsh);
  std::free(gh->cctk_lbnd);
  std::free(gh->cctk_ubnd);
  std::free(gh->cctk_tile_min);
  std::free(gh->cctk_tile_max);
  std::free(gh->cctk_ash);
  std::free(gh->cctk_to);
  std::free(gh->cctk_from);
  std::free(gh->cctk_levfac);
  std::free(gh->cctk_levoff);
  std::free(gh->cctk_levoffdenom);
  std::free(gh->cctk_bbox);
  std::free(gh->cctk_nghostzones);
  std::free(gh->cctk_delta_space);
  std::free(gh->cctk_origin_space);
  if (gh->data != nullptr) {
    int const nvars = CCTK_NumVars();
    for (int vi = 0; vi < nvars; ++vi)
      std::free(gh->data[vi]);
    std::free(gh->data);
  }
  std::free(gh->GroupData);
  std::free(gh);
}

void csamr_suppress_lib_timers(bool const on) {
  CarpetLib::Timer::suppress(on);
}

void csamr_require_thread_multiple() {
#ifdef CCTK_MPI
  int provided = MPI_THREAD_SINGLE;
  if (MPI_Query_thread(&provided) != MPI_SUCCESS)
    CCTK_ERROR("CSAMR could not query the MPI thread level.");
  if (provided < MPI_THREAD_MULTIPLE)
    CCTK_ERROR("CSAMR thread overlap needs MPI_THREAD_MULTIPLE.");
#else
  CCTK_ERROR("CSAMR thread overlap needs MPI.");
#endif
}

extern "C" CCTK_INT MoLNumEvolvedArrayVariables __attribute__((weak));

void csamr_reject_evolved_arrays() {
  if (&MoLNumEvolvedArrayVariables != nullptr and
      MoLNumEvolvedArrayVariables > 0) {
    CCTK_ERROR("CSAMR thread overlap has no private MoL array scratch.");
  }
}

void csamr_create_thread_comms() {
#ifdef CCTK_MPI
  if (fine_comm == MPI_COMM_NULL) {
    if (MPI_Comm_dup(MPI_COMM_WORLD, &fine_comm) != MPI_SUCCESS)
      CCTK_ERROR("CSAMR could not duplicate MPI_COMM_WORLD for the finest "
                 "update.");
  }
  if (coarse_comm == MPI_COMM_NULL) {
    if (MPI_Comm_dup(MPI_COMM_WORLD, &coarse_comm) != MPI_SUCCESS)
      CCTK_ERROR("CSAMR could not duplicate MPI_COMM_WORLD for the coarser "
                 "update.");
  }
  if (fine_comm == MPI_COMM_NULL or coarse_comm == MPI_COMM_NULL or
      fine_comm == coarse_comm) {
    CCTK_ERROR("CSAMR thread communicators were not created.");
  }
#else
  CCTK_ERROR("CSAMR thread overlap needs MPI.");
#endif
}

MPI_Comm csamr_fine_comm() {
#ifdef CCTK_MPI
  if (fine_comm == MPI_COMM_NULL)
    CCTK_ERROR("CSAMR finest communicator was used before it was created.");
  return fine_comm;
#else
  CCTK_ERROR("CSAMR thread overlap needs MPI.");
  return MPI_COMM_NULL;
#endif
}

MPI_Comm csamr_coarse_comm() {
#ifdef CCTK_MPI
  if (coarse_comm == MPI_COMM_NULL)
    CCTK_ERROR("CSAMR coarser communicator was used before it was created.");
  return coarse_comm;
#else
  CCTK_ERROR("CSAMR thread overlap needs MPI.");
  return MPI_COMM_NULL;
#endif
}

} // namespace Carpet
