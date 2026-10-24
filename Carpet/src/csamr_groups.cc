#include <cctk.h>
#include <cctk_Parameters.h>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <carpet.hh>
#include "csamr.hh"

namespace Carpet {

namespace {

int world_size() {
  int world = 0;
  if (MPI_Comm_size(MPI_COMM_WORLD, &world) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not read the MPI world size.");
  }
  return world;
}

int world_rank() {
  int rank = 0;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not read the MPI world rank.");
  }
  return rank;
}

} // namespace

void csamr_create_group_comm() {
  CCTK_ERROR("CSAMR no longer splits ranks across levels. A slot's two "
             "updates share every rank and run on separate OpenMP threads.");
}

MPI_Comm csamr_owner_reduction_comm() { return MPI_COMM_NULL; }

void csamr_broadcast_evolved_scalars(cGH *const cctkGH, int const source_rank) {
  if (cctkGH == NULL) {
    CCTK_ERROR("CSAMR scalar broadcast was called without a GH.");
  }
  int const world = world_size();
  if (source_rank < 0 or source_rank >= world) {
    CCTK_VERROR("CSAMR scalar broadcast source rank %d is outside 0..%d.",
                source_rank, world - 1);
  }
  // This runs outside level mode, where mglevel is -1. QueryGroupStorage
  // and ActiveTimeLevels index activetimelevels with that mglevel.
  // Scalars have one multigrid level and no refinement level.
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_GroupTypeI(group) != CCTK_SCALAR)
      continue;
    if (groupdata.AT(group).activetimelevels.empty() or
        groupdata.AT(group).activetimelevels.AT(0).empty()) {
      CCTK_VERROR("CSAMR scalar broadcast has no time-level count for "
                  "group %d.",
                  group);
    }
    int const ntl = groupdata.AT(group).activetimelevels.AT(0).AT(0);
    if (ntl <= 1)
      continue;
    int const nvars = CCTK_NumVarsInGroupI(group);
    int const first = CCTK_FirstVarIndexI(group);
    for (int var = 0; var < nvars; ++var) {
      int const vi = first + var;
      int const size = CCTK_VarTypeSize(CCTK_VarTypeI(vi));
      if (size <= 0) {
        CCTK_VERROR("CSAMR scalar broadcast does not know the size of "
                    "variable %d.",
                    vi);
      }
      ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
      if (ff->h.local_components(0) < 1) {
        CCTK_VERROR("CSAMR scalar broadcast has no local component for "
                    "variable %d.",
                    vi);
      }
      for (int tl = 0; tl < ntl; ++tl) {
        void *const ptr = ff->data_pointer(tl, 0, 0, 0)->storage();
        if (ptr == NULL) {
          CCTK_ERROR("CSAMR scalar broadcast found empty storage.");
        }
        int const status =
            MPI_Bcast(ptr, size, MPI_BYTE, source_rank, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) {
          CCTK_VERROR("CSAMR scalar broadcast failed for variable %d time "
                      "level %d (status %d).",
                      vi, tl, status);
        }
      }
    }
  }
}

namespace {

CCTK_INT *csamr_nan_counter_storage() {
  if (not CCTK_IsThornActive("NaNChecker"))
    return NULL;
  int const vi = CCTK_VarIndex("NaNChecker::NaNsFound");
  if (vi < 0) {
    CCTK_ERROR("NaNChecker is active but NaNChecker::NaNsFound has no "
               "variable index.");
  }
  if (CCTK_GroupTypeFromVarI(vi) != CCTK_SCALAR) {
    CCTK_ERROR("NaNChecker::NaNsFound is not a scalar.");
  }
  int const group = CCTK_GroupIndexFromVarI(vi);
  int const var = vi - CCTK_FirstVarIndexI(group);
  if (group < 0 or var < 0) {
    CCTK_ERROR("NaNChecker::NaNsFound has no group variable.");
  }
  ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
  if (ff->h.local_components(0) < 1) {
    CCTK_ERROR("NaNChecker::NaNsFound has no local component.");
  }
  void *const ptr = ff->data_pointer(0, 0, 0, 0)->storage();
  if (ptr == NULL) {
    CCTK_ERROR("NaNChecker::NaNsFound has empty storage.");
  }
  return static_cast<CCTK_INT *>(ptr);
}

} // namespace

void csamr_reset_nan_counter() {
  CCTK_INT *const counter = csamr_nan_counter_storage();
  if (counter != NULL)
    *counter = 0;
}

double csamr_nan_counter() {
  CCTK_INT const *const counter = csamr_nan_counter_storage();
  if (counter == NULL)
    return -1.0;
  return static_cast<double>(*counter);
}

void CarpetCsamrShareCoarseSpacing(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CarpetCsamrShareCoarseSpacing;
  DECLARE_CCTK_PARAMETERS;

  // Every rank is on this level together. Returning here stays collective.
  if (not use_csamr or reflevels < 2 or reflevel != 0)
    return;
  if (sizeof(CCTK_REAL) != sizeof(double)) {
    CCTK_ERROR("CSAMR coarse spacing broadcast requires CCTK_REAL to be "
               "double.");
  }
  if (cctkGH->cctk_dim <= 0) {
    CCTK_ERROR("CSAMR coarse spacing broadcast has no spatial dimensions.");
  }
  if (maps != 1) {
    CCTK_ERROR("CSAMR coarse spacing is published for one map.");
  }
  int root_owns = 0;
  if (world_rank() == 0) {
    root_owns = vhh.AT(0)->local_components(0) > 0 ? 1 : 0;
  }
  if (MPI_Bcast(&root_owns, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not broadcast the coarse-owner flag.");
  }
  if (root_owns != 1) {
    CCTK_ERROR("CSAMR rank 0 owns no component on the coarsest level, so it "
               "cannot publish the coarse spacing.");
  }

  vector<CCTK_REAL> spacing(cctkGH->cctk_dim);
  for (int d = 0; d < cctkGH->cctk_dim; ++d) {
    spacing.AT(d) = cctkGH->cctk_delta_space[d];
  }
  int const status = MPI_Bcast(&spacing.front(), cctkGH->cctk_dim, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not broadcast the coarse grid spacing.");
  }
  for (int d = 0; d < cctkGH->cctk_dim; ++d) {
    cctkGH->cctk_delta_space[d] = spacing.AT(d);
  }
  if (world_rank() == 0) {
    CCTK_VInfo(CCTK_THORNSTRING, "CSAMR coarse spacing is %g",
               static_cast<double>(spacing.AT(0)));
  }
}

} // namespace Carpet
