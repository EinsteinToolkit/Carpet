#ifndef CARPET_HH
#define CARPET_HH

#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Schedule.h>

#include <gh.hh>

#include "carpet_public.hh"

namespace Carpet {

using namespace std;

// Scheduled functions
extern "C" {
int CarpetStartup(void);
int CarpetMultiModelStartup(void);
void CarpetParamCheck(CCTK_ARGUMENTS);
void CarpetRefineTimeStep(CCTK_ARGUMENTS);
void CarpetCsamrShareCoarseSpacing(CCTK_ARGUMENTS);
void CarpetPoisonTaper(CCTK_ARGUMENTS);
void Carpet_HermitePreStep(CCTK_ARGUMENTS);
void CarpetUnusedMask(CCTK_ARGUMENTS);
}

void NoteTaperHaloRefilled();
// True while serial evolution should fill fine ghosts from Hermite output.
extern bool hermite_prolongate;
// Registered functions
void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cgh);

int Initialise(tFleshConfig *config);
int Evolve(tFleshConfig *config);
int Shutdown(tFleshConfig *config);
int OutputGH(const cGH *cgh);

int CallFunction(void *function, cFunctionData *attribute, void *data);

// Other functions
bool Regrid(cGH const *cctkGH, bool force_recompose, bool do_init);

void CycleTimeLevels(cGH *cgh, bool cycle_scalars = true);
// Rotate scalar and array time levels on this rank. Grid functions are
// left alone.
void CycleReplicatedLevels(cGH *cgh);
// Kept so an old call site fails instead of splitting ranks.
void csamr_create_group_comm();
// The coarse OpenMP thread of a CSAMR pair. It must not write the
// shared coordinate spacing.
bool csamr_overlap_worker();
// True on both threads of a CSAMR pair. Per-level timers are skipped.
bool csamr_overlap_active();
// Point replicated variables at the finest update's pre-cycle snapshot.
void csamr_retarget_replicated(cGH *cgh);
// Copy groupdata for the coarse thread. OpenMP copyin left that
// vector empty. Call this on the master thread before the pair.
void csamr_share_groupdata();
void csamr_release_shared_groupdata();
// Deep-copy the group-info arrays on this thread. The coarse thread
// assigns the shared snapshot first, then replaces the int arrays.
void csamr_privatize_group_info();
void csamr_release_group_info();
void csamr_snapshot_replicated();
void csamr_release_replicated();
void csamr_use_replicated_snapshot(bool on);
void csamr_set_overlap_worker(bool on);
void csamr_set_overlap_active(bool on);
cGH *csamr_clone_gh(cGH const *src);
void csamr_free_gh(cGH *gh);
// CarpetLib timers are process-wide. Two updates must not start one.
void csamr_suppress_lib_timers(bool on);
void csamr_require_thread_multiple();
void csamr_reject_evolved_arrays();
void csamr_create_thread_comms();
MPI_Comm csamr_fine_comm();
MPI_Comm csamr_coarse_comm();
// MPI_COMM_NULL. Reductions stay on dist::comm(), which is the
// per-update duplicate while a threaded pair is running.
MPI_Comm csamr_owner_reduction_comm();
// Copy multi-timelevel scalars from source_rank onto every rank.
void csamr_broadcast_evolved_scalars(cGH *cgh, int source_rank);
// The scheduled NaN reset is global, so the process group that skips
// PRESTEP would otherwise keep a stale count. Every rank calls this.
void csamr_reset_nan_counter();
// NaNChecker::NaNsFound on this rank, or -1 when that thorn is inactive.
double csamr_nan_counter();
void UncycleTimeLevels(cGH *cgh);
void FlipTimeLevels(cGH *cgh);
void FillTimeLevels(const cGH *cgh);
void SyncGroups(const cGH *cgh, const vector<int> &groups);
int SyncProlongateGroups(const cGH *cgh, const vector<int> &groups,
                         cFunctionData const *function_data = NULL);

// Sanity checks
enum checktimes {
  currenttime,
  currenttimebutnotifonly,
  previoustime,
  allbutlasttime,
  allbutcurrenttime,
  alltimes
};

int min_timelevel(checktimes where, int num_tl,
                  bool persistent) CCTK_ATTRIBUTE_CONST;
int max_timelevel(checktimes where, int num_tl,
                  bool persistent) CCTK_ATTRIBUTE_CONST;

void Poison(const cGH *cgh, checktimes where, int what = 0);
void PoisonGroup(const cGH *cgh, int group, checktimes where);
void PoisonCheck(const cGH *cgh, checktimes where);

void CalculateChecksums(const cGH *cgh, checktimes where);
void CheckChecksums(const cGH *cgh, checktimes where);

// Schedule
int CallBeforeRoutines(cGH const *cctkGH, void *const function,
                       cFunctionData *const attribute, void *const data);
int CallAfterRoutines(cGH const *cctkGH, void *const function,
                      cFunctionData *const attribute, void *const data);

// Debugging output
void Output(const char *fmt, ...);
void Waypoint(const char *fmt, ...);
void Checkpoint(const char *fmt, ...);

// Error output
void UnsupportedVarType(int vindex);

// Check for a map0group
bool IsMap0Group(int gindex);

} // namespace Carpet

#endif // !defined(CARPET_HH)
