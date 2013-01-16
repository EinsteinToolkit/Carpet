#include "loopcontrol.h"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <list>
#include <ostream>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#else

// Simple non-OpenMP implementations
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }

#endif

#ifdef HAVE_CAPABILITY_CYCLECLOCK
// We have a fast, accurate clock

#  include <cycleclock.h>

#else

#  ifdef _OPENMP
// We use the OpenMP clock

typedef double ticks;
static inline ticks getticks() { return omp_get_wtime(); }
static inline double elapsed(ticks t1, ticks t0) { return t1-t0; }
static inline double seconds_per_tick() { return 1.0; }

#  else
// We use gettimeofday as fallback

#    include <sys/time.h>
typedef timeval ticks;
static inline ticks getticks()
{
  timeval tp;
  gettimeofday(&tp, NULL);
  return tp;
}
static inline double elapsed(ticks t1, ticks t0)
{
  return 1.0e+6 * (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec);
}
static inline double seconds_per_tick() { return 1.0e-6; }

#  endif
#endif

using namespace std;



struct lc_thread_info_t {
  char padding1[128];      // pad to ensure cache lines are not shared
  int idx;                 // linear index of next coarse thread block
  char padding2[128];
};

struct lc_fine_thread_comm_t {
  char padding1[128];      // pad to ensure cache lines are not shared
  int state;               // waiting threads
  int value;               // broadcast value
  char padding2[128];
};



struct lc_stats_t {
  char const* name;
  char const* file;
  int line;
  int init_count;
  double points, threads;
  double count, sum, sum2, min, max;
  ticks start_time;
};

  

extern "C" CCTK_FCALL
void CCTK_FNAME(lc_get_fortran_type_sizes) (ptrdiff_t* type_sizes);
  


namespace {
  
  list<lc_stats_t*> all_stats;
  
  
  
  void check_fortran_type_sizes()
  {
    ptrdiff_t type_sizes[4];
    CCTK_FNAME(lc_get_fortran_type_sizes) (type_sizes);
    assert(type_sizes[0] == sizeof(lc_vec_t));
    assert(type_sizes[1] == sizeof(lc_ivec_t));
    assert(type_sizes[2] == sizeof(lc_space_t));
    assert(type_sizes[3] == sizeof(lc_control_t));
  }
  
  
  
  template<typename T>
  T divup(T const i, T const j)
  {
    assert(i >= 0 and j > 0);
    return (i + j - 1) / j;
  }
  
  template<typename T>
  T alignup(T const i, T const j)
  {
    return divup(i, j) * j;
  }
  
  template<typename T>
  T divdown(T const i, T const j)
  {
    assert(i >= 0 and j > 0);
    return i / j;
  }
  
  template<typename T>
  T aligndown(T const i, T const j)
  {
    return divdown(i, j) * j;
  }
  
  
  
  ostream& operator<<(ostream& os, lc_vec_t const& x)
  {
    os << "[";
    for (int d=0; d<LC_DIM; ++d) {
      if (d>0) os << ",";
      os << x.v[d];
    }
    os << "]";
    return os;
  }
  
  ostream& operator<<(ostream& os, lc_ivec_t const& x)
  {
    os << "[";
    for (int d=0; d<LC_DIM; ++d) {
      if (d>0) os << ",";
      os << x.v[d];
    }
    os << "]";
    return os;
  }
  
  ostream& operator<<(ostream& os, lc_space_t const& s)
  {
    os << "{"
       << "min:" << s.min << ","
       << "max:" << s.max << ","
       << "step:" << s.step << ","
       << "pos:" << s.pos << ","
       << "count:" << s.count << ","
       << "idx:" << s.idx
       << "}";
    return os;
  }
  
  ostream& operator<<(ostream& os, lc_control_t const& c)
  {
    os << "lc_control{\n"
       << "   ash:" << c.ash << ",\n"
       << "   loop:" << c.loop << ",\n"
       << "   thread:" << c.thread << ",\n"
       << "   coarse:" << c.coarse << ",\n"
       << "   fine:" << c.fine << "\n"
       << "   fine_thread:" << c.fine_thread << "\n"
       << "}\n";
    return os;
  }
  
  
  
  ptrdiff_t prod(lc_vec_t const& x)
  {
    ptrdiff_t r = 1;
    for (int d=0; d<LC_DIM; ++d) {
      assert(x.v[d] >= 0);
      r *= x.v[d];
    }
    return r;
  }
  
  ptrdiff_t ind(lc_vec_t const& shape, lc_vec_t const& pos)
  {
    ptrdiff_t r = 0;
    ptrdiff_t f = 1;
    for (int d=0; d<LC_DIM; ++d) {
      assert(pos.v[d] >= 0 and pos.v[d] < shape.v[d]);
      r += f * pos.v[d];
      assert(shape.v[d] >= 0);
      f *= shape.v[d];
    }
    return r;
  }
  
  ptrdiff_t ind(lc_vec_t const& shape,
                ptrdiff_t const i, ptrdiff_t const j, ptrdiff_t const k)
  {
    lc_vec_t const pos = {{ i, j, k }};
    return ind(shape, pos);
  }
  
  
  
  void space_set_count(lc_space_t& space)
  {
    for (int d=0; d<LC_DIM; ++d) {
      space.count.v[d] =
        divup(space.max.v[d] - space.min.v[d], space.step.v[d]);
    }
  }
  
  void space_idx2pos(lc_space_t& space)
  {
    for (int d=0; d<LC_DIM; ++d) {
      space.pos.v[d] = space.min.v[d] + space.idx.v[d] * space.step.v[d];
    }
  }
  
  bool space_global2local(lc_space_t& space, int gidx)
  {
    assert(gidx >= 0);
    for (int d=0; d<LC_DIM; ++d) {
      if (space.count.v[d] > 0) {
        space.idx.v[d] = gidx % space.count.v[d];
        gidx /= space.count.v[d];
      } else {
        space.idx.v[d] = 0;
      }
    }
    return gidx != 0;
  }
  
  int space_local2global(lc_space_t const& space)
  {
    int gidx = 0;
    int fact = 1;
    for (int d=0; d<LC_DIM; ++d) {
      assert(space.idx.v[d] >= 0 and space.idx.v[d] < space.count.v[d]);
      gidx += fact * space.idx.v[d];
      fact *= space.count.v[d];
    }
    return gidx;
  }
  
  
  
  int get_num_fine_threads()
  {
    DECLARE_CCTK_PARAMETERS;
    if (not use_smt_threads) return 1;
    if (omp_get_num_threads() == 1) return 1;
    static int num_smt_threads = -1;
    if (CCTK_BUILTIN_EXPECT(num_smt_threads<0, false)) {
#pragma omp barrier
      0;                        // PGI compiler needs this
#pragma omp master
      if (CCTK_IsFunctionAliased("GetNumSMTThreads")) {
        num_smt_threads = GetNumSMTThreads();
      } else {
        num_smt_threads = 1;
      }
#pragma omp barrier
    }
    return num_smt_threads;
  }
  
  int get_fine_thread_num()
  {
    DECLARE_CCTK_PARAMETERS;
    if (not use_smt_threads) return 0;
    if (omp_get_num_threads() == 1) return 0;
    int const thread_num = omp_get_thread_num();
    int const num_fine_threads = get_num_fine_threads();
    return thread_num % num_fine_threads;
  }
  
  int get_num_coarse_threads()
  {
    int const num_threads = omp_get_num_threads();
    int const num_fine_threads = get_num_fine_threads();
    assert(num_threads % num_fine_threads == 0);
    return num_threads / num_fine_threads;
  }
  
  int get_coarse_thread_num()
  {
    int const thread_num = omp_get_thread_num();
    int const num_fine_threads = get_num_fine_threads();
    return thread_num / num_fine_threads;
  }
  
  // Wait until *ptr is different from old_value
  void thread_wait(int const *const ptr, int const old_value)
  {
    while (*ptr == old_value) {
#pragma omp flush
      0;                        // PGI compiler needs this
    }
  }
  
  int fine_thread_broadcast(lc_fine_thread_comm_t* const comm, int value)
  {
    int const num_fine_threads = get_num_fine_threads();
    if (num_fine_threads == 1) return value;
    assert(num_fine_threads < 8 * int(sizeof comm->state));
    int const fine_thread_num = get_fine_thread_num();
    int const master_mask = 1;
    
    // Assume comm->count == 0 initially
    if (fine_thread_num == 0) { // if on master
      
      int const all_threads_mask = (1 << num_fine_threads) - 1;
      if (comm->state != 0) {
        // wait until everybody has acknowledged the previous value
#pragma omp flush
        for (;;) {
          int const state = comm->state;
          if (comm->state == all_threads_mask) break;
          thread_wait(&comm->state, state);
        }
        // mark the value as invalid
        comm->state = 0;
#pragma omp flush
      }
      // publish value
      comm->value = value;
#pragma omp flush
      // mark the value as valid
      comm->state = master_mask;
#pragma omp flush
      
    } else {                    // if not on master
      
      // wait until the value is valid, and it is a new value
      int const thread_mask = 1 << fine_thread_num;
#pragma omp flush
      for (;;) {
        int const state = comm->state;
        if ((comm->state & (master_mask | thread_mask)) == master_mask) break;
        thread_wait(&comm->state, state);
      }
      // read value
      value = comm->value;
#pragma omp flush
      0;                        // PGI compiler needs this
      // acknowledge the value
#pragma omp atomic
      comm->state |= thread_mask;
#pragma omp flush
      
    } // if not on master
    
    return value;
  }
  
} // namespace



void lc_stats_init(lc_stats_t** const stats_ptr,
                   char const* const name,
                   char const* const file,
                   int const line)
{
  if (*stats_ptr) return;
  
  lc_stats_t* stats;
#pragma omp single copyprivate(stats)
  {
    stats = new lc_stats_t;
    
    stats->name = strdup(name);
    stats->file = strdup(file);
    stats->line = line;
    stats->init_count = 0;
    stats->points = 0.0;
    stats->threads = 0.0;
    stats->count = 0.0;
    stats->sum = 0.0;
    stats->sum2 = 0.0;
    stats->min = numeric_limits<double>::max();
    stats->max = 0.0;
    
    all_stats.push_back(stats);
  }
#pragma omp single
  {
    *stats_ptr = stats;
  }
}



void lc_control_init(lc_control_t* restrict const control,
                     lc_stats_t* const stats,
                     ptrdiff_t imin, ptrdiff_t jmin, ptrdiff_t kmin,
                     ptrdiff_t imax, ptrdiff_t jmax, ptrdiff_t kmax,
                     ptrdiff_t iash, ptrdiff_t jash, ptrdiff_t kash,
                     ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dk)
{
  DECLARE_CCTK_PARAMETERS;
  
#pragma omp barrier
#pragma omp master
  {
    stats->start_time = getticks();
  }
  
  // Ensure thread counts are consistent
  assert(get_num_coarse_threads() * get_num_fine_threads() ==
         omp_get_num_threads());
  
  // Get cache line size
  static ptrdiff_t max_cache_linesize = -1;
  if (CCTK_BUILTIN_EXPECT(max_cache_linesize<0, false)) {
#pragma omp barrier
#pragma omp master
    if (CCTK_IsFunctionAliased("GetCacheInfo1")) {
      int const num_levels = GetCacheInfo1(NULL, NULL, 0);
      vector<int> linesizes(num_levels);
      vector<int> strides  (num_levels);
      GetCacheInfo1(&linesizes[0], &strides[0], num_levels);
      max_cache_linesize = 1;
      for (int level=0; level<num_levels; ++level) {
        max_cache_linesize =
          max(max_cache_linesize, ptrdiff_t(linesizes[level]));
      }
    } else {
      max_cache_linesize = 1;
    }
#pragma omp barrier
  }
  
  // Initialize everything with a large, bogus value
  memset(control, 123, sizeof *control);
  
  ptrdiff_t tilesize_alignment = 1;
  if (align_with_cachelines) {
    tilesize_alignment =
      divup(max_cache_linesize, ptrdiff_t(sizeof(CCTK_REAL)));
    tilesize_alignment = alignup(tilesize_alignment, di);
  }
  
  // Parameters (all in units of grid points)
  // ptrdiff_t const smt_size[LC_DIM] = { smtsize_i, smtsize_j, smtsize_k };
  // TODO: put fine threads into i direction, so that they share cache
  // lines
  ptrdiff_t smt_size[LC_DIM] = { 1, 1, 1 };
  {
    int const num_fine_threads = get_num_fine_threads();
    if (num_fine_threads <= loopsize_j) {
      smt_size[1] = num_fine_threads;
    } else if (num_fine_threads <= loopsize_k) {
      smt_size[2] = num_fine_threads;
    } else {
      smt_size[0] = num_fine_threads;
    }
  }
  ptrdiff_t const tile_size[LC_DIM] = {
    alignup(ptrdiff_t(tilesize_i), tilesize_alignment),
    tilesize_j,
    tilesize_k,
  };
  ptrdiff_t const loop_size[LC_DIM] = { loopsize_i, loopsize_j, loopsize_k };
  
  // Arguments
  ptrdiff_t const loop_min[LC_DIM] = { imin, jmin, kmin };
  ptrdiff_t const loop_max[LC_DIM] = { imax, jmax, kmax };
  ptrdiff_t const ash[LC_DIM] = { iash, jash, kash };
  ptrdiff_t const vect_size[LC_DIM] = { di, dj, dk };
  
  // Copy ash arguments
  for (int d=0; d<LC_DIM; ++d) {
    control->ash.v[d] = ash[d];
  }
  
  // Set up multithreading state
  lc_thread_info_t* thread_info_ptr;
#pragma omp single copyprivate(thread_info_ptr)
  {
    thread_info_ptr = new lc_thread_info_t;
  }
  control->thread_info_ptr = thread_info_ptr;
  
  {
    lc_fine_thread_comm_t** fine_thread_comm_ptrs;
#pragma omp single copyprivate(fine_thread_comm_ptrs)
    {
      fine_thread_comm_ptrs =
        new lc_fine_thread_comm_t*[get_num_coarse_threads()];
    }
    if (get_fine_thread_num() == 0) {
      lc_fine_thread_comm_t* fine_thread_comm_ptr;
      fine_thread_comm_ptr = new lc_fine_thread_comm_t;
      fine_thread_comm_ptr->state = 0;
      fine_thread_comm_ptrs[get_coarse_thread_num()] = fine_thread_comm_ptr;
    }
#pragma omp barrier
    control->fine_thread_comm_ptr =
      fine_thread_comm_ptrs[get_coarse_thread_num()];
#pragma omp barrier
#pragma omp single nowait
    {
      delete[] fine_thread_comm_ptrs;
    }
  }
  
  // Set loop sizes
  for(int d=0; d<LC_DIM; ++d) {
    // Overall loop: as specified
    control->loop.min.v[d] = loop_min[d];
    control->loop.max.v[d] = loop_max[d];
    // Thread loop
#if VECTORISE_ALIGNED_ARRAYS
    // Move start to be aligned with vector size
    control->thread.min.v[d] =
      aligndown(control->loop.min.v[d], vect_size[d]);
#else
    control->thread.min.v[d] = control->loop.min.v[d];
#endif
    control->thread.max.v[d] = loop_max[d];
    // Fine threads
    control->fine_thread.count.v[d] = smt_size[d];
  }
  {
    int const fine_thread_num = get_fine_thread_num();
    bool const outside =
      space_global2local(control->fine_thread, fine_thread_num);
    assert(not outside);
  }
  
  // Set loop step sizes
  if (CCTK_EQUALS(initial_setup, "legacy")) {
    // Like a non-LoopControl loop: no loop tiling (i.e. do not use
    // coarse loops), parallelise only in k direction (i.e. assign
    // equal k ranges to threads)
    for(int d=0; d<LC_DIM; ++d) {
      assert(smt_size[d] == 1); // TODO: implement this
      control->fine_thread.step.v[d] = vect_size[d];
      control->fine.step.v[d] = vect_size[d];
      ptrdiff_t const npoints = control->loop.max.v[d] - control->loop.min.v[d];
      ptrdiff_t const nthreads = d!=LC_DIM-1 ? 1 : get_num_coarse_threads();
      control->coarse.step.v[d] =
        alignup(divup(npoints, nthreads), control->fine.step.v[d]);
      control->thread.step.v[d] = alignup(npoints, control->coarse.step.v[d]);
    }
  } else if (CCTK_EQUALS(initial_setup, "tiled")) {
    // Basic LoopControl setup
    for(int d=0; d<LC_DIM; ++d) {
      control->fine_thread.step.v[d] = vect_size[d];
      control->fine.step.v[d] =
        alignup(smt_size[d], control->fine_thread.step.v[d]);
      control->coarse.step.v[d] =
        alignup(tile_size[d], control->fine.step.v[d]);
      control->thread.step.v[d] =
        alignup(loop_size[d], control->coarse.step.v[d]);
    }
  } else {
    CCTK_WARN(CCTK_WARN_ABORT, "internal error");
  }
  
  if (veryverbose) {
#pragma omp master
    CCTK_VInfo(CCTK_THORNSTRING,
               "Loop %s (%s:%d): imin=[%td,%td,%td] imax=[%td,%td,%td]\n"
               "   smt.step=[%td,%td,%td] fine.step=[%td,%td,%td] coarse.step=[%td,%td,%td] thread.step=[%td,%td,%td]",
               stats->name, stats->file, stats->line,
               control->loop.min.v[0], control->loop.min.v[1], control->loop.min.v[2],
               control->loop.max.v[0], control->loop.max.v[1], control->loop.max.v[2],
               control->fine_thread.step.v[0], control->fine_thread.step.v[1], control->fine_thread.step.v[2],
               control->fine.step.v[0], control->fine.step.v[1], control->fine.step.v[2],
               control->coarse.step.v[0], control->coarse.step.v[1], control->coarse.step.v[2],
               control->thread.step.v[0], control->thread.step.v[1], control->thread.step.v[2]);
  }
  
  // Initialise selftest
  if (selftest) {
    unsigned char* restrict selftest_array;
#pragma omp single copyprivate(selftest_array)
    {
      ptrdiff_t const npoints = prod(control->ash);
      selftest_array = new unsigned char[npoints];
      memset(selftest_array, 0, npoints * sizeof *selftest_array);
    }
    control->selftest_array = selftest_array;
  } else {
    control->selftest_array = NULL;
  }
}

void lc_control_finish(lc_control_t* restrict const control,
                       lc_stats_t* restrict const stats)
{
  DECLARE_CCTK_PARAMETERS;
  
#pragma omp barrier
  
  // Finish selftest
  if (selftest) {
    assert(control->selftest_array);
#pragma omp barrier
#pragma omp single nowait
    {
      ptrdiff_t nfailed = 0;
      for (ptrdiff_t k=0; k<control->ash.v[2]; ++k) {
        for (ptrdiff_t j=0; j<control->ash.v[1]; ++j) {
          for (ptrdiff_t i=0; i<control->ash.v[0]; ++i) {
            bool const inside =
              i >= control->loop.min.v[0] and i < control->loop.max.v[0] and
              j >= control->loop.min.v[1] and j < control->loop.max.v[1] and
              k >= control->loop.min.v[2] and k < control->loop.max.v[2];
            ptrdiff_t const ipos = ind(control->ash, i,j,k);
            nfailed += control->selftest_array[ipos] != inside;
          }
        }
      }
      if (nfailed > 0) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "LoopControl self-test failed");
      }
      delete[] control->selftest_array;
    }
    control->selftest_array = NULL;
  }
  
  // Collect statistics
#pragma omp barrier
#pragma omp master
  {
    ticks const end_time = getticks();
    double const elapsed_time =
      seconds_per_tick() * elapsed(end_time, stats->start_time);
    ptrdiff_t npoints = 1;
    for (int d=0; d<LC_DIM; ++d) {
      npoints *= control->loop.max.v[d] - control->loop.min.v[d];
    }
    if (stats->init_count < 1) {
      // Skip the first iteration
      ++stats->init_count;
      if (veryverbose) {
        double const time_point =
          elapsed_time * omp_get_num_threads() / npoints;
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Loop %s: time=%g, time/point=%g s",
                   stats->name, elapsed_time, time_point);
      }
    } else {
      stats->points += double(npoints);
      stats->threads += double(omp_get_num_threads());
      stats->count += 1.0;
      stats->sum += elapsed_time;
      stats->sum2 += pow(elapsed_time, 2.0);
      stats->min = fmin(stats->min, elapsed_time);
      stats->max = fmax(stats->max, elapsed_time);
      if (veryverbose) {
        double const avg_thread = stats->sum / stats->count;
        double const avg_point =
          stats->sum * stats->threads / (stats->count * stats->points);
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Loop %s: count=%g, avg/thread=%g s, avg/point=%g s",
                   stats->name, stats->count, avg_thread, avg_point);
      }
    }
  }
  
  // Tear down multithreading state
#pragma omp single nowait
  {
    delete control->thread_info_ptr;
  }
  control->thread_info_ptr = NULL;
  if (get_fine_thread_num() == 0) {
    delete control->fine_thread_comm_ptr;
  }
  control->fine_thread_comm_ptr = NULL;
}



void lc_thread_init(lc_control_t* restrict const control)
{
  space_set_count(control->thread);
#pragma omp single
  {
    control->thread_info_ptr->idx = get_num_coarse_threads();
  }
  control->thread_done =
    space_global2local(control->thread, get_coarse_thread_num());
  space_idx2pos(control->thread);
}

int lc_thread_done(lc_control_t const* restrict const control)
{
  return control->thread_done;
}

void lc_thread_step(lc_control_t* restrict const control)
{
  // Get next thread block
  int new_global_idx;
  if (get_fine_thread_num() == 0) {
#pragma omp critical(LoopControl_lc_thread_step)
    {
      new_global_idx = control->thread_info_ptr->idx++;
    }
  }
  new_global_idx =
    fine_thread_broadcast(control->fine_thread_comm_ptr, new_global_idx);
  control->thread_done = space_global2local(control->thread, new_global_idx);
  space_idx2pos(control->thread);
}



void lc_selftest_set(lc_control_t const* restrict control,
                     ptrdiff_t const lmin, ptrdiff_t const lmax,
                     ptrdiff_t const imin, ptrdiff_t const imax,
                     ptrdiff_t const di,
                     ptrdiff_t const i0, ptrdiff_t const j, ptrdiff_t const k)
{
  DECLARE_CCTK_PARAMETERS;
  assert(selftest);
  assert(imin>=0 and imin<imax and imax<=control->ash.v[0]);
  assert(di>0);
  assert(j>=0 and j<control->ash.v[1]);
  assert(k>=0 and k<control->ash.v[2]);
  assert(i0+di-1>=control->loop.min.v[0] and i0<control->loop.max.v[0]);
  if (imin>lmin) {
    ptrdiff_t const ipos_imin = ind(control->ash, imin,j,k);
    assert(ipos_imin % di == 0);
  }
  if (imax<lmax) {
    ptrdiff_t const ipos_imax = ind(control->ash, imax,j,k);
    assert(ipos_imax % di == 0);
  }
  assert(j>=control->loop.min.v[1] and j<control->loop.max.v[1]);
  assert(k>=control->loop.min.v[2] and k<control->loop.max.v[2]);
  for (ptrdiff_t i=i0; i<i0+di; ++i) {
    if (i>=imin and i<imax) {
      assert(i>=0 and i<control->ash.v[0]);
      assert(i>=control->loop.min.v[0] and i<control->loop.max.v[0]);
      ptrdiff_t const ipos = ind(control->ash, i,j,k);
      unsigned char& elt = control->selftest_array[ipos];
#pragma omp atomic
      ++elt;
      if (elt!=1) {
#pragma omp critical
        {
          fprintf(stderr,
                  "thread=%d/%d fine_thread=%d/%d ijk=[%td,%td,%td]\n",
                  get_coarse_thread_num(), get_num_coarse_threads(),
                  get_fine_thread_num(), get_num_fine_threads(),
                  i,j,k);
          assert(elt==1);
        }
      }
    }
  }
}



int lc_setup(void)
{
  check_fortran_type_sizes();
  return 0;
}



void lc_statistics(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INFO("LoopControl statistics:");
  for (list<lc_stats_t*>::const_iterator
         istats = all_stats.begin(); istats != all_stats.end(); ++istats)
  {
    lc_stats_t const* restrict const stats = *istats;
    if (stats->count == 0.0) {
      printf("   Loop %s (%s:%d):\n",
             stats->name, stats->file, stats->line);
    } else {
      double const avg_thread = stats->sum / stats->count;
      double const avg_point =
        stats->sum * stats->threads / (stats->count * stats->points);
      printf("   Loop %s (%s:%d): count=%g, avg/thread=%g s, avg/point=%g s\n",
             stats->name, stats->file, stats->line,
             stats->count, avg_thread, avg_point);
    }
  }
}

void lc_statistics_maybe(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (not verbose) return;
  
  lc_statistics(CCTK_PASS_CTOC);
}
