#ifndef LOOPCONTROL_H
#define LOOPCONTROL_H

/* This file uses the namespace LC_* for macros and lc_* for C
   identifiers. */

#define LC_DIM 3
  


#ifdef CCODE

#include <stddef.h>
#include <stdlib.h>

#include <cctk.h>

#ifdef __cplusplus
extern "C" {
#endif
  
  
  
  static inline ptrdiff_t lc_min(ptrdiff_t const i, ptrdiff_t const j)
  {
    return i < j ? i : j;
  }
  
  static inline ptrdiff_t lc_max(ptrdiff_t const i, ptrdiff_t const j)
  {
    return i > j ? i : j;
  }
  
  
  
  struct lc_thread_info_t;
  struct lc_fine_thread_comm_t;
  
  struct lc_stats_t;
  
  typedef struct {
    ptrdiff_t v[LC_DIM];
  } lc_vec_t;
  
  typedef struct {
    int v[LC_DIM];
  } lc_ivec_t;
  
  typedef struct {
    /* Traverse pos from min (inclusive) to max (exclusive) with a
       stride of step. Equivalently, traverse idx from 0 (inclusive)
       to count (exclusive). */
    lc_vec_t min, max, step, pos;
    lc_ivec_t count, idx;
  } lc_space_t;
  
  typedef struct {
    lc_vec_t ash;
    lc_space_t loop;
    lc_space_t thread;
    struct lc_thread_info_t* thread_info_ptr; /* shared between all
                                                 threads */
    int thread_done;
    lc_space_t coarse;          /* count, idx, pos are undefined */
    lc_space_t fine;            /* count, idx, pos are undefined */
    lc_space_t fine_thread;     /* min, max, pos are undefined */
    struct lc_fine_thread_comm_t* fine_thread_comm_ptr; /* shared
                                                           between SMT
                                                           threads */
    unsigned char* restrict selftest_array;
  } lc_control_t;
  
  
  
  void lc_stats_init(struct lc_stats_t** stats,
                     char const* name, char const* file, int line);
  void lc_control_init(lc_control_t* restrict control,
                       struct lc_stats_t *restrict stats,
                       ptrdiff_t imin, ptrdiff_t jmin, ptrdiff_t kmin,
                       ptrdiff_t imax, ptrdiff_t jmax, ptrdiff_t kmax,
                       ptrdiff_t iash, ptrdiff_t jash, ptrdiff_t kash,
                       ptrdiff_t di, ptrdiff_t dj, ptrdiff_t dk);
  void lc_control_finish(lc_control_t* restrict control,
                         struct lc_stats_t *restrict stats);
  
  void lc_thread_init(lc_control_t* restrict control);
  int lc_thread_done(lc_control_t const* restrict control);
  void lc_thread_step(lc_control_t* restrict control);
  
  void lc_selftest_set(lc_control_t const* restrict control,
                       ptrdiff_t lmin, ptrdiff_t lmax,
                       ptrdiff_t imin, ptrdiff_t imax, ptrdiff_t di,
                       ptrdiff_t i, ptrdiff_t j, ptrdiff_t k);
  
  
  
#define LC_COARSE_SETUP(D)                                              \
  lc_control.coarse.min.v[D] = lc_control.thread.pos.v[D];              \
  lc_control.coarse.max.v[D] =                                          \
    lc_min(lc_control.thread.max.v[D],                                  \
           lc_control.coarse.min.v[D] + lc_control.thread.step.v[D]);   \
  ptrdiff_t const lc_cmin##D = lc_control.coarse.min.v[D];              \
  ptrdiff_t const lc_cmax##D = lc_control.coarse.max.v[D];              \
  ptrdiff_t const lc_cstep##D = lc_control.coarse.step.v[D];
#define LC_COARSE_LOOP(D)                       \
  for (ptrdiff_t lc_cpos##D = lc_cmin##D;       \
       lc_cpos##D < lc_cmax##D;                 \
       lc_cpos##D += lc_cstep##D)
  
#define LC_FINE_SETUP(D)                                                \
  lc_control.fine.min.v[D] = lc_cpos##D;                                \
  lc_control.fine.max.v[D] =                                            \
    lc_min(lc_control.coarse.max.v[D],                                  \
           lc_control.fine.min.v[D] + lc_control.coarse.step.v[D]);     \
  ptrdiff_t /*const*/ lc_fmin##D = lc_control.fine.min.v[D];            \
  ptrdiff_t /*const*/ lc_fmax##D = lc_control.fine.max.v[D];            \
  ptrdiff_t const lc_fstep##D = lc_control.fine.step.v[D];              \
  ptrdiff_t const lc_ftoff##D =                                         \
    lc_control.fine_thread.idx.v[D] * lc_control.fine_thread.step.v[D];
#define LC_FINE_LOOP(I, NI, D)                          \
  for (ptrdiff_t I = lc_fmin##D + lc_ftoff##D;          \
       I < lc_fmax##D;                                  \
       I += lc_fstep##D)                                \
  {                                                     \
    ptrdiff_t const NI CCTK_ATTRIBUTE_UNUSED =          \
      CCTK_BUILTIN_EXPECT(lc_dir##D==0, 1) ? 0 :        \
       lc_dir##D<0 ? I+1 : lc_control.loop.max.v[D]-I;
  
#if VECTORISE_ALIGNED_ARRAYS
  /* Arrays are aligned: fmin0 is the aligned loop boundary; keep it,
     and set up imin to be the intended loop boundary */
#  define LC_ALIGN(i,j,k, imin,imax)                                    \
  ptrdiff_t const imin = lc_max(lc_control.loop.min.v[0], lc_fmin0);    \
  ptrdiff_t const imax = lc_fmax0;
#else
  /* Arrays are not aligned: fine.min[0] and fine.max[0] are the
     intended loop boundaries; override fmin0 and fmax0 to be aligned,
     and set imin and imax; this may move the fine loop boundaries
     except at outer boundaries to avoid partial stores */
#  define LC_ALIGN(i,j,k, imin,imax)                                    \
  lc_fmin0 = lc_control.fine.min.v[0];                                  \
  lc_fmax0 = lc_control.fine.max.v[0];                                  \
  ptrdiff_t imin = lc_fmin0;                                            \
  ptrdiff_t imax = lc_fmax0;                                            \
  int const lc_fmin0_is_outer = lc_fmin0 == lc_control.loop.min.v[0];   \
  int const lc_fmax0_is_outer = lc_fmax0 == lc_control.loop.max.v[0];   \
  ptrdiff_t const lc_iminpos = lc_fmin0 + lc_ash0 * (j + lc_ash1 * k);  \
  ptrdiff_t const lc_iminoffset = lc_iminpos % lc_align0;               \
  ptrdiff_t const lc_imaxpos = lc_fmax0 + lc_ash0 * (j + lc_ash1 * k);  \
  ptrdiff_t const lc_imaxoffset = lc_imaxpos % lc_align0;               \
  lc_fmin0 -= lc_iminoffset;                                            \
  if (!lc_fmax0_is_outer) lc_fmax0 -= lc_imaxoffset;                    \
  if (!lc_fmin0_is_outer) imin = lc_fmin0;                              \
  if (!lc_fmax0_is_outer) imax = lc_fmax0;
#endif
  
#define LC_SELFTEST(i,j,k, imin,imax)                                   \
  if (CCTK_BUILTIN_EXPECT(lc_control.selftest_array != NULL, 0)) {      \
    lc_selftest_set(&lc_control,                                        \
                    lc_control.loop.min.v[0], lc_control.loop.max.v[0], \
                    imin, imax, lc_align0, i, j, k);                    \
  }
  
  
  
#define LC_LOOP3STR_NORMAL(name, i,j,k, ni,nj,nk,               \
                           idir_, jdir_, kdir_,                 \
                           imin_,jmin_,kmin_,                   \
                           imax_,jmax_,kmax_,                   \
                           iash_,jash_,kash_,                   \
                           imin,imax, di_)                      \
  do {                                                          \
    typedef int lc_loop3vec_##name;                             \
                                                                \
    ptrdiff_t const lc_dir0 CCTK_ATTRIBUTE_UNUSED = (idir_);    \
    ptrdiff_t const lc_dir1 CCTK_ATTRIBUTE_UNUSED = (jdir_);    \
    ptrdiff_t const lc_dir2 CCTK_ATTRIBUTE_UNUSED = (kdir_);    \
                                                                \
    ptrdiff_t const lc_ash0 CCTK_ATTRIBUTE_UNUSED = (iash_);    \
    ptrdiff_t const lc_ash1 CCTK_ATTRIBUTE_UNUSED = (jash_);    \
    ptrdiff_t const lc_ash2 CCTK_ATTRIBUTE_UNUSED = (kash_);    \
                                                                \
    ptrdiff_t const lc_align0 CCTK_ATTRIBUTE_UNUSED = (di_);    \
    ptrdiff_t const lc_align1 CCTK_ATTRIBUTE_UNUSED = 1;        \
    ptrdiff_t const lc_align2 CCTK_ATTRIBUTE_UNUSED = 1;        \
                                                                \
    static struct lc_stats_t* lc_stats = NULL;                  \
    lc_stats_init(&lc_stats, #name, __FILE__, __LINE__);        \
                                                                \
    lc_control_t lc_control;                                    \
    lc_control_init(&lc_control, lc_stats,                      \
                    (imin_), (jmin_), (kmin_),                  \
                    (imax_), (jmax_), (kmax_),                  \
                    lc_ash0, lc_ash1, lc_ash2,                  \
                    lc_align0, lc_align1, lc_align2);           \
                                                                \
    /* Multithreading */                                        \
    for (lc_thread_init(&lc_control);                           \
         !lc_thread_done(&lc_control);                          \
         lc_thread_step(&lc_control))                           \
    {                                                           \
                                                                \
      /* Coarse loops */                                        \
      LC_COARSE_SETUP(2)                                        \
      LC_COARSE_SETUP(1)                                        \
      LC_COARSE_SETUP(0)                                        \
      LC_COARSE_LOOP(2) {                                       \
      LC_COARSE_LOOP(1) {                                       \
      LC_COARSE_LOOP(0) {                                       \
                                                                \
        /* Fine loops */                                        \
        LC_FINE_SETUP(2)                                        \
        LC_FINE_SETUP(1)                                        \
        LC_FINE_SETUP(0)                                        \
        LC_FINE_LOOP(k, nk, 2) {                                \
        LC_FINE_LOOP(j, nj, 1) {                                \
        LC_ALIGN(i,j,k, imin,imax)                              \
        LC_FINE_LOOP(i, ni, 0) {                                \
          LC_SELFTEST(i,j,k, imin,imax)                         \
          {
        
#define LC_ENDLOOP3STR_NORMAL(name)                             \
          }                     /* body */                      \
        }}}}}}                  /* fine */                      \
      }}}                       /* coarse */                    \
    }                           /* multithreading */            \
    lc_control_finish(&lc_control, lc_stats);                   \
    typedef lc_loop3vec_##name lc_ensure_proper_nesting;        \
  } while(0)
    
    
    
/* Definitions to ensure compatibility with earlier versions of
   LoopControl */
#define LC_LOOP3VEC(name, i,j,k,                        \
                    imin_,jmin,kmin, imax_,jmax,kmax,   \
                    iash,jash,kash, imin,imax, di)      \
  LC_LOOP3STR_NORMAL(name, i,j,k, lc_ni,lc_nj,lc_nk,    \
                     0,0,0,                             \
                     imin_,jmin,kmin, imax_,jmax,kmax,  \
                     iash,jash,kash, imin,imax, di)
#define LC_ENDLOOP3VEC(name)                    \
  LC_ENDLOOP3STR_NORMAL(name)

#define LC_LOOP3(name, i,j,k,                           \
                 imin,jmin,kmin, imax,jmax,kmax,        \
                 iash,jash,kash)                        \
  LC_LOOP3VEC(name, i,j,k,                              \
              imin,jmin,kmin, imax,jmax,kmax,           \
              iash,jash,kash, lc_imin_,lc_imax_, 1)
#define LC_ENDLOOP3(name)                       \
  LC_ENDLOOP3VEC(name)
  
  
  
/* Replace CCTK_LOOP macros */
#if !defined CCTK_LOOP3STR_NORMAL || !defined CCTK_ENDLOOP3STR_NORMAL
#  error "internal error"
#endif
#undef CCTK_LOOP3STR_NORMAL
#undef CCTK_ENDLOOP3STR_NORMAL
#define CCTK_LOOP3STR_NORMAL(name, i,j,k, ni,nj,nk,     \
                             idir, jdir, kdir,          \
                             imin_,jmin,kmin,           \
                             imax_,jmax,kmax,           \
                             iash,jash,kash,            \
                             imin,imax, di)             \
  LC_LOOP3STR_NORMAL(name, i,j,k, ni,nj,nk,             \
                     idir, jdir, kdir,                  \
                     imin_,jmin,kmin,                   \
                     imax_,jmax,kmax,                   \
                     iash,jash,kash,                    \
                     imin,imax, di)
#define CCTK_ENDLOOP3STR_NORMAL(name)           \
  LC_ENDLOOP3STR_NORMAL(name)
  
  
  
#ifdef __cplusplus 
    }
#endif

#endif /* #ifdef CCODE */

#ifdef FCODE
#  include "loopcontrol_fortran.h"
#endif

#endif  /* #ifndef LOOPCONTROL_H */
