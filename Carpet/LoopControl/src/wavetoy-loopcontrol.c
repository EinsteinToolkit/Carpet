// Evolve the scalar wave equation with Dirichlet boundaries
// (C) 2007-08-16 Erik Schnetter <schnetter@cct.lsu.edu>



// Redshift:
// ~/gcc/bin/gcc -std=gnu99 -fopenmp -Wall -g -O3 -fomit-frame-pointer -march=pentium4 -malign-double -o wavetoy-loopcontrol loopcontrol.c wavetoy-loopcontrol.c
// icc -std=c99 -openmp -Wall -g -fast -march=pentium4 -align -o wavetoy-loopcontrol loopcontrol.c wavetoy-loopcontrol.c

// Abe:
// icc -restrict -openmp -Wall -g -fast -fomit-frame-pointer -align -o wavetoy-loopcontrol loopcontrol.c wavetoy-loopcontrol.c

// Ducky:
// xlc_r -q64 -qlanglvl=stdc99 -qsmp=omp -g -O3 -qmaxmem=-1 -qhot -qarch=pwr5 -qtune=pwr5 -o wavetoy wavetoy.c

// Eric:
// icc -std=c99 -openmp -Wall -g -fast -fomit-frame-pointer -align -o wavetoy wavetoy.c
// LM_LICENSE_FILE=/usr/local/compilers/pgi/license.dat /usr/local/compilers/pgi/linux86/6.1/bin/pgcc -c9x -mp -g -fastsse -tp=piv -Mdalign -Mllalign -o wavetoy wavetoy.c

// Prism:
// /opt/intel/8.1.023/bin/icc -std=c99 -openmp -Wall -g -fast -o wavetoy wavetoy.c

// Santaka:
// icc -std=c99 -openmp -Wall -g -fast -fomit-frame-pointer -o wavetoy wavetoy.c



#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "loopcontrol.h"


/***
 define this if you want to build a standalone demo
 (but beware of hidden LoopControl dependencies on Cactus)
#define BUILD_STANDALONE
 ***/



#ifndef _OPENMP
/* Replacements for some OpenMP routines if OpenMP is not available */

static inline
int
omp_get_thread_num (void)
{
  return 0;
}

static inline
int
omp_get_num_threads (void)
{
  return 1;
}

static inline
double
omp_get_wtime (void)
{
  struct timeval tv;
  gettimeofday (& tv, NULL);
  return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

#endif



static int NI;
static int NJ;
static int NK;

static int NSTEPS;



#ifdef BUILD_STANDALONE
static
int
getint (char const * restrict const string)
{
  char * endptr;
  errno = 0;
  long const n = strtol (string, & endptr, 0);
  if (* string == '\0' || * endptr != '\0' || errno != 0) {
    fprintf (stderr, "Argument \"%s\" is not a legal number\n", string);
    if (errno != 0) {
      perror (NULL);
    }
    exit (2);
  }
  return n;
}
#else

#include "cctk_Parameters.h"

#endif



static inline
int
ind3d (int const i, int const j, int const k)
{
  return i + NI * (j + NJ * k);
}



#ifdef BUILD_STANDALONE
int main (int argc, char **argv)
{
#else
int lc_demo (void)
{
  DECLARE_CCTK_PARAMETERS
#endif
  printf ("WaveToy\n");
#if 0
  if (argc != 5) {
    fprintf (stderr,
             "Synopsis:\n"
             "   %s <NI> <NJ> <NK> <NSTEPS>\n"
             "      where <NI> <NJ> <NK>:       number of grid points\n"
             "            <NSTEPS>:             number of time steps\n",
             argv[0]);
    exit (1);
  }
  NI = getint (argv[1]);
  NJ = getint (argv[2]);
  NK = getint (argv[3]);
  NSTEPS = getint (argv[4]);
#else
  NI = nx;
  NJ = ny;
  NK = nz;
  NSTEPS = nsteps;
#endif
  printf ("   NI=%d NJ=%d NK=%d\n", NI, NJ, NK);
  printf ("   NSTEPS=%d\n", NSTEPS);
  
  //
  // Statistics
  //
#pragma omp parallel
  {
#pragma omp single
    {
      int const num_threads = omp_get_num_threads();
      printf ("Running with %d threads\n", num_threads);
    }
  }
  
  //
  // Setup
  //
  // double const time_setup_start = omp_get_wtime();
  double const alpha = 0.5;     // CFL factor
  double * restrict phi0, * restrict phi1, * restrict phi2;
  phi0 = malloc (NI*NJ*NK * sizeof *phi0);
  phi1 = malloc (NI*NJ*NK * sizeof *phi1);
  phi2 = malloc (NI*NJ*NK * sizeof *phi2);
  // double const time_setup_end = omp_get_wtime();
  // printf ("Setup time: %g\n", time_setup_end - time_setup_start);
  
  //
  // Initialise
  //
  printf ("Initialisation\n");
  double const time_init_start = omp_get_wtime();
  double const kx = M_PI/(NI-1);
  double const ky = M_PI/(NJ-1);
  double const kz = M_PI/(NK-1);
  double const kt = sqrt (pow(kx,2) + pow(ky,2) + pow(kz,2));
#pragma omp parallel
  {
    LC_LOOP3 (initialisation, i,j,k, 0,0,0, NI,NJ,NK, NI,NJ,NK) {
      int const ind = ind3d(i,j,k);
      // phi0[ind] = 0.0;
      // phi0[ind] = i==1 && j==1 && k==1 ? 1.0 : 0.0;
      phi0[ind] = 0.0;
      phi1[ind] = sin(kx*i) * sin(ky*j) * sin(kz*k) * sin(kt*(-  alpha));
      phi2[ind] = sin(kx*i) * sin(ky*j) * sin(kz*k) * sin(kt*(-2*alpha));
    } LC_ENDLOOP3 (initialisation);
  }
  // phi0[ind3d(NI/2,NJ/2,NK/2)] = 1.0;
  double const time_init_end = omp_get_wtime();
  printf ("Initialisation time: %g\n", time_init_end - time_init_start);
  
  //
  // Evolve
  //
  printf ("Evolution\n");
  double const time_evol_start = omp_get_wtime();
  for (int step=1; step<=NSTEPS; ++step) {
    // if (step % 10 == 0) {
    //   printf ("Step %d\n", step);
    // }
    
    //
    // Rotate
    //
    {
      double * const tmp=phi2; phi2=phi1; phi1=phi0; phi0=tmp;
    }
    
    //
    // Step
    //
    // double const time_step_start = omp_get_wtime();
    double const alpha2 = pow(alpha,2);
    double const factor = 2 * (1 - 3 * alpha2);
    int const di = ind3d(1,0,0), dj = ind3d(0,1,0), dk=ind3d(0,0,1);
#pragma omp parallel
    {
      LC_LOOP3 (evolution, i,j,k, 1,1,1, NI-1,NJ-1,NK-1, NI,NJ,NK) {
        int const ind = ind3d(i,j,k);
        phi0[ind] = (+ factor*phi1[ind] - phi2[ind]
                     + alpha2 * (+ phi1[ind-di] + phi1[ind+di]
                                 + phi1[ind-dj] + phi1[ind+dj]
                                 + phi1[ind-dk] + phi1[ind+dk]));
        // phi0[ind] = phi1[ind] + 1.0;
      } LC_ENDLOOP3 (evolution);
    } // omp parallel
    // double const time_step_end = omp_get_wtime();
    // printf ("Step time: %g\n", time_step_end - time_step_start);
    
  }
  double const time_evol_end = omp_get_wtime();
  printf ("Evolution time: %g\n", time_evol_end - time_evol_start);
  
  //
  // Analyse
  //
  printf ("Analysis\n");
  double const time_analysis_start = omp_get_wtime();
  double sum = 0.0;
#pragma omp parallel reduction(+: sum)
  {
    LC_LOOP3 (analysis, i,j,k, 1,1,1, NI-1,NJ-1,NK-1, NI,NJ,NK) {
      int const ind = ind3d(i,j,k);
      sum += phi0[ind];
    } LC_ENDLOOP3 (analysis);
  }
  double const avg = sum / ((NI-2)*(NJ-2)*(NK-2));
  double const time_analysis_end = omp_get_wtime();
  printf ("Analysis time: %g\n", time_analysis_end - time_analysis_start);
  printf ("Result: %.15g\n", avg);
  
  lc_printstats(0);
  
  //
  // Shutdown
  //
  return 0;
}
