/*@@
  @file      InitialData.c
  @date
  @author    Werner Benger
  @desc
             Initial data for the 3D Wave Equation
             Derived from Tom Goodale
  @enddesc
@@*/

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static CCTK_REAL sqr(CCTK_REAL val) { return val * val; }

void IDScalarWaveC_InitialData(CCTK_ARGUMENTS);

/*@@
  @routine    IDScalarWaveC_InitialData
  @date
  @author     Tom Goodale
  @desc
              Set up initial data for the wave equation
  @enddesc
  @calls
  @calledby
  @history
  @hdate Mon Oct 11 11:48:03 1999 @hauthor Werner Benger
  @hdesc  Converted to C++
  @hdate Mon Oct 11 11:48:20 1999 @hauthor Tom Goodale
  @hdesc Added the rest of the initial data.
  @hdate Thu Feb 17 09:22:01 2000 @hauthor Tom Goodale
  @hdesc Converted to C
  @endhistory

@@*/

void IDScalarWaveC_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_IDScalarWaveC_InitialData
  DECLARE_CCTK_PARAMETERS

  int i, j, k;

  CCTK_REAL dt;
  CCTK_REAL omega;
  int index;
  CCTK_REAL X, Y, Z, R;
  CCTK_REAL pi;

  dt = CCTK_DELTA_TIME;

  if (CCTK_Equals(initial_data, "plane")) {
    omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));

    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {
          index = CCTK_GFINDEX3D(cctkGH, i, j, k);

          phi[index] = amplitude * cos(kx * x[index] + ky * y[index] +
                                       kz * z[index] + omega * cctk_time);
          phi_p[index] =
              amplitude * cos(kx * x[index] + ky * y[index] + kz * z[index] +
                              omega * (cctk_time - dt));
        }
      }
    }
  } else if (CCTK_Equals(initial_data, "gaussian")) {
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {
          index = CCTK_GFINDEX3D(cctkGH, i, j, k);

          X = x[index];
          Y = y[index];
          Z = z[index];

          R = sqrt(X * X + Y * Y + Z * Z);

          phi[index] = amplitude * exp(-sqr((R - radius) / sigma));

          if (R == 0.0) {
            phi_p[index] = amplitude * (1.0 - 2.0 * dt * dt / sigma) *
                           exp(-dt * dt / sigma);
          } else {
            phi_p[index] = amplitude / 2.0 * (R - dt) / R *
                               exp(-sqr((R - radius - dt) / sigma)) +
                           amplitude / 2.0 * (R + dt) / R *
                               exp(-sqr((R - radius + dt) / sigma));
          }
        }
      }
    }
  } else if (CCTK_Equals(initial_data, "box")) {
    pi = 4.0 * atan(1.0);
    omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));

    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {
          index = CCTK_GFINDEX3D(cctkGH, i, j, k);

          phi[index] = amplitude * sin(kx * (x[index] - 0.5) * pi) *
                       sin(ky * (y[index] - 0.5) * pi) *
                       sin(kz * (z[index] - 0.5) * pi) *
                       cos(omega * cctk_time * pi);

          phi_p[index] = amplitude * sin(kx * (x[index] - 0.5) * pi) *
                         sin(ky * (y[index] - 0.5) * pi) *
                         sin(kz * (z[index] - 0.5) * pi) *
                         cos(omega * (cctk_time - dt) * pi);
        }
      }
    }
  } else if (CCTK_Equals(initial_data, "none")) {
    for (k = 0; k < cctk_lsh[2]; k++) {
      for (j = 0; j < cctk_lsh[1]; j++) {
        for (i = 0; i < cctk_lsh[0]; i++) {
          index = CCTK_GFINDEX3D(cctkGH, i, j, k);

          phi[index] = 0.0;

          phi_p[index] = 0.0;
        }
      }
    }
  }

  return;
}
