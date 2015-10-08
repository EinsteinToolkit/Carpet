#include <cmath>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "mask_boxexcluded.hh"

#include <CarpetReduce_bits.h>

namespace CarpetMask {

using namespace std;

// Number of excluded regions which are defined in param.ccl
int const num_excluded = 10;

/**
 * Set the weight in the excluded regions to zero.
 */

void CarpetBoxExcludedSetup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Some state verbose output
  static vector<int> last_output;
  if (last_output.empty()) {
    last_output.resize(num_excluded, -1);
  }

  CCTK_INT *restrict const iweight = static_cast<CCTK_INT *>(
      CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::iweight"));

  CCTK_REAL *restrict const excised_cells = static_cast<CCTK_REAL *>(
      CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::excised_cells"));

  if (not iweight or not excised_cells) {
    CCTK_WARN(CCTK_WARN_ABORT, "CarpetReduce is not active, or "
                               "CarpetReduce::iweight does not have storage");
  }

  // Volume of a grid cell on this level, in terms of coarse grid
  // cells
  CCTK_REAL const cell_volume =
      1.0 / (cctk_levfac[0] * cctk_levfac[1] * cctk_levfac[2]);

  unsigned const bits = BMSK(cctk_dim);
  CCTK_REAL const factor = 1.0 / bits;

  for (int n = 0; n < 10; ++n) {

    CCTK_REAL const r0_x = box_excluded_radius_x[n];
    CCTK_REAL const r0_y = box_excluded_radius_y[n];
    CCTK_REAL const r0_z = box_excluded_radius_z[n];

    if (r0_x >= 0.0 && r0_y >= 0.0 && r0_z >= 0.0) {

      CCTK_REAL const x0 = box_excluded_centre_x[n];
      CCTK_REAL const y0 = box_excluded_centre_y[n];
      CCTK_REAL const z0 = box_excluded_centre_z[n];

      bool const exterior = box_exclude_exterior[n];

      if (verbose) {
        if (cctk_iteration > last_output.at(n)) {
          last_output.at(n) = cctk_iteration;
          CCTK_VInfo(CCTK_THORNSTRING, "Setting mask from excluded boxes at "
                                       "[%g,%g,%g], with radii [%g,%g,%g], at "
                                       "iteration %d",
                     double(x0), double(y0), double(z0), double(r0_x),
                     double(r0_y), double(r0_z), int(cctk_iteration));
        }
      }

      CCTK_REAL local_excised = 0.0;
#pragma omp parallel reduction(+ : local_excised)
      CCTK_LOOP3_ALL(CarpetExcludedSetup, cctkGH, i, j, k) {
        int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL const dx = fabs(x[ind] - x0);
        CCTK_REAL const dy = fabs(y[ind] - y0);
        CCTK_REAL const dz = fabs(z[ind] - z0);

        if (exterior
                ?
                /* Exclude corners */
                (dx >= r0_x && dy >= r0_y && dz >= r0_z) ||
                    /* Exclude top and bottom edges */
                    (dx >= r0_x && dy <= r0_y && dz >= r0_z) ||
                    (dx <= r0_x && dy >= r0_y && dz >= r0_z) ||
                    /* Exclude top and bottom faces */
                    (dx <= r0_x && dy <= r0_y && dz >= r0_z) ||
                    /* Exclude edges along z */
                    (dx >= r0_x && dy >= r0_y && dz <= r0_z) ||
                    /* Exclude front and back faces */
                    (dx >= r0_x && dy <= r0_y && dz <= r0_z) ||
                    /* Exclude side faces */
                    (dx <= r0_x && dy >= r0_y && dz <= r0_z)
                :
                /* Exclude inside box */
                dx <= r0_x && dy <= r0_y && dz <= r0_z) {
          // Tally up the weight we are removing
          local_excised += cell_volume * factor * BCNT(iweight[ind]);
          iweight[ind] = 0;
        }
      }
      CCTK_ENDLOOP3_ALL(CarpetExcludedSetup);
      *excised_cells += local_excised;

    } // if r>=0
  }   // for n
}

} // namespace CarpetMask
