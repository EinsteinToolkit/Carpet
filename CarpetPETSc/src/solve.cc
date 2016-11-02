#include <cassert>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "TATelliptic.h"
#include "carpet.hh"

#include "dist.hh"

#include <petscvec.h>

#ifdef CCTK_CXX_RESTRICT
#define restrict CCTK_CXX_RESTRICT
#endif

namespace CarpetPETSc {

using namespace std;
using namespace Carpet;

char const *const solver_name = "CarpetPETSc";

struct options_t {

  // Width of outer boundary
  CCTK_INT bndwidth;

  // Symmetry information
  CCTK_INT symmetry_handle[2 * dim];
  CCTK_INT symmetry_zone_width[2 * dim];
};

void copy_in(cGH const *restrict const cctkGH, Vec dst,
             vector<CCTK_INT> const &src, options_t const &options);

void copy_out(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
              Vec const src, options_t const &options);

int solve(cGH const *const cctkGH, int const *const var, int const *const val,
          int const nvars, int const options_table,
          int (*const fun)(cGH const *cctkGH, int options_table, void *data),
          int (*const bnd)(cGH const *cctkGH, int options_table, void *data),
          void *const data) {
  DECLARE_CCTK_PARAMETERS;
}

void copy_in(cGH const *restrict const cctkGH, Vec dst,
             vector<CCTK_INT> const &src, options_t const &options) {
  double *restrict dstptr;
  VecGetArray(dst, &dstptr);
  PetscInt dstsize;
  VecGetSize(dst, &dstsize);

  int dstind = 0;

  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;

      vector<CCTK_REAL const *> srcptrs(src.size());
      for (int n = 0; n < src.size(); ++n) {
        srcptrs.at(n) = static_cast<CCTK_REAL const *>(
            CCTK_VarDataPtrI(cctkGH, 0, src.at(n)));
        assert(srcptrs.at(n));
      }

      int imin[3], imax[3];
      interior_shape(cctkGH, options, imin, imax);
      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
          for (int i = imin[0]; i < imax[0]; ++i) {
            int const srcind = CCTK_GFINDEX3D(cctkGH, i, j, k);

            for (int n = 0; n < src.size(); ++n) {
              assert(dstind < dstsize);
              dstptr[dstind] = srcptr[srcind];
              ++dstind;
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
}

void copy_out(cGH const *restrict const cctkGH, vector<CCTK_INT> const &dst,
              Vec const src, options_t const &options) {
  double const *restrict srcptr;
  VecGetArray(src, &srcptr);
  PetscInt srcsize;
  VecGetSize(dst, &srcsize);

  int srcind = 0;

  BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;

      vector<CCTK_REAL *> dstptrs(dst.size());
      for (int n = 0; n < dst.size(); ++n) {
        dstptrs.at(n) =
            static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, dst.at(n)));
        assert(dstptrs.at(n));
      }

      int imin[3], imax[3];
      interior_shape(cctkGH, options, imin, imax);
      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
          for (int i = imin[0]; i < imax[0]; ++i) {
            int const dstind = CCTK_GFINDEX3D(cctkGH, i, j, k);

            for (int n = 0; n < src.size(); ++n) {
              assert(srcind < srcsize);
              dstptr[dstind] = srcptr[srcind];
              ++srcind;
            }
          }
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_MAP_LOOP;
}

// Calculate the shape of the interior
void interior_shape(cGH const *restrict const cctkGH, options_t const &options,
                    int *restrict const imin, int *restrict const imax) {
  DECLARE_CCTK_ARGUMENTS;

  for (int d = 0; d < dim; ++d) {

    imin[d] = 0;
    if (cctk_bbox[2 * d]) {
      if (options.symmetry_handle[2 * d] >= 0) {
        imin[d] += options.symmetry_zone_width[2 * d];
      } else {
        imin[d] += options.bndwidth;
      }
    } else {
      imin[d] += cctk_nghostzones[d];
    }

    imax[d] = cctk_lsh[d];
    if (cctk_bbox[2 * d + 1]) {
      if (options.symmetry_handle[2 * d + 1] >= 0) {
        imax[d] -= options.symmetry_zone_width[2 * d + 1];
      } else {
        imax[d] -= options.bndwidth;
      }
    } else {
      imax[d] -= cctk_nghostzones[d];
    }

    assert(imin[d] <= imax[d]);

  } // for d
}

} // namespace CarpetPETSc
