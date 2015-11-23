#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>
#include <loopcontrol.h>

#include "iof5.hh"

namespace CarpetIOF5 {

using namespace std;
using namespace Carpet;

extern "C" void F5_Poison(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_global_mode());
  CCTK_VInfo(CCTK_THORNSTRING, "F5_Poison: iteration=%d",
             cctkGH->cctk_iteration);

  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
#pragma omp parallel
        CCTK_LOOP3_ALL(F5_Poison, cctkGH, i, j, k) {
          int const ind3d = CCTK_GFINDEX3D(cctkGH, i, j, k);
          x[ind3d] = nan;
          y[ind3d] = nan;
          z[ind3d] = nan;
          r[ind3d] = nan;
        }
        CCTK_ENDLOOP3_ALL(F5_Poison);
      }
      END_LOCAL_COMPONENT_LOOP;
    }
    END_LOCAL_MAP_LOOP;
  }
  END_REFLEVEL_LOOP;
}

extern "C" void F5_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_global_mode());
  CCTK_VInfo(CCTK_THORNSTRING, "F5_Check: iteration=%d",
             cctkGH->cctk_iteration);

  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
#pragma omp parallel
        CCTK_LOOP3_ALL(F5_Check, cctkGH, i, j, k) {
          int const ind3d = CCTK_GFINDEX3D(cctkGH, i, j, k);
          assert(not isnan(x[ind3d]));
          assert(not isnan(y[ind3d]));
          assert(not isnan(z[ind3d]));
          assert(not isnan(r[ind3d]));
        }
        CCTK_ENDLOOP3_ALL(F5_Check);
      }
      END_LOCAL_COMPONENT_LOOP;
    }
    END_LOCAL_MAP_LOOP;
  }
  END_REFLEVEL_LOOP;
}

} // end namespace CarpetIOF5
