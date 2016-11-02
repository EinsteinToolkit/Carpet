#include <cassert>
#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

using std::printf;
#define then /* empty */

// prototypes
namespace {
void do_singlemap_stuff(cGH *GH, int map_number, bool test_local_mode);
void do_local_stuff(cGH *GH, int map_number, int component_number);
};

//******************************************************************************

// this function is called in LEVEL mode by the Cactus scheduler
extern "C" void TestCarpetGridInfo_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  cGH *GH = cctkGH;

  assert(Carpet::is_level_mode());

  for (int map_number = 0; map_number < N_maps; ++map_number) {
    // switch Carpet to this map
    Carpet::enter_singlemap_mode(GH, map_number, CCTK_GF);

    do_singlemap_stuff(GH, map_number, (test_local_mode != 0));

    Carpet::leave_singlemap_mode(GH);
  }
}

//******************************************************************************

// this function is called in SINGLEMAP mode for each patch (= Carpet map)
namespace {
void do_singlemap_stuff(cGH *GH, const int map_number, bool test_local_mode) {
  assert(Carpet::is_singlemap_mode());

  cGH *cctkGH = GH;       // magic variable name for CCTK_*() macros
  DECLARE_CCTK_ARGUMENTS; // set up for CCTK_*() macros

  printf("in singlemap mode, map_number=%d\n", map_number);
  printf("   cctk_origin_space[]=(%g,%g,%g)\n",
         double(GH->cctk_origin_space[0]), double(GH->cctk_origin_space[1]),
         double(GH->cctk_origin_space[2]));
  printf("   cctk_delta_space[]=(%g,%g,%g)\n", double(GH->cctk_delta_space[0]),
         double(GH->cctk_delta_space[1]), double(GH->cctk_delta_space[2]));
  printf("   CCTK_ORIGIN_SPACE()=(%g,%g,%g)\n", double(CCTK_ORIGIN_SPACE(0)),
         double(CCTK_ORIGIN_SPACE(1)), double(CCTK_ORIGIN_SPACE(2)));
  printf("   CCTK_DELTA_SPACE()=(%g,%g,%g)\n", double(CCTK_DELTA_SPACE(0)),
         double(CCTK_DELTA_SPACE(1)), double(CCTK_DELTA_SPACE(2)));

  if (test_local_mode)
    then {
      using namespace Carpet; // needed for magic Carpet macros below
      BEGIN_LOCAL_COMPONENT_LOOP(GH, CCTK_GF) {
        assert(Carpet::is_local_mode());
        const int component_number = Carpet::component;

        do_local_stuff(GH, map_number, component_number);
      }
      END_LOCAL_COMPONENT_LOOP;
    }
}
}

//******************************************************************************

namespace {
void do_local_stuff(cGH *GH, int map_number, int component_number) {
  assert(Carpet::is_local_mode());

  cGH *cctkGH = GH;       // magic variable name for CCTK_*() macros
  DECLARE_CCTK_ARGUMENTS; // set up for CCTK_*() macros

  printf("in local mode, map_number=%d component_number=%d:\n", map_number,
         component_number);

  printf("   cctk_origin_space[]=(%g,%g,%g)\n",
         double(GH->cctk_origin_space[0]), double(GH->cctk_origin_space[1]),
         double(GH->cctk_origin_space[2]));
  printf("   cctk_delta_space[]=(%g,%g,%g)\n", double(GH->cctk_delta_space[0]),
         double(GH->cctk_delta_space[1]), double(GH->cctk_delta_space[2]));
  printf("   CCTK_ORIGIN_SPACE()=(%g,%g,%g)\n", double(CCTK_ORIGIN_SPACE(0)),
         double(CCTK_ORIGIN_SPACE(1)), double(CCTK_ORIGIN_SPACE(2)));
  printf("   CCTK_DELTA_SPACE()=(%g,%g,%g)\n", double(CCTK_DELTA_SPACE(0)),
         double(CCTK_DELTA_SPACE(1)), double(CCTK_DELTA_SPACE(2)));
}
}
