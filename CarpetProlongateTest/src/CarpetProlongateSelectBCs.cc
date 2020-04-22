#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <assert.h>

extern "C"
void CarpetProlongateTest_SelectBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  char const *const varnames[] = {"CARPETPROLONGATETEST::u",
                                  "CARPETPROLONGATETEST::u0",
                                  "CARPETPROLONGATETEST::du"};
  for (auto var : varnames) {
    int ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES,
                                       cctk_nghostzones[0], -1, var, "None");
    assert(ierr != -1);
  }
}
