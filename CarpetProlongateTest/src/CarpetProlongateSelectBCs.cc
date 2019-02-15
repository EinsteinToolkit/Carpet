#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <string>
#include <vector>
#include <assert.h>

void CarpetProlongateSelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  std::vector<std::string> v{
    "CARPETPROLONGATETEST::u",
    "CARPETPROLONGATETEST::u0",
    "CARPETPROLONGATETEST::du"
  };
  for(auto var : v) {
    int ierr = Boundary_SelectVarForBC(
      cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0],
      -1, var.c_str(), "None");
    assert(ierr != -1);
  }
}
