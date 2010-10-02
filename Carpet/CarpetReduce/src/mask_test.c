#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <assert.h>
#include <math.h>



void
MaskBase_TestMask (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose) {
    CCTK_INFO ("Testing weight");
  }
  
  
  
  int const sum = CCTK_ReductionHandle ("sum");
  assert (sum >= 0);
  
  int const proc = 0;
  
  int const weight_var = CCTK_VarIndex ("CarpetReduce::weight");
  int const one_var    = CCTK_VarIndex ("CarpetReduce::one");
  assert (weight_var >= 0);
  assert (one_var    >= 0);
  
  CCTK_REAL sum_weight;
  
  {
    int const ierr = CCTK_Reduce (cctkGH,
                                  proc,
                                  sum,
                                  1, CCTK_VARIABLE_REAL, &sum_weight,
                                  1, one_var);
    assert (ierr >= 0);
  }
  
  
  
  if (proc == -1 || CCTK_MyProc(cctkGH) == proc) {
    
    CCTK_REAL physical_min[cctk_dim];
    CCTK_REAL physical_max[cctk_dim];
    CCTK_REAL interior_min[cctk_dim];
    CCTK_REAL interior_max[cctk_dim];
    CCTK_REAL exterior_min[cctk_dim];
    CCTK_REAL exterior_max[cctk_dim];
    CCTK_REAL spacing     [cctk_dim];
    int const ierr = GetDomainSpecification (cctk_dim,
                                             physical_min,
                                             physical_max,
                                             interior_min,
                                             interior_max,
                                             exterior_min,
                                             exterior_max,
                                             spacing);
    assert (!ierr);
    
    CCTK_REAL domain_volume = 1.0;
    for (int d=0; d<cctk_dim; ++d) {
      domain_volume *= (physical_max[d] - physical_min[d]) / spacing[d];
    }
    
    
    
    int const there_is_a_problem =
      fabs(sum_weight - domain_volume) > 1.0e-12 * (sum_weight + domain_volume);
    
    
    
    if (verbose || there_is_a_problem) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Simulation domain volume: %.15g", (double)domain_volume);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Reduction weight sum:     %.15g", (double)sum_weight);
    }
    
    if (there_is_a_problem) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Simulation domain volume and reduction weight sum differ");
    }
    
  }
  
  
  
  CCTK_DisableGroupStorage (cctkGH, "CarpetReduce::iweight");
  CCTK_DisableGroupStorage (cctkGH, "CarpetReduce::one");
}
