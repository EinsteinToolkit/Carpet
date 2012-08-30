#ifndef REQUIREMENTS_HH
#define REQUIREMENTS_HH

#include <cctk.h>
#include <cctk_Schedule.h>

#include <vector>
  
namespace Carpet {
  namespace Requirements {
    
    using namespace std;
    
    namespace valid {
      enum valid_t { nowhere, interior, everywhere };
    }
    
    // Set up basic grid structure
    void Setup(int maps);
    // Change number of allocated time levels
    void ChangeStorage(vector<int> const& groups, vector<int> const& timelevels,
                       int reflevel);
    // Regrid, set new number of refinement levels, mark all levels as
    // invalid
    void Regrid(int reflevels);
    // Recompose, ensures valid data on one level, indicating whether
    // boundaries are valid (e.g. if recomposing was a no-op)
    void Recompose(int reflevel, valid::valid_t where);
    // Free data structures after regridding
    void RegridFree();
    // Cycle time levels
    void Cycle(int reflevel);
    // Before calling a routine: ensure all reads clauses are
    // satisfied
    // TODO: Either combine these "before" and "after" routines, or
    // split the other routines as well
    void BeforeRoutine(cFunctionData const* function_data,
                       int reflevel, int map, int timelevel);
    // After calling a routine: update according to writes clauses
    void AfterRoutine(cFunctionData const* function_data,
                      int reflevel, int map, int timelevel);
    // Synchronise and prolongate
    // TODO: This does not handle variables that are not prolongated
    // TODO: This does not handle buffer zones
    void Sync(cFunctionData const* function_data,
              vector<int> const& groups, int reflevel, int timelevel);
    // Restrict
    void Restrict(vector<int> const& groups, int reflevel);
    
  } // namespace Requirements
} // namespace Carpet

#endif // #ifndef REQUIREMENTS_HH
