#ifndef CARPET_HH
#define CARPET_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "gh.hh"

#include "carpet_public.hh"



namespace Carpet {
  
  using namespace std;
  
  // Scheduled functions
  extern "C" {
    int CarpetStartup (void);
    int CarpetMultiModelStartup (void);
    void CarpetParamCheck (CCTK_ARGUMENTS);
    void CarpetRefineTimeStep (CCTK_ARGUMENTS);
  }
  
  // Registered functions
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
  
  int Initialise (tFleshConfig* config);
  int Evolve (tFleshConfig* config);
  int Shutdown (tFleshConfig* config);
  int OutputGH (const cGH* cgh);
  
  int CallFunction (void* function, cFunctionData* attribute, void* data);
  
  // Other functions
  bool Regrid (cGH const * cctkGH, bool force_recompose);
  
  void CycleTimeLevels (const cGH* cgh);
  void FlipTimeLevels (const cGH* cgh);
  void FillTimeLevels (const cGH* cgh);
  void SyncGroups (const cGH* cgh, const vector<int>& groups);
  int  SyncProlongateGroups (const cGH* cgh, const vector<int>& groups);
 
  // Sanity checks
  enum checktimes { currenttime,
		    currenttimebutnotifonly,
                    previoustime,
		    allbutlasttime,
		    allbutcurrenttime,
		    alltimes };
  
  int min_timelevel (checktimes where, int num_tl);
  int max_timelevel (checktimes where, int num_tl);
  
  void Poison (const cGH* cgh, checktimes where, int what = 0);
  void PoisonGroup (const cGH* cgh, int group, checktimes where);
  void PoisonCheck (const cGH* cgh, checktimes where);
  
  void CalculateChecksums (const cGH* cgh, checktimes where);
  void CheckChecksums (const cGH* cgh, checktimes where);
  
  // Schedule
  void CallBeforeRoutines (cGH const * cctkGH,
                           void * const function,
                           cFunctionData * const attribute,
                           void * const data);
  void CallAfterRoutines (cGH const * cctkGH,
                          void * const function,
                          cFunctionData * const attribute,
                          void * const data);
  
  // Requirements
  namespace Requirements {
    void CheckRequirements (cGH const * cctkGH);
  }
  
  // Debugging output
  void Output (const char* fmt, ...);
  void Waypoint (const char* fmt, ...);
  void Checkpoint (const char* fmt, ...);
  
  // Error output
  void UnsupportedVarType (int vindex);
  
} // namespace Carpet

#endif // !defined(CARPET_HH)
