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
    void CarpetParamCheck (CCTK_ARGUMENTS);
    void CarpetStartup (void);
  }
  
  // Registered functions
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
  
  int Initialise (tFleshConfig* config);
  int Evolve (tFleshConfig* config);
  int Shutdown (tFleshConfig* config);
  int CallFunction (void* function, cFunctionData* attribute, void* data);
  
  // Other functions
  bool Regrid (const cGH* cgh, const bool force_recompose, const bool do_init);
  void CycleTimeLevels (const cGH* cgh);
  void FlipTimeLevels (const cGH* cgh);
  void Restrict (const cGH* cgh);
  
  // Sanity checks
  enum checktimes { currenttime,
		    currenttimebutnotifonly,
                    previoustime,
		    allbutlasttime,
		    allbutcurrenttime,
		    alltimes };
  
  int mintl (checktimes where, int num_tl);
  int maxtl (checktimes where, int num_tl);
  
  void Poison (const cGH* cgh, checktimes where);
  void PoisonGroup (const cGH* cgh, int group, checktimes where);
  void PoisonCheck (const cGH* cgh, checktimes where);
  
  void CalculateChecksums (const cGH* cgh, checktimes where);
  void CheckChecksums (const cGH* cgh, checktimes where);
  
  // Debugging output
  void Output (const char* fmt, ...);
  void Waypoint (const char* fmt, ...);
  void Checkpoint (const char* fmt, ...);
  
  // Error output
  void UnsupportedVarType (int vindex);
  
} // namespace Carpet

#endif // !defined(CARPET_HH)
