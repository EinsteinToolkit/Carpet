// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.14 2001/12/14 16:39:09 schnetter Exp $

#include "carpet_public.hh"

namespace Carpet {
  
  void Recompose (const cGH* cgh);
  void CycleTimeLevels (const cGH* cgh);
  void Restrict (const cGH* cgh);
  
  enum checktimes { currenttime,
		    currenttimebutnotifonly,
		    allbutlasttime,
		    allbutcurrenttime,
		    alltimes };
  
  int mintl (checktimes where, int num_tl);
  int maxtl (checktimes where, int num_tl);
  
  void Poison (cGH* cgh, checktimes where);
  void PoisonGroup (cGH* cgh, int group, checktimes where);
  void PoisonCheck (cGH* cgh, checktimes where);
  
  void CalculateChecksums (cGH* cgh, checktimes where);
  void CheckChecksums (cGH* cgh, checktimes where);
  
  // Debugging output
  void Waypoint (const char* fmt, ...);
  void Checkpoint (const char* fmt, ...);
  
  // Error output
  void UnsupportedVarType (int vindex);
  
} // namespace Carpet
