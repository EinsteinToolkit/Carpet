// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.12 2001/07/09 09:00:13 schnetter Exp $

#include "carpet_public.hh"

namespace Carpet {
  
  void Recompose (cGH* cgh);
  void CycleTimeLevels (cGH* cgh);
  void Restrict (cGH* cgh);
  
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
  void Checkpoint (const char* fmt, ...);
  
  // Error output
  void UnsupportedVarType (int vindex);
  
} // namespace Carpet
