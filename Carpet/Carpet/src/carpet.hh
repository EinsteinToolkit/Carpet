// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.13 2001/11/02 10:59:00 schnetter Exp $

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
  void Checkpoint (const char* fmt, ...);
  
  // Error output
  void UnsupportedVarType (int vindex);
  
} // namespace Carpet
