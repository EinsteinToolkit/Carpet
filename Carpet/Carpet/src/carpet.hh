// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.16 2002/03/11 13:17:10 schnetter Exp $

#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet_public.hh"

namespace Carpet {
  
  void Regrid (const cGH* cgh);
  void CycleTimeLevels (const cGH* cgh);
  void Restrict (const cGH* cgh);
  
  void Recompose (const cGH* cgh,
		  const gh<dim>::rexts& bbsss,
		  const gh<dim>::rbnds& obss,
		  const gh<dim>::rprocs& pss);
  
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
