// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.26 2003/11/05 16:18:37 schnetter Exp $

#ifndef CARPET_HH
#define CARPET_HH

#include "gh.hh"

#include "carpet_public.hh"



namespace Carpet {
  
  void Regrid (const cGH* cgh,
               const int initialise_from, const bool do_prolongate);
  void CycleTimeLevels (const cGH* cgh);
  void FlipTimeLevels (const cGH* cgh);
  void Restrict (const cGH* cgh);
  
  void Recompose (const cGH* cgh,
		  const gh<dim>::rexts& bbsss,
		  const gh<dim>::rbnds& obss,
		  const gh<dim>::rprocs& pss,
                  const int initialise_from,
                  const bool do_prolongate);
  
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
  void Waypoint (const char* fmt, ...);
  void Checkpoint (const char* fmt, ...);
  
  // Error output
  void UnsupportedVarType (int vindex);
  
} // namespace Carpet

#endif // !defined(CARPET_HH)
