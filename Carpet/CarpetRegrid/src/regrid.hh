// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.hh,v 1.1 2001/12/14 16:34:39 schnetter Exp $

#ifndef REGRID_HH
#define REGRID_HH

#include <list>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "Carpet/CarpetLib/src/gf.hh"

#include "carpet.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  // scheduled functions
  extern "C" {
    int CarpetRegridRegrid (CCTK_ARGUMENTS);
  }
  
  
  
  void MakeRegions_BaseLevel    (const cGH* cctkGH,
				 list<bbox<int,dim> >& bbl);
  
  void MakeRegions_RefineCentre (const cGH* cctkGH, const int reflevels,
				 list<bbox<int,dim> >& bbl);
  
  void MakeRegions_AsSpecified  (const cGH* cctkGH, const int reflevels,
				 const vector<vect<int,dim> > lower,
				 const vector<vect<int,dim> > upper,
				 list<bbox<int,dim> >& bbl);
  
  void MakeRegions_Adaptively   (const cGH* cctkGH,
				 const int minwidth, const double minfraction,
				 const CCTK_REAL maxerror,
				 const gf<CCTK_REAL,dim>& error,
				 list<bbox<int,dim> >& bbl);
  
} // namespace CarpetRegrid

#endif // ! defined(REGRID_HH)
