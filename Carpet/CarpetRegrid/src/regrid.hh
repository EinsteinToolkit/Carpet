// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.hh,v 1.4 2002/03/11 13:17:16 schnetter Exp $

#ifndef REGRID_HH
#define REGRID_HH

#include <list>

#include "cctk.h"

#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  // scheduled functions
  extern "C" {
    int CarpetRegridStartup ();
  }
  
  
  
  int CarpetRegridRegrid (const cGH * const cctkGH,
			  gh<dim>::rexts& bbsss,
			  gh<dim>::rbnds& obss,
			  gh<dim>::rprocs& pss);
  
  
  
  void MakeRegions_BaseLevel    (const cGH* cctkGH,
				 list<bbox<int,dim> >& bbl);
  
  void MakeRegions_RefineCentre (const cGH* cctkGH, const int reflevels,
				 list<bbox<int,dim> >& bbl);
  
  void MakeRegions_AsSpecified  (const cGH* cctkGH, const int reflevels,
				 const vector<vect<int,dim> > lower,
				 const vector<vect<int,dim> > upper,
				 list<bbox<int,dim> >& bbl);
  void MakeRegions_AsSpecified  (const cGH* cctkGH, const int reflevels,
				 const vector<vect<CCTK_REAL,dim> > lower,
				 const vector<vect<CCTK_REAL,dim> > upper,
				 list<bbox<int,dim> >& bbl);
  
  void MakeRegions_Adaptively   (const cGH* cctkGH,
				 const int minwidth,
				 const CCTK_REAL minfraction,
				 const CCTK_REAL maxerror,
				 const gf<CCTK_REAL,dim>& error,
				 list<bbox<int,dim> >& bbl);
  
} // namespace CarpetRegrid

#endif // ! defined(REGRID_HH)
