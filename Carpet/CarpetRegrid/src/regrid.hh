// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.hh,v 1.6 2002/05/16 23:25:54 schnetter Exp $

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
    int CarpetRegridParamcheck (CCTK_ARGUMENTS);
  }
  
  
  
  typedef vect<int,dim> ivect;
  typedef bbox<int,dim> ibbox;
  
  typedef vect<vect<bool,2>,dim> bvect;
  
  
  
  int CarpetRegridRegrid (const cGH * const cctkGH,
			  gh<dim>::rexts& bbsss,
			  gh<dim>::rbnds& obss,
			  gh<dim>::rprocs& pss);
  
  
  
  void MakeRegions_BaseLevel    (const cGH* cctkGH,
				 list<ibbox>& bbl, list<bvect>& obl);
  
  void MakeRegions_RefineCentre (const cGH* cctkGH, const int reflevels,
				 list<ibbox>& bbl, list<bvect>& obl);
  
  void MakeRegions_AsSpecified  (const cGH* cctkGH, const int reflevels,
				 const vector<ivect> lower,
				 const vector<ivect> upper,
				 list<ibbox>& bbl, list<bvect>& obl);
  void MakeRegions_AsSpecified  (const cGH* cctkGH, const int reflevels,
				 const vector<vect<CCTK_REAL,dim> > lower,
				 const vector<vect<CCTK_REAL,dim> > upper,
				 list<ibbox>& bbl, list<bvect>& obl);
  
  void MakeRegions_AsSpecified  (const cGH* cctkGH, const int reflevels,
				 const vector<vector<ibbox> > bbss,
				 const vector<vector<bvect> > obss,
				 list<ibbox>& bbl, list<bvect>& obl);
  
  void MakeRegions_Adaptively   (const cGH* cctkGH,
				 const int minwidth,
				 const CCTK_REAL minfraction,
				 const CCTK_REAL maxerror,
				 const gf<CCTK_REAL,dim>& error,
				 list<ibbox>& bbl, list<bvect>& obl);
  
} // namespace CarpetRegrid

#endif // ! defined(REGRID_HH)
