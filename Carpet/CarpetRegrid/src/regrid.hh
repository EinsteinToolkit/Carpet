// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.hh,v 1.8 2003/04/30 12:37:56 schnetter Exp $

#ifndef CARPETREGRID_HH
#define CARPETREGRID_HH

#include <list>

#include "cctk.h"

#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

#include "regrid.h"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  typedef vect<int,dim> ivect;
  typedef bbox<int,dim> ibbox;
  
  typedef vect<vect<bool,2>,dim> bvect;
  
  typedef vect<CCTK_REAL,dim> rvect;
  
  
  
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
				 const vector<rvect> lower,
				 const vector<rvect> upper,
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

#endif // !defined(CARPETREGRID_HH)
