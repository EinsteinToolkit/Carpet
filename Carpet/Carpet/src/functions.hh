// $Header:$

#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <vector>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "bbox.hh"
#include "gh.hh"
#include "vect.hh"

#include "defines.hh"


  
namespace Carpet {
  
  using namespace std;
  
  int SyncGroup (const cGH* cgh, const char* groupname);
  int EnableGroupStorage (const cGH* cgh, const char* groupname);
  int DisableGroupStorage (const cGH* cgh, const char* groupname); 
  int EnableGroupComm (const cGH* cgh, const char* groupname);
  int DisableGroupComm (const cGH* cgh, const char* groupname);
  int Barrier (const cGH* cgh);
  int Exit (cGH* cgh, int retval);
  int Abort (cGH* cgh, int retval);
  int MyProc (const cGH* cgh);
  int nProcs (const cGH* cgh);
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname);
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname);
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data);
   
   
   
  // Helpers for recomposing the grid hierarchy
  void CheckRegions (const gh<dim>::rexts & bbsss,
                     const gh<dim>::rbnds & obss,
                     const gh<dim>::rprocs& pss);
  
  void OutputGrids (const cGH* cgh, const int m, const gh<dim>& hh);
  
  void OutputGridStructure (const cGH *cgh,
                            const int m,
                            const gh<dim>::rexts & bbsss,
                            const gh<dim>::rbnds & obss,
                            const gh<dim>::rprocs& pss);
  
  
  
  // Functions for recomposing the grid hierarchy
  void SplitRegions (const cGH* cgh, vector<ibbox>& bbs,
		     vector<bbvect>& obs, vector<int>& ps);
  void SplitRegions_AlongZ (const cGH* cgh, vector<ibbox>& bbs,
			    vector<bbvect>& obs, vector<int>& ps);
  void SplitRegions_AlongDir (const cGH* cgh, vector<ibbox>& bbs,
                              vector<bbvect>& obs, vector<int>& ps,
                              const int dir);
  void SplitRegions_Automatic (const cGH* cgh, vector<ibbox>& bbs,
                               vector<bbvect>& obs, vector<int>& ps);
  
  void MakeMultigridBoxes (const cGH* cgh,
                           vector<ibbox> const & bbs,
                           vector<bbvect> const & obs,
                           vector<vector<ibbox> > & bbss);
  
} // namespace Carpet

#endif // !defined(FUNCTIONS_HH)
