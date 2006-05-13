#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "bbox.hh"
#include "gh.hh"
#include "vect.hh"

#include "defines.hh"


  
namespace Carpet {
  
  using namespace std;
  
  int SyncGroup (const cGH* cgh, const char* groupname);
  int EnableGroupComm (const cGH* cgh, const char* groupname);
  int DisableGroupComm (const cGH* cgh, const char* groupname);
  int EnableGroupStorage (const cGH* cgh, const char* groupname);
  int DisableGroupStorage (const cGH* cgh, const char* groupname); 
  int GroupStorageIncrease (const cGH* cgh, int n_groups, const int* groups,
                            const int* timelevels, int* status);
  int GroupStorageDecrease (const cGH* cgh, int n_groups, const int* groups,
                            const int* timelevels, int* status);
  int Barrier (const cGH* cgh);
  int Exit (cGH* cgh, int retval);
  int Abort (cGH* cgh, int retval);
  int MyProc (const cGH* cgh);
  int nProcs (const cGH* cgh);
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname);
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname);
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data);
   
   
   
  // Multi-Model
  void
  SplitUniverse (MPI_Comm const world, string const model, MPI_Comm & comm,
                 bool verbose);
  
  // Model id to model name
  vector <string> Models ();
  string Model (int id);
  
  // Model name to model id
  std::map <string, int> ModelMap ();
  int ModelMap (string name);
  
  // Processor to model id
  vector <int> ModelIds ();
  int ModelId (int proc);
  
  // Model id to processors
  vector <vector <int> > ModelProcs ();
  vector <int> ModelProcs (int proc);
  
  extern "C" {
    CCTK_POINTER_TO_CONST
    Carpet_GetMPICommUniverse (CCTK_POINTER_TO_CONST cctkGH);
    CCTK_POINTER_TO_CONST
    Carpet_GetMPICommWorld (CCTK_POINTER_TO_CONST cctkGH);
  }
  
  
  
  // Helpers for storage
  void GroupsStorageCheck (cGH const * const cctkGH);

  // Helpers for recomposing the grid hierarchy
  void Recompose (cGH const * const cctkGH,
                  int const m,
                  gh::mexts  const & bbsss,
                  gh::rbnds  const & obss,
                  gh::rprocs const & pss,
                  bool const do_init);
  
  void PostRecompose ();
  
  void CheckRegions (const gh::mexts & bbsss,
                     const gh::rbnds & obss,
                     const gh::rprocs& pss);
  
  void OutputGrids (const cGH* cgh, const int m, const gh& hh);
  
  void OutputGridStructure (const cGH *cgh,
                            const int m,
                            const gh::mexts & bbsss,
                            const gh::rbnds & obss,
                            const gh::rprocs& pss);
  
  
  
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
                           gh::rexts const & bbss,
                           gh::rbnds const & obss,
                           gh::mexts & bbsss);
  
} // namespace Carpet

#endif // !defined(FUNCTIONS_HH)
