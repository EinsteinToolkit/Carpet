#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "bbox.hh"
#include "dh.hh"
#include "gh.hh"
#include "vect.hh"


  
namespace Carpet {
  
  using namespace std;
  using namespace CarpetLib;
  
  int SyncGroupsByDirI (const cGH* cctkGH, int num_groups,
                        const int* groups, const int* directions);
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
  
  void Restrict (const cGH* cgh);
  
  
  
  // Strings
  vector <string>
  AllGatherString (MPI_Comm const world,
                   string const & data);
  
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
    CCTK_INT
    Carpet_GetCoordRange (CCTK_POINTER_TO_CONST const         cctkGH_,
                          CCTK_INT              const         m,
                          CCTK_INT              const         ml,
                          CCTK_INT              const         size,
                          CCTK_INT                    * const gsh,
                          CCTK_REAL                   * const lower,
                          CCTK_REAL                   * const upper,
                          CCTK_REAL                   * const delta);
  }
  
  
  
  void SetSystemLimits ();
  
  
  
  // Helpers for storage
  void GroupsStorageCheck (cGH const * const cctkGH);

  // Helpers for recomposing the grid hierarchy
  void
  RegridMap (cGH const * cctkGH,
             int m,
             gh::rregs const & supeerregss,
             gh::mregs const & regsss);
  void
  PostRegrid (cGH const * cctkGH);
  bool
  Recompose (cGH const * cctkGH,
             int rl,
             bool do_init);
  
  void
  CheckRegions (gh::mregs const & regsss);
  
  void
  OutputGrids (cGH const * cctkGH,
               int const m,
               gh const & hh,
               dh const & dd);
  
  void
  OutputGridStructure (cGH const * cctkGH,
                       int const m,
                       gh::mregs const & regsss);
  
  void
  OutputGridCoordinates (cGH const * cctkGH,
                         int const m,
                         gh::mregs const & regsss);
  
  void
  OutputGridStatistics (cGH const * cctkGH);
  
  
  
  // Functions for recomposing the grid hierarchy
  void
  SplitRegions (cGH const * cctkGH,
                vector<region_t> & superregs,
                vector<region_t> & regs);
  void
  SplitRegions_AlongZ (cGH const * cctkGH,
                       vector<region_t> & superregs,
                       vector<region_t> & regs);
  void
  SplitRegions_AlongDir (cGH const * cctkGH,
                         vector<region_t> & superregs,
                         vector<region_t> & regs,
                         int dir);
  void
  SplitRegions_Automatic (cGH const * cctkGH,
                          vector<region_t> & superregs,
                          vector<region_t> & regs);
  
  void
  SplitRegionsMaps (cGH const * cctkGH,
                    vector<vector<region_t> > & superregss,
                    vector<vector<region_t> > & regss);
  void
  SplitRegionsMaps_Automatic (cGH const * cctkGH,
                              vector<vector<region_t> > & superregss,
                              vector<vector<region_t> > & regss);
  
  void
  MakeMultigridBoxes (cGH const * cctkGH,
                      int m,
                      gh::rregs const & regss,
                      gh::mregs       & regsss);
  
  void
  MakeMultigridBoxesMaps (cGH const * cctkGH,
                          vector<gh::rregs> const & regsss,
                          vector<gh::mregs>       & regssss);
  
  
  
  // Timing statistics functions
  void InitTimingStats (cGH const * cctkGH);
  void BeginTiming (cGH const * cctkGH);
  void StepTiming (cGH const * cctkGH);
  void BeginTimingIO (cGH const * cctkGH);
  void EndTimingIO (cGH const * cctkGH,
                    CCTK_REAL files, CCTK_REAL bytes, bool is_binary);
  void BeginTimingCommunication (cGH const * cctkGH);
  void EndTimingCommunication (cGH const * cctkGH,
                               CCTK_REAL messages, CCTK_REAL bytes);
  void UpdateTimingStats (cGH const * cctkGH);
  void PrintTimingStats (cGH const * cctkGH);
  
} // namespace Carpet

#endif // !defined(FUNCTIONS_HH)
