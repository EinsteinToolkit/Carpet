#include <cassert>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "dist.hh"



namespace Carpet {
  
  using namespace std;
  
  int CarpetMultiModelStartup()
  {
    DECLARE_CCTK_PARAMETERS;
    
    comm_universe = MPI_COMM_WORLD;
    SplitUniverse (comm_universe, model, comm_world, true);
    dist::pseudoinit (comm_world);
    
    return 0;
  }
  
  int CarpetStartup()
  {
    CCTK_RegisterBanner ("AMR driver provided by Carpet");
    
    GHExtension = CCTK_RegisterGHExtension("Carpet");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    CCTK_OverloadInitialise (Initialise);
    CCTK_OverloadEvolve (Evolve);
    CCTK_OverloadShutdown (Shutdown);
    
    CCTK_OverloadOutputGH (OutputGH);
    
    CCTK_OverloadSyncGroupsByDirI (SyncGroupsByDirI);
    CCTK_OverloadEnableGroupStorage (EnableGroupStorage);
    CCTK_OverloadDisableGroupStorage (DisableGroupStorage); 
    CCTK_OverloadGroupStorageIncrease (GroupStorageIncrease);
    CCTK_OverloadGroupStorageDecrease (GroupStorageDecrease); 
    CCTK_OverloadEnableGroupComm (EnableGroupComm);
    CCTK_OverloadDisableGroupComm (DisableGroupComm);
    CCTK_OverloadBarrier (Barrier);
    CCTK_OverloadExit (Exit);
    CCTK_OverloadAbort (Abort);
    CCTK_OverloadMyProc (MyProc);
    CCTK_OverloadnProcs (nProcs);
    CCTK_OverloadArrayGroupSizeB (ArrayGroupSizeB);
    CCTK_OverloadQueryGroupStorageB (QueryGroupStorageB);
    CCTK_OverloadGroupDynamicData (GroupDynamicData);
    
    return 0;
  }
  
} // namespace Carpet
