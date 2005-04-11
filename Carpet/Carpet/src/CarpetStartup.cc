#include <cassert>
#include <cstdlib>

#include "cctk.h"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  void CarpetStartup()
  {
    CCTK_RegisterBanner ("AMR driver provided by Carpet");
    
    GHExtension = CCTK_RegisterGHExtension("Carpet");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    CCTK_OverloadInitialise (Initialise);
    CCTK_OverloadEvolve (Evolve);
    CCTK_OverloadShutdown (Shutdown);
    
    CCTK_OverloadSyncGroup (SyncGroup);
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

    if (CCTK_ParameterQueryTimesSet ("minimise_outstanding_communications",
                                     "CarpetLib") > 0) {
      CCTK_WARN (CCTK_WARN_COMPLAIN,
        "You set the parameter 'CarpetLib::minimise_outstanding_communications'"
        " in your parfile. This parameter is deprecated and should not be used"
        " anymore. Use 'CarpetLib::use_collective_communication_buffers'"
        " instead.");
    }
  }
  
} // namespace Carpet
