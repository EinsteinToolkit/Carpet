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

    // print warnings if the user set deprecated parameters in her/his parfile
    if (CCTK_ParameterQueryTimesSet ("minimise_outstanding_communications",
                                     "CarpetLib") > 0) {
      CCTK_WARN (CCTK_WARN_COMPLAIN,
        "You set the parameter 'CarpetLib::minimise_outstanding_communications'"
        " in your parfile. This parameter is deprecated and should not be used"
        " anymore. Use 'CarpetLib::use_collective_communication_buffers'"
        " instead.");
    }
    if (CCTK_ParameterQueryTimesSet ("combine_recv_send",
                                     "CarpetLib") > 0) {
      CCTK_WARN (CCTK_WARN_COMPLAIN,
        "You set the parameter 'CarpetLib::combine_recv_send'"
        " in your parfile. This parameter is deprecated and should not be used"
        " anymore.");
    }
    if (CCTK_ParameterQueryTimesSet ("use_lightweight_buffers",
                                     "CarpetLib") > 0) {
      CCTK_WARN (CCTK_WARN_COMPLAIN,
        "You set the parameter 'CarpetLib::use_lightweight_buffers'"
        " in your parfile. This parameter is deprecated and should not be used"
        " anymore.");
    }
    
    return 0;
  }
  
} // namespace Carpet
