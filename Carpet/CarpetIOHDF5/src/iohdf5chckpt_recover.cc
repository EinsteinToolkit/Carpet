#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"

extern "C" {
  static const char* rcsid = "$ $";
  CCTK_FILEVERSION(Carpet_CarpetIOHDF5_iohdf5chckpt_recover_cc);
}

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "iohdf5.hh"
#include "iohdf5GH.h"


namespace CarpetIOHDF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  int Checkpoint (const cGH* const cctkGH, int called_from);

  //  int RecoverParameters (IObase* reader);
  //int RecoverGHextensions (cGH* cgh, IObase* reader);
  //int RecoverVariables (cGH* cgh, IObase* reader, AmrGridReader* amrreader);

  void CarpetIOHDF5_EvolutionCheckpoint( const cGH* const cgh){
    
    DECLARE_CCTK_PARAMETERS

      if (checkpoint &&
      ((checkpoint_every > 0 && cgh->cctk_iteration % checkpoint_every == 0) ||
       checkpoint_next))
      {
	if (verbose)
	{
	  CCTK_INFO ("---------------------------------------------------------");
	  CCTK_VInfo (CCTK_THORNSTRING, "Dumping periodic checkpoint at "
		      "iteration %d", cgh->cctk_iteration);
	  CCTK_INFO ("---------------------------------------------------------");
	}

	Checkpoint (cgh, CP_EVOLUTION_DATA);

       	if (checkpoint_next)
	{
	  CCTK_ParameterSet ("checkpoint_next", CCTK_THORNSTRING, "no");
	}
      }
  } // CarpetIOHDF5_EvolutionCheckpoint

  

  int Checkpoint (const cGH* const cctkGH, int called_from)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    char cp_filename[1024], cp_tempname[1024];
    char *fullname;
    
    herr_t herr;
    
    int retval = 0;

    const ioGH *ioUtilGH;
    
    cGroup  gdata;
    ioRequest *request;

    CarpetIOHDF5GH *myGH;
    myGH = (CarpetIOHDF5GH *) CCTK_GHExtension (cctkGH, "CarpetIOHDF5");

    /* check if CarpetIOHDF5 was registered as I/O method */
    if (myGH == NULL) {
      CCTK_WARN (-1, "No CarpetIOHDF5 I/O methods registered");
      return (-1);
    }
  
    int myproc = CCTK_MyProc (cctkGH);

    // Invent a filename

    ioUtilGH = (const ioGH *) CCTK_GHExtension (cctkGH, "IO");
    IOUtil_PrepareFilename (cctkGH, NULL, cp_filename, called_from,
            myproc / ioUtilGH->ioproc_every, ioUtilGH->unchunked);

    /* ... and append the extension */
    sprintf (cp_tempname, "%s.tmp.h5", cp_filename);
    sprintf (cp_filename, "%s.h5",     cp_filename);

    hid_t writer = -1;

    if (myproc == 0) {

      if (verbose) {
	CCTK_VInfo (CCTK_THORNSTRING, "Creating temporary checkpoint file '%s'", cp_tempname);
      }
      writer = H5Fcreate (cp_tempname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      assert (writer>=0);
      herr = H5Fclose (writer);
      assert (!herr);
      writer = -1;

      // Now open the file

      writer = H5Fopen (cp_tempname, H5F_ACC_RDWR, H5P_DEFAULT);
      assert (writer>=0);

    } // myproc == 0 

    


    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,"This feature is not working yet. Sorry!");
    

    // Close the file
    H5Fclose(writer);

    return retval;

  } // Checkpoint

  
  
} // namespace CarpetIOHDF5
