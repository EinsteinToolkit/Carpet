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
#include "cctk_Version.h"

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

/* some macros for HDF5 group names */
#define PARAMETERS_GLOBAL_ATTRIBUTES_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"



namespace CarpetIOHDF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  int Checkpoint (const cGH* const cctkGH, int called_from);
  int DumpParametersGHExtentions (const cGH *cctkGH, int all, hid_t writer);  
  
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


  int CarpetIOHDF5_RecoverParameters(void){
    // Register with the Cactus Recovery Interface
    return (IOUtil_RecoverParameters (CarpetIOHDF5_Recover, ".h5", "HDF5"));
  }

  int CarpetIOHDF5_Recover (cGH* cctkGH, const char *basefilename, int called_from) {
    
    DECLARE_CCTK_PARAMETERS;

    int result,myproc;
    CarpetIOHDF5GH *myGH;
    char filename[1024];
    
    myGH = NULL;
    result = 0;

    myproc = CCTK_MyProc (cctkGH);
    
  
    if (called_from == CP_RECOVER_PARAMETERS) {
	CCTK_WARN (-1,"Sorry, this feature is not implemented yet.");
      }
    else {
	/* This is the case for CP_RECOVER_DATA.
	   CCTK_RECOVER_PARAMETERS must have been called before
	   and set up the file info structure. */
	if (myproc == 0) {
	  CCTK_WARN (-1,"Sorry, this feature is not implemented yet.");
	}
    }
  
  
    if (called_from == CP_RECOVER_DATA) {
      CCTK_VInfo (CCTK_THORNSTRING,
		  "Restarting simulation at iteration %d (physical time %g)",
		  cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
    }
  
    CCTK_WARN (-1,"STOPSTOPSTOP2");
  
    return (result);
  } // CarpetIOHDF5_Recover
  

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

    // Dump all parameters and GHExtentions
    retval += DumpParametersGHExtentions (cctkGH, 1, writer);
    assert(!retval);

    if (myproc==0) {
      // Close the file
      herr = H5Fclose(writer);
      assert(!herr);
    }
    
    if (retval == 0) {
      if (myproc==0) {
	if (rename (cp_tempname, cp_filename))
	  {
	    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
			"Could not rename temporary checkpoint file '%s' to '%s'",
			cp_tempname, cp_filename);
	    retval = -1;
	  }
	else {
	  if (myGH->cp_filename_list[myGH->cp_filename_index]) {
	    if (checkpoint_keep > 0) {
	      remove (myGH->cp_filename_list[myGH->cp_filename_index]);
	    }
	    free (myGH->cp_filename_list[myGH->cp_filename_index]);
	  }
	  myGH->cp_filename_list[myGH->cp_filename_index] = strdup (cp_filename);
	  myGH->cp_filename_index = (myGH->cp_filename_index+1) % abs (checkpoint_keep);
	} // else
      } // myproc == 0
    } // retval == 0

    
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,"This feature is not working yet. Sorry!");

    return retval;

  } // Checkpoint

  
  int DumpParametersGHExtentions (const cGH *cctkGH, int all, hid_t writer)
  {
    // large parts of this routine were taken from
    // Thomas Radke's IOHDF5Util. Thanks Thomas!

    DECLARE_CCTK_PARAMETERS;

    char *parameters;
    hid_t group, dataspace, dataset;
    hsize_t size;
    herr_t herr;
    
    CCTK_INT4 itmp;
    CCTK_REAL dtmp;
    const char *version;
    ioGH *ioUtilGH;
  
    if (verbose) {
      CCTK_INFO ("Dumping Parameters and GH Extentions...");
    }
  
    /* get the parameter string buffer */
    parameters = IOUtil_GetAllParameters (cctkGH, all);
  
    if (parameters)
    {
      size = strlen (parameters) + 1;
      group = H5Gcreate (writer, PARAMETERS_GLOBAL_ATTRIBUTES_GROUP, 0);
      assert(group>=0);
      dataspace = H5Screate_simple (1, &size, NULL);
      assert(dataspace>=0);
      dataset = H5Dcreate (group, ALL_PARAMETERS, H5T_NATIVE_UCHAR,
                                       dataspace, H5P_DEFAULT);
      assert(dataset>=0);
      herr = H5Dwrite (dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, parameters);
      assert(!herr);

      // now dump the GH Extentions

      /* get the handle for IOUtil extensions */
      ioUtilGH = (ioGH *) CCTK_GHExtension (cctkGH, "IO");
      
      itmp = CCTK_MainLoopIndex ();
      WriteAttribute(dataset,"main loop index",itmp);
      
      itmp = cctkGH->cctk_iteration;
      WriteAttribute(dataset,"GH$iteration",itmp);
      
      itmp = ioUtilGH->ioproc_every;
      WriteAttribute(dataset,"GH$ioproc_every",itmp);
      
      itmp = CCTK_nProcs (cctkGH);
      WriteAttribute(dataset,"GH$nprocs",itmp);
      
      dtmp = cctkGH->cctk_time;
      WriteAttribute(dataset,"GH$time", dtmp);

      version = CCTK_FullVersion();
      WriteAttribute(dataset,"Cactus version", version);

      /* finally, we need all the times on the individual refinement levels */
      const int numberoftimes=leveltimes[0].size();
      WriteAttribute(dataset,"numberoftimes",numberoftimes);
      for(int i=0;i < numberoftimes;i++) {
	char buffer[100];
	sprintf(buffer,"leveltime%d",i);
	WriteAttribute(dataset,buffer,leveltimes[0][i]);
      }


      herr = H5Dclose (dataset);
      assert(!herr);
      herr = H5Sclose (dataspace);
      assert(!herr);
      herr = H5Gclose (group);
      assert(!herr);
  
      free (parameters);
    }

    return 0;
  } // DumpParametersGHExtentions
  

    
  
} // namespace CarpetIOHDF5
