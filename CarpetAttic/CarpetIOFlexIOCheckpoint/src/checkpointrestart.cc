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



#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Version.h"

#include "AMRwriter.hh"
#include "AmrGridReader.hh"
#ifdef HDF4
#  include "HDFIO.hh"
#endif
#ifdef HDF5
#  include "H5IO.hh"
#endif
#include "IEEEIO.hh"
#include "IO.hh"

// Hack to stop FlexIO type clash
#undef BYTE
#undef CHAR

//#include "CactusBase/IOUtil/src/ioGH.h"
//#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
//#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

//#include "StoreNamedData.h"
#include "carpet.hh"

#include "ioflexio.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/checkpointrestart.cc,v 1.7 2003/09/23 12:34:43 cvs_anon Exp $";
  CCTK_FILEVERSION(Carpet_CarpetIOFlexIO_checkpointrestart_cc);
}




namespace CarpetCheckpointRestart {

  using namespace std;
  using namespace Carpet;
  using namespace CarpetIOFlexIO;
  using namespace CarpetIOFlexIOUtil;

  static int Checkpoint (const cGH* const cgh, int called_from);




  void CarpetIOFlexIO_EvolutionCheckpoint( const cGH* const cgh){
    
    DECLARE_CCTK_PARAMETERS

      CCTK_VInfo (CCTK_THORNSTRING, "CHECKPOINT? cgh->cctk_iteration %d, cgh->cctk_iteration mod checkpoint_every = %d, checkpoint = %d, checkpoint_next = %d", cgh->cctk_iteration, cgh->cctk_iteration % checkpoint_every,checkpoint,checkpoint_next);
      if (checkpoint &&
      ((checkpoint_every > 0 && cgh->cctk_iteration % checkpoint_every == 0) ||
       checkpoint_next))
      {
	if (! CCTK_Equals ("verbose", "none"))
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
    

  }

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/


  static int DumpParams (const cGH* const cgh, int all, IObase* writer){

    char *parameters;
    
    parameters = IOUtil_GetAllParameters(cgh,all);

    if(parameters)
    {
      writer->writeAttribute("global_parameters",IObase::Char,
			     strlen(parameters)+1,parameters);
      free(parameters);
    }
    return 0;
  }

  static int DumpGHExtensions (const cGH* const cgh, IObase* writer){

    CCTK_INT4 itmp;
    CCTK_REAL dtmp;
    const char *version;
    ioGH *ioUtilGH;


    /* get the handle for IOUtil extensions */
    ioUtilGH = (ioGH *) CCTK_GHExtension (cgh, "IO");

    itmp = CCTK_MainLoopIndex ();
    writer->writeAttribute("main loop index",FLEXIO_INT4,1,&itmp);

    itmp = cgh->cctk_iteration;
    writer->writeAttribute("GH$iteration",FLEXIO_INT4, 1, &itmp);

    itmp = ioUtilGH->ioproc_every;
    writer->writeAttribute("GH$ioproc_every",FLEXIO_INT4,1,&itmp);

    itmp = CCTK_nProcs (cgh);
    writer->writeAttribute("GH$nprocs",FLEXIO_INT4, 1, &itmp);

    dtmp = cgh->cctk_time;
    writer->writeAttribute("GH$time", FLEXIO_REAL, 1, &dtmp);

    version = CCTK_FullVersion ();
    writer->writeAttribute("Cactus version", FLEXIO_CHAR,
                                  strlen (version) + 1, version);

    return 0;
  }


  static int Checkpoint (const cGH* const cgh, int called_from)
  {
    char cp_filename[1024], cp_tempname[1024];
    int myproc, first_vindex, gindex, retval;
    char *fullname;
    const char *timer_descriptions[3] = {"Time to dump parameters: ",
					 "Time to dump datasets:   ",
					 "Total time to checkpoint:"};
    const ioGH *ioUtilGH;

    //    const int varindex = CCTK_VarIndex("ADMBASE:gxx");
    int varindex = 0;
    int group = 0;

    cGroup  gdata;
    IObase* writer = 0;
    AMRwriter* amrwriter = 0;
    ioRequest *request;

    DECLARE_CCTK_PARAMETERS
      //   CCTK_VInfo (CCTK_THORNSTRING, "boguscheck reflevel,component,mglevel %d,%d,%d",reflevel,component,mglevel);

  
    myproc = CCTK_MyProc (cgh);
    ioUtilGH = (const ioGH *) CCTK_GHExtension (cgh, "IO");
    IOUtil_PrepareFilename (cgh, NULL, cp_filename, called_from,
            myproc / ioUtilGH->ioproc_every, ioUtilGH->unchunked);

    // Invent a file name
    const char* extension = 0;
    if (CCTK_Equals(out3D_format, "IEEE")) {
      extension = ".raw";
#ifdef HDF4
    } else if (CCTK_Equals(out3D_format, "HDF4")) {
      extension = ".hdf";
#endif
#ifdef HDF5
    } else if (CCTK_Equals(out3D_format, "HDF5")) {
      extension = ".h5";
#endif
    } else {
      assert (0);
    }
  
    sprintf(cp_tempname,"%s.tmp.%s",cp_filename,extension);
    sprintf(cp_filename,"%s.%s",cp_filename,extension);
  

    if (CCTK_MyProc(cgh)==0)
      {
	if (CCTK_Equals ("verbose", "full"))
	  {
	    CCTK_VInfo (CCTK_THORNSTRING, "Creating temporary checkpoint file '%s'", cp_tempname);
	  }

	writer = new IEEEIO(cp_tempname, IObase::Create);

	if (! (writer->isValid()) )
	  {
	    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
			"Can't open checkpoint file '%s'. Checkpointing is skipped",
			cp_tempname);
	    return (-1);
	  }
	amrwriter = new AMRwriter(*writer);

	// dump parameters 
	DumpParams (cgh, 1, writer);
	// dump GH extentions
	DumpGHExtensions(cgh,writer);

      }

    
    // now dump the grid varibles for all reflevels and components, sorted by groups
    //   CCTK_VInfo (CCTK_THORNSTRING, "maxreflevelfact,reflevelfact,it %d,%d,%d",maxreflevelfact,reflevelfact,cgh->cctk_iteration);
    BEGIN_REFLEVEL_LOOP(cgh) {

     const int do_every = maxreflevelfact/reflevelfact;
      if (cgh->cctk_iteration % do_every == 0) {
	
	BEGIN_MGLEVEL_LOOP(cgh) {
	  const int do_every = mglevelfact * maxreflevelfact/reflevelfact;
          CCTK_VInfo (CCTK_THORNSTRING, "do_every %d",do_every);	 	  

	  if (cgh->cctk_iteration % do_every == 0) {
	      
	    if (CCTK_Equals ("verbose", "full"))
	    {
	      CCTK_INFO ("Dumping Grid Variables ...");
	    }
	    for (group = CCTK_NumGroups () - 1; group >= 0; group--)
	    {
	      /* only dump groups which have storage assigned */
	      if (CCTK_QueryGroupStorageI (cgh, group) <= 0)
	      {
		continue;
	      }
	       
	      /* get the number of allocated timelevels */
	      CCTK_GroupData (gindex, &gdata);
	      gdata.numtimelevels = 0;
	      gdata.numtimelevels = CCTK_GroupStorageIncrease (cgh, 1, &group,
							       &gdata.numtimelevels,NULL);
	      
	      /* dump all timelevels except the oldest (for multi-level groups) */
	      CCTK_GroupData (group, &gdata);
	      if (gdata.numtimelevels > 1)
	      {
		gdata.numtimelevels--;
	      }
	      
	      int first_vindex = CCTK_FirstVarIndexI (group);

	      const int grouptype = CCTK_GroupTypeI(group);

	      /* get the default I/O request for this group */
	      request = IOUtil_DefaultIORequest (cgh, first_vindex, 1);

	      /* disable checking for old data objects, disable datatype conversion
		 and downsampling */
	      request->check_exist = 0;
	      request->hdatatype = gdata.vartype;
	      for (request->hdim = 0; request->hdim < request->vdim; request->hdim++)
	      {
		request->downsample[request->hdim] = 1;
	      }

	      /* loop over all variables in this group */
	      for (request->vindex = first_vindex;
		   request->vindex < first_vindex + gdata.numvars;
		   request->vindex++)
	      {
		/* loop over all timelevels of this variable */
		for (request->timelevel = 0;
		     request->timelevel < gdata.numtimelevels;
		     request->timelevel++)
		{
		  
		  if (verbose)
		  {
		    fullname = CCTK_FullName (request->vindex);
		    CCTK_VInfo (CCTK_THORNSTRING, "  %s (timelevel %d)",
				fullname, request->timelevel);
		    free (fullname);
		  }
		  // write the var
	    
		  if (grouptype == CCTK_SCALAR)
		  {
		    retval += WriteGS(cgh,writer,request);
		  }
		  //		  else if (grouptype == CCTK_ARRAY || grouptype == CCTK_GF)
		  else if (grouptype == CCTK_GF)
		  {
		    retval += WriteGF(cgh,writer,amrwriter,request);
		  }
		  else
		  {
		    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
				"Invalid group type %d for variable '%s'", grouptype, fullname);
		    retval = -1;
		  }

		}
	      } /* end of loop over all variables */
	
	    } /* end of loop over all groups */
	  }
	} END_MGLEVEL_LOOP;
      }
    } END_REFLEVEL_LOOP;




    // Close the file
    if (CCTK_MyProc(cgh)==0) {
      delete amrwriter;
      amrwriter = 0;
      delete writer;
      writer = 0;
    }


    return 0;
      }


} // namespace CarpetCheckpointRestart














