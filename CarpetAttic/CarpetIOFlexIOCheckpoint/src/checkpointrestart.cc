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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/checkpointrestart.cc,v 1.18 2004/01/07 13:14:04 cott Exp $";
  CCTK_FILEVERSION(Carpet_CarpetIOFlexIO_checkpointrestart_cc);
}


namespace CarpetCheckpointRestart {

  using namespace std;
  using namespace Carpet;
  using namespace CarpetIOFlexIO;
  using namespace CarpetIOFlexIOUtil;

  int Checkpoint (const cGH* const cgh, int called_from);

  int RecoverParameters (IObase* reader);
  int RecoverGHextensions (cGH* cgh, IObase* reader);
  int RecoverVariables (cGH* cgh, IObase* reader, AmrGridReader* amrreader);

  void CarpetIOFlexIO_EvolutionCheckpoint( const cGH* const cgh){
    
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
    

  }

/*@@
   @routine    CarpetIOFlexIO_RecoverParameters
   @date       Fri Oct 10 2003
   @author     Christian Ott, Thomas Radke
   @desc
   @desc
               Recovers the parameters from an HDF5 checkpoint file.
               This routine is scheduled at CCTK_RECOVER_PARAMETERS.

               Note that it cannot be registered with IOUtil to be scheduled
               from there (as done with the CarpetIOFlexIO_Recover routine) because
               the registration mechanism isn't activated yet
               at CCTK_RECOVER_PARAMETERS.
               Instead we call the generic parameter recovery routine
               from IOUtil here, and just pass the necessary callback function
               and its arguments.

               Note also that this routine doesn't get passed any parameters,
               not even a GH, because this doesn't exist yet at the time it is
               being called.
   @enddesc

   @calls      IOUtil_RecoverParameters

   @returntype int
   @returndesc
               return code of @seeroutine IOUtil_RecoverParameters, ie.
               positive for successful parameter recovery, or<BR>
               0 if recovery wasn't requested, or<BR>
               negative if parameter recovery failed
   @endreturndesc
@@*/


  int CarpetIOFlexIO_RecoverParameters(void){
    return (IOUtil_RecoverParameters (CarpetIOFlexIO_Recover, ".hdf5", "HDF5"));
  }

/*@@
   @routine    CarpetIOFlexIO_Recover
   @date       Fri Oct 10 2003
   @author     Christian Ott (Tom Goodale IOFlexIO version)
   @desc
               Recovers a Carpet GH from an HDF5 file.
               This routine is registered with IOUtil as CarpetIOFlexIO's recovery
               routine.
   @enddesc

   @calls      OpenFile
               RecoverParameters
               RecoverGHextensions
               CarpetIOFlexIOi_RecoverVariables
               IOUtil_PrintTimings

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar
   @var        basefilename
   @vdesc      the basefilename of the file to recover from
               The file suffix is appended by the routine.
   @vtype      const char *
   @vio        in
   @endvar
   @var        called_from
   @vdesc      flag indicating where this routine was called from
               (either CP_RECOVER_DATA or FILEREADER_DATA)
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               >0 = success
               -1 = recovery failed
   @endreturndesc
@@*/


int CarpetIOFlexIO_Recover (cGH* cgh, const char *basefilename, int called_from)
{
  int result,myproc;
  CarpetIOFlexIOGH *myGH;
  char filename[1024];
  
  static IObase* reader = NULL;
  
  DECLARE_CCTK_PARAMETERS

  /* to make the compiler happy */
  myGH = NULL;
  result = 0;

  myproc = CCTK_MyProc (cgh);

  fprintf(stderr,"\n reflevel: %d\n",reflevel);


  if (called_from == CP_RECOVER_PARAMETERS)
  {
    CCTK_VInfo (CCTK_THORNSTRING, "got this far... '%s'", basefilename);  
    if (myproc == 0){
      reader = new H5IO(basefilename,IObase::Read);
      CCTK_VInfo (CCTK_THORNSTRING, "blah '%s'", basefilename);  
      if ( ! reader->isValid() )
	{
	  CCTK_VInfo(CCTK_THORNSTRING,"file is not open");
	  return (-1);
	}
      CCTK_VInfo(CCTK_THORNSTRING,"file is open");
    }
  }
  else
  {
    /* This is the case for CP_RECOVER_DATA.
       CCTK_RECOVER_PARAMETERS must have been called before
       and set up the file info structure. */
    if (myproc == 0){
      if (! reader->isValid() )
	{
	  CCTK_VInfo(CCTK_THORNSTRING,"file is not open2");
	  return (-1);
	}
    }
  }

  /* Recover parameters */

  if (called_from == CP_RECOVER_PARAMETERS)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"called from recover parameters");
    return (RecoverParameters (reader));
  }


  if (called_from == CP_RECOVER_DATA) {
    if(myproc ==0){
      AmrGridReader* amrreader = 0;    
      amrreader = new AmrGridReader(*reader);

      BEGIN_REFLEVEL_LOOP(cgh) {
	BEGIN_MGLEVEL_LOOP(cgh) {
	  
	/* Recover GH extentions */
	  CCTK_INFO ("Recovering GH extensions");
	  result = RecoverGHextensions (cgh, reader);
	  
	if (! result)
	  {
	    /* Recover variables */
	    CCTK_VInfo (CCTK_THORNSTRING, "Recovering data! ");
	    result = RecoverVariables (cgh, reader,amrreader);
	  }
	
	} END_MGLEVEL_LOOP;
      } END_REFLEVEL_LOOP;

      delete reader;
      delete amrreader;
    }
  }


  
  if (called_from == CP_RECOVER_DATA)
    {
	  CCTK_VInfo (CCTK_THORNSTRING,
		      "Restarting simulation at iteration %d (physical time %g)",
		      cgh->cctk_iteration, (double) cgh->cctk_time);
    }
  
  //  CCTK_WARN (-1,"STOPSTOPSTOP2");
  
  return (result);
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

  int RecoverParameters(IObase* reader){

    int myproc, retval;
    int i, asize;
    char *parameters;
    CCTK_INT4 parameterSize;

    IObase::DataType datatype;
    CCTK_REAL bogusdata;

    int dims[3];
    int rank=0;
    int maxdims=3;

    DECLARE_CCTK_PARAMETERS 

    myproc = CCTK_MyProc (NULL);

    if (myproc == 0){
      CCTK_VInfo (CCTK_THORNSTRING, "Recovering parameters from checkpoint ");
      
      /* read the first (bogus) dataset to which the parameters an GHExtensions
         are attached */
      
      reader->readInfo(datatype,rank,dims,maxdims);
      CCTK_VInfo (CCTK_THORNSTRING, "blahbalh: datatype %d FLEXIO_REAL: %d",datatype,FLEXIO_REAL);
      
      if(datatype != FLEXIO_REAL || rank !=1 )
	CCTK_WARN (-1,"Wrong recover file format! First dataset type mismatch!");
   
      reader->read(&bogusdata);

      /* get the parameters attribute. */
      i = reader->readAttributeInfo ("global_parameters", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_CHAR && asize > 0)
      {
	parameterSize = (CCTK_INT4) asize;
	parameters = (char*) malloc (parameterSize + 1);
	reader->readAttribute (i, parameters);
      }
      else
	{
	  CCTK_WARN (1, "Can't read global parameters. "
		     "Is this really a Cactus IEEEIO checkpoint file ?");
	}
      CCTK_VInfo (CCTK_THORNSTRING, "\n%s\n",parameters);
    }

#ifdef CCTK_MPI
  /* Broadcast the parameter buffer size to all processors */
  /* NOTE: We have to use MPI_COMM_WORLD here
     because CARPET_COMM_WORLD is not yet set up at parameter recovery time.
     We also assume that CARPET_MPI_INT4 is a compile-time defined datatype. */
    CACTUS_MPI_ERROR (MPI_Bcast (&parameterSize, 1, CARPET_MPI_INT4, 0,
				 MPI_COMM_WORLD));
#endif

    if (parameterSize > 0)
      {
#ifdef CCTK_MPI
	if (myproc)
	  {
	    parameters = (char*) malloc (parameterSize + 1);
	  }
	
    CACTUS_MPI_ERROR (MPI_Bcast (parameters, parameterSize + 1, CARPET_MPI_CHAR,
				 0, MPI_COMM_WORLD));
#endif

    IOUtil_SetAllParameters (parameters);

    free (parameters);
  }

    /* return positive value for success otherwise negative */
    retval = (parameterSize > 0 ? 1 : -1);

    return (retval);


    //    CCTK_WARN (-1,"STOPSTOPSTOP");


  }


  static int RecoverGHextensions (cGH *GH, IObase* reader)
  {
    int i, type,dim;
    CCTK_REAL realBuffer;
    CCTK_INT4 int4Buffer[2];

    IObase::DataType datatype;
    
    if (CCTK_MyProc (GH) == 0)
      {
		
	/* get the iteration number */
	i = reader->readAttributeInfo ("GH$iteration", datatype, dim);


	if (i >= 0 && datatype == FLEXIO_INT4 && dim == 1)
	  {
	    reader->readAttribute (i, &int4Buffer[0]);
	  }
	else
	  {
	    CCTK_WARN (1, "Unable to restore GH->cctk_iteration, defaulting to 0");
	    int4Buffer[0] = 0;
	  }

	/* get the main loop index */
	i = reader->readAttributeInfo ( "main loop index", datatype, dim);
	if (i >= 0 && datatype == FLEXIO_INT4 && dim == 1)
	  {
	    reader->readAttribute (i, &int4Buffer[1]);
	  }
	else
	  {
	    CCTK_WARN (1, "Unable to restore main loop index, defaulting to 0");
	    int4Buffer[1] = 0;
	  }
	
	/* get cctk_time */
	i = reader->readAttributeInfo ("GH$time", datatype, dim);
	if (i >= 0 && datatype == FLEXIO_REAL && dim == 1)
	  {
	    reader->readAttribute (i, &realBuffer);
	  }
	else
	  {
	    CCTK_WARN (1, "Unable to restore GH->cctk_time, defaulting to 0.0");
	    realBuffer = 0.0;
	  }
      }
    
#ifdef CCTK_MPI
  /* Broadcast the GH extensions to all processors */
  /* NOTE: We have to use MPI_COMM_WORLD here
     because PUGH_COMM_WORLD is not yet set up at parameter recovery time.
     We also assume that PUGH_MPI_INT4 is a compile-time defined datatype. */
    CACTUS_MPI_ERROR (MPI_Bcast (int4Buffer, 2, CARPET_MPI_INT4, 0,MPI_COMM_WORLD));
    CACTUS_MPI_ERROR (MPI_Bcast (&realBuffer, 1, CARPET_MPI_REAL,0,MPI_COMM_WORLD));
#endif

    GH->cctk_time = realBuffer;
    GH->cctk_iteration = (int) int4Buffer[0];
    CCTK_SetMainLoopIndex ((int) int4Buffer[1]);

    return (0);
}




  int DumpParams (const cGH* const cgh, int all, IObase* writer){

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


  int DumpGHExtensions (const cGH* const cgh, IObase* writer){

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


  int Checkpoint (const cGH* const cgh, int called_from)
  {
    char cp_filename[1024], cp_tempname[1024];
    int myproc, first_vindex, gindex;
    char *fullname;
    const char *timer_descriptions[3] = {"Time to dump parameters: ",
					 "Time to dump datasets:   ",
					 "Total time to checkpoint:"};
    const ioGH *ioUtilGH;
    

    //    const int varindex = CCTK_VarIndex("ADMBASE:gxx");
    int varindex = 0;
    int group = 0;
    int retval = 0;

    cGroup  gdata;
    IObase* writer = 0;
    AMRwriter* amrwriter = 0;
    ioRequest *request;

    DECLARE_CCTK_PARAMETERS

    CarpetIOFlexIOGH *myGH;
    myGH = (CarpetIOFlexIOGH *) CCTK_GHExtension (cgh, "CarpetIOFlexIO");

    /* check if CarpetIOFlexIO was registered as I/O method */
    if (myGH == NULL)
    {
      CCTK_WARN (-1, "No CarpetIOFlexIO I/O methods registered");
      return (-1);
    }

  
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
  
    sprintf(cp_tempname,"%s.tmp%s",cp_filename,extension);


    sprintf(cp_filename,"%s%s",cp_filename,extension);
  

    if (CCTK_MyProc(cgh)==0)
      {
	if (verbose)
	  {
	    CCTK_VInfo (CCTK_THORNSTRING, "Creating temporary checkpoint file '%s'", cp_tempname);
	  }

	//	writer = new IEEEIO(cp_tempname, IObase::Create);

	if (CCTK_Equals(out3D_format, "IEEE")) {
	  writer = new IEEEIO(cp_tempname, IObase::Create);
#ifdef HDF4
	} else if (CCTK_Equals(out3D_format, "HDF4")) {
	    writer = new HDFIO(cp_tempname, IObase::Create);
#endif
#ifdef HDF5
	} else if (CCTK_Equals(out3D_format, "HDF5")) {
	  writer = new H5IO(cp_tempname, IObase::Create);
#endif
	} else {
	  assert (0);
	}

	if (! (writer->isValid()) )
	  {
	    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
			"Can't open checkpoint file '%s'. Checkpointing is skipped",
			cp_tempname);
	    return (-1);
	  }

	amrwriter = new AMRwriter(*writer);

	/* now we are writing a first (bogus) dataset to which
	   we will attach all parameters and GHextensions as Attributes
	*/

	
       	CCTK_REAL startdata = 666.66;
	int rank=1;
	int dim[1]={1};
	writer->write(FLEXIO_REAL,rank,dim,&startdata);
	
	/* now dump parameters */ 
	if (verbose)
	  {
	    CCTK_VInfo (CCTK_THORNSTRING, "Dumping Parameters'");
	  }
	DumpParams (cgh, 1, writer);


	/* and now dump GH extentions */
	if (verbose)
	  {
	    CCTK_VInfo (CCTK_THORNSTRING, "Dumping GHExtensions");
	  }

	DumpGHExtensions(cgh,writer);

      }

    
    // now dump the grid varibles for all reflevels and components, sorted by groups
    BEGIN_REFLEVEL_LOOP(cgh) {

  	BEGIN_MGLEVEL_LOOP(cgh) {
	      
	    if (verbose)
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

	      const int grouptype = CCTK_GroupTypeI(group);

	      /* scalars and grid arrays only have 1 reflevel: */
	      if ( (grouptype != CCTK_GF) && (reflevel != 0) )
		continue;

	      /* now check if there is any memory allocated
                 for GFs and GAs. GSs should always have 
		 memory allocated and there is at this point
                 no CCTK function to check this :/
	      */

	      if ( (grouptype == CCTK_GF) || (grouptype == CCTK_ARRAY)){
		const int gpdim = CCTK_GroupDimI(group);
		int gtotalsize=1;
		for(int d=0;d<gpdim;d++){
		  const int* gpsize= CCTK_ArrayGroupSizeI(cgh,d,group);
		  assert(gpsize != NULL);		  
		  gtotalsize*=gpsize[d];		  
		}
		if(gtotalsize == 0){
		  if (verbose) CCTK_VInfo(CCTK_THORNSTRING, "Group %s is zero-sized. No checkpoint info written",CCTK_GroupName(group));
		  continue;
		}
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

		  //#if 1	    
		  if (grouptype == CCTK_SCALAR)
		  {
		    //		    retval += WriteGS(cgh,writer,request);
		    retval += WriteGF(cgh,writer,amrwriter,request);
		  }
		  else 
		    //#endif
		  if (grouptype == CCTK_ARRAY || grouptype == CCTK_GF)
		    //else if (grouptype == CCTK_GF)
		  {
		    char* fullname = CCTK_FullName (request->vindex);
		    if (verbose)
		      CCTK_VInfo (CCTK_THORNSTRING,"%s:: reflevel: %d component: %d grouptype: %d ",fullname,reflevel,component,grouptype);
		    free(fullname);
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
	
	} END_MGLEVEL_LOOP;
    
    } END_REFLEVEL_LOOP;


    // Close the temporary file
    if (CCTK_MyProc(cgh)==0) {
      delete amrwriter;
      amrwriter = 0;
      delete writer;
      writer = 0;
    }

    CCTK_VInfo(CCTK_THORNSTRING,"retval: %d",retval);

    if (retval == 0)
    {
      if (rename (cp_tempname, cp_filename))
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not rename temporary checkpoint file '%s' to '%s'",
                    cp_tempname, cp_filename);
        retval = -1;
      }
      else
      {
        if (myGH->cp_filename_list[myGH->cp_filename_index])
        {
          if (checkpoint_keep > 0)
          {
            remove (myGH->cp_filename_list[myGH->cp_filename_index]);
          }
          free (myGH->cp_filename_list[myGH->cp_filename_index]);
        }
        myGH->cp_filename_list[myGH->cp_filename_index] = strdup (cp_filename);
        myGH->cp_filename_index = (myGH->cp_filename_index+1) % abs (checkpoint_keep);
      }
    }
    

    return 0;
  }

  int RecoverVariables (cGH* cgh, IObase* reader, AmrGridReader* amrreader){

    CCTK_VInfo(CCTK_THORNSTRING,"Starting to recover data for refinement level %d",reflevel);




    return 0;
  }


} // namespace CarpetCheckpointRestart














