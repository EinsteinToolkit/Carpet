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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5chckpt_recover.cc,v 1.18 2004/03/23 19:30:14 schnetter Exp $";
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

  // linked list for reading in the checkpoint file

  list<string> datasetnamelist;
  
  
  int Checkpoint (const cGH* const cctkGH, int called_from);
  int DumpParametersGHExtentions (const cGH *cctkGH, int all, hid_t writer);  
  
  int RecoverParameters (hid_t reader);
  int RecoverGHextensions (cGH* cctkGH, hid_t reader);
  int RecoverVariables (cGH* cctkGH, hid_t reader);

  void CarpetIOHDF5_InitialDataCheckpoint( const cGH* const cgh){
    
    DECLARE_CCTK_PARAMETERS;
    
    if (checkpoint && checkpoint_ID)
      {
	CCTK_INFO ("---------------------------------------------------------");
	CCTK_INFO ("Dumping initial data checkpoint");
	CCTK_INFO ("---------------------------------------------------------");
	
	Checkpoint (cgh, CP_INITIAL_DATA);
      }
  } // CarpetIOHDF5_InitialDataCheckpoint


  void CarpetIOHDF5_EvolutionCheckpoint( const cGH* const cgh){
    
    DECLARE_CCTK_PARAMETERS

      if (checkpoint &&
      ((checkpoint_every > 0 && cgh->cctk_iteration % checkpoint_every == 0) ||
       checkpoint_next))
      {
	CCTK_INFO ("---------------------------------------------------------");
	CCTK_VInfo (CCTK_THORNSTRING, "Dumping periodic checkpoint at "
		    "iteration %d", cgh->cctk_iteration);
	CCTK_INFO ("---------------------------------------------------------");

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

    
    static hid_t reader; //this thing absolutely needs to be static!!!

    myGH = NULL;
    result = 0;

    myproc = CCTK_MyProc (cctkGH);
    
  
    if (called_from == CP_RECOVER_PARAMETERS) {
      // Okay, let's see what we can do about the parameters
    
       // Invent a file name
      ostringstream filenamebuf;

      if(CCTK_nProcs(cctkGH) == 1)
	filenamebuf << recover_dir << "/" << basefilename << ".h5";
      else
	filenamebuf << recover_dir << "/" << basefilename << ".file_0.h5";

      string filenamestr = filenamebuf.str();
      const char * const filename = filenamestr.c_str();
 
      if (myproc == 0) {
	// First, open the file
	if (h5verbose) 
	  CCTK_VInfo(CCTK_THORNSTRING, "Opening Checkpoint file %s for recovery",filename);
	reader = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (reader<0) {
	  CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "Could not recover from \"%s\"", filename);
	}
	assert (reader>=0);
      } // myproc == 0
    }
    else {
	/* This is the case for CP_RECOVER_DATA.
	   CCTK_RECOVER_PARAMETERS must have been called before
	   and set up the file info structure. */
	if (myproc == 0) {
	  assert(reader>=0);
	}
    }
  
    if (called_from == CP_RECOVER_PARAMETERS)
      {
	return (RecoverParameters (reader));
      }
  
    if (called_from == CP_RECOVER_DATA) {
      CCTK_INT4 numberofmgtimes;

      CCTK_VInfo(CCTK_THORNSTRING,"Starting to recover data on reflevel %d!!!",reflevel);

      if (myproc == 0) {

	/* we need all the times on the individual levels */
	// these are a bit messy to extract

	// Actually, we do have to do this only once

	  hid_t group = H5Gopen (reader, PARAMETERS_GLOBAL_ATTRIBUTES_GROUP);
	  assert(group>=0);
	  hid_t dataset = H5Dopen (group, ALL_PARAMETERS);
	  assert(dataset >= 0);
	  hid_t attr = H5Aopen_name (dataset, "numberofmgtimes");
	  assert(attr >= 0);
	  hid_t atype = H5Aget_type (attr);
	  if(H5Tequal(atype, H5T_NATIVE_INT)) {
	    herr_t herr = H5Aread(attr, atype, &numberofmgtimes);
	    assert(!herr);
	    herr = H5Aclose(attr);
	    assert(numberofmgtimes==mglevels);
	    char buffer[100];
	    for(int lcv=0;lcv<numberofmgtimes;lcv++) {
	      sprintf(buffer,"mgleveltimes %d",lcv);
	      attr = H5Aopen_name(dataset, buffer);
	      assert (attr>=0);
	      atype = H5Aget_type (attr);
	      assert (atype>=0);
	      herr = H5Aread (attr, atype, &leveltimes.at(lcv).at(0));
	      assert(!herr);
	      herr = H5Aclose(attr);
	      assert(!herr);
	    }  
	    herr = H5Dclose(dataset);
	    assert(!herr);
	    herr = H5Gclose(group);
	    assert(!herr);
	  } else {
	    CCTK_WARN(0,"BAD BAD BAD! Can't read leveltimes!!");
	  }

      } // myproc == 0  	


      // communicate the time on all the mglevels and reflevels

      int mpierr = MPI_Bcast (&numberofmgtimes, 1, CARPET_MPI_INT4, 0,MPI_COMM_WORLD);
      assert(!mpierr);


      for(int i=0;i<numberofmgtimes;i++) {
	mpierr = MPI_Bcast (&(leveltimes.at(i).at(0)), reflevels, CARPET_MPI_REAL, 0, MPI_COMM_WORLD);
	assert(!mpierr);
      }

      if (h5verbose) cout << "leveltimes: " << leveltimes << endl;

      cctkGH->cctk_time = leveltimes.at(mglevel).at(reflevel);

      result += RecoverGHextensions(cctkGH,reader);

	  
      if (h5verbose) cout << "reflevel: " << reflevel << endl;
      result += RecoverVariables (cctkGH,reader);


      CCTK_VInfo (CCTK_THORNSTRING,
		  "Restarting simulation at iteration %d (physical time %g)",
		  cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
    } // called_from == CP_RECOVER_DATA
  
    if (myproc == 0 && reflevel==maxreflevels) 	
      H5Fclose(reader);

    return (result);
  } // CarpetIOHDF5_Recover

  
  int RecoverVariables (cGH* cctkGH, hid_t reader) {

    DECLARE_CCTK_PARAMETERS;

    int retval = 0;
    int myproc = CCTK_MyProc (cctkGH);
    char * name;

    int ndatasets;

    int varindex; 

    double datasettime;
    double leveltime;
    static double totaltime;

    hid_t dataset;
    herr_t herr;

    list<string> refleveldatasetnamelist;

    if (reflevel==0) {
      totaltime = 0;
    }

    leveltime = MPI_Wtime();


    if(myproc==0) {
      ndatasets=GetnDatasets(reader);
      assert (ndatasets>=0);
    }

    // Broadcast number of datasets
    MPI_Bcast (&ndatasets, 1, MPI_INT, 0, dist::comm);

    assert ((ndatasets)>=0);
    
    //if (h5verbose && reflevel==0) cout << "ndatasets: " << ndatasets << endl;
    if ( reflevel == 0) cout << "ndatasets: " << ndatasets << endl;

    if (reflevel==0) {
      for (int currdataset=0;currdataset<ndatasets+1;currdataset++){
	char datasetname[256];
	if (myproc==0) {
	  GetDatasetName(reader,currdataset,datasetname);
	  datasetnamelist.push_back(datasetname);
	} //myproc = 0
	else {
	  datasetnamelist.push_back("blah");
	}
      }
    } 

    cout << "I have " << datasetnamelist.size() << endl;

    double comparetime = MPI_Wtime();
      
    if(myproc==0) {
      //      for(currdataset=datasetnamelist.begin();
      //	  currdataset!=datasetnamelist.end();
      //	  ++currdataset) {

      list<string>::iterator currdataset;
      currdataset=datasetnamelist.begin();
       while(currdataset!=datasetnamelist.end()) {
	char tempstr[10];
	sprintf(tempstr,"rl=%d m=",reflevel);
	if ( (*currdataset).find(tempstr) < (*currdataset).size() ) {
	  refleveldatasetnamelist.push_back((*currdataset).c_str());
	  currdataset = datasetnamelist.erase(currdataset);
	} else {
	  ++currdataset;
	} // if ...
      } // while ...
    } // if(myproc==0)




    long reflevelnamenum;
    if(myproc==0) reflevelnamenum=refleveldatasetnamelist.size();
    MPI_Bcast (&reflevelnamenum, 1, MPI_LONG, 0, dist::comm);
    assert ((reflevelnamenum)>=0);
    
    // fill bogus namelist for non-IO cpus
    if(myproc !=0) {
      for(long i=0;i < reflevelnamenum;i++) {
	refleveldatasetnamelist.push_back("blah");
      }
    }

    comparetime = MPI_Wtime() - comparetime;
    cout << "Time for string comparison: " << comparetime << endl;
    cout << "I have for this reflevel " << refleveldatasetnamelist.size() << endl;

    list<string>::iterator currdataset;

      currdataset=refleveldatasetnamelist.begin();
 
      while(currdataset!=refleveldatasetnamelist.end()) {

	//	cout << "name: " << (*currdataset).c_str() << endl;

      if (myproc==0) {
	dataset = H5Dopen (reader, (*currdataset).c_str());
	assert(dataset);
	// Read data
	ReadAttribute (dataset, "name", name);
	varindex = CCTK_VarIndex(name);
      }

      MPI_Bcast (&varindex, 1, MPI_INT, 0, dist::comm);

      name = CCTK_FullName(varindex);

      if (h5verbose) cout << name << "  rl: " << reflevel << endl;
      vector<ibset> regions_read(Carpet::maps);

      assert (varindex>=0 && varindex<CCTK_NumVars());
      const int group = CCTK_GroupIndexFromVarI (varindex);
      const int grouptype = CCTK_GroupTypeI(group);

      int did_read_something = ReadVar(cctkGH,reader,name,dataset,regions_read,1);

      MPI_Bcast (&did_read_something, 1, MPI_INT, 0, dist::comm);

      if (did_read_something) {
      	currdataset = refleveldatasetnamelist.erase(currdataset);
      } else {
	++currdataset;
      }

      if(myproc==0) {
	herr = H5Dclose(dataset);
	assert(!herr);
      }
      free(name);

    }  // for (currdataset ... )
    leveltime = MPI_Wtime() - leveltime;
    totaltime = totaltime + leveltime;

    cout << "Timers: leveltime: " << leveltime << " totaltime: " << totaltime << endl; 

    return retval;
  }
   


  
  int RecoverGHextensions (cGH *cctkGH, hid_t reader)
  {
    const int myproc = CCTK_MyProc(cctkGH);
    CCTK_INT4 int4Buffer[3];
    CCTK_REAL realBuffer;
    CCTK_REAL realBuffer2;
    CCTK_INT4 intbuffer;

    int mpierr = 0;

    if (myproc==0)
      {
		
	// First open group and dataset
	hid_t group = H5Gopen (reader, PARAMETERS_GLOBAL_ATTRIBUTES_GROUP);
	assert(group>=0);
	hid_t dataset = H5Dopen (group, ALL_PARAMETERS);
	assert(dataset >= 0);

	ReadAttribute(dataset,"GH$iteration",int4Buffer[0]);
	ReadAttribute(dataset,"main loop index",int4Buffer[1]);
	ReadAttribute(dataset,"carpet_global_time",realBuffer);
	//	ReadAttribute(dataset,"carpet_reflevels",int4Buffer[2]);
	ReadAttribute(dataset,"carpet_delta_time",realBuffer2);

	herr_t herr = H5Dclose(dataset);
	assert(!herr);
	herr = H5Gclose(group);
	assert(!herr);
	
      }
  /* Broadcast the GH extensions to all processors */
  /* NOTE: We have to use MPI_COMM_WORLD here
     because PUGH_COMM_WORLD is not yet set up at parameter recovery time.
     We also assume that PUGH_MPI_INT4 is a compile-time defined datatype. */

    mpierr = MPI_Bcast (int4Buffer, 3, CARPET_MPI_INT4, 0,MPI_COMM_WORLD);
    assert(!mpierr);
    mpierr = MPI_Bcast (int4Buffer, 3, CARPET_MPI_INT4, 0,MPI_COMM_WORLD);
    assert(!mpierr);
    mpierr = MPI_Bcast (&realBuffer, 1, CARPET_MPI_REAL,0,MPI_COMM_WORLD);
    assert(!mpierr);
    mpierr = MPI_Bcast (&realBuffer2, 1, CARPET_MPI_REAL,0,MPI_COMM_WORLD);
    assert(!mpierr);

    global_time = (CCTK_REAL) realBuffer;
    delta_time = (CCTK_REAL) realBuffer2;
    //    reflevels = (int) int4Buffer[2];
    cctkGH->cctk_iteration = (int) int4Buffer[0];
    CCTK_SetMainLoopIndex ((int) int4Buffer[1]);


    return (0);

  } // RecoverGHExtensions

  int RecoverParameters(hid_t reader){

    DECLARE_CCTK_PARAMETERS; 
    
    int myproc, retval;
    char *parameters;
    CCTK_INT4 parameterSize;

    hid_t group,dataset;
    herr_t herr;

    int mpierr;

    myproc = CCTK_MyProc (NULL);

    if (myproc == 0){
      CCTK_VInfo (CCTK_THORNSTRING, "Recovering parameters from checkpoint ");
      
      group = H5Gopen (reader, PARAMETERS_GLOBAL_ATTRIBUTES_GROUP);
      assert(group >= 0);
      dataset = H5Dopen (group, ALL_PARAMETERS);
      assert(dataset>= 0);

      parameterSize = H5Dget_storage_size (dataset);
      parameters = (char *) malloc (parameterSize);
      herr = H5Dread (dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, parameters);
      assert(!herr);
      herr = H5Dclose(dataset);
      assert(!herr);
      herr = H5Gclose(group);
      assert(!herr);

      if(h5verbose) 
	CCTK_VInfo (CCTK_THORNSTRING, "\n%s\n",parameters);
      
      CCTK_VInfo(CCTK_THORNSTRING, "Successfully recovered parameters!");
    } // myproc == 0

    /* Broadcast the parameter buffer size to all processors */
    /* NOTE: We have to use MPI_COMM_WORLD here
       because CARPET_COMM_WORLD is not yet set up at parameter recovery time.
       We also assume that CARPET_MPI_INT4 is a compile-time defined datatype. */
    mpierr = MPI_Bcast (&parameterSize, 1, CARPET_MPI_INT4, 0,
			MPI_COMM_WORLD);
    assert(!mpierr);

    if (parameterSize > 0) {
      if (myproc) {
	parameters = (char*) malloc (parameterSize + 1);
      }
	
      mpierr = MPI_Bcast (parameters, parameterSize + 1, CARPET_MPI_CHAR,
			  0,MPI_COMM_WORLD);
      assert(!mpierr);
      
      IOUtil_SetAllParameters (parameters);
    
      free (parameters);
    }
  
    /* return positive value for success otherwise negative */
    retval = (parameterSize > 0 ? 1 : -1);

    return (retval);

  } // RecoverParameters



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

    // I didn't like what the flesh provides:
    
    IOUtil_PrepareFilename (cctkGH, NULL, cp_filename, called_from,
           myproc / ioUtilGH->ioproc_every, ioUtilGH->unchunked);


    /* ... and append the extension */
    sprintf (cp_tempname, "%s.tmp.h5", cp_filename);
    sprintf (cp_filename, "%s.h5",     cp_filename);

    hid_t writer = -1;

    if (myproc == 0) {

      if (h5verbose) {
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

      // Dump all parameters and GHExtentions
      retval += DumpParametersGHExtentions (cctkGH, 1, writer);
      assert(!retval);

    } // myproc == 0 


    // now dump the grid variables on all mglevels, reflevels, maps and components
    BEGIN_MGLEVEL_LOOP(cctkGH) {

      BEGIN_REFLEVEL_LOOP(cctkGH) {

	if (h5verbose)
	  {
	    CCTK_INFO ("Dumping Grid Variables ...");
	  }
	for (int group = CCTK_NumGroups () - 1; group >= 0; group--)
	  {
	    /* only dump groups which have storage assigned */

	    if (CCTK_QueryGroupStorageI (cctkGH, group) <= 0)
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
	      int tlsh[gpdim];
	      assert(!CCTK_GrouplshGI(cctkGH,gpdim,tlsh,group));
	      for(int i=0;i<gpdim;i++) {
		gtotalsize=tlsh[i];		  
	      }
	      if(gtotalsize == 0){
		if (h5verbose) CCTK_VInfo(CCTK_THORNSTRING, 
		    "Group %s is zero-sized. No checkpoint info written",CCTK_GroupName(group));
		continue;
	      }
	    }
	    
	    /* get the number of allocated timelevels */
	    CCTK_GroupData (group, &gdata);
	    gdata.numtimelevels = 0;
	    gdata.numtimelevels = CCTK_GroupStorageIncrease (cctkGH, 1, &group,
							     &gdata.numtimelevels,NULL);
	    
	    
	    
	    int first_vindex = CCTK_FirstVarIndexI (group);

	    /* get the default I/O request for this group */
	    request = IOUtil_DefaultIORequest (cctkGH, first_vindex, 1);
	    
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
		    if (h5verbose)
		      {
			fullname = CCTK_FullName (request->vindex);
			CCTK_VInfo (CCTK_THORNSTRING, "  %s (timelevel %d)",
				    fullname, request->timelevel);
			free (fullname);
		      }
		    // write the var
		    
		    if (grouptype == CCTK_ARRAY || grouptype == CCTK_GF || grouptype == CCTK_SCALAR)
		      {
			char* fullname = CCTK_FullName (request->vindex);
			if (h5verbose)
			  CCTK_VInfo (CCTK_THORNSTRING,"%s: reflevel: %d map: %d component: %d grouptype: %d ",
				      fullname,reflevel,Carpet::map,component,grouptype);
			free(fullname);
		    retval += WriteVar(cctkGH,writer,request,1);
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
      } END_REFLEVEL_LOOP;
      
    } END_MGLEVEL_LOOP;
  

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
  
    if (h5verbose) {
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

      dtmp = global_time;
      WriteAttribute(dataset,"carpet_global_time", dtmp);

      itmp = reflevels;
      WriteAttribute(dataset,"carpet_reflevels", itmp);

      dtmp = delta_time;
      WriteAttribute(dataset,"carpet_delta_time", dtmp);

      version = CCTK_FullVersion();
      WriteAttribute(dataset,"Cactus version", version);

      /* finally, we need all the times on the individual refinement levels */

      const int numberofmgtimes=mglevels;
      WriteAttribute(dataset,"numberofmgtimes",numberofmgtimes);
      for(int i=0;i < numberofmgtimes;i++) {
	char buffer[100];
	sprintf(buffer,"mgleveltimes %d",i);
	WriteAttribute(dataset,buffer,(double *) &leveltimes.at(i).at(0), reflevels);
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
