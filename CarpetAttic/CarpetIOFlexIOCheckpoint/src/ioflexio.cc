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

// Hack to stop FlexIO data type clash with LAM MPI
#undef BYTE
#undef CHAR

#include "CactusBase/IOUtil/src/ioGH.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "ioflexio.hh"


extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/ioflexio.cc,v 1.20 2004/01/13 15:46:52 cott Exp $";
  CCTK_FILEVERSION(Carpet_CarpetIOFlexIO_ioflexio_cc);
}



namespace CarpetIOFlexIO {
  
  using namespace std;
  using namespace Carpet;
  using namespace CarpetIOFlexIOUtil;  
  using namespace CarpetCheckpointRestart;
  
  // Variable definitions
  //  int GHExtension;
  int IOMethod;
  vector<bool> do_truncate;
  vector<vector<int> > last_output;
  
  
  
  static const char* GetStringParameter (const char* const parametername,
					 const char* const fallback);
  static int GetIntParameter (const char* const parametername, int fallback);
  static bool CheckForVariable (const cGH* const cgh,
				const char* const varlist, const int vindex);
  static void SetFlag (int index, const char* optstring, void* arg);
  
  
  
  int CarpetIOFlexIO_Startup ()
  {
    CCTK_RegisterBanner ("AMR 3D FlexIO I/O provided by CarpetIOFlexIO");
    
    //    GHExtension = CCTK_RegisterGHExtension("CarpetIOFlexIO");
    CCTK_RegisterGHExtensionSetupGH (CCTK_RegisterGHExtension("CarpetIOFlexIO"),SetupGH);
    
    IOMethod = CCTK_RegisterIOMethod ("CarpetIOFlexIO");
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
  
    /* register the CarpetIOFlexIO recovery routine to thorn IOUtil */
    if (IOUtil_RegisterRecover ("CarpetIOFlexIO recovery", CarpetIOFlexIO_Recover) < 0)
      {
	CCTK_WARN (1, "Failed to register IOFlexIO recovery routine");
      }
    
  
    return 0;
  }
  
  
  
  void* SetupGH (tFleshConfig* const fc, const int convLevel, cGH* const cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    CarpetIOFlexIOGH* myGH;
    CCTK_INT i;

    // Truncate all files if this is not a restart
    do_truncate.resize(CCTK_NumVars(), true);
    
    // No iterations have yet been output
    last_output.resize(maxreflevels);
    for (int rl=0; rl<maxreflevels; ++rl) {
      last_output[rl].resize(CCTK_NumVars(), INT_MIN);
    }
    
    // We register only once, ergo we get only one handle.  We store
    // that statically, so there is no need to pass anything to
    // Cactus.

    /* allocate a new GH extension structure */


    CCTK_INT numvars = CCTK_NumVars ();
    myGH            = (CarpetIOFlexIOGH*) malloc (sizeof (CarpetIOFlexIOGH));
    myGH->out_last  = (int *) malloc (numvars * sizeof (int));
    myGH->requests  = (ioRequest **) calloc (numvars, sizeof (ioRequest *));
    myGH->cp_filename_list = (char **) calloc (abs (checkpoint_keep), sizeof (char *));
    myGH->cp_filename_index = 0;
    myGH->out_vars = strdup ("");
    myGH->out_every_default = out_every - 1;

    for (i = 0; i < numvars; i++)
    {
      myGH->out_last [i] = -1;
    }

    myGH->open_output_files = NULL;

  
    return (myGH);

    return 0;
  }
  
  
  int OutputGH (const cGH* const cgh) {
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cgh, vindex)) {
	TriggerOutput(cgh, vindex);
      }
    }
    return 0;
  }
  
  

  int WriteGF (const cGH* const cgh, IObase* writer, AMRwriter* amrwriter, ioRequest* request, const int called_from_checkpoint)
  {

    DECLARE_CCTK_PARAMETERS;

    const int varindex  = request->vindex;

    const int group = CCTK_GroupIndexFromVarI (varindex);
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = varindex - n0;
    assert (var>=0 && var<CCTK_NumVars());
    int tl = 0;
    const int grouptype = CCTK_GroupTypeI(group);

    // let's get the correct Carpet time level (which is the (-1) * Cactus timelevel):
    if (request->timelevel==0)
      tl = 0;
    else
      tl = - request->timelevel;

    assert (! ( (grouptype != CCTK_GF) && reflevel>0));

    if (CCTK_MyProc(cgh)==0) {

      amrwriter->setType (FlexIODataType(CCTK_VarTypeI(varindex)));
      
      int gpdim = CCTK_GroupDimI(group);

      // need gpdim=1 if scalar (flexio wants this)
      if(gpdim == 0)
	if(grouptype == CCTK_SCALAR)
	  gpdim = 1;
        else
	  CCTK_WARN(0,"Non-scalar variable with dimension 0!!!");
	

      // Set coordinate information
      CCTK_REAL lower[dim], upper[dim];
      double origin[dim], delta[dim], timestep;
      for (int d=0; d<dim; ++d) {
	const int ierr = CCTK_CoordRange
          (cgh, &lower[d], &upper[d], d+1, 0, "cart3d");
	assert (ierr==0);
	origin[d] = lower[d];
	delta[d] = cgh->cctk_delta_space[d];
      }
      timestep = cgh->cctk_delta_time;
      amrwriter->setTopLevelParameters
	(gpdim, origin, delta, timestep, maxreflevels);
      
      // Set refinement information
      int interlevel_timerefinement;
      int interlevel_spacerefinement[dim];
      int initial_gridplacementrefinement[dim];
      interlevel_timerefinement = hh->reffact;
      for (int d=0; d<dim; ++d) {
	interlevel_spacerefinement[d] = hh->reffact;
	initial_gridplacementrefinement[d] = 1;
      }
      amrwriter->setRefinement
	(interlevel_timerefinement, interlevel_spacerefinement,
	 initial_gridplacementrefinement);
      
      // Set level
      amrwriter->setLevel (reflevel);
      
      // Set current time
      amrwriter->setTime (cgh->cctk_iteration);
    }


    
    // Traverse all components on this refinement and multigrid level
    BEGIN_COMPONENT_LOOP(cgh, grouptype) {

      const ggf<dim>* ff = 0;

      assert (var < (int)arrdata[group].data.size());
      ff = (ggf<dim>*)arrdata[group].data[var];
      
      const gdata<dim>* const data = (*ff) (tl, reflevel, component, mglevel);

      // get some more group information
      cGroupDynamicData gdyndata;

      int ierr = CCTK_GroupDynamicData(cgh,group,&gdyndata);
      assert(ierr==0);

      cGroup cgdata;
      ierr = CCTK_GroupData(group,&cgdata);
      assert(ierr==0);

      /* handle CCTK_DISTRIB_CONSTANT scalar and arrays */
#if 0
      if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
	assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
	if (hh->processors[reflevel][component] == 0) {
	  if (grouptype == CCTK_SCALAR) {
	    CCTK_VInfo (CCTK_THORNSTRING, "dumping SCALAR distrib const");
	    int rank=1;
	    int dim[1]={1};
	    writer -> write(FlexIODataType(CCTK_VarTypeI(varindex)),rank,dim,CCTK_VarDataPtrI(cgh,tl,varindex));
	    DumpCommonAttributes(cgh,writer,request);
	    continue;
	  }
	  else {
	    writer -> write(FlexIODataType(CCTK_VarTypeI(varindex)),cgdata.dim,gdyndata.lsh,CCTK_VarDataPtrI(cgh,tl,varindex));
	    DumpCommonAttributes(cgh,writer,request);
	    continue;
	  }
	}
	else {
	  continue;
	}

      }
#endif

      // Make temporary copy on processor 0
      bbox<int,dim> ext = data->extent();
      vect<int,dim> lo = ext.lower();
      vect<int,dim> hi = ext.upper();
      vect<int,dim> str = ext.stride();

      // Ignore ghost zones if desired

      const int out3D_output_outer_boundary_var = (called_from_checkpoint) ? -1 : out3D_output_outer_boundary;
      const int out3D_max_num_lower_ghosts_var = (called_from_checkpoint) ? -1 : out3D_max_num_lower_ghosts;
      const int out3D_max_num_upper_ghosts_var = (called_from_checkpoint) ? -1 : out3D_max_num_upper_ghosts;
      

      for (int d=0; d<dim; ++d) {
	const int max_lower_ghosts = (gdyndata.bbox[2*d  ] && out3D_output_outer_boundary_var) ? -1 : out3D_max_num_lower_ghosts_var;
	const int max_upper_ghosts = (gdyndata.bbox[2*d+1] && out3D_output_outer_boundary_var) ? -1 : out3D_max_num_upper_ghosts_var;
	
	const int num_lower_ghosts = max_lower_ghosts == -1 ? gdyndata.nghostzones[d] : min(out3D_max_num_lower_ghosts_var, gdyndata.nghostzones[d]);
	const int num_upper_ghosts = max_upper_ghosts == -1 ? gdyndata.nghostzones[d] : min(out3D_max_num_upper_ghosts_var, gdyndata.nghostzones[d]);
	
	lo[d] += (gdyndata.nghostzones[d] - num_lower_ghosts) * str[d];
	hi[d] -= (gdyndata.nghostzones[d] - num_upper_ghosts) * str[d];
      }
      
      ext = bbox<int,dim>(lo,hi,str);
      
      gdata<dim>* const tmp = data->make_typed (varindex);
      tmp->allocate (ext, 0);
      //fprintf(stderr,"\n writing1 %d\n",CCTK_MyProc(cgh));
      if ( !((cgdata.disttype == CCTK_DISTRIB_CONSTANT) && (hh->processors[reflevel][component]!=0))) {

	if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
	  assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
	  //fprintf(stderr,"\n scalar %d %d comp: %d\n",CCTK_MyProc(cgh),varindex,component);
	  int origin[dim], dims[dim];
	  for (int d=0; d<dim; ++d) {
	    origin[d] = (ext.lower() / ext.stride())[d];
	    dims[d]   = (ext.shape() / ext.stride())[d];
	  }
	  if (CCTK_MyProc(cgh)==0) {
	    amrwriter->write (origin, dims, (void*)data->storage());
	    DumpCommonAttributes(cgh,writer,request);
	  }
	  delete tmp;
	  continue;
	} else {

	  for (comm_state<dim> state; !state.done(); state.step()) {
	    tmp->copy_from (state, data, ext);
	  }

	  //fprintf(stderr,"\n writing2 %d component: %d varindex: %d distrib_const: %d\n",CCTK_MyProc(cgh),component,varindex,(cgdata.disttype == CCTK_DISTRIB_CONSTANT));
	  // Write data
	  if (CCTK_MyProc(cgh)==0) {
	    int origin[dim], dims[dim];
	    for (int d=0; d<dim; ++d) {
	      origin[d] = (ext.lower() / ext.stride())[d];
	      dims[d]   = (ext.shape() / ext.stride())[d];
	    }
	    
	    amrwriter->write (origin, dims, (void*)tmp->storage());
	    
	    // dump attributes
	    DumpCommonAttributes(cgh,writer,request);
	    
	  }
	  // Delete temporary copy

	  delete tmp;

	}
      }
    } END_COMPONENT_LOOP;


    return 0;
  }


  int OutputVarAs (const cGH* const cgh, const char* const varname,
		   const char* const alias) {
    DECLARE_CCTK_PARAMETERS;
    const int n = CCTK_VarIndex(varname);
    assert (n>=0 && n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVars());
    const int tl = 0;
    
    // Check for storage
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot output variable \"%s\" because it has no storage",
		  varname);
      return 0;
    }
    
    const int grouptype = CCTK_GroupTypeI(group);
    if (grouptype != CCTK_GF && reflevel>0) return 0;

    int first_vindex = CCTK_FirstVarIndexI (group);
    /* get the default I/O request for this group */
    ioRequest* request = IOUtil_DefaultIORequest (cgh, first_vindex, 1);

    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cgh, "IO");
    assert (iogh);
    
    // Create the output directory
    const char* myoutdir = GetStringParameter("out3D_dir", out_dir);
    if (CCTK_MyProc(cgh)==0) {
      CCTK_CreateDirectory (0755, myoutdir);
    }
    
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
    extension = GetStringParameter ("out3D_extension", extension);
    
    ostringstream filenamebuf;
    filenamebuf << myoutdir << "/" << alias << extension;
    string filenamestr = filenamebuf.str();
    const char * const filename = filenamestr.c_str();
    
    IObase* writer = 0;
    AMRwriter* amrwriter = 0;
    
    // Write the file only on the root processor
    if (CCTK_MyProc(cgh)==0) {
      
      // If this is the first time, then create and truncate the file
      if (do_truncate[n]) {
	struct stat fileinfo;
	if (! iogh->recovered
	    || stat(filename, &fileinfo)!=0) {
	  writer = 0;
	  if (CCTK_Equals(out3D_format, "IEEE")) {
	    writer = new IEEEIO(filename, IObase::Create);
#ifdef HDF4
	  } else if (CCTK_Equals(out3D_format, "HDF4")) {
	    writer = new HDFIO(filename, IObase::Create);
#endif
#ifdef HDF5
	  } else if (CCTK_Equals(out3D_format, "HDF5")) {
	    writer = new H5IO(filename, IObase::Create);
#endif
	  } else {
	    assert (0);
	  }
	  delete writer;
	  writer = 0;
	}
      }
      
      // Open the file 
      if (CCTK_Equals(out3D_format, "IEEE")) {
	writer = new IEEEIO(filename, IObase::Append);
#ifdef HDF4
      } else if (CCTK_Equals(out3D_format, "HDF4")) {
	writer = new HDFIO(filename, IObase::Append);
#endif
#ifdef HDF5
      } else if (CCTK_Equals(out3D_format, "HDF5")) {
	writer = new H5IO(filename, IObase::Append);
#endif
      } else {
	assert (0);
      }
      assert (writer->isValid());
    
    
      amrwriter = new AMRwriter(*writer);
    }


      WriteGF(cgh,writer,amrwriter,request,0);
    
    // Close the file
    if (CCTK_MyProc(cgh)==0) {
      delete amrwriter;
      amrwriter = 0;
      delete writer;
      writer = 0;
    }
    
    // Don't truncate again
    do_truncate[n] = false;
    
    return 0;
  }
  
    
  
  int TimeToOutput (const cGH* const cgh, const int vindex) {
    DECLARE_CCTK_PARAMETERS;
    
    assert (vindex>=0 && vindex<(int)last_output[reflevel].size());
    
    const int myoutevery = GetIntParameter("out3D_every", out_every);
    
    if (myoutevery < 0) {
      // Nothing should be output at all
      return 0;
    }
    
    if (cgh->cctk_iteration % myoutevery != 0) {
      // Nothing should be output during this iteration
      return 0;
    }
    
    if (! CheckForVariable(cgh, GetStringParameter("out3D_vars",""), vindex)) {
      // This variable should not be output
      return 0;
    }
    
    if (last_output[reflevel][vindex] == cgh->cctk_iteration) {
      // Has already been output during this iteration
      char* varname = CCTK_FullName(vindex);
      CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Skipping output for variable \"%s\", because this variable "
		  "has already been output during the current iteration -- "
		  "probably via a trigger during the analysis stage",
		  varname);
      free (varname);
      return 0;
    }
    
    assert (last_output[reflevel][vindex] < cgh->cctk_iteration);
    
    // Should be output during this iteration
    return 1;
  }
  
  
  
  int TriggerOutput (const cGH* const cgh, const int vindex) {
    assert (vindex>=0 && vindex<CCTK_NumVars());
    
    char* varname = CCTK_FullName(vindex);
    const int retval = OutputVarAs (cgh, varname, CCTK_VarName(vindex));
    free (varname);
    
    last_output[reflevel][vindex] = cgh->cctk_iteration;
    
    return retval;
  }
  


  int ReadGF (const cGH* const cgh, IObase* reader, AmrGridReader* amrreader,int currdataset) {

    /* this functions reads in a variable on the current reflevel from an already open file. At
       some point it should be called from InputVarAs */


    DECLARE_CCTK_PARAMETERS;
    
    int tl = -1;
    int mglevel = -1;
    int rl = -1;
    int comp = -1;
    int myproc = CCTK_MyProc (cgh);
    int rank;
    int dims[dim];
    int nbytes;
    
    char* varname;
    char warnstring[256];
    int asize,i;
    IObase::DataType datatype;
    int group,varindex;
    CCTK_REAL cctk_time;

    if(myproc==0) {
      // read the name of the variable
      i = reader->readAttributeInfo ("name", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_CHAR && asize > 0)
	{
	  varname = (char*) malloc(sizeof(char)*asize+1);
	  reader->readAttribute (i, varname);
	}
      else
	{
	  CCTK_WARN (0, "Something is wrong! Can't read dataset names!!!"); 
	}
      
      varindex = CCTK_VarIndex(varname);

      assert(varindex > -1);
      
      group = CCTK_GroupIndexFromVarI(varindex);
      assert(group > -1);

      // Check for storage
      if (! CCTK_QueryGroupStorageI(cgh, group)) {
	CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "Cannot recover variable \"%s\" because it has no storage",
		    varname);
	return 0;
      }


      // read reflevel
      i = reader->readAttributeInfo ("reflevel", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_INT && asize > 0)
	{
	  reader->readAttribute (i, &rl);
	}
      else
	{
	  CCTK_WARN (0, "Something is wrong! Can't read refinement level!!!"); 
	}

      i = reader->readAttributeInfo ("component", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_INT && asize > 0)
	{
	  reader->readAttribute (i, &comp);
	}
      else
	{
	  CCTK_WARN (0, "Something is wrong! Can't read component!!!"); 
	}

      i = reader->readAttributeInfo ("timelevel", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_INT && asize > 0)
	{
	  reader->readAttribute (i, &tl);
	}
      else
	{
	  CCTK_WARN (0, "Something is wrong! Can't read timelevel!!!"); 
	}
      
      i = reader->readAttributeInfo ("mglevel", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_INT && asize > 0)
	{
	  reader->readAttribute (i, &mglevel);
	}
      else
	{
	  CCTK_WARN (0, "Something is wrong! Can't read multi group level!!!"); 
	}

      i = reader->readAttributeInfo ("cctk_time", datatype, asize);
      if (i >= 0 && datatype == FLEXIO_REAL && asize > 0)
	{
	  reader->readAttribute (i, &cctk_time);
	}
      else
	{
	  CCTK_WARN (0, "Something is wrong! Can't read coordinate time!!!"); 
	}
      


      // Read information about dataset
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Reading dataset info");
      reader->readInfo (datatype, rank, dims);
      nbytes = IObase::nBytes(datatype,rank,dims);
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "type=%d rank=%d dims=[%d,%d,%d] nbytes=%d", (int)datatype, rank, dims[0], dims[1], dims[2], nbytes);
      
      int gpdim = CCTK_GroupDimI(group);
      const int grouptype = CCTK_GroupTypeI(group);

      // need gpdim=1 if scalar (flexio wants this)
      if(gpdim == 0)
	if(grouptype == CCTK_SCALAR)
	  gpdim = 1;
      else
	CCTK_WARN(0,"Non-scalar variable with dimension 0!!!");


      CCTK_VInfo(CCTK_THORNSTRING,"Recovering varindex: %d grouptype: %d varname: %s tl: %d, rl: %d, c: %d",varindex,grouptype,varname,tl,rl,comp);

      free(varname);	
      
      // Check rank
      assert (rank==gpdim);
    }


    // Broadcast varindex,group,rank, dimensions, and nbytes,rl,tl,mglevel
    MPI_Bcast (&varindex, 1, MPI_INT, 0, dist::comm);
    assert (varindex>=0);
    MPI_Bcast (&group, 1, MPI_INT, 0, dist::comm);
    assert (group>=0);
    MPI_Bcast (&rank, 1, MPI_INT, 0, dist::comm);
    assert (rank>=1);
    MPI_Bcast (&dims, rank, MPI_INT, 0, dist::comm);
    for (int d=0; d<rank; ++d) assert (dims[d]>=0);
    MPI_Bcast (&nbytes, 1, MPI_INT, 0, dist::comm);
    assert (nbytes>=0);

    MPI_Bcast (&rl, 1, MPI_INT, 0, dist::comm);
    MPI_Bcast (&tl, 1, MPI_INT, 0, dist::comm);
    MPI_Bcast (&mglevel, 1, MPI_INT, 0, dist::comm);
    MPI_Bcast (&comp, 1, MPI_INT, 0, dist::comm);


    int gpdim = CCTK_GroupDimI(group);
    const int grouptype = CCTK_GroupTypeI(group);

    cGroup cgdata;
    int ierr = CCTK_GroupData(group,&cgdata);
    assert(ierr==0);

    // Read grid
    AmrGrid* amrgrid = 0;
    int amr_origin[dim];
    int amr_dims[dim];
    if (myproc==0) {
        
        // Read data
        if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Reading AMR data");

        amrgrid = amrreader->getGrid(currdataset);
        assert (amrgrid!=0);
        assert (amrgrid->data!=0);

	IObase::DataType atype;
        int alength;
        // If iorigin attribute is absent, assume file has unigrid
        // data.  Initialize iorigin to 0.
        if (reader->readAttributeInfo("iorigin", atype, alength) < 0) {
          for (int d=0; d<gpdim; ++d) {
            amrgrid->iorigin[d] = 0;
          }
        }
        for (int d=0; d<gpdim; ++d) {
          amr_origin[d] = amrgrid->iorigin[d];
	  //	  fprintf(stderr,"\n amr_origin[%d]=%d",d,amr_origin[d]);
          amr_dims[d] = amrgrid->dims[d];
	  //fprintf(stderr,"\n amr_dims[%d]=%d",d,amr_dims[d]);
        }
        for (int d=gpdim; d<dim; ++d) {
          amr_origin[d] = 0;
          amr_dims[d] = 1;
        }
        
    } // MyProc == 0


    MPI_Bcast (amr_origin, dim, MPI_INT, 0, dist::comm);
    MPI_Bcast (amr_dims, dim, MPI_INT, 0, dist::comm);

    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = varindex - n0;
    assert (var>=0 && var<CCTK_NumVars());


    // Traverse all components on this refinement and multigrid
    // level


    //    fprintf(stderr,"\n bogus! reflevel:%d mglevel:%d\n",reflevel,mglevel);        
    //fprintf(stderr,"blahblah: rank: %d dims[0,1,2]: %d,%d,%d\n",rank,dims[0],dims[1],dims[2]);

    //    cout << "var " << varindex << " has " << CCTK_NumTimeLevelsFromVarI(varindex) << " timelevels" << endl;

    //    BEGIN_COMPONENT_LOOP(cgh, grouptype) {

    //cout << "compontents " << hh->components(rl) << endl;

    //cout << "myproc: " <<  CCTK_MyProc(cgh) << endl;
    //    fprintf(stderr,"%d amr_dims: %d,%d,%d\n",CCTK_MyProc(cgh),amr_dims[0],amr_dims[1],amr_dims[2]);
    //fprintf(stderr,"%d amr_origin: %d,%d,%d\n",CCTK_MyProc(cgh),amr_origin[0],amr_origin[1],amr_origin[2]);

    for(int c=0;c<hh->components(rl);c++) {

      ggf<dim>* ff = 0;

      assert (var < (int)arrdata[group].data.size());
      ff = (ggf<dim>*)arrdata[group].data[var];

      gdata<dim>* const data = (*ff) (tl, rl, c, mglevel);

        
      // Create temporary data storage on processor 0
      const int reflevelfact_local=ipow(reffact,rl);
      vect<int,dim> str = vect<int,dim>(maxreflevelfact/reflevelfact_local);

      if(grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY)
	str = vect<int,dim> (1);


      vect<int,dim> lb = vect<int,dim>(amr_origin) * str;
      vect<int,dim> ub = lb + (vect<int,dim>(amr_dims) - 1) * str;



      gdata<dim>* const tmp = data->make_typed (varindex);

        
      // Copy into grid function

      if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
	assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
	if (grouptype == CCTK_SCALAR) {
	  lb[0] = hh->processors.at(rl).at(c);
	  ub[0] = hh->processors.at(rl).at(c);
	} else {
	  lb[dim-1] = lb[dim-1] + (ub[dim-1]-lb[dim-1]+1)*hh->processors.at(rl).at(c);
	  ub[dim-1] = ub[dim-1] + (ub[dim-1]-lb[dim-1]+1)*hh->processors.at(rl).at(c);
	}
      }
      
      const bbox<int,dim> ext(lb,ub,str);


      if (myproc==0) {
	tmp->allocate (ext, 0, amrgrid->data);
      } else {
	tmp->allocate (ext, 0, 0);
      }

      for (comm_state<dim> state; !state.done(); state.step()) {
	data->copy_from (state, tmp, ext & data->extent() );
      }

       
      // Delete temporary copy
      delete tmp;
        
    } // manual component loop
      
    if (myproc==0) {
        free (amrgrid->data);
        free (amrgrid);
        amrgrid = 0;
    }
  

    return 0;
  }
  
  
  
  int InputGH (const cGH* const cgh) {
    int retval = 0;
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (CheckForVariable(cgh, GetStringParameter("in3D_vars",""), vindex)) {
	char* varname = CCTK_FullName(vindex);
	retval = InputVarAs (cgh, varname, CCTK_VarName(vindex));
	free (varname);
	if (retval != 0) return retval;
      }
    }
    return retval;
  }
  
  
  
  int InputVarAs (const cGH* const cgh, const char* const varname,
		  const char* const alias) {
    DECLARE_CCTK_PARAMETERS;
    
    const int n = CCTK_VarIndex(varname);
    assert (n>=0 && n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVars());
    const int tl = 0;      //   CCTK_VInfo (CCTK_THORNSTRING, "boguscheck reflevel,component,mglevel %d,%d,%d",reflevel,component,mglevel);

    
    // Check for storage
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot input variable \"%s\" because it has no storage",
                  varname);
      return 0;
    }
    
    const int grouptype = CCTK_GroupTypeI(group);
    const int rl = grouptype==CCTK_GF ? reflevel : 0;
    
    // Find the input directory
    const char* myindir = GetStringParameter("indir3D", "");
    
    // Invent a file name
    const char* extension = 0;
    if (CCTK_Equals(in3D_format, "IEEE")) {
      extension = ".raw";
#ifdef HDF4
    } else if (CCTK_Equals(in3D_format, "HDF4")) {
      extension = ".hdf";
#endif
#ifdef HDF5
    } else if (CCTK_Equals(in3D_format, "HDF5")) {
      extension = ".h5";
#endif
    } else {
      assert (0);
    }
    extension = GetStringParameter ("in3D_extension", extension);
    
    ostringstream filenamebuf;
    filenamebuf << myindir << "/" << alias << extension;
    string filenamestr = filenamebuf.str();
    const char * const filename = filenamestr.c_str();
    
    IObase* reader = 0;
    AmrGridReader* amrreader = 0;
    int ndatasets = -1;
    
    const int gpdim = CCTK_GroupDimI(group);
    
    int rank;
    int dims[dim];
    int nbytes;
    
    // Read the file only on the root processor
    if (CCTK_MyProc(cgh)==0) {
      
      // Open the file 
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Opening file \"%s\"", filename);
      if (CCTK_Equals(in3D_format, "IEEE")) {
        reader = new IEEEIO(filename, IObase::Read);
#ifdef HDF4
      } else if (CCTK_Equals(in3D_format, "HDF4")) {
        reader = new HDFIO(filename, IObase::Read);
#endif
#ifdef HDF5
      } else if (CCTK_Equals(in3D_format, "HDF5")) {
        reader = new H5IO(filename, IObase::Read);
#endif
      } else {
        assert (0);
      }
      if (!reader->isValid()) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not open file \"%s\" for reading", filename);
      }
      assert (reader->isValid());
      
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Reading AMR info");
      amrreader = new AmrGridReader(*reader);
      
      // Read information about dataset
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Reading dataset info");
      IObase::DataType numbertype;
      reader->readInfo (numbertype, rank, dims);
      nbytes = IObase::nBytes(numbertype,rank,dims);
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "type=%d rank=%d dims=[%d,%d,%d] nbytes=%d", (int)numbertype, rank, dims[0], dims[1], dims[2], nbytes);
      
      // Check rank
      assert (rank==gpdim);
      
      // Check datatype
      // TODO: Check datatype correctly
      assert (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL8
              || (sizeof(CCTK_REAL) == sizeof(CCTK_REAL8)
                  && CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL));
      
      // TODO: check grid spacing
      
      // Number of datasets
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Reading number of datasets");
      ndatasets = reader->nDatasets();
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "ndatasets=%d", ndatasets);
      assert (ndatasets>=0);
    }
    
    // Broadcast rank, dimensions, and nbytes
    MPI_Bcast (&rank, 1, MPI_INT, 0, dist::comm);
    assert (rank>=1);
    MPI_Bcast (&dims, rank, MPI_INT, 0, dist::comm);
    for (int d=0; d<rank; ++d) assert (dims[d]>=0);
    MPI_Bcast (&nbytes, 1, MPI_INT, 0, dist::comm);
    assert (nbytes>=0);
    
    // Broadcast number of datasets
    MPI_Bcast (&ndatasets, 1, MPI_INT, 0, dist::comm);
    assert (ndatasets>=0);
    
    // Read all datasets
    // TODO: read only some datasets
    for (int dataset=0; dataset<ndatasets; ++dataset) {
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Handling dataset #%d", dataset);
      
      // Read grid
      AmrGrid* amrgrid = 0;
      int amr_origin[dim];
      int amr_dims[dim];
      
      if (CCTK_MyProc(cgh)==0) {
        
        // Read data
        if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Reading AMR data");
        amrgrid = amrreader->getGrid(dataset);
        assert (amrgrid!=0);
        assert (amrgrid->data!=0);
        
        // If iorigin attribute is absent, assume file has unigrid
        // data.  Initialize iorigin to 0.
        IObase::DataType atype;
        int alength;
        if (reader->readAttributeInfo("iorigin", atype, alength) < 0) {
          for (int d=0; d<gpdim; ++d) {
            amrgrid->iorigin[d] = 0;
          }
        }
        
        for (int d=0; d<gpdim; ++d) {
          amr_origin[d] = amrgrid->iorigin[d];
          amr_dims[d] = amrgrid->dims[d];
        }
        for (int d=gpdim; d<dim; ++d) {
          amr_origin[d] = 0;
          amr_dims[d] = 1;
        }
        
      } // MyProc == 0
      MPI_Bcast (amr_origin, dim, MPI_INT, 0, dist::comm);
      MPI_Bcast (amr_dims, dim, MPI_INT, 0, dist::comm);
      
      // Traverse all components on this refinement and multigrid
      // level
      BEGIN_COMPONENT_LOOP(cgh, grouptype) {
        
        ggf<dim>* ff = 0;
        
        assert (var < (int)arrdata[group].data.size());
        ff = (ggf<dim>*)arrdata[group].data[var];
        
        gdata<dim>* const data = (*ff) (tl, rl, component, mglevel);
        
        // Create temporary data storage on processor 0
        const vect<int,dim> str = vect<int,dim>(reflevelfact);
        const vect<int,dim> lb = vect<int,dim>(amr_origin) * str;
        const vect<int,dim> ub
          = lb + (vect<int,dim>(amr_dims) - 1) * str;
        const bbox<int,dim> ext(lb,ub,str);
        gdata<dim>* const tmp = data->make_typed (n);
        
        if (CCTK_MyProc(cgh)==0) {
          tmp->allocate (ext, 0, amrgrid->data);
        } else {
          tmp->allocate (ext, 0, 0);
        }
        
        // Copy into grid function
	for (comm_state<dim> state; !state.done(); state.step()) {
           data->copy_from (state, tmp, ext);
         }
        
        // Delete temporary copy
        delete tmp;
        
      } END_COMPONENT_LOOP;
      
      if (CCTK_MyProc(cgh)==0) {
        free (amrgrid->data);
        free (amrgrid);
        amrgrid = 0;
      }
      
    } // loop over datasets
    
    // Close the file
    if (CCTK_MyProc(cgh)==0) {
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Deleting AMR info");
      delete amrreader;
      amrreader = 0;
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Closing file");
      delete reader;
      reader = 0;
    }
    
    return 0;
  }
  
  
  
  int CarpetIOFlexIO_ReadData (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    return InputGH(cctkGH);
  }
  
  
  
  const char* GetStringParameter (const char* const parametername,
				  const char* const fallback)
  {
    if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
      int ptype;
      const char* const* const ppval = (const char* const*)CCTK_ParameterGet
	(parametername, CCTK_THORNSTRING, &ptype);
      assert (ppval);
      const char* const pval = *ppval;
      assert (ptype == PARAMETER_STRING);
      return pval;
    }
    
    return fallback;
  }
  
  
  
  int GetIntParameter (const char* const parametername, int fallback)
  {
    if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
      int ptype;
      const int* const ppval = (const int*)CCTK_ParameterGet
	(parametername, CCTK_THORNSTRING, &ptype);
      assert (ppval);
      const int pval = *ppval;
      assert (ptype == PARAMETER_INT);
      return pval;
    }
    
    return fallback;
  }
  
  
  
  bool CheckForVariable (const cGH* const cgh,
			 const char* const varlist, const int vindex)
  {
    const int numvars = CCTK_NumVars();
    assert (vindex>=0 && vindex<numvars);
    
    vector<bool> flags(numvars);
    
    CCTK_TraverseString (varlist, SetFlag, &flags, CCTK_GROUP_OR_VAR);
    
    return flags[vindex];
  }
  
  void SetFlag (int index, const char* optstring, void* arg)
  {
    vector<bool>& flags = *(vector<bool>*)arg;
    flags[index] = true;
  }
  

  
} // namespace CarpetIOFlexIO
