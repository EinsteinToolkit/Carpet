#include <alloca.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <vector>

#include <AMRwriter.hh>
#include <AmrGridReader.hh>
#include <H5IO.hh>
#include <HDFIO.hh>
#include <IEEEIO.hh>
#include <IO.hh>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/data.hh"
#include "Carpet/CarpetLib/src/gdata.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "ioflexio.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIO/src/ioflexio.cc,v 1.10 2001/07/04 12:29:50 schnetter Exp $";



namespace CarpetIOFlexIO {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Variable definitions
  int GHExtension;
  int IOMethod;
  vector<bool> do_truncate;
  vector<vector<int> > last_output;
  
  
  
  static const char* GetStringParameter (const char* const parametername,
					 const char* const fallback);
  static int GetIntParameter (const char* const parametername, int fallback);
  static bool CheckForVariable (cGH* const cgh,
				const char* const varlist, const int vindex);
  static void SetFlag (int index, const char* optstring, void* arg);
  
  
  
  int CarpetIOFlexIOStartup ()
  {
    CCTK_RegisterBanner ("AMR 3D FlexIO I/O provided by CarpetIOFlexIO");
    
    GHExtension = CCTK_RegisterGHExtension("CarpetIOFlexIO");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    IOMethod = CCTK_RegisterIOMethod ("IOFlexIO");
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
    
    return 0;
  }
  
  
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
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
    return 0;
  }
  
  
  
  int OutputGH (cGH* const cgh) {
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cgh, vindex)) {
	TriggerOutput(cgh, vindex);
      }
    }
    return 0;
  }
  
  
  
  int OutputVarAs (cGH* const cgh, const char* const varname,
		   const char* const alias) {
    DECLARE_CCTK_PARAMETERS;
    
    const int n = CCTK_VarIndex(varname);
    assert (n>=0 && n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::gfdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVars());
    const int tl = 0;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR: {
      // Don't output scalars
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannout output variable \"%s\" because it is a scalar",
		  varname);
      return 0;
    }
    
    case CCTK_ARRAY:
    case CCTK_GF: {
      
      // Check for storage
      if (! CCTK_QueryGroupStorageI(cgh, group)) {
	CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "Cannot output variable \"%s\" because it has no storage",
		    varname);
	return 0;
      }
      
      // Create the output directory
      const char* myoutdir = GetStringParameter("outdir3D", outdir);
      if (CCTK_MyProc(cgh)==0) {
	CCTK_CreateDirectory (0755, myoutdir);
      }
      
      // Invent a file name
      const char* extension = 0;
      if (CCTK_Equals(out3D_format, "IEEE")) {
	extension = ".raw";
#ifdef HDF5
      } else if (CCTK_Equals(out3D_format, "HDF4")) {
	extension = ".hdf";
      } else if (CCTK_Equals(out3D_format, "HDF5")) {
	extension = ".h5";
#endif
      } else {
	abort();
      }
      extension = GetStringParameter ("out3D_extension", extension);
      
      char* const filename = (char*)alloca(strlen(myoutdir)
					   +strlen(alias)
					   +strlen(extension)+100);
      sprintf (filename, "%s/%s%s", myoutdir, alias, extension);
      
      IObase* writer = 0;
      AMRwriter* amrwriter = 0;
      
      // Write the file only on the root processor
      if (CCTK_MyProc(cgh)==0) {
	
	// If this is the first time, then create and truncate the
	// file
	if (do_truncate[n]) {
	  struct stat fileinfo;
	  if (! IOUtil_RestartFromRecovery(cgh)
	      || stat(filename, &fileinfo)!=0) {
	    writer = 0;
	    if (CCTK_Equals(out3D_format, "IEEE")) {
	      writer = new IEEEIO(filename, IObase::Create);
#ifdef HDF5
	    } else if (CCTK_Equals(out3D_format, "HDF4")) {
	      writer = new HDFIO(filename, IObase::Create);
	    } else if (CCTK_Equals(out3D_format, "HDF5")) {
	      writer = new H5IO(filename, IObase::Create);
#endif
	    } else {
	      abort();
	    }
	    delete writer;
	    writer = 0;
	  }
	}
	
	// Open the file 
	if (CCTK_Equals(out3D_format, "IEEE")) {
	  writer = new IEEEIO(filename, IObase::Append);
#ifdef HDF5
	} else if (CCTK_Equals(out3D_format, "HDF4")) {
	  writer = new HDFIO(filename, IObase::Append);
	} else if (CCTK_Equals(out3D_format, "HDF5")) {
	  writer = new H5IO(filename, IObase::Append);
#endif
	} else {
	  abort();
	}
	assert (writer->isValid());
	amrwriter = new AMRwriter(*writer);
	
	// Set datatype
	assert (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL8
		|| (sizeof(CCTK_REAL) == sizeof(CCTK_REAL8)
		    && CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL));
	// TODO: Set datatype correctly
	amrwriter->setType (IObase::Float64);
	
	const int gpdim = CCTK_GroupDimI(group);
	
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
	timestep = base_delta_time;
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
      
      // Traverse all components on this refinement and multigrid
      // level
      BEGIN_COMPONENT_LOOP(cgh) {
	
	const generic_gf<dim>* ff = 0;
	
	switch (CCTK_GroupTypeI(group)) {
	  
	case CCTK_ARRAY:
	  assert (var < (int)arrdata[group].data.size());
	  ff = (generic_gf<dim>*)arrdata[group].data[var];
	  break;
	  
	case CCTK_GF:
	  assert (var < (int)gfdata[group].data.size());
	  ff = (generic_gf<dim>*)gfdata[group].data[var];
	  break;
	  
	default:
	  abort();
	}
	
	const generic_data<dim>* const data
	  = (*ff) (tl, reflevel, component, mglevel);
	
	/* Make temporary copy on processor 0 */
	const bbox<int,dim> ext = data->extent();
	generic_data<dim>* const tmp = data->make_typed ();
	tmp->allocate (ext, 0);
	tmp->copy_from (data, ext);
	
	/* Write data */
	if (CCTK_MyProc(cgh)==0) {
	  int origin[dim], dims[dim];
	  for (int d=0; d<dim; ++d) {
	    origin[d] = (ext.lower() / ext.stride())[d];
	    dims[d]   = (ext.shape() / ext.stride())[d];
	  }
	  amrwriter->write (origin, dims, (void*)tmp->storage());
	}
	
	/* Delete temporary copy */
	delete tmp;
	
      } END_COMPONENT_LOOP(cgh);
      
      // Close the file
      if (CCTK_MyProc(cgh)==0) {
	delete amrwriter;
	amrwriter = 0;
	delete writer;
	writer = 0;
      }
      
      break;
    } // ARRAY or GROUP
      
    default:
      abort();
    }
    
    // Don't truncate again
    do_truncate[n] = false;
    
    return 0;
  }
  
  
  
  int TimeToOutput (cGH* const cgh, const int vindex) {
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
  
  
  
  int TriggerOutput (cGH* const cgh, const int vindex) {
    assert (vindex>=0 && vindex<CCTK_NumVars());
    
    char* varname = CCTK_FullName(vindex);
    const int retval = OutputVarAs (cgh, varname, CCTK_VarName(vindex));
    free (varname);
    
    last_output[reflevel][vindex] = cgh->cctk_iteration;
    
    return retval;
  }
  
  
  
  int InputGH (cGH* const cgh) {
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
  
  
  
  int InputVarAs (cGH* const cgh, const char* const varname,
		   const char* const alias) {
    DECLARE_CCTK_PARAMETERS;
    
    const int n = CCTK_VarIndex(varname);
    assert (n>=0 && n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::gfdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVars());
    const int tl = 0;
    
    switch (CCTK_GroupTypeI(group)) {
      
    case CCTK_SCALAR: {
      // Don't input scalars
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannout input variable \"%s\" because it is a scalar",
		  varname);
      return 0;
    }
    
    case CCTK_ARRAY:
    case CCTK_GF: {
      
      // Check for storage
      if (! CCTK_QueryGroupStorageI(cgh, group)) {
	CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "Cannot input variable \"%s\" because it has no storage",
		    varname);
	return 0;
      }
      
      // Find the input directory
      const char* myindir = GetStringParameter("indir3D", "");
      
      // Invent a file name
      const char* extension = 0;
      if (CCTK_Equals(in3D_format, "IEEE")) {
	extension = ".raw";
#ifdef HDF5
      } else if (CCTK_Equals(in3D_format, "HDF4")) {
	extension = ".hdf";
      } else if (CCTK_Equals(in3D_format, "HDF5")) {
	extension = ".h5";
#endif
      } else {
	abort();
      }
      extension = GetStringParameter ("in3D_extension", extension);
      
      char* const filename = (char*)alloca(strlen(myindir)
					   +strlen(alias)
					   +strlen(extension)+100);
      sprintf (filename, "%s/%s.%s", myindir, alias, extension);
      
      IObase* reader = 0;
      AmrGridReader* amrreader = 0;
      int ndatasets = -1;
      
      const int gpdim = CCTK_GroupDimI(group);
      
      // Read the file only on the root processor
      if (CCTK_MyProc(cgh)==0) {
	
	// Open the file 
	if (CCTK_Equals(in3D_format, "IEEE")) {
	  reader = new IEEEIO(filename, IObase::Read);
#ifdef HDF5
	} else if (CCTK_Equals(in3D_format, "HDF4")) {
	  reader = new HDFIO(filename, IObase::Read);
	} else if (CCTK_Equals(in3D_format, "HDF5")) {
	  reader = new H5IO(filename, IObase::Read);
#endif
	} else {
	  abort();
	}
	assert (reader->isValid());
	amrreader = new AmrGridReader(*reader);
	
	// Read information about dataset
	IObase::DataType numbertype;
	int rank;
	int dims[dim];
	reader->readInfo (numbertype, rank, dims);
	
	// Check rank
	assert (rank==gpdim);
	
	// Check datatype
	// TODO: Check datatype correctly
	assert (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL8
		|| (sizeof(CCTK_REAL) == sizeof(CCTK_REAL8)
		    && CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL));
	
	// TODO: check grid spacing
	
	// Number of datasets
	ndatasets = reader->nDatasets();
	assert (ndatasets>=0);
      }
      
      // Broadcast number of datasets
      MPI_Bcast (&ndatasets, 1, MPI_INT, 0, dist::comm);
      assert (ndatasets>=0);
      
      // Read all datasets
      for (int dataset=0; dataset<ndatasets; ++dataset) {
	
	// Read grid
	AmrGrid* amrgrid = 0;
	int amr_origin[dim];
	int amr_dims[dim];
	
	if (CCTK_MyProc(cgh)==0) {
	  
	  // Read data
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
	}
	MPI_Bcast (amr_origin, dim, MPI_INT, 0, dist::comm);
	MPI_Bcast (amr_dims, dim, MPI_INT, 0, dist::comm);
	
	// Traverse all components on this refinement and multigrid
	// level
	BEGIN_COMPONENT_LOOP(cgh) {
	  
	  generic_gf<dim>* ff = 0;
	  
	  switch (CCTK_GroupTypeI(group)) {
	    
	  case CCTK_ARRAY:
	    assert (var < (int)arrdata[group].data.size());
	    ff = (generic_gf<dim>*)arrdata[group].data[var];
	    break;
	    
	  case CCTK_GF:
	    assert (var < (int)gfdata[group].data.size());
	    ff = (generic_gf<dim>*)gfdata[group].data[var];
	    break;
	    
	  default:
	    abort();
	  }
	  
	  generic_data<dim>* const data
	    = (*ff) (tl, reflevel, component, mglevel);
	  
	  /* Create temporary data storage on processor 0 */
	  const vect<int,dim> str = vect<int,dim>(reflevelfact);
	  const vect<int,dim> lb = vect<int,dim>(amr_origin) * str;
	  const vect<int,dim> ub
	    = lb + (vect<int,dim>(amr_dims) - 1) * str;
	  const bbox<int,dim> ext(lb,ub,str);
	  generic_data<dim>* const tmp = data->make_typed ();
	  
	  if (CCTK_MyProc(cgh)==0) {
	    tmp->allocate (ext, 0, amrgrid->data);
	  } else {
	    tmp->allocate (ext, 0, 0);
	  }
	  
	  /* Copy into grid function */
	  data->copy_from (tmp, ext);
	  
	  /* Delete temporary copy */
	  delete tmp;
	  
	} END_COMPONENT_LOOP(cgh);
	
	if (CCTK_MyProc(cgh)==0) {
	  free (amrgrid->data);
	  free (amrgrid);
	}
	
      }	// loop over datasets
      
      // Close the file
      if (CCTK_MyProc(cgh)==0) {
	delete amrreader;
	amrreader = 0;
	delete reader;
	reader = 0;
      }
      
      break;
    } // ARRAY or GROUP
      
    default:
      abort();
    }
    
    return 0;
  }
  
  
  
  int CarpetIOFlexIOReadData (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    return InputGH(cctkGH);
  }
  
  
  
  const char* GetStringParameter (const char* const parametername,
				  const char* const fallback)
  {
    if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
      int ptype;
      const char** const ppval = (const char**)CCTK_ParameterGet
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
      const int* const ppval
	= (int*)CCTK_ParameterGet (parametername, CCTK_THORNSTRING, &ptype);
      assert (ppval);
      const int pval = *ppval;
      assert (ptype == PARAMETER_INT);
      return pval;
    }
    
    return fallback;
  }
  
  
  
  bool CheckForVariable (cGH* const cgh,
			 const char* const varlist, const int vindex)
  {
    const int numvars = CCTK_NumVars();
    assert (vindex>=0 && vindex<numvars);
    
    bool* const flags = (bool*)alloca(numvars * sizeof(bool));
    
    for (int i=0; i<numvars; ++i) {
      flags[i] = false;
    }
    
    CCTK_TraverseString (varlist, SetFlag, flags, CCTK_GROUP_OR_VAR);
    
    return flags[vindex];
  }
  
  void SetFlag (int index, const char* optstring, void* arg)
  {
    ((bool*)arg)[index] = true;
  }
  
  
  
} // namespace CarpetIOFlexIO
