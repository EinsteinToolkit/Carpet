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

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIO/src/ioflexio.cc,v 1.40 2004/02/07 16:21:56 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetIOFlexIO_ioflexio_cc);
}

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
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "ioflexio.hh"



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
  static bool CheckForVariable (const cGH* const cgh,
				const char* const varlist, const int vindex);
  static void SetFlag (int index, const char* optstring, void* arg);
  
  static void WriteAttribute (IObase* writer, const char* name,
                              int value);
  static void WriteAttribute (IObase* writer, const char* name,
                              const int* values, int nvalues);
  static void WriteAttribute (IObase* writer, const char* name,
                              CCTK_REAL value);
  static void WriteAttribute (IObase* writer, const char* name,
                              const CCTK_REAL* values, int nvalues);
  static void WriteAttribute (IObase* writer, const char* name,
                              char value);
  static void WriteAttribute (IObase* writer, const char* name,
                              const char* values);
  static void WriteAttribute (IObase* writer, const char* name,
                              const char* values, int nvalues);

  
  
  int CarpetIOFlexIOStartup ()
  {
    int ierr;
    
    CCTK_RegisterBanner ("AMR 3D FlexIO I/O provided by CarpetIOFlexIO");
    
    GHExtension = CCTK_RegisterGHExtension ("CarpetIOFlexIO");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    IOMethod = CCTK_RegisterIOMethod ("CarpetIOFlexIO");
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
    
    ierr = IOUtil_RegisterRecover ("CarpetIOFlexIO", Recover);
    assert (! ierr);
    
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
  
  
  
  int OutputGH (const cGH* const cgh) {
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cgh, vindex)) {
	TriggerOutput(cgh, vindex);
      }
    }
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
    const int rl = grouptype==CCTK_GF ? reflevel : 0;
    
    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cgh, "IO");
    assert (iogh);
    
    // Create the output directory
    const char* const myoutdir = GetStringParameter("out3D_dir", out_dir);
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
      timestep = cgh->cctk_delta_time;
      amrwriter->setTopLevelParameters
	(gpdim, origin, delta, timestep, maxreflevels);
      
      // Set refinement information
      int interlevel_timerefinement;
      int interlevel_spacerefinement[dim];
      int initial_gridplacementrefinement[dim];
      interlevel_timerefinement = vhh[0]->reffact;
      for (int d=0; d<dim; ++d) {
	interlevel_spacerefinement[d] = vhh[0]->reffact;
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
    
    // Traverse all components
    BEGIN_MAP_LOOP(cgh, grouptype) {
      BEGIN_COMPONENT_LOOP(cgh, grouptype) {
        
        const ggf<dim>* ff = 0;
        
        assert (var < (int)arrdata[group][Carpet::map].data.size());
        ff = (ggf<dim>*)arrdata[group][Carpet::map].data[var];
        
        const gdata<dim>* const data
          = (*ff) (tl, rl, component, mglevel);
        
        // Make temporary copy on processor 0
        bbox<int,dim> ext = data->extent();
        vect<int,dim> lo = ext.lower();
        vect<int,dim> hi = ext.upper();
        vect<int,dim> str = ext.stride();
        
        // Ignore ghost zones if desired
        for (int d=0; d<dim; ++d) {
          const int max_lower_ghosts = (cgh->cctk_bbox[2*d  ] && out3D_output_outer_boundary) ? -1 : out3D_max_num_lower_ghosts;
          const int max_upper_ghosts = (cgh->cctk_bbox[2*d+1] && out3D_output_outer_boundary) ? -1 : out3D_max_num_upper_ghosts;
          
          const int num_lower_ghosts = max_lower_ghosts == -1 ? cgh->cctk_nghostzones[d] : min(out3D_max_num_lower_ghosts, cgh->cctk_nghostzones[d]);
          const int num_upper_ghosts = max_upper_ghosts == -1 ? cgh->cctk_nghostzones[d] : min(out3D_max_num_upper_ghosts, cgh->cctk_nghostzones[d]);
          
          lo[d] += (cgh->cctk_nghostzones[d] - num_lower_ghosts) * str[d];
          hi[d] -= (cgh->cctk_nghostzones[d] - num_upper_ghosts) * str[d];
        }
        
        ext = bbox<int,dim>(lo,hi,str);
        
        gdata<dim>* const tmp = data->make_typed (n);
        tmp->allocate (ext, 0);
        for (comm_state<dim> state; !state.done(); state.step()) {
          tmp->copy_from (state, data, ext);
        }
        
        // Write data
        if (CCTK_MyProc(cgh)==0) {
          int origin[dim], dims[dim];
          for (int d=0; d<dim; ++d) {
            origin[d] = (ext.lower() / ext.stride())[d];
            dims[d]   = (ext.shape() / ext.stride())[d];
          }
          amrwriter->write (origin, dims, (void*)tmp->storage());
          
          // Write some additional attributes
          WriteAttribute (writer, "group_version", 1);
          {
            char * fullname = CCTK_FullName(n);
            assert (fullname);
            WriteAttribute (writer, "group_fullname", fullname);
            free (fullname);
          }
          WriteAttribute (writer, "group_varname", CCTK_VarName(n));
          {
            char * groupname = CCTK_GroupName(group);
            assert (groupname);
            WriteAttribute (writer, "group_groupname", groupname);
            free (groupname);
          }
          WriteAttribute (writer, "group_grouptype", grouptype);
          WriteAttribute (writer, "group_dim", CCTK_GroupDimI(group));
          WriteAttribute (writer, "group_timelevel", tl);
          WriteAttribute (writer, "group_numtimelevels", CCTK_NumTimeLevelsI(group));
          
          WriteAttribute (writer, "cctk_version", 1);
          WriteAttribute (writer, "cctk_dim", cgh->cctk_dim);
          WriteAttribute (writer, "cctk_iteration", cgh->cctk_iteration);
// TODO: disable temporarily
//           WriteAttribute (writer, "cctk_nmaps", cgh->cctk_nmaps);
//           WriteAttribute (writer, "cctk_map", cgh->cctk_map);
          WriteAttribute (writer, "cctk_gsh", cgh->cctk_gsh, dim);
          WriteAttribute (writer, "cctk_lsh", cgh->cctk_lsh, dim);
          WriteAttribute (writer, "cctk_lbnd", cgh->cctk_lbnd, dim);
          WriteAttribute (writer, "cctk_delta_time", cgh->cctk_delta_time);
          WriteAttribute (writer, "cctk_delta_space", cgh->cctk_delta_space, dim);
          WriteAttribute (writer, "cctk_origin_space", cgh->cctk_origin_space, dim);
          WriteAttribute (writer, "cctk_bbox", cgh->cctk_bbox, 2*dim);
          WriteAttribute (writer, "cctk_levfac", cgh->cctk_levfac, dim);
          WriteAttribute (writer, "cctk_levoff", cgh->cctk_levoff, dim);
          WriteAttribute (writer, "cctk_levoffdenom", cgh->cctk_levoffdenom, dim);
          WriteAttribute (writer, "cctk_timefac", cgh->cctk_timefac);
          WriteAttribute (writer, "cctk_convlevel", cgh->cctk_convlevel);
          WriteAttribute (writer, "cctk_convfac", cgh->cctk_convfac);
          WriteAttribute (writer, "cctk_nghostzones", cgh->cctk_nghostzones, dim);
          WriteAttribute (writer, "cctk_time", cgh->cctk_time);
          
          WriteAttribute (writer, "carpet_version", 1);
          WriteAttribute (writer, "carpet_dim", dim);
          WriteAttribute (writer, "carpet_basemglevel", basemglevel);
          WriteAttribute (writer, "carpet_mglevel", mglevel);
          WriteAttribute (writer, "carpet_mglevels", mglevels);
          WriteAttribute (writer, "carpet_mgface", mgfact);
          WriteAttribute (writer, "carpet_reflevel", reflevel);
          WriteAttribute (writer, "carpet_reflevels", reflevels);
          WriteAttribute (writer, "carpet_reffact", reffact);
          WriteAttribute (writer, "carpet_map", Carpet::map);
          WriteAttribute (writer, "carpet_maps", maps);
          WriteAttribute (writer, "carpet_component", component);
          WriteAttribute (writer, "carpet_components", vhh[Carpet::map]->components(reflevel));
        }
        
        // Delete temporary copy
        delete tmp;
        
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
    
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
    const int tl = 0;
    
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
    const char* myindir = GetStringParameter("in3D_dir", ".");
    
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
      
      // Traverse all components on all levels
      BEGIN_MAP_LOOP(cgh, grouptype) {
        BEGIN_COMPONENT_LOOP(cgh, grouptype) {
          
          ggf<dim>* ff = 0;
          
          assert (var < (int)arrdata[group][Carpet::map].data.size());
          ff = (ggf<dim>*)arrdata[group][Carpet::map].data[var];
          
          gdata<dim>* const data = (*ff) (tl, rl, component, mglevel);
          
          // Create temporary data storage on processor 0
          const vect<int,dim> str
            = vect<int,dim>(maxreflevelfact/reflevelfact);
          const vect<int,dim> lb = vect<int,dim>::ref(amr_origin) * str;
          const vect<int,dim> ub
            = lb + (vect<int,dim>::ref(amr_dims) - 1) * str;
          const bbox<int,dim> ext(lb,ub,str);
          
          if (data->extent() != ext) {
            CCTK_WARN (0, "The stored data have a different extent than the grid function in memory.  This is not yet supported.");
          }
          
          gdata<dim>* const tmp = data->make_typed (n);
          
          if (CCTK_MyProc(cgh)==0) {
            tmp->allocate (ext, 0, amrgrid->data);
          } else {
            tmp->allocate (ext, 0);
          }
          
          // Copy into grid function
          for (comm_state<dim> state; !state.done(); state.step()) {
            data->copy_from (state, tmp, ext);
          }
          
          // Delete temporary copy
          delete tmp;
          
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
      
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
  
  
  
  /** returns the number of recovered variables */
  int Recover (cGH* const cgh, const char *basefilename,
               const int called_from)
  {
    assert (cgh);
    assert (basefilename);
    assert (called_from == CP_INITIAL_DATA
            || called_from == CP_EVOLUTION_DATA
            || called_from == CP_RECOVER_PARAMETERS
            || called_from == CP_RECOVER_DATA
            || called_from == FILEREADER_DATA);
    
    // the other modes are not supported yet
    assert (called_from == FILEREADER_DATA);
    
    const ioGH * const iogh = (const ioGH *) CCTK_GHExtension (cgh, "IO");
    assert (iogh);
    
    int num_vars_read = 0;
    assert (iogh->do_inVars);
    for (int n=0; n<CCTK_NumVars(); ++n) {
      if (iogh->do_inVars[n]) {
        const char * p = strrchr (basefilename, '/');
        if (p) {
          ++ p;
        } else {
          p = basefilename;
        }
        if (strcmp (p, CCTK_VarName(n)) == 0) {
          char * const fullname = CCTK_FullName(n);
          assert (fullname);
          cout << "Recover: Reading variable \"" << fullname << "\" from file \"" << basefilename << "\"" << endl;
          const int ierr = InputVarAs (cgh, fullname, basefilename);
          assert (! ierr);
          free (fullname);
          ++ num_vars_read;
        }
      }
    }
    
    return num_vars_read;
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
  
  
  
  void WriteAttribute (IObase* writer, const char* name, int value)
  {
    WriteAttribute (writer, name, &value, 1);
  }
  
  void WriteAttribute (IObase* writer, const char* name,
                       const int* values, int nvalues)
  {
    assert (writer);
    assert (name);
    assert (values);
    vector<CCTK_INT4> values1(nvalues);
    for (int i=0; i<nvalues; ++i) {
      values1[i] = values[i];
    }
    writer->writeAttribute (name, IObase::Int32, nvalues, &values1[0]);
  }
  
  void WriteAttribute (IObase* writer, const char* name, CCTK_REAL value)
  {
    WriteAttribute (writer, name, &value, 1);
  }
  
  void WriteAttribute (IObase* writer, const char* name,
                       const CCTK_REAL* values, int nvalues)
  {
    assert (writer);
    assert (name);
    assert (values);
    vector<CCTK_REAL8> values1(nvalues);
    for (int i=0; i<nvalues; ++i) {
      values1[i] = values[i];
    }
    writer->writeAttribute (name, IObase::Float64, nvalues, &values1[0]);
  }
  
  void WriteAttribute (IObase* writer, const char* name, char value)
  {
    WriteAttribute (writer, name, &value, 1);
  }
  
  void WriteAttribute (IObase* writer, const char* name, const char * values)
  {
    WriteAttribute (writer, name, values, strlen(values));
  }
  
  void WriteAttribute (IObase* writer, const char* name,
                       const char * values, int nvalues)
  {
    assert (writer);
    assert (name);
    assert (values);
    vector<char> values1(nvalues);
    for (int i=0; i<nvalues; ++i) {
      values1[i] = values[i];
    }
    writer->writeAttribute (name, IObase::Char8, nvalues, &values1[0]);
  }
  
  
  
} // namespace CarpetIOFlexIO
