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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5.cc,v 1.30 2004/05/21 18:11:34 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetIOHDF5_iohdf5_cc);
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
  
  
  // Variable definitions
  int GHExtension;
  int IOMethod;
  vector<bool> do_truncate;     // [var]
  vector<vector<vector<int> > > last_output; // [ml][rl][var]
  
  
  int CarpetIOHDF5Startup ()
  {
    DECLARE_CCTK_PARAMETERS;

    int ierr;
    
    CCTK_RegisterBanner ("AMR 3D HDF5 I/O provided by CarpetIOHDF5");
    
    GHExtension = CCTK_RegisterGHExtension ("CarpetIOHDF5");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    IOMethod = CCTK_RegisterIOMethod ("CarpetIOHDF5");
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
    


    // Christian's Recovery routine
    if ( !(CCTK_Equals(recover,"no")) ) {
      ierr = IOUtil_RegisterRecover ("CarpetIOHDF5 recovery", CarpetIOHDF5_Recover);
      assert (! ierr);
    } else {
      // Erik's Recovery routine
      ierr = IOUtil_RegisterRecover ("CarpetIOHDF5", Recover);
      assert (! ierr);
    }


    return 0;
  }
  
  
  
  void CarpetIOHDF5Init (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    *this_iteration = -1;
    *next_output_iteration = 0;
    *next_output_time = cctk_time;
  }
  
  
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    CarpetIOHDF5GH* myGH;

    // Truncate all files if this is not a restart
    do_truncate.resize(CCTK_NumVars(), true);
    
    // No iterations have yet been output
    last_output.resize(mglevels);
    for (int ml=0; ml<mglevels; ++ml) {
      last_output.at(ml).resize(maxreflevels);
      for (int rl=0; rl<maxreflevels; ++rl) {
        last_output.at(ml).at(rl).resize(CCTK_NumVars(), INT_MIN);
      }
    }
    
    // We register only once, ergo we get only one handle.  We store
    // that statically, so there is no need to pass anything to
    // Cactus.

    /* allocate a new GH extension structure */

    CCTK_INT numvars = CCTK_NumVars ();

    myGH            = (CarpetIOHDF5GH*) malloc (sizeof (CarpetIOHDF5GH));
    myGH->out_last  = (int *) malloc (numvars * sizeof (int));
    myGH->requests  = (ioRequest **) calloc (numvars, sizeof (ioRequest *));
    myGH->cp_filename_list = (char **) calloc (abs (checkpoint_keep), sizeof (char *));
    myGH->cp_filename_index = 0;
    myGH->out_vars = strdup ("");
    myGH->out_every_default = out_every - 1;

    for (int i = 0; i < numvars; i++)
    {
      myGH->out_last [i] = -1;
    }

    myGH->open_output_files = NULL;

    // Now set hdf5 complex datatypes
    // Stolen from Thomas Radke
    HDF5_ERROR (myGH->HDF5_COMPLEX =
              H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX)));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX, "real",
                         offsetof (CCTK_COMPLEX, Re), HDF5_REAL));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX, "imag",
                         offsetof (CCTK_COMPLEX, Im), HDF5_REAL));
#ifdef CCTK_REAL4
    HDF5_ERROR (myGH->HDF5_COMPLEX8 =
		H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX8)));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX8, "real",
			   offsetof (CCTK_COMPLEX8, Re), HDF5_REAL4));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX8, "imag",
			   offsetof (CCTK_COMPLEX8, Im), HDF5_REAL4));
#endif
#ifdef CCTK_REAL8
    HDF5_ERROR (myGH->HDF5_COMPLEX16 =
		H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX16)));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX16, "real",
			   offsetof (CCTK_COMPLEX16, Re), HDF5_REAL8));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX16, "imag",
			   offsetof (CCTK_COMPLEX16, Im), HDF5_REAL8));
#endif
#ifdef CCTK_REAL16
    HDF5_ERROR (myGH->HDF5_COMPLEX32 =
		H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX32)));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX32, "real",
			   offsetof (CCTK_COMPLEX32, Re), HDF5_REAL16));
    HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX32, "imag",
			   offsetof (CCTK_COMPLEX32, Im), HDF5_REAL16));
#endif
    return (myGH);
  }
  
  int OutputGH (const cGH* const cctkGH) {
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cctkGH, vindex)) {
	TriggerOutput(cctkGH, vindex);
      }
    }
    return 0;
  }
  
  
  
  int OutputVarAs (const cGH* const cctkGH, const char* const varname,
		   const char* const alias) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    herr_t herr;
    
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
    if (! CCTK_QueryGroupStorageI(cctkGH, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot output variable \"%s\" because it has no storage",
		  varname);
      return 0;
    }
    
    const int grouptype = CCTK_GroupTypeI(group);
    switch (grouptype) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      assert (do_global_mode);
      break;
    case CCTK_GF:
      /* do nothing */
      break;
    default:
      assert (0);
    }
    const int rl = grouptype==CCTK_GF ? reflevel : 0;
    
    /* get the default I/O request for this variable */
    ioRequest* request = IOUtil_DefaultIORequest (cctkGH, n, 1);
 
    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cctkGH, "IO");
    assert (iogh);
    
    // Create the output directory
    const char* myoutdir = out3D_dir;
    if (CCTK_EQUALS(myoutdir, "")) {
      myoutdir = out_dir;
    }
    if (CCTK_MyProc(cctkGH)==0) {
      CCTK_CreateDirectory (0755, myoutdir);
    }
    
    // Invent a file name
    ostringstream filenamebuf;
    filenamebuf << myoutdir << "/" << alias << out3D_extension;
    string filenamestr = filenamebuf.str();
    const char * const filename = filenamestr.c_str();
    
    hid_t writer = -1;
    
    // Write the file only on the root processor
    if (CCTK_MyProc(cctkGH)==0) {
      
      // If this is the first time, then create and truncate the file
      if (do_truncate.at(n)) {
	struct stat fileinfo;
	if (! iogh->recovered || stat(filename, &fileinfo)!=0) {
	  writer = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          assert (writer>=0);
          herr = H5Fclose (writer);
          assert (!herr);
	  writer = -1;
	}
      }
      
      // Open the file 
      writer = H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT);
      assert (writer>=0);
      
    }

    if(verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, 
		  "Writing variable %s on refinement level %d",varname,reflevel);
    }
    
    WriteVar(cctkGH,writer,request,0);

    // Close the file
    if (CCTK_MyProc(cctkGH)==0) {
      herr = H5Fclose (writer);
      assert (!herr);
      writer = -1;
    }
    
    // Don't truncate again
    do_truncate.at(n) = false;
    
    return 0;
  }

  int WriteVar (const cGH* const cctkGH, const hid_t writer, const ioRequest* request,
		   const int called_from_checkpoint) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    herr_t herr=0;
    
    void * h5data=NULL;

    const int n = request->vindex;
    assert (n>=0 && n<CCTK_NumVars());
    const char * varname = CCTK_FullName(n);
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVars());
    int tl = 0;

    const int gpdim = CCTK_GroupDimI(group);
    
    // Check for storage
    if (! CCTK_QueryGroupStorageI(cctkGH, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot output variable \"%s\" because it has no storage",
		  varname);
      return 0;
    }

    const int grouptype = CCTK_GroupTypeI(group);
    switch (grouptype) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      assert (do_global_mode);
      break;
    case CCTK_GF:
      /* do nothing */
      break;
    default:
      assert (0);
    }
    const int rl = grouptype==CCTK_GF ? reflevel : 0;

    cGroup cgdata;
    int ierr = CCTK_GroupData(group,&cgdata);
    assert(ierr==0);

    // let's get the correct Carpet time level (which is the (-1) * Cactus timelevel):
    tl = - request->timelevel;
    
    // Traverse all components
    BEGIN_MAP_LOOP(cctkGH, grouptype) {
      BEGIN_COMPONENT_LOOP(cctkGH, grouptype) {
        
        const ggf<dim>* ff = 0;
        
        ff = (ggf<dim>*)arrdata.at(group).at(Carpet::map).data.at(var);
        
        const gdata<dim>* const data = (*ff) (tl, rl, component, mglevel);
        
        // Make temporary copy on processor 0
        const ibbox ext = data->extent();
//         vect<int,dim> lo = ext.lower();
//         vect<int,dim> hi = ext.upper();
//         vect<int,dim> str = ext.stride();
        
        gdata<dim>* const tmp = data->make_typed (n);
        tmp->allocate (ext, 0);

	if ( !((cgdata.disttype == CCTK_DISTRIB_CONSTANT) && 
	       (arrdata.at(group).at(Carpet::map).hh->processors.at(reflevel).at(component)!=0))) {

	  if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
	    assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
	    if(component!=0) continue;
	    h5data = CCTK_VarDataPtrI(cctkGH,tl,n);
	  } else {
	    for (comm_state<dim> state; !state.done(); state.step()) {
	      tmp->copy_from (state, data, ext);
	    }
	  }
	    // Write data
	    if (CCTK_MyProc(cctkGH)==0) {
          
	      int ldim=0;
	      if ( grouptype==CCTK_SCALAR ) {
		ldim = 1;
	      } else {
		ldim = gpdim;
	      }
	      
	      //	      hsize_t shape[ldim];

	      vector<hsize_t> shape(ldim);

	      for (int d=0; d<ldim; ++d) {
		shape[ldim-1-d] = (ext.shape() / ext.stride())[d];
	      }
	      const hid_t dataspace = H5Screate_simple (ldim, &shape[0], NULL);
	      assert (dataspace>=0);


	      // Select datatype

	      const hid_t datatype = h5DataType(cctkGH,CCTK_VarTypeI(n));
          
	      ostringstream datasetnamebuf;
	      datasetnamebuf << varname
			     << " it=" << cctk_iteration
			     << " tl=" << tl
			     << " ml=" << mglevel
			     << " rl=" << rl
			     << " m=" << Carpet::map
			     << " c=" << component;
	      string datasetnamestr = datasetnamebuf.str();
	      assert (datasetnamestr.size() <= 256); // limit dataset name size
	      const char * const datasetname = datasetnamestr.c_str();
	      const hid_t dataset = H5Dcreate (writer, datasetname, datatype, dataspace, H5P_DEFAULT);
	      
	      if (dataset>=0) {
          
		  if (cgdata.disttype != CCTK_DISTRIB_CONSTANT) {
		      h5data = (void*)tmp->storage();
		  }

		  herr = H5Dwrite (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, h5data);
		  assert (!herr);
          
		  // Write FlexIO attributes
		  WriteAttribute (dataset, "level", rl);
		  {
		      CCTK_REAL origin[dim], delta[dim];
		      CCTK_REAL min_ext[dim], max_ext[dim];
		      for (int d=0; d<dim; ++d) {
			  origin[d] = CCTK_ORIGIN_SPACE(d) + cctk_lbnd[d] * delta[d];
			  delta[d] = CCTK_DELTA_SPACE(d);
			  min_ext[d] = origin[d];
			  max_ext[d] = origin[d] + cctk_lsh[d] * delta[d];
		      }
		      WriteAttribute (dataset, "origin", origin, dim);
		      WriteAttribute (dataset, "delta", delta, dim);
		      WriteAttribute (dataset, "min_ext", min_ext, dim);
		      WriteAttribute (dataset, "max_ext", max_ext, dim);
		  }
		  WriteAttribute (dataset, "time", cctk_time);
		  WriteAttribute (dataset, "timestep", cctk_iteration);
		  WriteAttribute (dataset, "level_timestep", cctk_iteration / reflevelfact);
		  WriteAttribute (dataset, "persistence", maxreflevelfact / reflevelfact);
		  {
		      int time_refinement=0;
		      int spatial_refinement[dim];
		      int grid_placement_refinement[dim];
		      time_refinement = reflevelfact;
		      for (int d=0; d<dim; ++d) {
			  spatial_refinement[d] = reflevelfact;
			  grid_placement_refinement[d] = reflevelfact;
		      }
		      WriteAttribute (dataset, "time_refinement", time_refinement);
		      WriteAttribute (dataset, "spatial_refinement", spatial_refinement, dim);
		      WriteAttribute (dataset, "grid_placement_refinement", grid_placement_refinement, dim);
		  }
		  {
		      int iorigin[dim];
		      for (int d=0; d<dim; ++d) {
			  iorigin[d] = (ext.lower() / ext.stride())[d];
		      }
		      WriteAttribute (dataset, "iorigin", iorigin, dim);
		  }
	      
		  // Write some additional attributes
	      
		  // Legacy arguments
		  {
		      char * fullname = CCTK_FullName(n);
		      assert (fullname);
		      WriteAttribute (dataset, "name", fullname);
		      free (fullname);
		  }
	      
		  // Group arguments
		  WriteAttribute (dataset, "group_version", 1);
		  {
		      char * fullname = CCTK_FullName(n);
		      assert (fullname);
		      WriteAttribute (dataset, "group_fullname", fullname);
		      free (fullname);
		  }
		  WriteAttribute (dataset, "group_varname", CCTK_VarName(n));
		  {
		      char * groupname = CCTK_GroupName(group);
		      assert (groupname);
		      WriteAttribute (dataset, "group_groupname", groupname);
		      free (groupname);
		  }
		  switch (grouptype) {
		      case CCTK_GF:
			  WriteAttribute (dataset, "group_grouptype", "CCTK_GF");
			  break;
		      case CCTK_ARRAY:
			  WriteAttribute (dataset, "group_grouptype", "CCTK_ARRAY");
			  break;
		      case CCTK_SCALAR:
			  WriteAttribute (dataset, "group_grouptype", "CCTK_SCALAR");
			  break;
		      default:
			  assert (0);
		  }
		  WriteAttribute (dataset, "group_dim", CCTK_GroupDimI(group));
		  WriteAttribute (dataset, "group_timelevel", tl);
		  WriteAttribute (dataset, "group_numtimelevels", CCTK_NumTimeLevelsI(group));
	      
		  // Cactus arguments
		  WriteAttribute (dataset, "cctk_version", 1);
		  WriteAttribute (dataset, "cctk_dim", cctk_dim);
		  WriteAttribute (dataset, "cctk_iteration", cctk_iteration);
		  WriteAttribute (dataset, "cctk_gsh", cctk_gsh, dim);
		  WriteAttribute (dataset, "cctk_lsh", cctk_lsh, dim);
		  WriteAttribute (dataset, "cctk_lbnd", cctk_lbnd, dim);
		  WriteAttribute (dataset, "cctk_delta_time", cctk_delta_time);
		  WriteAttribute (dataset, "cctk_delta_space", cctk_delta_space, dim);
		  WriteAttribute (dataset, "cctk_origin_space", cctk_origin_space, dim);
		  WriteAttribute (dataset, "cctk_bbox", cctk_bbox, 2*dim);
		  WriteAttribute (dataset, "cctk_levfac", cctk_levfac, dim);
		  WriteAttribute (dataset, "cctk_levoff", cctk_levoff, dim);
		  WriteAttribute (dataset, "cctk_levoffdenom", cctk_levoffdenom, dim);
		  WriteAttribute (dataset, "cctk_timefac", cctk_timefac);
		  WriteAttribute (dataset, "cctk_convlevel", cctk_convlevel);
		  WriteAttribute (dataset, "cctk_convfac", cctk_convfac);
		  WriteAttribute (dataset, "cctk_nghostzones", cctk_nghostzones, dim);
		  WriteAttribute (dataset, "cctk_time", cctk_time);
	      
		  // Carpet arguments
		  WriteAttribute (dataset, "carpet_version", 1);
		  WriteAttribute (dataset, "carpet_dim", dim);
		  WriteAttribute (dataset, "carpet_basemglevel", basemglevel);
		  WriteAttribute (dataset, "carpet_mglevel", mglevel);
		  WriteAttribute (dataset, "carpet_mglevels", mglevels);
		  WriteAttribute (dataset, "carpet_mgface", mgfact);
		  WriteAttribute (dataset, "carpet_reflevel", reflevel);
		  WriteAttribute (dataset, "carpet_reflevels", reflevels);
		  WriteAttribute (dataset, "carpet_reffact", reffact);
		  WriteAttribute (dataset, "carpet_map", Carpet::map);
		  WriteAttribute (dataset, "carpet_maps", maps);
		  WriteAttribute (dataset, "carpet_component", component);
		  WriteAttribute (dataset, "carpet_components", vhh.at(Carpet::map)->components(reflevel));
	      	    
		  herr = H5Dclose (dataset);
		  assert (!herr);

	      }
	      
	      herr = H5Sclose (dataspace);
	      assert (!herr);
	      
	    } // if on root processor
	    
	    // Delete temporary copy
	    delete tmp;
	} // if ! CCTK_DISTRIB_BLAH
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
    
    return 0;
  }
  
  
  
  int TimeToOutput (const cGH* const cctkGH, const int vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (vindex>=0 && vindex<CCTK_NumVars());
    
    
    
    const int grouptype = CCTK_GroupTypeFromVarI(vindex);
    switch (grouptype) {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      if (! do_global_mode) return 0;
      break;
    case CCTK_GF:
      // do nothing
      break;
    default:
      assert (0);
    }
    
    
    
    // check whether to output at this iteration
    bool output_this_iteration;
    
    const char* myoutcriterion = out3D_criterion;
    if (CCTK_EQUALS(myoutcriterion, "default")) {
      myoutcriterion = out_criterion;
    }
    
    if (CCTK_EQUALS (myoutcriterion, "never")) {
      
      // Never output
      output_this_iteration = false;
      
    } else if (CCTK_EQUALS (myoutcriterion, "iteration")) {
      
      int myoutevery = out3D_every;
      if (myoutevery == -2) {
        myoutevery = out_every;
      }
      if (myoutevery <= 0) {
        // output is disabled
        output_this_iteration = false;
      } else if (cctk_iteration == *this_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >= *next_output_iteration) {
        // it is time for the next output
        output_this_iteration = true;
        *next_output_iteration = cctk_iteration + myoutevery;
        *this_iteration = cctk_iteration;
      } else {
        // we want no output at this iteration
        output_this_iteration = false;
      }
      
    } else if (CCTK_EQUALS (myoutcriterion, "time")) {
      
      CCTK_REAL myoutdt = out3D_dt;
      if (myoutdt == -2) {
        myoutdt = out_dt;
      }
      if (myoutdt < 0) {
        // output is disabled
        output_this_iteration = false;
      } else if (myoutdt == 0) {
        // output all iterations
        output_this_iteration = true;
      } else if (cctk_iteration == *this_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_time / cctk_delta_time
                 >= *next_output_time / cctk_delta_time - 1.0e-12) {
        // it is time for the next output
        output_this_iteration = true;
        *next_output_time = cctk_time + myoutdt;
        *this_iteration = cctk_iteration;
      } else {
        // we want no output at this iteration
        output_this_iteration = false;
      }
      
    } else {
      
      assert (0);
      
    } // select output criterion
    
    if (! output_this_iteration) return 0;
    
    
    
    if (! CheckForVariable(cctkGH, out3D_vars, vindex)) {
      // This variable should not be output
      return 0;
    }
    
    
    
    if (last_output.at(mglevel).at(reflevel).at(vindex) == cctk_iteration) {
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
    
    assert (last_output.at(mglevel).at(reflevel).at(vindex) < cctk_iteration);
    
    // Should be output during this iteration
    return 1;
  }
  
  
  
  int TriggerOutput (const cGH* const cctkGH, const int vindex) {
    assert (vindex>=0 && vindex<CCTK_NumVars());
    
    char* varname = CCTK_FullName(vindex);
    const int retval = OutputVarAs (cctkGH, varname, CCTK_VarName(vindex));
    free (varname);
    
    last_output.at(mglevel).at(reflevel).at(vindex) = cctkGH->cctk_iteration;
    
    return retval;
  }
  
  
  

  int InputGH (const cGH* const cctkGH) {
    DECLARE_CCTK_PARAMETERS;
    int retval = 0;
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (CheckForVariable(cctkGH, in3D_vars, vindex)) {
	char* varname = CCTK_FullName(vindex);
	retval = InputVarAs (cctkGH, varname, CCTK_VarName(vindex));
	free (varname);
	if (retval != 0) return retval;
      }
    }
    return retval;
  }
  
  

  int ReadVar (const cGH* const cctkGH, const hid_t reader, const char* const varname, 
	       const hid_t dataset, vector<ibset> &regions_read, 
	       const int called_from_recovery) {

    DECLARE_CCTK_PARAMETERS;

    const int n = CCTK_VarIndex(varname);
    assert (n>=0 && n<CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI (n);
    assert (group>=0 && group<(int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0 && n0<CCTK_NumVars());
    const int var = n - n0;
    assert (var>=0 && var<CCTK_NumVars());
    int tl = 0;

    herr_t herr = 1;
    bool did_read_something = false;

    // Stuff needed for Recovery 
    int recovery_tl = -1;
    int recovery_mglevel = -1;
    int recovery_rl = -1;
    int recovery_comp = -1;

    void * h5data;

    // Check for storage
    if (! CCTK_QueryGroupStorageI(cctkGH, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot input variable \"%s\" because it has no storage",
                  varname);
      return 0;
    }
    
    const int grouptype = CCTK_GroupTypeI(group);
    const int rl = grouptype==CCTK_GF ? reflevel : 0;

    const int gpdim = CCTK_GroupDimI(group);
       
    
    int amr_level;
    int amr_origin[dim];
    int amr_dims[dim];

    // We need to initialize amr_dims to 0
    for(int i=0;i<dim;i++) amr_dims[i] = 0;

    if (CCTK_MyProc(cctkGH)==0) {
	 
	 // get dataset dimensions
	 const hid_t dataspace = H5Dget_space(dataset);
	 assert (dataspace>=0);
	 hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
	 vector<hsize_t> shape(rank);
	 int rank2 = H5Sget_simple_extent_dims (dataspace, &shape[0], NULL);
	 herr = H5Sclose(dataspace);
	 assert(!herr);
	 assert (rank2 == rank);
	 
	 if(grouptype == CCTK_SCALAR) {
	   assert (gpdim+1 == rank);
	 } else {
	   assert (gpdim == rank);
	 }
 
	 int datalength=1;
	 
	 for(int i=0;i<rank;i++) {
	   datalength=datalength*shape[i];
	   amr_dims[i]=shape[rank-i-1];
	 }

	 if(grouptype == CCTK_ARRAY) {

	   for(int i=rank;i<dim;i++) {
	     amr_dims[i]=1;
	   }

	 }
         const int cctkDataType = CCTK_VarTypeI(n);
	 const hid_t datatype = h5DataType(cctkGH,cctkDataType);
	 
         //cout << "datalength: " << datalength << " rank: " << rank << "\n";
         //cout << shape[0] << " " << shape[1] << " " << shape[2] << "\n";
	 
	 // to do: read in an allocate with correct datatype

	 h5data = (void*) malloc( CCTK_VarTypeSize(cctkDataType) *datalength);

	 herr = H5Dread(dataset,datatype,H5S_ALL, H5S_ALL, H5P_DEFAULT,(void*)h5data);
	 assert(!herr);
	 

	 herr = ReadAttribute(dataset,"level",amr_level);
	 assert(herr>=0);
       	 ReadAttribute(dataset,"iorigin",amr_origin,dim);
	 assert(herr>=0);

	 if(called_from_recovery) {
	   recovery_rl = amr_level;
	   ReadAttribute(dataset,"carpet_component",recovery_comp);
	   ReadAttribute(dataset,"group_timelevel",recovery_tl);
	   ReadAttribute(dataset,"carpet_mglevel",recovery_mglevel);
	 }

     } // MyProc == 0



    if(called_from_recovery) {
      MPI_Bcast (&recovery_rl, 1, MPI_INT, 0, dist::comm);
      MPI_Bcast (&recovery_tl, 1, MPI_INT, 0, dist::comm);
      MPI_Bcast (&recovery_mglevel, 1, MPI_INT, 0, dist::comm);
      MPI_Bcast (&recovery_comp, 1, MPI_INT, 0, dist::comm);
    }
	 
    MPI_Bcast (&amr_level, 1, MPI_INT, 0, dist::comm);
    MPI_Bcast (amr_origin, dim, MPI_INT, 0, dist::comm);
    MPI_Bcast (amr_dims, dim, MPI_INT, 0, dist::comm);

    if ((grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY) && reflevel != 0) {
      return 0;
    }  

    if (verbose) cout << "amr_level: " << amr_level << " reflevel: " << reflevel << endl;
    
    if (amr_level == rl) {
	  
       // Traverse all components on all levels
       BEGIN_MAP_LOOP(cctkGH, grouptype) {
	 BEGIN_COMPONENT_LOOP(cctkGH, grouptype) {

	   did_read_something = true;
	   
	   ggf<dim>* ff = 0;
	     
	   assert (var < (int)arrdata.at(group).at(Carpet::map).data.size());
	   ff = (ggf<dim>*)arrdata.at(group).at(Carpet::map).data.at(var);
	   
	   if(called_from_recovery) tl = recovery_tl;
	     
	     
	   gdata<dim>* const data = (*ff) (tl, rl, component, mglevel);
	     
	   // Create temporary data storage on processor 0
	   vect<int,dim> str
	     = vect<int,dim>(maxreflevelfact/reflevelfact);

	   if(grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY)
	     str = vect<int,dim> (1);

	   vect<int,dim> lb = vect<int,dim>::ref(amr_origin) * str;
	   vect<int,dim> ub
	     = lb + (vect<int,dim>::ref(amr_dims) - 1) * str;
	     
	   gdata<dim>* const tmp = data->make_typed (n);
	   

	   cGroup cgdata;
	   int ierr = CCTK_GroupData(group,&cgdata);
	   assert(ierr==0);
	   //cout << "lb_before: " << lb << endl;
	   //cout << "ub_before: " << ub << endl;
	   if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
	     if (verbose) cout << "CCTK_DISTRIB_CONSTANT: " << varname << endl;
	     assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
	     if (grouptype == CCTK_SCALAR) {
	       lb[0] = arrdata.at(group).at(Carpet::map).hh->processors.at(rl).at(component);
	       ub[0] = arrdata.at(group).at(Carpet::map).hh->processors.at(rl).at(component);
	       for(int i=1;i<dim;i++) {
		 lb[i]=0;
		 ub[i]=0;
	       }
	     } else {
	       const int newlb = lb[gpdim-1] + 
		 (ub[gpdim-1]-lb[gpdim-1]+1)*
		 (arrdata.at(group).at(Carpet::map).hh->processors.at(rl).at(component));
	       const int newub = ub[gpdim-1] + 
		 (ub[gpdim-1]-lb[gpdim-1]+1)*
		 (arrdata.at(group).at(Carpet::map).hh->processors.at(rl).at(component));
	       lb[gpdim-1] = newlb;
	       ub[gpdim-1] = newub;
	     }
	      if (verbose)  cout << "lb: " << lb << endl;
	      if (verbose)  cout << "ub: " << ub << endl;
	   }
	   const bbox<int,dim> ext(lb,ub,str);
	   
	   if (verbose) cout << "ext: " << ext << endl;
	     
	   if (CCTK_MyProc(cctkGH)==0) {
	     tmp->allocate (ext, 0, h5data);
	   } else {
	     tmp->allocate (ext, 0);
	   }

	   // Initialise with what is found in the file -- this does
	   // not guarantee that everything is initialised.
	   const bbox<int,dim> overlap = tmp->extent() & data->extent();
	   regions_read.at(Carpet::map) |= overlap;

	   if (verbose) {
	     cout << "working on component: " << component << endl;
	     cout << "tmp->extent " << tmp->extent() << endl;
	     cout << "data->extent " << data->extent() << endl;
	     cout << "overlap " << overlap << endl;
	     cout << "-----------------------------------------------------" << endl;
	   }
	   MPI_Barrier(MPI_COMM_WORLD);
	   
	   // Copy into grid function
	   for (comm_state<dim> state; !state.done(); state.step()) {
	     data->copy_from (state, tmp, overlap);
	   }

	 	   
	   // Delete temporary copy
	   delete tmp;
	     
	 	 
	 } END_COMPONENT_LOOP;

	 if (called_from_recovery) {
	   arrdata.at(group).at(Carpet::map).tt->set_time(reflevel,mglevel,
		   (CCTK_REAL) ((cctkGH->cctk_time - cctk_initial_time)
	                        / (delta_time * mglevelfact)) );
	 }
	 

       } END_MAP_LOOP;

    } // if amr_level == rl       
     if (CCTK_MyProc(cctkGH)==0) {
       free (h5data);
     }
  
    return did_read_something;

  }



  int InputVarAs (const cGH* const cctkGH, const char* const varname, 
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

    herr_t herr = 1;
    int want_dataset = 0;
    bool did_read_something = false;
    int ndatasets = 0;
    hid_t dataset = 0;

    char datasetname[1024];

    CCTK_REAL *h5data;

    // Check for storage
    if (! CCTK_QueryGroupStorageI(cctkGH, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot input variable \"%s\" because it has no storage",
                  varname);
      return 0;
    }
    
    const int grouptype = CCTK_GroupTypeI(group);
    const int rl = grouptype==CCTK_GF ? reflevel : 0;
    //cout << "want level " << rl << " reflevel " << reflevel << endl;
    
    // Find the input directory
    const char* const myindir = in3D_dir;
    
    // Invent a file name
    ostringstream filenamebuf;
    filenamebuf << myindir << "/" << alias << in3D_extension;
    string filenamestr = filenamebuf.str();
    const char * const filename = filenamestr.c_str();
    
    hid_t reader = -1;
    
    const int gpdim = CCTK_GroupDimI(group);
    
    // Read the file only on the root processor
    if (CCTK_MyProc(cctkGH)==0) {
      
      // Open the file 
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Opening file \"%s\"", filename);
      reader = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (reader<0) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not open file \"%s\" for reading", filename);
      }
      assert (reader>=0);
      // get the number of datasets in the file 
      ndatasets=GetnDatasets(reader);
    }
    
    vector<ibset> regions_read(Carpet::maps);
    
    // Broadcast number of datasets
    MPI_Bcast (&ndatasets, 1, MPI_INT, 0, dist::comm);
    assert (ndatasets>=0);
      

    for (int datasetid=0; datasetid<ndatasets; ++datasetid) {
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Handling dataset #%d", datasetid);

        
      // Read data
      if (CCTK_MyProc(cctkGH)==0) {
	GetDatasetName(reader,datasetid,datasetname);
	// 	cout << datasetname << "\n";
  
	dataset = H5Dopen (reader, datasetname);
	assert(dataset);
      }
    
      
      int amr_level;
      int amr_origin[dim];
      int amr_dims[dim];

      if (CCTK_MyProc(cctkGH)==0) {
	  
       // Read data
       char * name;
       ReadAttribute (dataset, "name", name);
       //        cout << "dataset name is " << name << endl;
       if (verbose) {
	 if (name) {
		CCTK_VInfo (CCTK_THORNSTRING, "Dataset name is \"%s\"", name);
	 }
       }
       want_dataset = name && CCTK_EQUALS(name, varname);
       free (name);
      } // myproc == 0 

      MPI_Bcast (&want_dataset, 1, MPI_INT, 0, dist::comm);
     

       if(want_dataset) {

	 did_read_something = ReadVar(cctkGH,reader,varname,dataset,regions_read,0);

       } // want_dataset

    } // loop over datasets
      
    // Close the file
    if (CCTK_MyProc(cctkGH)==0) {
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Closing file");
      herr = H5Fclose(reader);
      //	  cout << "blah! " << reader << "\n";
      // cout << "closing file " << herr << "\n";
      assert(!herr);
      reader=-1;
    }
  
    // Was everything initialised?
    if (did_read_something) {
      for (int m=0; m<Carpet::maps; ++m) {
	dh<dim>& thedd = *arrdata.at(group).at(m).dd;
	ibset all_exterior;
	for (size_t c=0; c<thedd.boxes.at(rl).size(); ++c) {
	  all_exterior |= thedd.boxes.at(rl).at(c).at(mglevel).exterior;
        }
	if (regions_read.at(m) != all_exterior) {
          cout << "read: " << regions_read.at(m) << endl
               << "want: " << all_exterior << endl;
	  CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		      "Variable \"%s\" could not be initialised from file -- the file may be missing data",
		      varname);
	}
      }
    } // if did_read_something
    //	CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,"stop!");    

    return did_read_something ? 0 : -1;

  }

  
  /** returns the number of recovered variables */
  int Recover (cGH* const cctkGH, const char *basefilename,
               const int called_from)
  {
    assert (cctkGH);
    assert (basefilename);
    assert (called_from == CP_INITIAL_DATA
            || called_from == CP_EVOLUTION_DATA
            || called_from == CP_RECOVER_PARAMETERS
            || called_from == CP_RECOVER_DATA
            || called_from == FILEREADER_DATA);
    
    // the other modes are not supported yet
    assert (called_from == FILEREADER_DATA);
    
    const ioGH * const iogh = (const ioGH *) CCTK_GHExtension (cctkGH, "IO");
    assert (iogh);
    
    int num_vars_read = 0;
    assert (iogh->do_inVars);
    for (int n=0; n<CCTK_NumVars(); ++n) {
      if (iogh->do_inVars[n]) {
        char * const fullname = CCTK_FullName(n);
        assert (fullname);
        const int ierr = InputVarAs (cctkGH, fullname, basefilename);
        if (! ierr) {
          ++ num_vars_read;
        }
        free (fullname);
      }
    }
    
    return num_vars_read;
  }

  
  
  void CarpetIOHDF5ReadData (CCTK_ARGUMENTS)
  {
    InputGH (cctkGH);
  }

  
  
} // namespace CarpetIOHDF5
