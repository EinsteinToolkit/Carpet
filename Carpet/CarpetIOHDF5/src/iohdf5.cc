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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5.cc,v 1.1 2004/03/03 09:44:26 schnetter Exp $";
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



namespace CarpetIOHDF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Variable definitions
  int GHExtension;
  int IOMethod;
  vector<bool> do_truncate;     // [var]
  vector<vector<vector<int> > > last_output; // [ml][rl][var]
  
  
  
  const char* GetStringParameter (const char* const parametername,
                                  const char* const fallback);
  int GetIntParameter (const char* const parametername, int fallback);
  
  bool CheckForVariable (const cGH* const cctkGH,
                         const char* const varlist, const int vindex);
  void SetFlag (int index, const char* optstring, void* arg);
  
  void WriteAttribute (const hid_t dataset, const char* name, int value);
  void WriteAttribute (const hid_t dataset, const char* name, const int* values, int nvalues);
  void WriteAttribute (const hid_t dataset, const char* name, double value);
  void WriteAttribute (const hid_t dataset, const char* name, const double* values, int nvalues);
  void WriteAttribute (const hid_t dataset, const char* name, char value);
  void WriteAttribute (const hid_t dataset, const char* name, const char* values);
  void WriteAttribute (const hid_t dataset, const char* name, const char* values, int nvalues);
  
  int ReadAttribute (const hid_t dataset, const char* name, int& value);
  int ReadAttribute (const hid_t dataset, const char* name, int* values, int nvalues);
  int ReadAttribute (const hid_t dataset, const char* name, double& value);
  int ReadAttribute (const hid_t dataset, const char* name, double* values, int nvalues);
  int ReadAttribute (const hid_t dataset, const char* name, char& value);
  int ReadAttribute (const hid_t dataset, const char* name, char*& values);
  int ReadAttribute (const hid_t dataset, const char* name, char* values, int nvalues);
  
  
  
  int CarpetIOHDF5Startup ()
  {
    int ierr;
    
    CCTK_RegisterBanner ("AMR 3D HDF5 I/O provided by CarpetIOHDF5");
    
    GHExtension = CCTK_RegisterGHExtension ("CarpetIOHDF5");
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    IOMethod = CCTK_RegisterIOMethod ("CarpetIOHDF5");
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
    
#if 0
    ierr = IOUtil_RegisterRecover ("CarpetIOHDF5", Recover);
    assert (! ierr);
#endif
    
    return 0;
  }
  
  
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
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
    return 0;
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
    
    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cctkGH, "IO");
    assert (iogh);
    
    // Create the output directory
    const char* const myoutdir = GetStringParameter("out3D_dir", out_dir);
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
      
//       const int gpdim = CCTK_GroupDimI(group);
      
//       // Set coordinate information
//       double origin[dim], delta[dim], timestep;
//       for (int d=0; d<dim; ++d) {
// 	origin[d] = cctkGH->cctk_origin_space[d];
// 	delta[d] = cctkGH->cctk_delta_space[d];
//       }
//       timestep = cctkGH->cctk_delta_time;
//       amrwriter->setTopLevelParameters
// 	(gpdim, origin, delta, timestep, maxreflevels);
      
//       // Set refinement information
//       int interlevel_timerefinement;
//       int interlevel_spacerefinement[dim];
//       int initial_gridplacementrefinement[dim];
//       interlevel_timerefinement = reffact;
//       for (int d=0; d<dim; ++d) {
// 	interlevel_spacerefinement[d] = reffact;
// 	initial_gridplacementrefinement[d] = 1;
//       }
//       amrwriter->setRefinement
// 	(interlevel_timerefinement, interlevel_spacerefinement,
// 	 initial_gridplacementrefinement);
      
//       // Set level
//       amrwriter->setLevel (rl);
      
//       // Set current time
//       amrwriter->setTime (cctk_iteration);
    }
    
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
        for (comm_state<dim> state; !state.done(); state.step()) {
          tmp->copy_from (state, data, ext);
        }
        
        // Write data
        if (CCTK_MyProc(cctkGH)==0) {
          
          hsize_t shape[dim];
          for (int d=0; d<dim; ++d) {
            shape[dim-1-d] = (ext.shape() / ext.stride())[d];
          }
          const hid_t dataspace = H5Screate_simple (dim, shape, NULL);
          assert (dataspace>=0);
          
          // Select datatype
          assert (true
                  || (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL8
                      && sizeof(CCTK_REAL8) == sizeof(double))
                  || (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL
                      && sizeof(CCTK_REAL) == sizeof(double)));
          // TODO: Set datatype correctly
          const hid_t datatype = H5T_NATIVE_DOUBLE;
          
          ostringstream datasetnamebuf;
          datasetnamebuf << varname
                         << " it=" << cctk_iteration
                         << " ml=" << mglevel
                         << " rl=" << rl
                         << " m=" << Carpet::map
                         << " c=" << component;
          string datasetnamestr = datasetnamebuf.str();
          const char * const datasetname = datasetnamestr.c_str();
          const hid_t dataset = H5Dcreate (writer, datasetname, datatype, dataspace, H5P_DEFAULT);
          assert (dataset>=0);
          
          const void * const data = (void*)tmp->storage();
          herr = H5Dwrite (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
          assert (!herr);
          
          // Write FlexIO attributes
          WriteAttribute (dataset, "level", rl);
          {
            CCTK_REAL origin[dim], delta[dim];
            CCTK_REAL min_ext[dim], max_ext[dim];
            for (int d=0; d<dim; ++d) {
              origin[d] = CCTK_ORIGIN_SPACE(d);
              delta[d] = CCTK_DELTA_SPACE(d);
              min_ext[d] = origin[d];
              max_ext[d] = origin[d] + (cctk_lsh[d] - 1) * delta[d];
            }
            WriteAttribute (dataset, "origin", origin, dim);
            WriteAttribute (dataset, "delta", delta, dim);
            WriteAttribute (dataset, "min_ext", min_ext, dim);
            WriteAttribute (dataset, "max_ext", max_ext, dim);
          }
          WriteAttribute (dataset, "time", cctk_time);
          WriteAttribute (dataset, "timestep", cctk_iteration);
          WriteAttribute (dataset, "level_timestep", cctk_iteration / reflevelfact);
          WriteAttribute (dataset, "persistence", reflevelfact);
          {
            int time_refinement;
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
// TODO: disable temporarily
//           WriteAttribute (dataset, "cctk_nmaps", cctk_nmaps);
//           WriteAttribute (dataset, "cctk_map", cctk_map);
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
          
          herr = H5Sclose (dataspace);
          assert (!herr);
          
        } // if on root processor
        
        // Delete temporary copy
        delete tmp;
        
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
    
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
  
  
  
  int TimeToOutput (const cGH* const cctkGH, const int vindex) {
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
      /* do nothing */
      break;
    default:
      assert (0);
    }
    
    const int myoutevery = GetIntParameter("out3D_every", out_every);
    
    if (myoutevery < 0) {
      // Nothing should be output at all
      return 0;
    }
    
    if (cctk_iteration % myoutevery != 0) {
      // Nothing should be output during this iteration
      return 0;
    }
    
    if (! CheckForVariable(cctkGH, GetStringParameter("out3D_vars",""), vindex)) {
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
  
  
  
#if 0
  int InputGH (const cGH* const cctkGH) {
    int retval = 0;
    for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
      if (CheckForVariable(cctkGH, GetStringParameter("in3D_vars",""), vindex)) {
	char* varname = CCTK_FullName(vindex);
	retval = InputVarAs (cctkGH, varname, CCTK_VarName(vindex));
	free (varname);
	if (retval != 0) return retval;
      }
    }
    return retval;
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
    
    // Check for storage
    if (! CCTK_QueryGroupStorageI(cctkGH, group)) {
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
    ostringstream filenamebuf;
    filenamebuf << myindir << "/" << alias << in3D_extension;
    string filenamestr = filenamebuf.str();
    const char * const filename = filenamestr.c_str();
    
    hid_t reader = -1;
    
//     const int gpdim = CCTK_GroupDimI(group);
    
//     int have_dataset;
    
//     int rank;
//     int dims[dim];
//     int nbytes;
    
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
      
    }
    
    vector<ibset> regions_read(Carpet::maps);
    
    // Traverse all components on all levels
    BEGIN_MAP_LOOP(cctkGH, grouptype) {
      BEGIN_COMPONENT_LOOP(cctkGH, grouptype) {
        
        // Read data
        if (CCTK_MyProc(cctkGH)==0) {
          
          ostringstream datasetnamebuf;
          datasetnamebuf << varname
                         << " it=" << cctk_iteration
                         << " ml=" << mglevel
                         << " rl=" << rl
                         << " m=" << Carpet::map
                         << " c=" << component;
          string datasetnamestr = datasetnamebufs.str();
          const char * const datasetname = datasetnamestr.c_str();
          const hid_t dataset = H5Dopen (reader, datasetname);
          have_dataset = dataset>=0;
          
//       // Check rank
//       assert (rank==gpdim);
      
//       // Check datatype
//           assert (true
//                   || (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL8
//                       && sizeof(CCTK_REAL8) == sizeof(double))
//                   || (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL
//                       && sizeof(CCTK_REAL) == sizeof(double)));
//       // TODO: Check datatype correctly
//       assert (CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL8
//               || (sizeof(CCTK_REAL) == sizeof(double)
//                   && CCTK_VarTypeI(n) == CCTK_VARIABLE_REAL));
      
//       // TODO: check grid spacing
        }
    
//     // Broadcast rank, dimensions, and nbytes
//     MPI_Bcast (&rank, 1, MPI_INT, 0, dist::comm);
//     assert (rank>=1);
//     MPI_Bcast (&dims, rank, MPI_INT, 0, dist::comm);
//     for (int d=0; d<rank; ++d) assert (dims[d]>=0);
//     MPI_Bcast (&nbytes, 1, MPI_INT, 0, dist::comm);
//     assert (nbytes>=0);
    
//     // Broadcast number of datasets
//     MPI_Bcast (&ndatasets, 1, MPI_INT, 0, dist::comm);
//     assert (ndatasets>=0);
        
        // Broadcast dataset
        MPI_Bcast (&have_dataset, 1, MPI_INT, 0, dist::comm);
        
        if (have_dataset) {
          
          ###
      // Read grid
      int amr_level;
      int amr_origin[dim];
      int amr_dims[dim];
      
      if (CCTK_MyProc(cctkGH)==0) {
        
        // Read data
        const hid_t dataset
        
        {
          char * name;
          ReadAttribute (reader, "name", name);
          if (verbose) {
            if (name) {
              CCTK_VInfo (CCTK_THORNSTRING, "Dataset name is \"%s\"", name);
            }
          }
          want_dataset = name && CCTK_EQUALS(name, varname);
          free (name);
        }
        
        // If iorigin attribute is absent, assume file has unigrid
        // data.
        {
          IObase::DataType atype;
          int alength;
          if (reader->readAttributeInfo("iorigin", atype, alength) < 0) {
            amrgrid->level = 0;
            for (int d=0; d<gpdim; ++d) {
              amrgrid->iorigin[d] = 0;
            }
          }
        }
        
        amr_level = amrgrid->level;
        for (int d=0; d<gpdim; ++d) {
          amr_origin[d] = amrgrid->iorigin[d];
          amr_dims[d] = amrgrid->dims[d];
        }
        for (int d=gpdim; d<dim; ++d) {
          amr_origin[d] = 0;
          amr_dims[d] = 1;
        }
        
      } // MyProc == 0
      
      MPI_Bcast (&want_dataset, 1, MPI_INT, 0, dist::comm);
      MPI_Bcast (&amr_level, 1, MPI_INT, 0, dist::comm);
      MPI_Bcast (amr_origin, dim, MPI_INT, 0, dist::comm);
      MPI_Bcast (amr_dims, dim, MPI_INT, 0, dist::comm);
      
      if (want_dataset && amr_level == reflevel) {
        did_read_something = true;
        
        // Traverse all components on all levels
        BEGIN_MAP_LOOP(cctkGH, grouptype) {
          BEGIN_COMPONENT_LOOP(cctkGH, grouptype) {
            
            ggf<dim>* ff = 0;
            
            assert (var < (int)arrdata.at(group).at(Carpet::map).data.size());
            ff = (ggf<dim>*)arrdata.at(group).at(Carpet::map).data.at(var);
            
            gdata<dim>* const data = (*ff) (tl, rl, component, mglevel);
            
            // Create temporary data storage on processor 0
            const vect<int,dim> str
              = vect<int,dim>(maxreflevelfact/reflevelfact);
            const vect<int,dim> lb = vect<int,dim>::ref(amr_origin) * str;
            const vect<int,dim> ub
              = lb + (vect<int,dim>::ref(amr_dims) - 1) * str;
            const bbox<int,dim> ext(lb,ub,str);
            
            gdata<dim>* const tmp = data->make_typed (n);
            
            if (CCTK_MyProc(cctkGH)==0) {
              tmp->allocate (ext, 0, amrgrid->data);
            } else {
              tmp->allocate (ext, 0);
            }
            
            // Initialise with what is found in the file -- this does
            // not guarantee that everything is initialised.
            const bbox<int,dim> overlap = tmp->extent() & data->extent();
            regions_read.at(Carpet::map) |= overlap;
            
            // Copy into grid function
            for (comm_state<dim> state; !state.done(); state.step()) {
              data->copy_from (state, tmp, overlap);
            }
            
            // Delete temporary copy
            delete tmp;
            
          } END_COMPONENT_LOOP;
        } END_MAP_LOOP;
        
      } // if want_dataset && level == reflevel
      
      if (CCTK_MyProc(cctkGH)==0) {
        free (amrgrid->data);
        free (amrgrid);
        amrgrid = 0;
      }
      
    } // loop over datasets
    
    // Close the file
    if (CCTK_MyProc(cctkGH)==0) {
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Deleting AMR info");
      delete amrreader;

      amrreader = 0;
      if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Closing file");
      delete reader;
      reader = 0;
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
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Variable \"%s\" could not be initialised from file -- the file may be missing data",
                      varname);
        }
      }
    } // if did_read_something
    
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
  
  
  
  int CarpetIOHDF5ReadData (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    return InputGH(cctkGH);
  }
#endif
  
  
  
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
  
  
  
  bool CheckForVariable (const cGH* const cctkGH,
			 const char* const varlist, const int vindex)
  {
    const int numvars = CCTK_NumVars();
    assert (vindex>=0 && vindex<numvars);
    
    vector<bool> flags(numvars);
    
    CCTK_TraverseString (varlist, SetFlag, &flags, CCTK_GROUP_OR_VAR);
    
    return flags.at(vindex);
  }
  
  void SetFlag (int index, const char* optstring, void* arg)
  {
    vector<bool>& flags = *(vector<bool>*)arg;
    flags.at(index) = true;
  }
  
  
  
  void WriteAttribute (const hid_t dataset, const char* const name, const int value)
  {
    WriteAttribute (dataset, name, &value, 1);
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const int* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    hsize_t shape[1];
    shape[0] = nvalues;
    const hid_t dataspace = nvalues==1 ? H5Screate (H5S_SCALAR) : H5Screate_simple (1, shape, NULL);
    assert (dataspace>=0);
    
    const hid_t datatype = H5T_NATIVE_INT;
    
    const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
    assert (attribute>=0);
    herr = H5Awrite (attribute, datatype, values);
    assert (!herr);
    herr = H5Aclose (attribute);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
  }
  
  
  
  void WriteAttribute (const hid_t dataset, const char* const name, const double value)
  {
    WriteAttribute (dataset, name, &value, 1);
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const double* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    hsize_t shape[1];
    shape[0] = nvalues;
    const hid_t dataspace = nvalues==1 ? H5Screate (H5S_SCALAR) : H5Screate_simple (1, shape, NULL);
    assert (dataspace>=0);
    
    const hid_t datatype = H5T_NATIVE_DOUBLE;
    
    const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
    assert (attribute>=0);
    herr = H5Awrite (attribute, datatype, values);
    assert (!herr);
    herr = H5Aclose (attribute);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
  }
  
  
  
  void WriteAttribute (const hid_t dataset, const char* const name, const char value)
  {
    WriteAttribute (dataset, name, &value, 1);
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const char* const values)
  {
    WriteAttribute (dataset, name, values, strlen(values));
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const char* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t dataspace = H5Screate (H5S_SCALAR);
    assert (dataspace>=0);
    
    const hid_t datatype = H5Tcopy (H5T_C_S1);
    assert (datatype>=0);
    herr = H5Tset_size (datatype, nvalues);
    assert (!herr);
    
    const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
    assert (attribute>=0);
    herr = H5Awrite (attribute, datatype, values);
    assert (!herr);
    herr = H5Aclose (attribute);
    assert (!herr);
    
    herr = H5Tclose (datatype);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
  }
  
  
  
  int ReadAttribute (const hid_t dataset, const char* const name, int& value)
  {
    return ReadAttribute (dataset, name, &value, 1);
  }
  
  int ReadAttribute (const hid_t dataset, const char* name, int* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    hsize_t shape[1];
    if (rank==0) {
      shape[0] = 1;
    } else if (rank==1) {
      herr = H5Sget_simple_extent_dims (dataspace, shape, NULL);
      assert (!herr);
    } else {
      assert (0);
    }
    const int length = shape[0];
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (datatype != H5T_NATIVE_INT) return -100;
    
    vector<int> values1(length);
    
    herr = H5Aread (attribute, datatype, &values1.at(0));
    assert (!herr);
    
    for (int i=0; i<min(length, nvalues); ++i) {
      values[i] = values1[i];
    }
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  
  
  int ReadAttribute (const hid_t dataset, const char* const name, double& value)
  {
    return ReadAttribute (dataset, name, &value, 1);
  }
  
  int ReadAttribute (const hid_t dataset, const char* const name, double* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    hsize_t shape[1];
    if (rank==0) {
      shape[0] = 1;
    } else if (rank==1) {
      herr = H5Sget_simple_extent_dims (dataspace, shape, NULL);
      assert (!herr);
    } else {
      assert (0);
    }
    const int length = shape[0];
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (datatype != H5T_NATIVE_DOUBLE) return -100;
    
    vector<double> values1(length);
    
    herr = H5Aread (attribute, datatype, &values1.at(0));
    assert (!herr);
    
    for (int i=0; i<min(length, nvalues); ++i) {
      values[i] = values1[i];
    }
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  
  
  int ReadAttribute (const hid_t dataset, const char* const name, char& value)
  {
    return ReadAttribute (dataset, name, &value, 1);
  }
  
  int ReadAttribute (const hid_t dataset, const char* const name, char*& values)
  {
    assert (dataset>=0);
    assert (name);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    assert (rank==0);
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (H5Tget_class (datatype) != H5T_STRING) return -100;
    const int length = H5Tget_size (datatype);
    assert (length>=0);
    
    values = (char*) malloc (length+1);
    assert (values);
    
    herr = H5Aread (attribute, datatype, values);
    assert (!herr);
    values[length] = '\0';
    
    return length;
  }
  
  int ReadAttribute (const hid_t dataset, const char* const name, char* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    assert (rank==0);
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (H5Tget_class (datatype) != H5T_STRING) return -100;
    const int length = H5Tget_size (datatype);
    assert (length>=0);
    
    vector<char> values1(length);
    
    herr = H5Aread (attribute, datatype, &values1.at(0));
    assert (!herr);
    
    for (int i=0; i<min(length, nvalues); ++i) {
      values[i] = values1[i];
    }
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  
  
} // namespace CarpetIOHDF5
