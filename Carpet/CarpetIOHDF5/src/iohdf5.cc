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
  static const char* rcsid = "$Header: /home/cvs/carpet/Carpet/CarpetIOHDF5/src/iohdf5.cc,v 1.42 2004/10/08 13:30:18 cott Exp $";
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
vector<bool> do_truncate;     // [var]
vector<vector<vector<int> > > last_output; // [ml][rl][var]


static void AddAttributes (const cGH *const cctkGH, const char *fullname,
                           int vdim, int refinementlevel, int timelevel,
                           ibset::const_iterator bbox, hid_t dataset);

void CarpetIOHDF5Startup (void)
{
  DECLARE_CCTK_PARAMETERS


  CCTK_RegisterBanner ("AMR 3D HDF5 I/O provided by CarpetIOHDF5");

  int GHExtension = CCTK_RegisterGHExtension ("CarpetIOHDF5");
  CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);

  int IOMethod = CCTK_RegisterIOMethod ("IOHDF5");
  CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);

  /* initial I/O parameter check */
  int numvars = CCTK_NumVars ();
  vector<bool> flags(numvars);

  if (CCTK_TraverseString (out3D_vars, SetFlag, &flags,CCTK_GROUP_OR_VAR) < 0)
  {
    CCTK_VWarn (strict_io_parameter_check ? 0 : 1,
                __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing parameter 'IOHDF5::out3D_vars'");
  }

  if (IOUtil_RegisterRecover ("CarpetIOHDF5 recovery", Recover) < 0)
  {
    CCTK_WARN (1, "Failed to register " CCTK_THORNSTRING " recovery routine");
  }
}



int CarpetIOHDF5Init (const cGH* const cctkGH)
{
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = -1;
  *next_output_iteration = 0;
  *next_output_time = cctk_time;

  return (0);
}



void* SetupGH (tFleshConfig* const fc,
               const int convLevel, cGH* const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  CarpetIOHDF5GH* myGH;

  const void *dummy;
  dummy = &fc;
  dummy = &convLevel;
  dummy = &cctkGH;
  dummy = &dummy;

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


int OutputGH (const cGH* const cctkGH)
{
  for (int vindex=0; vindex<CCTK_NumVars(); ++vindex)
  {
    if (TimeToOutput(cctkGH, vindex))
    {
      TriggerOutput(cctkGH, vindex);
    }
  }
  return 0;
}


int OutputVarAs (const cGH* const cctkGH, const char* const varname,
                 const char* const alias)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int numvars = CCTK_NumVars ();
  vector<bool> flags (numvars);

  if (CCTK_TraverseString (varname, SetFlag, &flags, CCTK_VAR) < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing variable name '%s' (alias name '%s')",
                varname, alias);
    return (-1);
  }

  int vindex = 0;
  while (! flags.at (vindex) && vindex < numvars) vindex++;
  if (vindex >= numvars)
  {
    return (-1);
  }

  const int group = CCTK_GroupIndexFromVarI (vindex);
  assert (group>=0 && group<(int)Carpet::arrdata.size());

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, group))
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot output variable '%s' because it has no storage",
                varname);
    return (0);
  }

  const int grouptype = CCTK_GroupTypeI(group);
  if (grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY)
  {
    assert (do_global_mode);
  }

  const int myproc = CCTK_MyProc (cctkGH);

  /* get the default I/O request for this variable */
  ioRequest* request = IOUtil_DefaultIORequest (cctkGH, vindex, 1);

  // Get grid hierarchy extentsion from IOUtil
  const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cctkGH, "IO");
  assert (iogh);

  // Create the output directory
  const char* myoutdir = *out3D_dir ? out3D_dir : out_dir;
  if (myproc == 0)
  {
    CCTK_CreateDirectory (0755, myoutdir);
  }

  // Invent a file name
  ostringstream filenamebuf;
  filenamebuf << myoutdir << "/" << alias << out3D_extension;
  string filenamestr = filenamebuf.str();
  const char * const filename = filenamestr.c_str();

  hid_t writer = -1;

  // Write the file only on the root processor
  if (myproc == 0)
  {
    // If this is the first time, then create and truncate the file
    if (do_truncate.at(vindex))
    {
      struct stat fileinfo;
      if (IO_TruncateOutputFiles (cctkGH) || stat(filename, &fileinfo)!=0)
      {
        HDF5_ERROR (writer = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT,
                                        H5P_DEFAULT));
        assert (writer>=0);
        HDF5_ERROR (H5Fclose (writer));
        writer = -1;
      }
    }

    // Open the file
    HDF5_ERROR (writer = H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT));
  }

  if (verbose)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "Writing variable '%s' on mglevel %d reflevel %d",
                varname, mglevel, reflevel);
  }
  WriteVar (cctkGH, writer, request, 0);

  // Close the file
  if (writer >= 0)
  {
    HDF5_ERROR (H5Fclose (writer));
  }

  // Don't truncate again
  do_truncate.at(vindex) = false;

  return (0);
}


int WriteVar (const cGH* const cctkGH, const hid_t writer,
              const ioRequest* request, const int called_from_checkpoint)
{
  DECLARE_CCTK_PARAMETERS

  char *fullname = CCTK_FullName(request->vindex);
  const int gindex = CCTK_GroupIndexFromVarI (request->vindex);
  assert (gindex >= 0 && gindex < (int) Carpet::arrdata.size ());
  const int var = request->vindex - CCTK_FirstVarIndexI (gindex);
  assert (var >= 0 && var < CCTK_NumVars ());
  cGroup group;
  CCTK_GroupData (gindex, &group);


  // Scalars and arrays have only one refinement level 0,
  // regardless of what the current refinement level is.
  // Output for them must be called in global mode.
  int refinementlevel = reflevel;
  if (group.grouptype == CCTK_SCALAR || group.grouptype == CCTK_ARRAY)
  {
    assert (do_global_mode);
    refinementlevel = 0;
  }

  // HDF5 doesn't like 0-dimensional arrays
  if (group.grouptype == CCTK_SCALAR)
  {
    group.dim = 1;
  }

  // If the user requested so, select single precision data output
  hid_t memdatatype, filedatatype;
  HDF5_ERROR (memdatatype = filedatatype = h5DataType (cctkGH, group.vartype));
  if (out_single_precision && ! called_from_checkpoint)
  {
    if (group.vartype == CCTK_VARIABLE_REAL)
    {
      HDF5_ERROR (filedatatype = h5DataType (cctkGH, CCTK_VARIABLE_REAL4));
    }
    else if (group.vartype == CCTK_VARIABLE_COMPLEX)
    {
      HDF5_ERROR (filedatatype = h5DataType (cctkGH, CCTK_VARIABLE_COMPLEX8));
    }
#ifdef CCTK_INT2
    else if (group.vartype == CCTK_VARIABLE_INT)
    {
      HDF5_ERROR (filedatatype = h5DataType (cctkGH, CCTK_VARIABLE_INT2));
    }
#endif
  }

  const int myproc = CCTK_MyProc (cctkGH);

  // create a file access property list to use the CORE virtual file driver
  hid_t plist;
  HDF5_ERROR (plist = H5Pcreate (H5P_FILE_ACCESS));
  HDF5_ERROR (H5Pset_fapl_core (plist, 0, 0));

  // Traverse all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype)
  {
    // Collect the set of all components' bboxes
    ibset bboxes;
    BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype)
    {
      // Per Cactus convention, scalars are assumed the same values on each
      // processor, so only component 0 is output.
      if (component == 0 || group.grouptype != CCTK_SCALAR)
      {
        // Using "interior" removes ghost zones and refinement boundaries.
        bboxes += arrdata.at(gindex).at(Carpet::map).dd->
                  boxes.at(refinementlevel).at(component).at(mglevel).interior;
      }
    } END_COMPONENT_LOOP;

    // Normalise the set, i.e., try to represent the set with fewer bboxes.
    bboxes.normalize();

    // Loop over all components in the bbox set
    int bbox_id = 0;
    for (ibset::const_iterator bbox  = bboxes.begin();
                               bbox != bboxes.end();
                               bbox++, bbox_id++)
    {
      // Get the shape of the HDF5 dataset (in Fortran index order)
      hsize_t shape[dim];
      hsize_t num_elems = 1;
      for (int d = 0; d < group.dim; ++d)
      {
        shape[group.dim-1-d] = (bbox->shape() / bbox->stride())[d];
        num_elems *= shape[group.dim-1-d];
      }

      // Don't create zero-sized components
      if (num_elems <= 0)
      {
        continue;
      }

      hid_t memfile = -1, memdataset = -1;
      hid_t dataspace = -1, dataset = -1;
      if (myproc == 0)
      {
        // Construct a file-wide unique HDF5 dataset name
        // (only add parts with varying metadata)
        ostringstream datasetname;
        datasetname << fullname
                    << " it=" << cctkGH->cctk_iteration
                    << " tl=" << request->timelevel;
        if (mglevels > 1)
        {
          datasetname << " ml=" << mglevel;
        }
        if (Carpet::maps > 1)
        {
          datasetname << " m="  << Carpet::map;
        }
        if (group.grouptype == CCTK_GF)
        {
          datasetname << " rl=" << refinementlevel;
        }
        if (bboxes.setsize () > 1)
        {
          datasetname << " c=" << bbox_id;
        }

        // We create a temporary HDF5 file in memory (using the core VFD)
        // in order to do the recombination of individual processor components
        // for a single dataset.
        // Although this requires some more memory, it should be much faster
        // than recombining an HDF5 dataset on a disk file.
        HDF5_ERROR (memfile = H5Fcreate ("tempfile", H5F_ACC_EXCL, H5P_DEFAULT,
                                         plist));
        HDF5_ERROR (dataspace = H5Screate_simple (group.dim, shape, NULL));
        HDF5_ERROR (memdataset = H5Dcreate (memfile, datasetname.str().c_str(),
                                         filedatatype, dataspace, H5P_DEFAULT));
        HDF5_ERROR (dataset = H5Dcreate (writer, datasetname.str().c_str(),
                                         filedatatype, dataspace, H5P_DEFAULT));
      }

      // Loop over all components
      int first_time = 1;
      BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype)
      {
        // Get the intersection of the current component with this combination
        // (use either the interior or exterior here, as we did above)
        ibbox const overlap = *bbox &
          arrdata.at(gindex).at(Carpet::map).dd->
          boxes.at(refinementlevel).at(component).at(mglevel).interior;

        // Continue if this component is not part of this combination
        if (overlap.empty())
        {
          continue;
        }

        // Copy the overlap to the local processor
        const ggf<dim>* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
        const gdata<dim>* const data = (*ff) (-request->timelevel,
                                              refinementlevel,
                                              component, mglevel);
        gdata<dim>* const processor_component =
          data->make_typed (request->vindex);

        void *h5data;
        if (group.grouptype == CCTK_SCALAR)
        {
          h5data = CCTK_VarDataPtrI (cctkGH, request->timelevel,
                                     request->vindex);
        }
        else
        {
          processor_component->allocate (overlap, 0);
          for (comm_state<dim> state; !state.done(); state.step())
          {
            processor_component->copy_from (state, data, overlap);
          }
          h5data = (void *) processor_component->storage();
        }

        // Write data
        if (myproc == 0)
        {
          hsize_t overlapshape[dim];
          hssize_t overlaporigin[dim];
          for (int d = 0; d < group.dim; ++d)
          {
            overlaporigin[group.dim-1-d] =
             ((overlap.lower() - bbox->lower()) / overlap.stride())[d];
            overlapshape[group.dim-1-d]  = (overlap.shape() / overlap.stride())[d];
          }

          // Write the processor component as a hyperslab into the recombined
          // dataset.
          hid_t overlap_dataspace;
          HDF5_ERROR (overlap_dataspace =
                      H5Screate_simple (group.dim, overlapshape, NULL));
          HDF5_ERROR (H5Sselect_hyperslab (dataspace, H5S_SELECT_SET,
                                           overlaporigin, NULL,
                                           overlapshape, NULL));
          HDF5_ERROR (H5Dwrite (memdataset, memdatatype, overlap_dataspace,
                                dataspace, H5P_DEFAULT, h5data));
          HDF5_ERROR (H5Sclose (overlap_dataspace));

          // Add metadata information on the first time through
          // (have to do it inside of the COMPONENT_LOOP so that we have
          //  access to the cGH elements)
          if (first_time)
          {
            AddAttributes (cctkGH, fullname, group.dim, refinementlevel,
                           request->timelevel, bbox, dataset);
            first_time = 0;
          }
        }

        // Delete temporary copy of this component
        delete processor_component;

      } END_COMPONENT_LOOP;

      // Finally create the recombined dataset in the real HDF5 file on disk
      if (myproc == 0)
      {
        void *h5data = malloc (H5Dget_storage_size (memdataset));
        HDF5_ERROR (H5Dread (memdataset, filedatatype, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, h5data));
        HDF5_ERROR (H5Dclose (memdataset));
        HDF5_ERROR (H5Fclose (memfile));
        HDF5_ERROR (H5Dwrite (dataset, filedatatype, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, h5data));
        free (h5data);
        HDF5_ERROR (H5Sclose (dataspace));
        HDF5_ERROR (H5Dclose (dataset));
      }

    } // for bboxes
  } END_MAP_LOOP;

  HDF5_ERROR (H5Pclose (plist));

  free (fullname);

  return 0;
}


static void AddAttributes (const cGH *const cctkGH, const char *fullname,
                           int vdim, int refinementlevel, int timelevel,
                           ibset::const_iterator bbox, hid_t dataset)
{
  DECLARE_CCTK_ARGUMENTS


  // Write FlexIO attributes
  WriteAttribute (dataset, "level", refinementlevel);

  CCTK_REAL origin[dim], delta[dim];
  for (int d = 0; d < vdim; ++d)
  {
    delta[d]  = CCTK_DELTA_SPACE(d);
    origin[d] = CCTK_ORIGIN_SPACE(d) + cctk_lbnd[d] * delta[d];

    // it's the interior bbox so we have to substract the ghostzones
    if (! cctk_bbox[2*d])
    {
      origin[d] += cctk_nghostzones[d] * delta[d];
    }
  }
  WriteAttribute (dataset, "origin", origin, vdim);
  WriteAttribute (dataset, "delta", delta, vdim);
#ifdef WRITE_ADDITIONAL_ATTRIBUTES
  WriteAttribute (dataset, "min_ext", min_ext, vdim);
  WriteAttribute (dataset, "max_ext", max_ext, vdim);
  WriteAttribute (dataset, "level_timestep", cctk_iteration / reflevelfact);
  WriteAttribute (dataset, "persistence", maxreflevelfact / reflevelfact);
#endif

  WriteAttribute (dataset, "time", cctk_time);
  WriteAttribute (dataset, "timestep", cctk_iteration);

#ifdef WRITE_ADDITIONAL_ATTRIBUTES
  const int time_refinement = reflevelfact;
  const vect<int, dim> spatial_refinement = reflevelfact;
  WriteAttribute (dataset, "time_refinement", time_refinement);
  WriteAttribute (dataset, "spatial_refinement", spatial_refinement, vdim);
  WriteAttribute (dataset, "grid_placement_refinement", spatial_refinement, vdim);
#endif

  vect<int, dim> iorigin = bbox->lower() / bbox->stride();
  WriteAttribute (dataset, "iorigin", &iorigin[0], vdim);

  // Legacy arguments
  WriteAttribute (dataset, "name", fullname);

#ifdef WRITE_ADDITIONAL_ATTRIBUTES
  // Group arguments
  const int gindex = CCTK_GroupIndexFromVar (fullname);

  WriteAttribute (dataset, "group_version", 1);
  WriteAttribute (dataset, "group_fullname", fullname);
  WriteAttribute (dataset, "group_varname",
                  CCTK_VarName (CCTK_VarIndex (fullname)));
  {
    char * groupname = CCTK_GroupName (gindex);
    assert (groupname);
    WriteAttribute (dataset, "group_groupname", groupname);
    free (groupname);
  }
  const char *group_grouptype = NULL;
  switch (CCTK_GroupTypeI (gindex))
  {
    case CCTK_GF:     group_grouptype = "CCTK_GF"; break;
    case CCTK_ARRAY:  group_grouptype = "CCTK_ARRAY"; break;
    case CCTK_SCALAR: group_grouptype = "CCTK_SCALAR"; break;
  }
  assert (group_grouptype);
  WriteAttribute (dataset, "group_grouptype", group_grouptype);
  WriteAttribute (dataset, "group_dim", vdim);
#endif
  WriteAttribute (dataset, "group_timelevel", -timelevel);
#ifdef WRITE_ADDITIONAL_ATTRIBUTES
  WriteAttribute (dataset, "group_numtimelevels", CCTK_NumTimeLevelsI(gindex));

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
  WriteAttribute (dataset, "cctk_levfac", cctk_levfac, dim);
  WriteAttribute (dataset, "cctk_levoff", cctk_levoff, dim);
  WriteAttribute (dataset, "cctk_levoffdenom", cctk_levoffdenom, dim);
  WriteAttribute (dataset, "cctk_timefac", cctk_timefac);
  WriteAttribute (dataset, "cctk_convlevel", cctk_convlevel);
  WriteAttribute (dataset, "cctk_convfac", cctk_convfac);
  WriteAttribute (dataset, "cctk_time", cctk_time);
#endif
  WriteAttribute (dataset, "cctk_bbox", cctk_bbox, 2*vdim);
  WriteAttribute (dataset, "cctk_nghostzones", cctk_nghostzones, vdim);

  // Carpet arguments
#ifdef WRITE_ADDITIONAL_ATTRIBUTES
  WriteAttribute (dataset, "carpet_version", 1);
  WriteAttribute (dataset, "carpet_dim", dim);
  WriteAttribute (dataset, "carpet_basemglevel", basemglevel);
  WriteAttribute (dataset, "carpet_mglevels", mglevels);
  WriteAttribute (dataset, "carpet_mgface", mgfact);
  WriteAttribute (dataset, "carpet_reflevels", reflevels);
  WriteAttribute (dataset, "carpet_reffact", reffact);
  WriteAttribute (dataset, "carpet_map", Carpet::map);
  WriteAttribute (dataset, "carpet_maps", maps);
  WriteAttribute (dataset, "carpet_component", component);
  WriteAttribute (dataset, "carpet_components",
                  vhh.at(Carpet::map)->components(refinementlevel));
#endif
  WriteAttribute (dataset, "carpet_mglevel", mglevel);
  WriteAttribute (dataset, "carpet_reflevel", refinementlevel);
}


int TimeToOutput (const cGH* const cctkGH, const int vindex)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert (vindex>=0 && vindex<CCTK_NumVars());

  switch (CCTK_GroupTypeFromVarI (vindex))
  {
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      if (! do_global_mode)
      {
        return 0;
      }
      break;
    case CCTK_GF:
      // do nothing
      break;
    default:
      assert (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration;

  const char *myoutcriterion = out3D_criterion;
  if (CCTK_EQUALS(myoutcriterion, "default"))
  {
    myoutcriterion = out_criterion;
  }

  if (CCTK_EQUALS (myoutcriterion, "never"))
  {
    // Never output
    output_this_iteration = false;
  }
  else if (CCTK_EQUALS (myoutcriterion, "iteration"))
  {
    int myoutevery = out3D_every;
    if (myoutevery == -2)
    {
      myoutevery = out_every;
    }
    if (myoutevery <= 0)
    {
      // output is disabled
      output_this_iteration = false;
    }
    else if (cctk_iteration == *this_iteration)
    {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
    else if (cctk_iteration >= *next_output_iteration)
    {
      // it is time for the next output
      output_this_iteration = true;
      *next_output_iteration = cctk_iteration + myoutevery;
      *this_iteration = cctk_iteration;
    }
    else
    {
      // we want no output at this iteration
      output_this_iteration = false;
    }
  }
  else if (CCTK_EQUALS (myoutcriterion, "divisor"))
  {
    int myoutevery = out3D_every;
    if (myoutevery == -2)
    {
      myoutevery = out_every;
    }
    if (myoutevery <= 0)
    {
      // output is disabled
      output_this_iteration = false;
    }
    else if ((cctk_iteration % myoutevery) == 0)
    {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
    else
    {
      // we want no output at this iteration
      output_this_iteration = false;
    }
  }
  else if (CCTK_EQUALS (myoutcriterion, "time"))
  {
    CCTK_REAL myoutdt = out3D_dt;
    if (myoutdt == -2)
    {
      myoutdt = out_dt;
    }
    if (myoutdt < 0)
    {
      // output is disabled
      output_this_iteration = false;
    }
    else if (myoutdt == 0)
    {
      // output all iterations
      output_this_iteration = true;
    }
    else if (cctk_iteration == *this_iteration)
    {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
    else if (cctk_time / cctk_delta_time
               >= *next_output_time / cctk_delta_time - 1.0e-12)
    {
      // it is time for the next output
      output_this_iteration = true;
      *next_output_time = cctk_time + myoutdt;
      *this_iteration = cctk_iteration;
    }
    else
    {
      // we want no output at this iteration
      output_this_iteration = false;
    }
  }
  else
  {
    assert (0);
  } // select output criterion

  if (! output_this_iteration)
  {
    return 0;
  }

  if (! CheckForVariable(out3D_vars, vindex))
  {
    // This variable should not be output
    return 0;
  }

  if (last_output.at(mglevel).at(reflevel).at(vindex) == cctk_iteration)
  {
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


int TriggerOutput (const cGH* const cctkGH, const int vindex)
{
  assert (vindex>=0 && vindex<CCTK_NumVars());

  char *fullname = CCTK_FullName (vindex);
  const char *varname = CCTK_VarName (vindex);
  const int retval = OutputVarAs (cctkGH, fullname, varname);
  free (fullname);

  last_output.at(mglevel).at(reflevel).at(vindex) = cctkGH->cctk_iteration;

  return retval;
}


int ReadVar (const cGH* const cctkGH, const int vindex,
             const hid_t dataset, vector<ibset> &regions_read,
             const int called_from_recovery)
{
  DECLARE_CCTK_PARAMETERS;

  const int group = CCTK_GroupIndexFromVarI (vindex);
  assert (group>=0 && group<(int)Carpet::arrdata.size());
  char *fullname = CCTK_FullName (vindex);
  const int n0 = CCTK_FirstVarIndexI(group);
  assert (n0>=0 && n0<CCTK_NumVars());
  const int var = vindex - n0;
  assert (var>=0 && var<CCTK_NumVars());
  int tl = 0;

  bool did_read_something = false;

  // Stuff needed for Recovery

  void *h5data = NULL;

  if (verbose)
  {
    CCTK_VInfo (CCTK_THORNSTRING, "  reading '%s'", fullname);
  }

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, group))
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable \"%s\" because it has no storage",
                fullname);
    free (fullname);
    return 0;
  }

  const int grouptype = CCTK_GroupTypeI(group);
  if ((grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY) && reflevel > 0)
  {
    free (fullname);
    return 0;
  }

  const int gpdim = CCTK_GroupDimI(group);


  int intbuffer[2 + 2*dim];
  int &group_timelevel = intbuffer[0],
      &amr_level       = intbuffer[1];
  int *amr_origin      = &intbuffer[2],
      *amr_dims        = &intbuffer[2+dim];

  if (CCTK_MyProc(cctkGH)==0)
  {
    // get dataset dimensions
    hid_t dataspace;
    HDF5_ERROR (dataspace = H5Dget_space (dataset));
    int rank = (int) H5Sget_simple_extent_ndims (dataspace);
    assert (0 < rank && rank <= dim);
    assert ((grouptype == CCTK_SCALAR ? gpdim+1 : gpdim) == rank);
    vector<hsize_t> shape(rank);
    HDF5_ERROR (H5Sget_simple_extent_dims (dataspace, &shape[0], NULL));
    HDF5_ERROR (H5Sclose (dataspace));

    for (int i = 0; i < dim; i++)
    {
      amr_dims[i]   = 1;
      amr_origin[i] = 0;
    }
    int datalength = 1;
    for (int i = 0; i < rank; i++)
    {
      datalength *= shape[i];
      amr_dims[i] = shape[rank-i-1];
    }

    const int cctkDataType = CCTK_VarTypeI(vindex);
    const hid_t datatype = h5DataType(cctkGH,cctkDataType);

    //cout << "datalength: " << datalength << " rank: " << rank << "\n";
    //cout << shape[0] << " " << shape[1] << " " << shape[2] << "\n";

    // to do: read in an allocate with correct datatype

    h5data = malloc (CCTK_VarTypeSize (cctkDataType) * datalength);
    HDF5_ERROR (H5Dread (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         h5data));

    ReadAttribute (dataset, "level", amr_level);
    ReadAttribute (dataset, "iorigin", amr_origin, rank);

    if(called_from_recovery)
    {
      ReadAttribute(dataset,"group_timelevel", group_timelevel);
    }
  } // MyProc == 0

  MPI_Bcast (intbuffer, sizeof (intbuffer) / sizeof (*intbuffer), MPI_INT, 0, dist::comm);

#ifdef CARPETIOHDF5_DEBUG
  cout << "amr_level: " << amr_level << " reflevel: " << reflevel << endl;
#endif

  if (amr_level == reflevel)
  {
    // Traverse all components on all levels
    BEGIN_MAP_LOOP (cctkGH, grouptype)
    {
      BEGIN_COMPONENT_LOOP (cctkGH, grouptype)
      {
        did_read_something = true;

        ggf<dim>* ff = 0;

        assert (var < (int)arrdata.at(group).at(Carpet::map).data.size());
        ff = (ggf<dim>*)arrdata.at(group).at(Carpet::map).data.at(var);

        if(called_from_recovery) tl = group_timelevel;

        gdata<dim>* const data = (*ff) (tl, reflevel, component, mglevel);

        // Create temporary data storage on processor 0
        vect<int,dim> str = vect<int,dim>(maxreflevelfact/reflevelfact);

        if(grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY)
          str = vect<int,dim> (1);

        vect<int,dim> lb = vect<int,dim>::ref(amr_origin) * str;
        vect<int,dim> ub
          = lb + (vect<int,dim>::ref(amr_dims) - 1) * str;

        gdata<dim>* const tmp = data->make_typed (vindex);


        cGroup cgdata;
        int ierr = CCTK_GroupData(group,&cgdata);
        assert(ierr==0);
        //cout << "lb_before: " << lb << endl;
        //cout << "ub_before: " << ub << endl;
        if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
#ifdef CARPETIOHDF5_DEBUG
          cout << "CCTK_DISTRIB_CONSTANT: " << fullname << endl;
#endif
          assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
          if (grouptype == CCTK_SCALAR)
          {
            lb[0] = arrdata.at(group).at(Carpet::map).hh->processors.at(reflevel).at(component);
            ub[0] = arrdata.at(group).at(Carpet::map).hh->processors.at(reflevel).at(component);
            for(int i=1;i<dim;i++)
            {
              lb[i]=0;
              ub[i]=0;
            }
          }
          else
          {
            const int newlb = lb[gpdim-1] +
              (ub[gpdim-1]-lb[gpdim-1]+1)*
              (arrdata.at(group).at(Carpet::map).hh->processors.at(reflevel).at(component));
            const int newub = ub[gpdim-1] +
              (ub[gpdim-1]-lb[gpdim-1]+1)*
              (arrdata.at(group).at(Carpet::map).hh->processors.at(reflevel).at(component));
            lb[gpdim-1] = newlb;
            ub[gpdim-1] = newub;
          }
#ifdef CARPETIOHDF5_DEBUG
          cout << "lb: " << lb << endl;
          cout << "ub: " << ub << endl;
#endif
        }
        const bbox<int,dim> ext(lb,ub,str);

#ifdef CARPETIOHDF5_DEBUG
        cout << "ext: " << ext << endl;
#endif

        if (CCTK_MyProc(cctkGH)==0)
        {
          tmp->allocate (ext, 0, h5data);
        }
        else
        {
          tmp->allocate (ext, 0);
        }

        // Initialise with what is found in the file -- this does
        // not guarantee that everything is initialised.
        const bbox<int,dim> overlap = tmp->extent() & data->extent();
        regions_read.at(Carpet::map) |= overlap;

#ifdef CARPETIOHDF5_DEBUG
        cout << "working on component: " << component << endl;
        cout << "tmp->extent " << tmp->extent() << endl;
        cout << "data->extent " << data->extent() << endl;
        cout << "overlap " << overlap << endl;
        cout << "-----------------------------------------------------" << endl;
#endif

        // FIXME: is this barrier really necessary ??
        MPI_Barrier(MPI_COMM_WORLD);

        // Copy into grid function
        for (comm_state<dim> state; !state.done(); state.step())
        {
          data->copy_from (state, tmp, overlap);
        }


        // Delete temporary copy
        delete tmp;


      } END_COMPONENT_LOOP;

      if (called_from_recovery)
      {
        arrdata.at(group).at(Carpet::map).tt->set_time(reflevel,mglevel,
                (CCTK_REAL) ((cctkGH->cctk_time - cctk_initial_time)
                             / (delta_time * mglevelfact)) );
      }


    } END_MAP_LOOP;

  } // if amr_level == reflevel

  free (h5data);
  free (fullname);

  return did_read_something;
}


static int InputVarAs (const cGH* const cctkGH, const int vindex,
                       const char* const alias)
{
  DECLARE_CCTK_PARAMETERS;

  char *fullname = CCTK_FullName (vindex);
  const int group = CCTK_GroupIndexFromVarI (vindex);
  assert (group>=0 && group<(int)Carpet::arrdata.size());

  int want_dataset = 0;
  bool did_read_something = false;
  int ndatasets = 0;
  hid_t dataset = 0;

  char datasetname[1024];

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, group))
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable \"%s\" because it has no storage",
                fullname);
    free (fullname);
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

  // Read the file only on the root processor
  if (CCTK_MyProc(cctkGH)==0)
  {
    // Open the file
    if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Opening file \"%s\"", filename);
    reader = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (reader<0)
    {
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


  for (int datasetid=0; datasetid<ndatasets; ++datasetid)
  {
    if (verbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Handling dataset #%d", datasetid);
    }

    // Read data
    if (CCTK_MyProc(cctkGH)==0)
    {
      GetDatasetName(reader,datasetid,datasetname);
      //         cout << datasetname << "\n";

      HDF5_ERROR (dataset = H5Dopen (reader, datasetname));
    }


#if 0
    int amr_level;
    int amr_origin[dim];
    int amr_dims[dim];
#endif

    if (CCTK_MyProc(cctkGH)==0)
    {
      // Read data
      char * name;
      ReadAttribute (dataset, "name", name);
      //        cout << "dataset name is " << name << endl;
      if (verbose && name)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Dataset name is \"%s\"", name);
      }
      want_dataset = name && CCTK_EQUALS(name, fullname);
      free (name);
    } // myproc == 0

    MPI_Bcast (&want_dataset, 1, MPI_INT, 0, dist::comm);

    if(want_dataset)
    {
      did_read_something = ReadVar(cctkGH,vindex,dataset,regions_read,0);
    } // want_dataset

  } // loop over datasets

  // Close the file
  if (CCTK_MyProc(cctkGH)==0)
  {
    if (verbose) CCTK_VInfo (CCTK_THORNSTRING, "Closing file");
    HDF5_ERROR (H5Fclose(reader));
    reader=-1;
  }

  // Was everything initialised?
  if (did_read_something)
  {
    for (int m=0; m<Carpet::maps; ++m)
    {
      dh<dim>& thedd = *arrdata.at(group).at(m).dd;
      ibset all_exterior;
      for (size_t c=0; c<thedd.boxes.at(rl).size(); ++c)
      {
        all_exterior |= thedd.boxes.at(rl).at(c).at(mglevel).exterior;
      }
      if (regions_read.at(m) != all_exterior)
      {
#ifdef CARPETIOHDF5_DEBUG
        cout << "read: " << regions_read.at(m) << endl
             << "want: " << all_exterior << endl;
#endif
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Variable \"%s\" could not be initialised from file -- the file may be missing data",
                    fullname);
      }
    }
  } // if did_read_something
  //        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,"stop!");

  return did_read_something ? 0 : -1;
}


int CarpetIOHDF5ReadData (const cGH* const cctkGH)
{
  int retval = 0;
  DECLARE_CCTK_PARAMETERS


  int numvars = CCTK_NumVars ();
  vector<bool> flags (numvars);

  if (CCTK_TraverseString (in3D_vars, SetFlag, &flags, CCTK_GROUP_OR_VAR) < 0)
  {
    CCTK_VWarn (strict_io_parameter_check ? 0 : 1,
                __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing parameter 'IOHDF5::in3D_vars'");
  }

  for (int vindex = 0; vindex < numvars && retval == 0; vindex++)
  {
    if (flags.at (vindex))
    {
      retval = InputVarAs (cctkGH, vindex, CCTK_VarName (vindex));
    }
  }

  return (retval);
}



#if 0
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
      const int ierr = InputVarAs (cctkGH, vindex, basefilename);
      if (! ierr) {
        ++ num_vars_read;
      }
    }
  }

  return num_vars_read;
}
#endif

} // namespace CarpetIOHDF5
