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
  static const char* rcsid = "$Header:$";
  CCTK_FILEVERSION(Carpet_CarpetIOHDF5_Output_cc);
}

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "CarpetIOHDF5.hh"


namespace CarpetIOHDF5
{

using namespace std;
using namespace Carpet;


// Variable definitions
vector<bool> do_truncate;     // [var]
vector<vector<vector<int> > > last_output; // [ml][rl][var]


// registered GH extension setup routine
static void* SetupGH (tFleshConfig* const fleshconfig,
                      const int convLevel, cGH* const cctkGH);

// callbacks for CarpetIOHDF5's I/O method
static int OutputGH (const cGH* const cctkGH);
static int OutputVarAs (const cGH* const cctkGH, const char* const varname,
                        const char* const alias);
static int TimeToOutput (const cGH* const cctkGH, const int vindex);
static int TriggerOutput (const cGH* const cctkGH, const int vindex);

// add attributes to an HDF5 dataset
static void AddAttributes (const cGH *const cctkGH, const char *fullname,
                           int vdim, int refinementlevel, int timelevel,
                           ibset::const_iterator bbox, hid_t dataset);

// callback for I/O parameter parsing routine
static void SetFlag (int vindex, const char* optstring, void* arg);


int CarpetIOHDF5_Startup (void)
{
  DECLARE_CCTK_PARAMETERS


  CCTK_RegisterBanner ("AMR HDF5 I/O provided by CarpetIOHDF5");

  // initial I/O parameter check
  vector<bool> flags (CCTK_NumVars ());

  if (CCTK_TraverseString (out3D_vars, SetFlag, &flags, CCTK_GROUP_OR_VAR) < 0)
  {
    CCTK_VWarn (strict_io_parameter_check ? 0 : 1,
                __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing parameter 'IOHDF5::out3D_vars'");
  }

  const int GHExtension = CCTK_RegisterGHExtension ("CarpetIOHDF5");
  CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);

  const int IOMethod = CCTK_RegisterIOMethod ("IOHDF5");
  CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);

  // register CarpetIOHDF5's recovery routine
  if (IOUtil_RegisterRecover ("CarpetIOHDF5 recovery", Recover) < 0)
  {
    CCTK_WARN (1, "Failed to register " CCTK_THORNSTRING " recovery routine");
  }

  return (0);
}


static void* SetupGH (tFleshConfig* const fleshconfig,
                      const int convLevel, cGH* const cctkGH)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CarpetIOHDF5GH* myGH;

  const void *dummy;
  dummy = &fleshconfig;
  dummy = &convLevel;
  dummy = &cctkGH;
  dummy = &dummy;

  const int numvars = CCTK_NumVars ();

  // Truncate all files if this is not a restart
  do_truncate.resize(numvars, true);

  // No iterations have yet been output
  last_output.resize(mglevels);
  for (int ml=0; ml<mglevels; ++ml)
  {
    last_output.at(ml).resize(maxreflevels);
    for (int rl=0; rl<maxreflevels; ++rl)
    {
      last_output.at(ml).at(rl).resize(numvars, INT_MIN);
    }
  }

  // allocate a new GH extension structure
  myGH            = (CarpetIOHDF5GH*) malloc (sizeof (CarpetIOHDF5GH));
  myGH->out_last  = (int *) malloc (numvars * sizeof (int));
  myGH->requests  = (ioRequest **) calloc (numvars, sizeof (ioRequest *));
  myGH->cp_filename_list = (char **) calloc (abs (checkpoint_keep), sizeof (char *));
  myGH->cp_filename_index = 0;
  myGH->out_vars = strdup ("");
  myGH->out_every_default = out_every - 1;

  for (int i = 0; i < numvars; i++)
  {
    myGH->out_last[i] = -1;
  }

  // Now set hdf5 complex datatypes
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


int CarpetIOHDF5_Init (const cGH* const cctkGH)
{
  DECLARE_CCTK_ARGUMENTS

  *this_iteration = -1;
  *next_output_iteration = 0;
  *next_output_time = cctk_time;

  return (0);
}


static int OutputGH (const cGH* const cctkGH)
{
  for (int vindex = CCTK_NumVars () - 1; vindex >= 0; vindex--)
  {
    if (TimeToOutput (cctkGH, vindex))
    {
      TriggerOutput (cctkGH, vindex);
    }
  }

  return (0);
}


static int TimeToOutput (const cGH* const cctkGH, const int vindex)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const int numvars = CCTK_NumVars();
  assert (vindex>=0 && vindex<numvars);

  if (CCTK_GroupTypeFromVarI (vindex) != CCTK_GF && ! do_global_mode)
  {
    return 0;
  }

  const char *myoutcriterion = CCTK_EQUALS (out3D_criterion, "default") ?
                               out_criterion : out3D_criterion;

  if (CCTK_EQUALS (myoutcriterion, "never"))
  {
    return (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration;

  if (CCTK_EQUALS (myoutcriterion, "iteration"))
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

  vector<bool> flags(numvars, false);

  CCTK_TraverseString (out3D_vars, SetFlag, &flags, CCTK_GROUP_OR_VAR);

  if (! flags.at(vindex))
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


static int TriggerOutput (const cGH* const cctkGH, const int vindex)
{
  char *fullname = CCTK_FullName (vindex);
  const char *varname = CCTK_VarName (vindex);
  const int retval = OutputVarAs (cctkGH, fullname, varname);
  free (fullname);

  last_output.at(mglevel).at(reflevel).at(vindex) = cctkGH->cctk_iteration;

  return (retval);
}


static int OutputVarAs (const cGH* const cctkGH, const char* const fullname,
                        const char* const alias)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int numvars = CCTK_NumVars ();
  vector<bool> flags (numvars, false);

  if (CCTK_TraverseString (fullname, SetFlag, &flags, CCTK_VAR) < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing variable name '%s' (alias name '%s')",
                fullname, alias);
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
                fullname);
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
                fullname, mglevel, reflevel);
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
  HDF5_ERROR (memdatatype = h5DataType (cctkGH, group.vartype, 0));
  HDF5_ERROR (filedatatype =
              h5DataType (cctkGH, group.vartype,
                          out_single_precision && ! called_from_checkpoint));

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

  vect<int, dim> iorigin = bbox->lower() / bbox->stride();
  WriteAttribute (dataset, "iorigin", &iorigin[0], vdim);

  WriteAttribute (dataset, "time", cctk_time);
  WriteAttribute (dataset, "timestep", cctk_iteration);

  // Legacy arguments
  WriteAttribute (dataset, "name", fullname);
  WriteAttribute (dataset, "group_timelevel", -timelevel);

  // Cactus arguments
  WriteAttribute (dataset, "cctk_bbox", cctk_bbox, 2*vdim);
  WriteAttribute (dataset, "cctk_nghostzones", cctk_nghostzones, vdim);

  // Carpet arguments
  WriteAttribute (dataset, "carpet_mglevel", mglevel);
  WriteAttribute (dataset, "carpet_reflevel", refinementlevel);
}


static void SetFlag (int vindex, const char* optstring, void* arg)
{
  if (optstring)
  {
    char *fullname = CCTK_FullName (vindex);
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Option string '%s' will be ignored for HDF5 output of "
                "variable '%s'", optstring, fullname);
    free (fullname);
  }
  vector<bool>& flags = *(vector<bool>*)arg;
  flags.at(vindex) = true;
}


} // namespace CarpetIOHDF5
