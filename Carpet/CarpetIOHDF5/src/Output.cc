#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "util_Table.h"

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
                           int vdim, int refinementlevel,
                           const ioRequest* request,
                           const ibbox& bbox, hid_t dataset);

// callback for I/O parameter parsing routine
static void GetVarIndex (int vindex, const char* optstring, void* arg);

static void CheckSteerableParameters (const cGH *const cctkGH,
                                      CarpetIOHDF5GH *myGH);
static int WarnAboutDeprecatedParameters (void);

//////////////////////////////////////////////////////////////////////////////
// public routines
//////////////////////////////////////////////////////////////////////////////

int CarpetIOHDF5_Startup (void)
{
  CCTK_RegisterBanner ("AMR HDF5 I/O provided by CarpetIOHDF5");

  const int GHExtension = CCTK_RegisterGHExtension ("CarpetIOHDF5");
  CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);

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

  // register CarpetIOHDF5's routines as a new I/O method
  const int IOMethod = CCTK_RegisterIOMethod ("IOHDF5");
  CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);

  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("I/O Method 'IOHDF5' registered: AMR output of grid variables "
               "to HDF5 files");
  }

  // register CarpetIOHDF5's recovery routine with IOUtil
  if (IOUtil_RegisterRecover ("CarpetIOHDF5 recovery", Recover) < 0)
  {
    CCTK_WARN (1, "Failed to register " CCTK_THORNSTRING " recovery routine");
  }

  const int numvars = CCTK_NumVars ();

  // allocate a new GH extension structure
  myGH            = (CarpetIOHDF5GH*) malloc (sizeof (CarpetIOHDF5GH));
  myGH->out_last  = (int *) malloc (numvars * sizeof (int));
  myGH->requests  = (ioRequest **) calloc (numvars, sizeof (ioRequest *));
  myGH->cp_filename_list = (char **) calloc (abs (checkpoint_keep), sizeof (char *));
  myGH->cp_filename_index = 0;
  myGH->out_vars = strdup ("");
  myGH->out_every_default = out_every - 1;

  // initial I/O parameter check
  myGH->stop_on_parse_errors = strict_io_parameter_check;
  CheckSteerableParameters (cctkGH, myGH);
  myGH->stop_on_parse_errors = 0;

  for (int i = 0; i < numvars; i++)
  {
    myGH->out_last[i] = -1;
  }

  // create the output directory (if it doesn't match ".")
  const char *my_out_dir = *out_dir ? out_dir : io_out_dir;
  myGH->out_dir = (char *) calloc (1, strlen (my_out_dir) + 2);
  if (strcmp (myGH->out_dir, "."))
  {
    sprintf (myGH->out_dir, "%s/", my_out_dir);
    if (CCTK_MyProc (cctkGH) == 0)
    {
      int result = CCTK_CreateDirectory (0755, myGH->out_dir);
      if (result < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Problem creating CarpetIOHDF5 output directory '%s'",
                    myGH->out_dir);
      }
      else if (result > 0 && CCTK_Equals (verbose, "full"))
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "CarpetIOHDF5 output directory '%s' already exists",
                    myGH->out_dir);
      }
    }
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

  return (myGH);
}


int CarpetIOHDF5_Init (const cGH* const cctkGH)
{
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = -1;
  *next_output_iteration = 0;
  *next_output_time = cctk_time;

  return (0);
}


int WriteVarChunked (const cGH* const cctkGH,
                     const hid_t writer,
                     const ioRequest* request,
                     const int called_from_checkpoint)
{
  DECLARE_CCTK_PARAMETERS;

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

  // Traverse all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype)
  {
    BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype)
    {
      // Using "exterior" includes ghost zones and refinement boundaries.
      ibbox& bbox = arrdata.at(gindex).at(Carpet::map).dd->
                    boxes.at(mglevel).at(refinementlevel).at(component).exterior;

      // Get the shape of the HDF5 dataset (in Fortran index order)
      hsize_t shape[dim];
      hsize_t num_elems = 1;
      for (int d = 0; d < group.dim; ++d)
      {
        shape[group.dim-1-d] = (bbox.shape() / bbox.stride())[d];
        num_elems *= shape[group.dim-1-d];
      }

      // Don't create zero-sized components
      if (num_elems == 0)
      {
        continue;
      }

      // Copy the overlap to the local processor
      const ggf* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
      const gdata* const data = (*ff) (request->timelevel, refinementlevel,
                                       component, mglevel);
      gdata* const processor_component = data->make_typed (request->vindex);

      processor_component->allocate (bbox, 0);
      for (comm_state state(group.vartype); !state.done(); state.step())
      {
        processor_component->copy_from (state, data, bbox);
      }

      // Write data on I/O processor 0
      if (myproc == 0)
      {
        const void *data = (const void *) processor_component->storage();

        // As per Cactus convention, DISTRIB=CONSTANT arrays
        // (including grid scalars) are assumed to be the same on
        // all processors and are therefore stored only by processor 0.
        //
        // Warn the user if this convention is violated.
        if (component > 0 && group.disttype == CCTK_DISTRIB_CONSTANT)
        {
          const void *proc0 = CCTK_VarDataPtrI (cctkGH, request->timelevel,
                                                request->vindex);

          if (memcmp (proc0, data, num_elems*CCTK_VarTypeSize(group.vartype)))
          {
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "values for DISTRIB=CONSTANT grid variable '%s' "
                        "(timelevel %d) differ between processors 0 and %d; "
                        "only the array from processor 0 will be stored",
                        fullname, request->timelevel, component);
          }
        }
        else
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
          if (group.grouptype == CCTK_GF)
          {
            if (Carpet::maps > 1)
            {
              datasetname << " m="  << Carpet::map;
            }
            datasetname << " rl=" << refinementlevel;
          }
          if (arrdata.at(gindex).at(Carpet::map).dd->
              boxes.at(mglevel).at(refinementlevel).size () > 1 &&
              group.disttype != CCTK_DISTRIB_CONSTANT)
          {
            datasetname << " c=" << component;
          }

          // remove an already existing dataset of the same name
          if (request->check_exist)
          {
            H5E_BEGIN_TRY
            {
              H5Gunlink (writer, datasetname.str().c_str());
            } H5E_END_TRY;
          }

          hsize_t shape[dim];
          hssize_t origin[dim];
          for (int d = 0; d < group.dim; ++d)
          {
            origin[group.dim-1-d] = (bbox.lower() / bbox.stride())[d];
            shape[group.dim-1-d]  = (bbox.shape() / bbox.stride())[d];
          }

          // Write the component as an individual dataset
          hid_t dataspace, dataset;
          HDF5_ERROR (dataspace = H5Screate_simple (group.dim, shape, NULL));
          HDF5_ERROR (dataset = H5Dcreate (writer, datasetname.str().c_str(),
                                           filedatatype, dataspace, H5P_DEFAULT));
          HDF5_ERROR (H5Sclose (dataspace));
          HDF5_ERROR (H5Dwrite (dataset, memdatatype, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, data));
          AddAttributes (cctkGH, fullname, group.dim, refinementlevel,
                         request, bbox, dataset);
          HDF5_ERROR (H5Dclose (dataset));
        }

      } // if myproc == 0

      // Delete temporary copy of this component
      delete processor_component;

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  free (fullname);

  return 0;
}


int WriteVarUnchunked (const cGH* const cctkGH,
                       const hid_t writer,
                       const ioRequest* request,
                       const int called_from_checkpoint)
{
  DECLARE_CCTK_PARAMETERS;

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
      // Using "interior" removes ghost zones and refinement boundaries.
      bboxes += arrdata.at(gindex).at(Carpet::map).dd->
                boxes.at(mglevel).at(refinementlevel).at(component).interior;
    } END_COMPONENT_LOOP;

    // Normalise the set, i.e., try to represent the set with fewer bboxes.
    //
    // According to Cactus conventions, DISTRIB=CONSTANT arrays
    // (including grid scalars) are assumed to be the same on all
    // processors and are therefore stored only by processor 0.
    if (group.disttype != CCTK_DISTRIB_CONSTANT)
    {
      bboxes.normalize();
    }

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
      if (num_elems == 0)
      {
        continue;
      }

      // create the dataset on the I/O processor
      // skip DISTRIB=CONSTANT components from processors other than 0
      hid_t memfile = -1, memdataset = -1;
      hid_t dataspace = -1, dataset = -1;
      if (myproc == 0 &&
          (bbox_id == 0 || group.disttype != CCTK_DISTRIB_CONSTANT))
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
        if (group.grouptype == CCTK_GF)
        {
          if (Carpet::maps > 1)
          {
            datasetname << " m="  << Carpet::map;
          }
          datasetname << " rl=" << refinementlevel;
        }
        if (bboxes.setsize () > 1 && group.disttype != CCTK_DISTRIB_CONSTANT)
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

        // remove an already existing dataset of the same name
        if (request->check_exist)
        {
          H5E_BEGIN_TRY
          {
            H5Gunlink (writer, datasetname.str().c_str());
          } H5E_END_TRY;
        }
        HDF5_ERROR (dataset = H5Dcreate (writer, datasetname.str().c_str(),
                                         filedatatype, dataspace, H5P_DEFAULT));
      }

      // Loop over all components
      bool first_time = true;
      BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype)
      {
        // Get the intersection of the current component with this combination
        // (use either the interior or exterior here, as we did above)
        ibbox const overlap = *bbox &
          arrdata.at(gindex).at(Carpet::map).dd->
          boxes.at(mglevel).at(refinementlevel).at(component).interior;

        // Continue if this component is not part of this combination
        if (overlap.empty())
        {
          continue;
        }

        // Copy the overlap to the local processor
        const ggf* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
        const gdata* const data = (*ff) (request->timelevel,
                                         refinementlevel,
                                         component, mglevel);
        gdata* const processor_component = data->make_typed (request->vindex);

        processor_component->allocate (overlap, 0);
        for (comm_state state(group.vartype); !state.done(); state.step())
        {
          processor_component->copy_from (state, data, overlap);
        }

        // Write data
        if (myproc == 0)
        {
          const void *data = (const void *) processor_component->storage();

          // As per Cactus convention, DISTRIB=CONSTANT arrays
          // (including grid scalars) are assumed to be the same on
          // all processors and are therefore stored only by processor
          // 0.
          //
          // Warn the user if this convention is violated.
          if (bbox_id > 0 && group.disttype == CCTK_DISTRIB_CONSTANT)
          {
            const void *proc0 = CCTK_VarDataPtrI (cctkGH, request->timelevel,
                                                  request->vindex);

            if (memcmp (proc0, data, num_elems*CCTK_VarTypeSize(group.vartype)))
            {
              CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "values for DISTRIB=CONSTANT grid variable '%s' "
                          "(timelevel %d) differ between processors 0 and %d; "
                          "only the array from processor 0 will be stored",
                          fullname, request->timelevel, bbox_id);
            }
          }
          else
          {
            hsize_t overlapshape[dim];
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 6 && H5_VERS_RELEASE == 4)
            hsize_t overlaporigin[dim];
#else
            hssize_t overlaporigin[dim];
#endif
            for (int d = 0; d < group.dim; ++d)
            {
              overlaporigin[group.dim-1-d] =
                ((overlap.lower() - bbox->lower()) / overlap.stride())[d];
              overlapshape[group.dim-1-d]  =
                (overlap.shape() / overlap.stride())[d];
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
                                  dataspace, H5P_DEFAULT, data));
            HDF5_ERROR (H5Sclose (overlap_dataspace));

            // Add metadata information on the first time through
            // (have to do it inside of the COMPONENT_LOOP so that we have
            //  access to the cGH elements)
            if (first_time)
            {
              AddAttributes (cctkGH, fullname, group.dim, refinementlevel,
                             request, *bbox, dataset);
              first_time = false;
            }
          }
        }

        // Delete temporary copy of this component
        delete processor_component;

      } END_COMPONENT_LOOP;

      // Finally create the recombined dataset in the real HDF5 file on disk
      // (skip DISTRIB=CONSTANT components from processors other than 0)
      if (myproc == 0 &&
          (bbox_id == 0 || group.disttype != CCTK_DISTRIB_CONSTANT))
      {
        void *data = malloc (H5Dget_storage_size (memdataset));
        HDF5_ERROR (H5Dread (memdataset, filedatatype, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, data));
        HDF5_ERROR (H5Dclose (memdataset));
        HDF5_ERROR (H5Fclose (memfile));
        HDF5_ERROR (H5Dwrite (dataset, filedatatype, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, data));
        free (data);
        HDF5_ERROR (H5Sclose (dataspace));
        HDF5_ERROR (H5Dclose (dataset));
      }

    } // for bboxes
  } END_MAP_LOOP;

  HDF5_ERROR (H5Pclose (plist));

  free (fullname);

  return 0;
}


//////////////////////////////////////////////////////////////////////////////
// private routines
//////////////////////////////////////////////////////////////////////////////
static void CheckSteerableParameters (const cGH *const cctkGH,
                                      CarpetIOHDF5GH *myGH)
{
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOHDF5::out_vars' parameter if it has changed
  if (strcmp (out_vars, myGH->out_vars))
  {
    IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING, "IOHDF5::out_vars",
                               myGH->stop_on_parse_errors, out_vars,
                               -1, myGH->requests);

    // notify the user about the new setting
    if (! CCTK_Equals (verbose, "none"))
    {
      char *msg = NULL;
      for (int i = CCTK_NumVars () - 1; i >= 0; i--)
      {
        if (myGH->requests[i])
        {
          char *fullname = CCTK_FullName (i);
          if (! msg)
          {
            Util_asprintf (&msg, "Periodic HDF5 output requested for '%s'",
                           fullname);
          }
          else
          {
            Util_asprintf (&msg, "%s, '%s'", msg, fullname);
          }
          free (fullname);
        }
      }
      if (msg)
      {
        CCTK_INFO (msg);
        free (msg);
      }
    }

    // save the last setting of 'IOHDF5::out_vars' parameter
    free (myGH->out_vars);
    myGH->out_vars = strdup (out_vars);
  }
}


static int OutputGH (const cGH* const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  static bool first_time = true;

  // check if any deprecated parameters have been set in the parameter file
  //  (don't check after recovery though)
  if (first_time)
  {

    if (CCTK_Equals (recover, "no") || ! *recover_file)
    {
      if (WarnAboutDeprecatedParameters ())
      {
#if 0
        // annoy the user if (s)he still used deprecated parameters
        CCTK_WARN (1, "Now waiting 5 seconds to let your notice the "
                      "warning message(s) above...");
        sleep (5);
#endif
      }
    }
    first_time = false;
  }

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
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int numvars = CCTK_NumVars();
  assert (vindex>=0 && vindex<numvars);

  if (CCTK_GroupTypeFromVarI (vindex) != CCTK_GF && ! do_global_mode)
  {
    return 0;
  }

  CarpetIOHDF5GH *myGH =
    (CarpetIOHDF5GH *) CCTK_GHExtension (cctkGH, CCTK_THORNSTRING);
  CheckSteerableParameters (cctkGH, myGH);

  // check if output for this variable was requested
  if (! myGH->requests[vindex])
  {
    return (0);
  }

  // check whether this refinement level should be output
  if (! (myGH->requests[vindex]->refinement_levels & (1 << reflevel)))
  {
    return (0);
  }

  // check if output for this variable was requested individually
  // by a "<varname>{ out_every = <number> }" option string
  // this will overwrite the output criterion setting
  const char *myoutcriterion = CCTK_EQUALS (out_criterion, "default") ?
                               io_out_criterion : out_criterion;
  if (myGH->requests[vindex]->out_every >= 0)
  {
    myoutcriterion = "divisor";
  }

  if (CCTK_EQUALS (myoutcriterion, "never"))
  {
    return (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration = false;

  if (CCTK_EQUALS (myoutcriterion, "iteration"))
  {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myoutevery > 0)
    {
      if (*this_iteration == cctk_iteration)
      {
        // we already decided to output this iteration
        output_this_iteration = true;
      }
      else if (cctk_iteration >= *next_output_iteration)
      {
        // it is time for the next output
        output_this_iteration = true;
        *this_iteration = cctk_iteration;
        *next_output_iteration = cctk_iteration + myoutevery;
      }
    }
  }
  else if (CCTK_EQUALS (myoutcriterion, "divisor"))
  {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myGH->requests[vindex]->out_every >= 0)
    {
      myoutevery = myGH->requests[vindex]->out_every;
    }
    if (myoutevery > 0 && (cctk_iteration % myoutevery) == 0)
    {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
  }
  else if (CCTK_EQUALS (myoutcriterion, "time"))
  {
    CCTK_REAL myoutdt = out_dt == -2 ? io_out_dt : out_dt;
    if (myoutdt == 0 || *this_iteration == cctk_iteration)
    {
      output_this_iteration = true;
    }
    else if (myoutdt > 0 && (cctk_time / cctk_delta_time
                             >= *next_output_time / cctk_delta_time - 1.0e-12))
    {
      // it is time for the next output
      output_this_iteration = true;
      *this_iteration = cctk_iteration;
      *next_output_time = cctk_time + myoutdt;
    }
  }

  if (! output_this_iteration)
  {
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


static void GetVarIndex (int vindex, const char* optstring, void* arg)
{
  if (optstring)
  {
    char *fullname = CCTK_FullName (vindex);
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Option string '%s' will be ignored for HDF5 output of "
                "variable '%s'", optstring, fullname);
    free (fullname);
  }

  *((int *) arg) = vindex;
}


static int OutputVarAs (const cGH* const cctkGH, const char* const fullname,
                        const char* const alias)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vindex = -1;

  if (CCTK_TraverseString (fullname, GetVarIndex, &vindex, CCTK_VAR) < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing variable name '%s' (alias name '%s')",
                fullname, alias);
    return (-1);
  }

  if (vindex < 0)
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

  // Invent a file name
  const CarpetIOHDF5GH *myGH =
    (CarpetIOHDF5GH *) CCTK_GHExtension (cctkGH, CCTK_THORNSTRING);
  ostringstream filenamebuf;
  filenamebuf << myGH->out_dir << alias << out_extension;
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

  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "Writing variable '%s' on mglevel %d reflevel %d",
                fullname, mglevel, reflevel);
  }
  if (out_unchunked)
  {
    WriteVarUnchunked (cctkGH, writer, request, 0);
  }
  else
  {
    WriteVarChunked (cctkGH, writer, request, 0);
  }

  // Close the file
  if (writer >= 0)
  {
    HDF5_ERROR (H5Fclose (writer));
  }

  // Don't truncate again
  do_truncate.at(vindex) = false;

  return (0);
}


static void AddAttributes (const cGH *const cctkGH, const char *fullname,
                           int vdim, int refinementlevel,
                           const ioRequest* request,
                           const ibbox& bbox, hid_t dataset)
{
  DECLARE_CCTK_ARGUMENTS;


  // write bbox attributes if we have coordinate system info
  CCTK_REAL origin[dim], delta[dim];
  char *groupname = CCTK_GroupNameFromVarI (request->vindex);
  int coord_system_handle = Coord_GroupSystem (cctkGH, groupname);
  free (groupname);

  CCTK_INT coord_handles[dim];
  if (coord_system_handle >= 0 &&
      Util_TableGetIntArray (coord_system_handle, vdim,
                             coord_handles, "COORDINATES") >= 0)
  {
    const ibbox& baseext =
    vdd.at(Carpet::map)->bases.at(mglevel).at(reflevel).exterior;

    const ivect pos = (bbox.lower() - baseext.lower()) / bbox.stride();

    for (int d = 0; d < vdim; d++)
    {
      Util_TableGetReal (coord_handles[d], &origin[d], "COMPMIN");
      Util_TableGetReal (coord_handles[d], &delta[d], "DELTA");
      delta[d]  /= cctk_levfac[d];
      origin[d] += delta[d] * (cctk_levoff[d] / cctk_levoffdenom[d] + pos[d]);
    }

    WriteAttribute (dataset, "origin", origin, vdim);
    WriteAttribute (dataset, "delta", delta, vdim);
  }

  // Write FlexIO attributes
  WriteAttribute (dataset, "level", refinementlevel);

  vect<int, dim> iorigin = bbox.lower() / bbox.stride();
  WriteAttribute (dataset, "iorigin", &iorigin[0], vdim);

  WriteAttribute (dataset, "time", cctk_time);
  WriteAttribute (dataset, "timestep", cctk_iteration);

  // Legacy arguments
  WriteAttribute (dataset, "name", fullname);
  WriteAttribute (dataset, "group_timelevel", request->timelevel);

#if 0
  // FIXME TR: output bbox and nghostzones again for chunked output
  // Cactus arguments
  WriteAttribute (dataset, "cctk_bbox", cctk_bbox, 2*vdim);
  WriteAttribute (dataset, "cctk_nghostzones", cctk_nghostzones, vdim);
#endif

  // Carpet arguments
  WriteAttribute (dataset, "carpet_version", CARPET_VERSION);
  WriteAttribute (dataset, "carpet_mglevel", mglevel);
  WriteAttribute (dataset, "carpet_reflevel", refinementlevel);

  // Simulation arguments
  if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
    char const * const job_id
      = static_cast<char const *> (UniqueSimulationID (cctkGH));
    WriteAttribute (dataset, "simulation id", job_id);
  }
}


static int WarnAboutDeprecatedParameters (void)
{
  DECLARE_CCTK_PARAMETERS;
  int warnings = 0;
  char buffer[20];

  if (CCTK_ParameterQueryTimesSet ("out3D_dir", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("out_dir", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::out3D_dir' is deprecated, please use "
                  "'IOHDF5::out_dir' instead");
    CCTK_ParameterSet ("out_dir", CCTK_THORNSTRING, out3D_dir);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("out3D_vars", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("out_vars", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::out3D_vars' is deprecated, please use "
                  "'IOHDF5::out_vars' instead");
    CCTK_ParameterSet ("out_vars", CCTK_THORNSTRING, out3D_vars);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("out3D_extension", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("out_extension", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::out3D_extension' is deprecated, please use "
                  "'IOHDF5::out_extension' instead");
    CCTK_ParameterSet ("out_extension", CCTK_THORNSTRING, out3D_extension);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("out3D_criterion", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("out_criterion", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::out3D_criterion' is deprecated, please use "
                  "'IOHDF5::out_criterion' instead");
    CCTK_ParameterSet ("out_criterion", CCTK_THORNSTRING, out3D_criterion);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("out3D_every", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("out_every", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::out3D_every' is deprecated, please use "
                  "'IOHDF5::out_every' instead");
    snprintf (buffer, sizeof (buffer), "%d", out3D_every);
    CCTK_ParameterSet ("out_every", CCTK_THORNSTRING, buffer);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("out3D_dt", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("out_dt", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::out3D_dt' is deprecated, please use "
                  "'IOHDF5::out_dt' instead");
    snprintf (buffer, sizeof (buffer), "%f", (double)out3D_dt);
    CCTK_ParameterSet ("out_dt", CCTK_THORNSTRING, buffer);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("in3D_dir", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("in_dir", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::in3D_dir' is deprecated, please use "
                  "'IOHDF5::in_dir' instead");
    CCTK_ParameterSet ("in_dir", CCTK_THORNSTRING, in3D_dir);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("in3D_vars", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("in_vars", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::in3D_vars' is deprecated, please use "
                  "'IOHDF5::in_vars' instead");
    CCTK_ParameterSet ("in_vars", CCTK_THORNSTRING, in3D_vars);
    warnings++;
  }
  if (CCTK_ParameterQueryTimesSet ("in3D_extension", CCTK_THORNSTRING) >
      CCTK_ParameterQueryTimesSet ("in_extension", CCTK_THORNSTRING))
  {
    CCTK_WARN (2, "Parameter 'IOHDF5::in3D_extension' is deprecated, please use "
                  "'IOHDF5::in_extension' instead");
    CCTK_ParameterSet ("in_extension", CCTK_THORNSTRING, in3D_extension);
    warnings++;
  }

  return (warnings);
}


} // namespace CarpetIOHDF5
