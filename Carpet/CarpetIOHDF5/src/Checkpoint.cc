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
#include <string>
#include <vector>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Version.h"

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

// when was the last checkpoint written ?
static int last_checkpoint_iteration = -1;


static int Checkpoint (const cGH* const cctkGH, int called_from);
static int DumpParametersGHExtentions (const cGH *cctkGH, int all, hid_t writer);


int CarpetIOHDF5_InitialDataCheckpoint (const cGH* const cctkGH)
{
  DECLARE_CCTK_PARAMETERS


  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("---------------------------------------------------------");
    CCTK_INFO ("Dumping initial data checkpoint");
    CCTK_INFO ("---------------------------------------------------------");
  }
  int retval = Checkpoint (cctkGH, CP_INITIAL_DATA);

  return (retval);
}


int CarpetIOHDF5_EvolutionCheckpoint (const cGH* const cctkGH)
{
  int retval = 0;
  DECLARE_CCTK_PARAMETERS


  if (checkpoint &&
    ((checkpoint_every > 0 && cctkGH->cctk_iteration % checkpoint_every == 0) ||
     checkpoint_next))
  {
    if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("---------------------------------------------------------");
      CCTK_VInfo (CCTK_THORNSTRING, "Dumping periodic checkpoint at "
                  "iteration %d", cctkGH->cctk_iteration);
      CCTK_INFO ("---------------------------------------------------------");
    }

    retval = Checkpoint (cctkGH, CP_EVOLUTION_DATA);

    if (checkpoint_next)
    {
      CCTK_ParameterSet ("checkpoint_next", CCTK_THORNSTRING, "no");
    }
  }

  return (retval);
}


int CarpetIOHDF5_TerminationCheckpoint (const cGH *const GH)
{
  int retval = 0;
  DECLARE_CCTK_PARAMETERS


  if (checkpoint && checkpoint_on_terminate)
  {
    if (last_checkpoint_iteration < GH->cctk_iteration)
    {
      if (! CCTK_Equals (verbose, "none"))
      {
        CCTK_INFO ("---------------------------------------------------------");
        CCTK_VInfo (CCTK_THORNSTRING, "Dumping termination checkpoint at "
                    "iteration %d", GH->cctk_iteration);
        CCTK_INFO ("---------------------------------------------------------");
      }

      retval = Checkpoint (GH, CP_EVOLUTION_DATA);
    }
    else if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("---------------------------------------------------------");
      CCTK_VInfo (CCTK_THORNSTRING, "Termination checkpoint already dumped "
                  "as last evolution checkpoint at iteration %d",
                  last_checkpoint_iteration);
      CCTK_INFO ("---------------------------------------------------------");
    }
  }

  return (retval);
}


static int Checkpoint (const cGH* const cctkGH, int called_from)
{
  int retval = 0;
  DECLARE_CCTK_PARAMETERS


  CarpetIOHDF5GH *myGH = (CarpetIOHDF5GH *) CCTK_GHExtension (cctkGH,
                                                              "CarpetIOHDF5");
  // check if CarpetIOHDF5 was registered as I/O method
  if (myGH == NULL)
  {
    CCTK_WARN (0, "No CarpetIOHDF5 I/O methods registered");
    return -1;
  }

  int const myproc = CCTK_MyProc (cctkGH);


  /* get the filenames for both the temporary and real checkpoint file */
  char *filename = IOUtil_AssembleFilename (cctkGH, NULL, "", ".h5",
                                            called_from, 0, 1);
  char *tempname = IOUtil_AssembleFilename (cctkGH, NULL, ".tmp", ".h5",
                                            called_from, 0, 1);

  hid_t writer = -1;

  if (myproc == 0)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Creating temporary checkpoint file '%s'",
                  tempname);
    }

    HDF5_ERROR (writer = H5Fcreate (tempname, H5F_ACC_TRUNC, H5P_DEFAULT,
                                    H5P_DEFAULT));

    // Dump all parameters and GHExtentions
    retval = DumpParametersGHExtentions (cctkGH, 1, writer);

  }

  // now dump the grid variables on all mglevels, reflevels, maps and components
  BEGIN_MGLEVEL_LOOP (cctkGH)
  {
    BEGIN_REFLEVEL_LOOP (cctkGH)
    {
      if (CCTK_Equals (verbose, "full"))
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Dumping grid variables on mglevel %d reflevel %d ...",
                    mglevel, reflevel);
      }

      for (int group = CCTK_NumGroups () - 1; group >= 0; group--)
      {
        /* only dump groups which have storage assigned */
        if (CCTK_QueryGroupStorageI (cctkGH, group) <= 0)
        {
          continue;
        }

        cGroup gdata;
        CCTK_GroupData (group, &gdata);
        assert (gdata.grouptype == CCTK_ARRAY || gdata.grouptype == CCTK_GF ||
                gdata.grouptype == CCTK_SCALAR);

        // scalars and grid arrays only have one reflevel
        if (gdata.grouptype != CCTK_GF && reflevel > 0)
        {
          continue;
        }

        /* get the number of allocated timelevels */
        gdata.numtimelevels = 0;
        gdata.numtimelevels = CCTK_GroupStorageIncrease (cctkGH, 1, &group,
                                                         &gdata.numtimelevels,
                                                         NULL);

        int first_vindex = CCTK_FirstVarIndexI (group);

        /* get the default I/O request for this group */
        ioRequest *request = IOUtil_DefaultIORequest (cctkGH, first_vindex, 1);

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
          char *fullname = CCTK_FullName (request->vindex);
          assert (fullname);

          /* loop over all timelevels of this variable */
          for (request->timelevel = 0;
               request->timelevel < gdata.numtimelevels;
               request->timelevel++)
          {
            if (CCTK_Equals (verbose, "full"))
            {
              CCTK_VInfo (CCTK_THORNSTRING, "  %s (timelevel %d)",
                          fullname, request->timelevel);
            }

            // write the var
            retval += WriteVar (cctkGH, writer, request, 1);
          }
          free (fullname);

        } /* end of loop over all variables */

        // free I/O request structure
        IOUtil_FreeIORequest (&request);

      } /* end of loop over all groups */
    } END_REFLEVEL_LOOP;

  } END_MGLEVEL_LOOP;


  // Close the file
  if (writer >= 0)
  {
    HDF5_ERROR (H5Fclose(writer));
  }

  if (retval == 0 && myproc == 0)
  {
    if (rename (tempname, filename))
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not rename temporary checkpoint file '%s' to '%s'",
                  tempname, filename);
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
      myGH->cp_filename_list[myGH->cp_filename_index] = strdup (filename);
      myGH->cp_filename_index = (myGH->cp_filename_index+1) % abs (checkpoint_keep);
    }
  }

  // save the iteration number of this checkpoint
  last_checkpoint_iteration = cctkGH->cctk_iteration;

  // free allocated resources
  free (tempname);
  free (filename);

  return retval;

} // Checkpoint


static int DumpParametersGHExtentions (const cGH *cctkGH, int all, hid_t writer)
{
  // large parts of this routine were taken from
  // Thomas Radke's IOHDF5Util. Thanks Thomas!

  DECLARE_CCTK_PARAMETERS;

  char *parameters;
  hid_t group, dspace, dset;
  hsize_t size;

  int itmp;
  double dtmp;
  const char *version;
  const ioGH *ioUtilGH;

  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_INFO ("Dumping Parameters and GH Extentions...");
  }

  /* get the parameter string buffer */
  parameters = IOUtil_GetAllParameters (cctkGH, all);
  if (parameters)
  {
    size = strlen (parameters) + 1;
    HDF5_ERROR (group = H5Gcreate (writer, METADATA_GROUP, 0));
    HDF5_ERROR (dspace = H5Screate_simple (1, &size, NULL));
    HDF5_ERROR (dset = H5Dcreate (group, ALL_PARAMETERS, H5T_NATIVE_UCHAR,
                                  dspace, H5P_DEFAULT));
    HDF5_ERROR (H5Dwrite (dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, parameters));
    free (parameters);

    // now dump the GH Extentions

    /* get the handle for IOUtil extensions */
    ioUtilGH = (const ioGH *) CCTK_GHExtension (cctkGH, "IO");

    itmp = CCTK_MainLoopIndex ();
    WriteAttribute(dset,"main loop index",itmp);

    itmp = cctkGH->cctk_iteration;
    WriteAttribute(dset,"GH$iteration",itmp);

    itmp = ioUtilGH->ioproc_every;
    WriteAttribute(dset,"GH$ioproc_every",itmp);

    itmp = CCTK_nProcs (cctkGH);
    WriteAttribute(dset,"GH$nprocs",itmp);

    dtmp = cctkGH->cctk_time;
    WriteAttribute(dset,"GH$time", dtmp);

    dtmp = global_time;
    WriteAttribute(dset,"carpet_global_time", dtmp);

    itmp = reflevels;
    WriteAttribute(dset,"carpet_reflevels", itmp);

    dtmp = delta_time;
    WriteAttribute(dset,"carpet_delta_time", dtmp);

    version = CCTK_FullVersion();
    WriteAttribute(dset,"Cactus version", version);

    /* finally, we need all the times on the individual refinement levels */

    const int numberofmgtimes=mglevels;
    WriteAttribute(dset,"numberofmgtimes",numberofmgtimes);
    for(int i=0;i < numberofmgtimes;i++)
    {
      char buffer[100];
      snprintf (buffer, sizeof (buffer), "mgleveltimes %d", i);
      WriteAttribute(dset,buffer,(double *) &leveltimes.at(i).at(0), reflevels);
    }

    HDF5_ERROR (H5Dclose (dset));
    HDF5_ERROR (H5Sclose (dspace));
    HDF5_ERROR (H5Gclose (group));
  }

  return 0;
} // DumpParametersGHExtentions

} // namespace CarpetIOHDF5
