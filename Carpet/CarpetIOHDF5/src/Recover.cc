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

#include "CarpetIOHDF5.hh"

/* some macros for HDF5 group names */
#define METADATA_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"


namespace CarpetIOHDF5
{

using namespace std;
using namespace Carpet;

// structure describing a single dataset of an HDF5 file to read from
typedef struct
{
  char *datasetname;

  int vindex;
  int mglevel;
  int reflevel;
} dataset_t;

// structure describing the contents of an HDF5 file to read from
typedef struct
{
  int num_mglevels;
  int num_reflevels;
  int parameter_len;
  int cctk_iteration;
  int main_loop_index;
  CCTK_REAL global_time;
  CCTK_REAL delta_time;
  CCTK_REAL *mgleveltimes;  // [num_mglevels*num_reflevels]

  char *filename;
  hid_t file;
  int num_datasets;
  list<dataset_t> datasets;  // [num_datasets]
} file_t;

static file_t infile = {0, 0, 0, 0, 0, 0, 0, NULL, NULL, -1, -1,
                        list<dataset_t> ()};


static int OpenFile (const char *basefilename, file_t *file, int called_from);
static int RecoverVariables (cGH* cctkGH, file_t *file);
static herr_t ReadMetadata (hid_t group, const char *objectname, void *arg);

// callback for I/O parameter parsing routine
static void SetFlag (int vindex, const char* optstring, void* arg);


// Register with the Cactus Recovery Interface
int CarpetIOHDF5_RecoverParameters (void)
{
  int retval = IOUtil_RecoverParameters (Recover, ".h5", "HDF5");

  return (retval);
}


// close a checkpoint/filereader file after recovering grid variables
int CarpetIOHDF5_CloseFile (void)
{
  DECLARE_CCTK_PARAMETERS


  if (infile.num_datasets <= 0) {
    return (0);
  }
  infile.num_datasets = -1;

  if (CCTK_Equals (verbose, "full")) {
    CCTK_VInfo (CCTK_THORNSTRING, "closing file '%s' after recovery",
                infile.filename);
  }

  if (infile.file >= 0) {
    HDF5_ERROR (H5Fclose (infile.file));
    infile.file = -1;
  }
  free (infile.filename);
  delete[] infile.mgleveltimes;

  for (list<dataset_t>::iterator dataset = infile.datasets.begin ();
       dataset != infile.datasets.end ();
       dataset++) {
    free (dataset->datasetname);
  }
  infile.datasets.clear ();

  return (0);
}

static int OpenFile (const char *basefilename, file_t *file, int called_from)
{
  hid_t dset = -1;
  DECLARE_CCTK_PARAMETERS


  // generate filename for an unchunked checkpoint file */
  file->filename = IOUtil_AssembleFilename (NULL, basefilename, "", ".h5",
                                            called_from, 0, 1);
  if (CCTK_Equals (verbose, "full")) {
    CCTK_VInfo (CCTK_THORNSTRING, "opening %s file '%s'",
                called_from == CP_RECOVER_PARAMETERS ? "checkpoint" : "input",
                file->filename);
  }

  // try to open the file (prevent HDF5 error messages if it fails)
  H5E_BEGIN_TRY {
    file->file = H5Fopen (file->filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  } H5E_END_TRY;

  if (file->file >= 0) {
    if (called_from == CP_RECOVER_PARAMETERS) {
      HDF5_ERROR (dset = H5Dopen (file->file,
                                  METADATA_GROUP"/"ALL_PARAMETERS));
      ReadAttribute (dset, "carpet_reflevels", file->num_reflevels);
      ReadAttribute (dset, "numberofmgtimes", file->num_mglevels);
      ReadAttribute (dset, "GH$iteration", file->cctk_iteration);
      ReadAttribute (dset, "main loop index", file->main_loop_index);
      ReadAttribute (dset, "carpet_global_time", file->global_time);
      ReadAttribute (dset, "carpet_delta_time", file->delta_time);
      file->parameter_len = H5Dget_storage_size (dset) + 1;
      assert (file->parameter_len > 1);
    }

    HDF5_ERROR (H5Giterate (file->file, "/", NULL, ReadMetadata, file));
  }
  file->num_datasets = file->datasets.size ();
  if (file->num_datasets <= 0) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No valid HDF5 file '%s' found", file->filename);
    free (file->filename);
  }

  // return if no valid checkpoint could be found
  if (file->num_datasets <= 0) {
    return (-1);
  }

  if (called_from == FILEREADER_DATA) {
    return (0);
  }

  // leave space at the end for global_time and delta_time
  // so that all double variables can be broadcasted in one go
  int num_times = file->num_mglevels*file->num_reflevels + 2;
  file->mgleveltimes = new CCTK_REAL[num_times];

  // FIXME: should store all mgleveltimes in a single contiguous array
  //        to get rid of this loop and save some attributes
  for (int i = 0; i < file->num_mglevels; i++) {
    char buffer[32];

    snprintf (buffer, sizeof (buffer), "mgleveltimes %d", i);
    ReadAttribute (dset, buffer, file->mgleveltimes + i*file->num_reflevels,
                   file->num_reflevels);
  }

  // recover parameters
  char *parameters = new char[file->parameter_len];
  CCTK_VInfo (CCTK_THORNSTRING, "Recovering parameters from checkpoint");

  HDF5_ERROR (H5Dread (dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       parameters));
  HDF5_ERROR (H5Dclose (dset));
  IOUtil_SetAllParameters (parameters);
  delete[] parameters;

  return (0);
}


int Recover (cGH* cctkGH, const char *basefilename, int called_from)
{
  DECLARE_CCTK_PARAMETERS
  int retval = 0;


  assert (called_from == CP_RECOVER_PARAMETERS ||
          called_from == CP_RECOVER_DATA ||
          called_from == FILEREADER_DATA);

  if (called_from == CP_RECOVER_PARAMETERS ||
      called_from == FILEREADER_DATA)
  {
    // open the file, read and broadcast its metadata information
    // for CP_RECOVER_PARAMETERS: also recover all parameters
    retval = OpenFile (basefilename, &infile, called_from);

    if (called_from == CP_RECOVER_PARAMETERS || retval) {
      return (retval == 0 ? 1 : -1);
    }
  }

  // can only proceed with a valid checkpoint file from here on
  assert (infile.num_datasets > 0);

  // set global Cactus/Carpet variables
  if (called_from == CP_RECOVER_DATA) {
    global_time = infile.global_time;
    delta_time = infile.delta_time;
    CCTK_SetMainLoopIndex (infile.main_loop_index);

    cctkGH->cctk_iteration = infile.cctk_iteration;
    cctkGH->cctk_time = infile.mgleveltimes[mglevel*infile.num_reflevels +
                                            reflevel];
  }

  // now recover all grid variables on current mglevel and reflevel
  retval = RecoverVariables (cctkGH, &infile);

  if (called_from == CP_RECOVER_DATA) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "restarting simulation at iteration %d (physical time %g)",
                cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
  } else {
    // FIXME: keep filereader input file open for all mglevels, reflevels
    //        this doesn't work now because we don't know the number of
    //        mglevels and reflevels
    //        also: there may be multiple files to be read in, and they
    //              must not share the same data structures
    CarpetIOHDF5_CloseFile ();
  }

  return (retval);
}


int ReadVar (const cGH* const cctkGH, const int vindex,
             const hid_t dataset, vector<ibset> &regions_read,
             const int called_from_recovery)
{
  DECLARE_CCTK_PARAMETERS;

  const int gindex = CCTK_GroupIndexFromVarI (vindex);
  assert (gindex >= 0 && gindex < Carpet::arrdata.size());
  const int var = vindex - CCTK_FirstVarIndexI(gindex);
  cGroup group;
  CCTK_GroupData(gindex,&group);


  if (group.grouptype != CCTK_GF && reflevel > 0) {
    return 0;
  }

  int refinement_level;
  ReadAttribute (dataset, "level", refinement_level);
  if (refinement_level != reflevel) {
    return 0;
  }

  if (CCTK_Equals (verbose, "full")) {
    char *fullname = CCTK_FullName (vindex);
    CCTK_VInfo (CCTK_THORNSTRING, "  reading '%s'", fullname);
    free (fullname);
  }

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, gindex)) {
    char *fullname = CCTK_FullName (vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable '%s' (no storage)", fullname);
    free (fullname);
    return 0;
  }

  const hid_t datatype = h5DataType (cctkGH, group.vartype, 0);

  int timelevel = 0;
  if(called_from_recovery) {
    ReadAttribute(dataset,"group_timelevel", timelevel);
    // in old days (before February 2005)
    // timelevels were stored as negative numbers
    timelevel = abs (timelevel);
  }

  // get dataset dimensions
  hid_t filespace;
  HDF5_ERROR (filespace = H5Dget_space (dataset));
  int rank = (int) H5Sget_simple_extent_ndims (filespace);
  assert (group.dim == (group.grouptype == CCTK_SCALAR) ? rank-1 : rank);
  vector<hsize_t> h5shape(rank);
  HDF5_ERROR (H5Sget_simple_extent_dims (filespace, &h5shape[0], NULL));

  ivect shape(1);
  int elems = 1;
  for (int i = 0; i < rank; i++) {
    shape[i] = h5shape[rank-i-1];
    elems   *= shape[i];
  }

  ivect iorigin(0);
  ReadAttribute (dataset, "iorigin", &iorigin[0], rank);

  const ivect stride = group.grouptype == CCTK_GF ?
                       maxspacereflevelfact/spacereflevelfact : 1;
  ivect lower = iorigin * stride;
  ivect upper = lower + (shape - 1) * stride;

  // Traverse all local components on all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype) {
    struct arrdesc& data = arrdata.at(gindex).at(Carpet::map);

    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, group.grouptype) {
      ibbox& interior_membox =
        data.dd->boxes.at(mglevel).at(reflevel).at(component).interior;

      // reset bounds for DISTRIB=CONST variables
      if (group.disttype == CCTK_DISTRIB_CONSTANT) {
        assert (group.grouptype==CCTK_SCALAR || group.grouptype==CCTK_ARRAY);
        const int newlower = lower[rank-1] +
          (upper[rank-1]-lower[rank-1]+1)*
          (data.hh->processors().at(reflevel).at(component));
        const int newupper = upper[rank-1] +
          (upper[rank-1]-lower[rank-1]+1)*
          (data.hh->processors().at(reflevel).at(component));
        lower[rank-1] = newlower;
        upper[rank-1] = newupper;
      }

      const ibbox filebox(lower, upper, stride);
      const ibbox interior_overlap = interior_membox & filebox;

      // skip this dataset if it doesn't overlap with this component's interior
      if (interior_overlap.empty()) {
        continue;
      }

      // get the overlap with this component's exterior
      ibbox& membox =
        data.dd->boxes.at(mglevel).at(reflevel).at(component).exterior;
      const ibbox overlap = membox & filebox;

      // calculate hyperslab selection parameters
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 6 && H5_VERS_RELEASE == 4)
      hsize_t memorigin[dim], fileorigin[dim];
#else
      hssize_t memorigin[dim], fileorigin[dim];
#endif
      hsize_t memdims[dim], count[dim];
      for (int i = 0; i < rank; i++) {
        memdims[rank-i-1] = (membox.upper()-membox.lower())[i] / stride[i] + 1;
        count[rank-i-1] = (overlap.upper()-overlap.lower())[i] / stride[i] + 1;
        memorigin[rank-i-1] = (overlap.lower()-membox.lower())[i] / stride[i];
        fileorigin[rank-i-1] = group.disttype == CCTK_DISTRIB_CONSTANT ?
                               0 : (overlap.lower()-lower)[i] / stride[i];
      }
      // read the hyperslab
      gdata* const gv = (*data.data.at(var)) (timelevel, reflevel,
                                              component, mglevel);
      hid_t memspace;
      HDF5_ERROR (memspace  = H5Screate_simple (rank, memdims, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, fileorigin,
                                       NULL, count, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (memspace, H5S_SELECT_SET, memorigin,
                                       NULL, count, NULL));
      HDF5_ERROR (H5Dread (dataset, datatype, memspace, filespace,
                           H5P_DEFAULT, gv->storage()));
      HDF5_ERROR (H5Sclose (memspace));

      // Initialise with what is found in the file -- this does
      // not guarantee that everything is initialised.
      regions_read.at(Carpet::map) |= overlap;

    } END_LOCAL_COMPONENT_LOOP;

    if (called_from_recovery) {
      data.tt->set_time(reflevel, mglevel,
                        ((cctkGH->cctk_time - cctk_initial_time)
                         / (delta_time * mglevelfact)) );
    }

  } END_MAP_LOOP;

  HDF5_ERROR (H5Sclose (filespace));

  return 1;
}


static int InputVarAs (const cGH* const cctkGH, const int vindex,
                       const char* const alias)
{
  DECLARE_CCTK_PARAMETERS;

  char *fullname = CCTK_FullName (vindex);
  const int gindex = CCTK_GroupIndexFromVarI (vindex);
  assert (gindex>=0 && gindex<(int)Carpet::arrdata.size());

  int want_dataset = 0;
  bool did_read_something = false;
  int ndatasets = 0;
  hid_t dataset = 0;

  char datasetname[1024];

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, gindex))
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable \"%s\" because it has no storage",
                fullname);
    free (fullname);
    return 0;
  }

  const int grouptype = CCTK_GroupTypeI(gindex);
  const int rl = grouptype==CCTK_GF ? reflevel : 0;
  //cout << "want level " << rl << " reflevel " << reflevel << endl;

  // Find the input directory
  const char* const myindir = in_dir;

  // Invent a file name
  ostringstream filenamebuf;
  filenamebuf << myindir << "/" << alias << in_extension;
  string filenamestr = filenamebuf.str();
  const char * const filename = filenamestr.c_str();

  hid_t reader = -1;

  // Open the file
  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Opening file \"%s\"", filename);
  }
  reader = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (reader<0)
  {
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not open file \"%s\" for reading", filename);
  }
  assert (reader>=0);
  // get the number of datasets in the file
  ndatasets=GetnDatasets(reader);
  assert (ndatasets>=0);

  vector<ibset> regions_read(Carpet::maps);
  for (int datasetid=0; datasetid<ndatasets; ++datasetid)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Handling dataset #%d", datasetid);
    }

    // Read data
    GetDatasetName(reader,datasetid,datasetname);
    //         cout << datasetname << "\n";
    HDF5_ERROR (dataset = H5Dopen (reader, datasetname));

    // Read data
    char * name;
    ReadAttribute (dataset, "name", name);
    //        cout << "dataset name is " << name << endl;
    if (CCTK_Equals (verbose, "full") && name)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Dataset name is \"%s\"", name);
    }
    want_dataset = name && CCTK_EQUALS(name, fullname);
    free (name);

    if(want_dataset)
    {
      did_read_something = ReadVar(cctkGH,vindex,dataset,regions_read,0);
    } // want_dataset

  } // loop over datasets

  // Close the file
  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Closing file");
  }
  HDF5_ERROR (H5Fclose(reader));
  reader=-1;

  // Was everything initialised?
  if (did_read_something)
  {
    for (int m=0; m<Carpet::maps; ++m)
    {
      dh& thedd = *arrdata.at(gindex).at(m).dd;
      ibset all_exterior;
      for (size_t c=0; c<thedd.boxes.at(rl).size(); ++c)
      {
        all_exterior |= thedd.boxes.at(mglevel).at(rl).at(c).exterior;
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

  return did_read_something ? 0 : -1;
}


int CarpetIOHDF5_ReadData (const cGH* const cctkGH)
{
  int retval = 0;
  DECLARE_CCTK_PARAMETERS


  int numvars = CCTK_NumVars ();
  vector<bool> flags (numvars);

  if (CCTK_TraverseString (in_vars, SetFlag, &flags, CCTK_GROUP_OR_VAR) < 0)
  {
    CCTK_VWarn (strict_io_parameter_check ? 0 : 1,
                __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing parameter 'IOHDF5::in_vars'");
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


static herr_t ReadMetadata (hid_t group, const char *objectname, void *arg)
{
  file_t *file = (file_t *) arg;
  dataset_t dataset;
  hid_t dset, dspace;
  H5G_stat_t object_info;


  // we are interested only in datasets
  HDF5_ERROR (H5Gget_objinfo (group, objectname, 0, &object_info));
  if (object_info.type != H5G_DATASET)
  {
    return (0);
  }

  dataset.datasetname = strdup (objectname);
  assert (dataset.datasetname);

  HDF5_ERROR (dset = H5Dopen (group, objectname));
  char *varname = NULL;
  ReadAttribute (dset, "name", varname);
  dataset.vindex = CCTK_VarIndex (varname);
  free (varname);
  ReadAttribute (dset, "carpet_mglevel", dataset.mglevel);
  ReadAttribute (dset, "carpet_reflevel", dataset.reflevel);
  HDF5_ERROR (dspace = H5Dget_space (dset));
  HDF5_ERROR (H5Sclose (dspace));
  HDF5_ERROR (H5Dclose (dset));

  // add this dataset to our list and count the number of int elements
  file->datasets.push_back (dataset);

  return (0);
}


static int RecoverVariables (cGH* cctkGH, file_t *file)
{
  DECLARE_CCTK_PARAMETERS;


  // Use refinement levels parameter from checkpointing file ?
  if (use_reflevels_from_checkpoint)
  {
    char buffer[32];

    snprintf (buffer, sizeof (buffer), "%d", file->num_reflevels);
    CCTK_ParameterSet ("refinement_levels", "CarpetRegrid", buffer);

    CCTK_VInfo (CCTK_THORNSTRING, "Using %i reflevels read from checkpoint "
                "file. Ignoring value in parameter file.", file->num_reflevels);
  }

  CCTK_VInfo (CCTK_THORNSTRING,
              "reading grid variables on mglevel %d reflevel %d",
              mglevel, reflevel);

  int num_vars = CCTK_NumVars ();
  for (list<dataset_t>::iterator dataset = file->datasets.begin ();
       dataset != file->datasets.end ();
       dataset++)
  {
    // only recover grid variables for the current mglevel/reflevel
    if (dataset->mglevel != mglevel || dataset->reflevel != reflevel)
    {
      continue;
    }

    if (dataset->vindex <  0 || dataset->vindex >= num_vars)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Ignoring dataset '%s' (invalid variable name)",
                  dataset->datasetname);
      continue;
    }

    hid_t dset;
    HDF5_ERROR (dset = H5Dopen (file->file, dataset->datasetname));
    assert (dset >= 0);

    vector<ibset> regions_read (Carpet::maps);
    ReadVar (cctkGH, dataset->vindex, dset, regions_read, 1);
    HDF5_ERROR (H5Dclose (dset));
  }

  return (0);
}


static void SetFlag (int vindex, const char* optstring, void* arg)
{
  if (optstring)
  {
    char *fullname = CCTK_FullName (vindex);
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Option string '%s' will be ignored for HDF5 input of "
                "variable '%s'", optstring, fullname);
    free (fullname);
  }
  vector<bool>& flags = *(vector<bool>*)arg;
  flags.at(vindex) = true;
}


} // namespace CarpetIOHDF5
