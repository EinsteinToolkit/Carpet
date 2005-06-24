#include <assert.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CarpetIOHDF5.hh"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"


namespace CarpetIOHDF5
{

using namespace std;
using namespace Carpet;

// structure describing a patch as a single dataset of an HDF5 file to read from
typedef struct {
  string patchname;
  int vindex;
  int mglevel;
  int reflevel;
  int timestep;
  int timelevel;
  int rank;
  ivect iorigin;
} patch_t;

// structure describing the contents of an HDF5 file to read from
typedef struct {
  char *filename;
  hid_t file;
  list<patch_t> patches;
} file_t;

// structure describing a set of HDF5 files
typedef struct {
  string setname;
  string basefilename;
  int first_ioproc;
  vector<file_t> files;            // [nioprocs]

  int nioprocs;
  int num_mglevels;
  int num_reflevels;
  int cctk_iteration;
  int main_loop_index;
  CCTK_REAL global_time;
  CCTK_REAL delta_time;
  vector<CCTK_REAL> mgleveltimes;  // [num_mglevels*num_reflevels]
} fileset_t;

static list<fileset_t> filesets;


static list<fileset_t>::iterator OpenFileSet (const cGH* const cctkGH,
                                              const string setname,
                                              const char *basefilename,
                                              int called_from);
static void ReadMetadata (fileset_t& fileset, hid_t file);
static herr_t BrowseDatasets (hid_t group, const char *objectname, void *arg);
static int ReadVar (const cGH* const cctkGH,
                    hid_t file,
                    list<patch_t>::const_iterator patch,
                    vector<ibset> &bboxes_read,
                    bool in_recovery);


//////////////////////////////////////////////////////////////////////////////
// Register with the Cactus Recovery Interface
//////////////////////////////////////////////////////////////////////////////
int CarpetIOHDF5_RecoverParameters (void)
{
  int retval = IOUtil_RecoverParameters (Recover, ".h5", "HDF5");

  return (retval);
}


//////////////////////////////////////////////////////////////////////////////
// close all open checkpoint/filereader files after recovering grid variables
//////////////////////////////////////////////////////////////////////////////
int CarpetIOHDF5_CloseFiles (void)
{
  DECLARE_CCTK_PARAMETERS;


  for (list<fileset_t>::const_iterator set = filesets.begin();
       set != filesets.end(); set++) {
    for (int i = 0; i < set->files.size(); i++) {
      if (set->files[i].file >= 0) {

        if (CCTK_Equals (verbose, "full")) {
          CCTK_VInfo (CCTK_THORNSTRING, "closing file '%s' after recovery",
                      set->files[i].filename);
        }

        free (set->files[i].filename);
        HDF5_ERROR (H5Fclose (set->files[i].file));
      }
    }
  }

  filesets.clear();

  return (0);
}


//////////////////////////////////////////////////////////////////////////////
// Top-level recovery/filereader routine
// Opens the checkpoint/data file on the first time through,
// recovers parameters and variables
//////////////////////////////////////////////////////////////////////////////
int Recover (cGH* cctkGH, const char *basefilename, int called_from)
{
  DECLARE_CCTK_PARAMETERS;


  assert (called_from == CP_RECOVER_PARAMETERS or
          called_from == CP_RECOVER_DATA or
          called_from == FILEREADER_DATA);
  const bool in_recovery = called_from == CP_RECOVER_PARAMETERS or
                           called_from == CP_RECOVER_DATA;

  const char* const dir = in_recovery ? recover_dir : filereader_ID_dir;
  string setname (dir);
  if (called_from == FILEREADER_DATA) {
    setname.append (basefilename);
  }

  list<fileset_t>::iterator fileset = filesets.begin();
  while (fileset != filesets.end()) {
    if (fileset->setname == setname) {
      break;
    }
    fileset++;
  }
  if (fileset == filesets.end()) {
    assert (called_from == CP_RECOVER_PARAMETERS or
            called_from == FILEREADER_DATA);

    fileset = OpenFileSet (cctkGH, setname, basefilename, called_from);
    if (fileset == filesets.end()) {
      return (-1);
    }

    if (called_from == CP_RECOVER_PARAMETERS) {
      // return here if only parameters should be recovered
      // (come back later with called_from == CP_RECOVER_DATA)
      return (1);
    }
  }

  // set global Cactus/Carpet variables
  if (in_recovery) {
    global_time = fileset->global_time;
    delta_time = fileset->delta_time;
    CCTK_SetMainLoopIndex (fileset->main_loop_index);

    cctkGH->cctk_iteration = fileset->cctk_iteration;
    cctkGH->cctk_time =
      fileset->mgleveltimes.at(mglevel*fileset->num_reflevels + reflevel);
  }

  if (not CCTK_Equals (verbose, "none")) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "reading grid variables on mglevel %d reflevel %d",
                mglevel, reflevel);
  }

  // create a bbox set for each variable to mark how much has been read already
  const int numvars = CCTK_NumVars ();
  vector<bool> read_completely(numvars, false);
  vector<vector<ibset> > bboxes_read (numvars);
  for (int i = 0; i < bboxes_read.size(); i++) {
    bboxes_read[i].resize (Carpet::maps);
  }

  // create a bbox set for each group to list how much needs to be read
  const int numgroups = CCTK_NumGroups ();
  vector<vector<ibset> > group_bboxes (numgroups);
  for (int gindex = 0; gindex < group_bboxes.size(); gindex++) {
    group_bboxes[gindex].resize (Carpet::maps);
    const int grouptype = CCTK_GroupTypeI (gindex);
    BEGIN_MAP_LOOP (cctkGH, grouptype) {
      struct arrdesc& data = arrdata.at(gindex).at(Carpet::map);

      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, grouptype) {
        if (grouptype == CCTK_GF or (mglevel == 0 and reflevel == 0)) {
          group_bboxes[gindex][Carpet::map] |=
            data.dd->boxes.at(mglevel).at(reflevel).at(component).exterior;
        }
      } END_LOCAL_COMPONENT_LOOP;
    } END_MAP_LOOP;
  }

  const ioGH* ioUtilGH = (const ioGH*) CCTK_GHExtension (cctkGH, "IO");

  // loop over all input files of this set
  for (int i = 0; i < fileset->files.size(); i++) {
    const int file_idx = (i + fileset->first_ioproc) % fileset->nioprocs;
    file_t& file = fileset->files[file_idx];

    // open the file (if it hasn't been already) and read its metadata
    if (file.file < 0) {
      file.filename =
        IOUtil_AssembleFilename (NULL, fileset->basefilename.c_str(),
                                 "", ".h5", called_from, file_idx, 0);
      assert (file.filename);
      HDF5_ERROR (file.file = H5Fopen (file.filename, H5F_ACC_RDONLY,
                                       H5P_DEFAULT));

      // browse through all datasets contained in this file
      HDF5_ERROR (H5Giterate (file.file, "/", NULL, BrowseDatasets, &file));
    }
    assert (file.patches.size() > 0);

    // now loop over all patches contained in this file
    for (list<patch_t>::const_iterator patch = file.patches.begin();
         patch != file.patches.end();
         patch++) {

      // only recover grid variables for the current mglevel/reflevel
      if (patch->mglevel != mglevel or patch->reflevel != reflevel) {
        continue;
      }

      // When called by the filereader, IOUtil will set the 'do_inVars' vector
      // if the parameter 'IO::filereader_ID_vars' has been set.
      // In this parameter the user can specify individual variables to be read
      // as well as (optionally) a specific timestep.
      //
      // If 'IO::filereader_ID_vars' is left empty, 'do_inVars' will be NULL
      // which implicitely selects all variables for input.
      // If 'do_inVars[vindex]' is 0, variable 'vindex' was not selected.
      // If 'do_inVars[vindex]' is positive, a specific timestep of 'vindex'
      // was selected.
      // If 'do_inVars[vindex]' is negative, 'vindex' was selected with no
      // specific timestep; all timesteps will be read, the last one wins.
      if (ioUtilGH->do_inVars and
          ioUtilGH->do_inVars[patch->vindex] >= 0 and
          ioUtilGH->do_inVars[patch->vindex] != patch->timestep + 1) {
        continue;
      }

      // actually read the patch
      if (not read_completely[patch->vindex]) {
        ReadVar (cctkGH, file.file, patch, bboxes_read[patch->vindex],
                 in_recovery);
      }
    }

    // check if all variables have been read completely already
    int all_read_completely = true;
    for (int vindex = 0; vindex < read_completely.size(); vindex++) {
      if (not read_completely[vindex]) {
        const int gindex = CCTK_GroupIndexFromVarI (vindex);
        const int gtype = CCTK_GroupTypeI (gindex);
        if (gtype == CCTK_GF or (reflevel == 0 and mglevel == 0)) {
          read_completely[vindex] = bboxes_read[vindex] == group_bboxes[gindex];
          all_read_completely &= read_completely[vindex];
        }
      }
    }
    if (all_read_completely) {
      break;
    }
  }

  // check that all variables have been read completely on this mglevel/reflevel
  int num_incomplete = 0;
  for (int vindex = 0; vindex < read_completely.size(); vindex++) {
    if (not read_completely[vindex]) {
      // check if the variable has been read partially
      size_t size = 0;
      for (int map = 0; map < Carpet::maps; map++) {
        size += bboxes_read[vindex][map].size();
      }
      char* fullname = CCTK_FullName (vindex);
      if (size == 0) {
        CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "variable '%s' has not been read", fullname);
      } else {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "variable '%s' has been read only partially", fullname);
        num_incomplete++;
      }
      free (fullname);
    }
  }
  if (num_incomplete) {
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%d variables on mglevel %d reflevel %d have been read "
                "only partially", num_incomplete, mglevel, reflevel);
  }

  if (in_recovery and not CCTK_Equals (verbose, "none")) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "restarting simulation on mglevel %d reflevel %d "
                "at iteration %d (physical time %g)", mglevel, reflevel,
                cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
  }

  return (0);
}


//////////////////////////////////////////////////////////////////////////////
// Open a set of checkpoint/filereader files
//////////////////////////////////////////////////////////////////////////////
static list<fileset_t>::iterator OpenFileSet (const cGH* const cctkGH,
                                              const string setname,
                                              const char *basefilename,
                                              int called_from)
{
  file_t file;
  fileset_t fileset;
  DECLARE_CCTK_PARAMETERS;

  // first try to open a chunked file written on this processor
  // (note that dist::rank() cannot be called yet during RECOVER_PARAMETERS)
  fileset.first_ioproc = CCTK_MyProc (cctkGH);
  file.filename = IOUtil_AssembleFilename (NULL, basefilename, "", ".h5",
                                           called_from, fileset.first_ioproc, 0);

  // try to open the file (prevent HDF5 error messages if it fails)
  H5E_BEGIN_TRY {
    file.file = H5Fopen (file.filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  } H5E_END_TRY;

  // if that failed, try a chunked file written on processor 0
  // (which always is an I/O proc)
  if (file.file < 0) {
    free (file.filename);
    fileset.first_ioproc = 0;
    file.filename = IOUtil_AssembleFilename (NULL, basefilename, "", ".h5",
                                             called_from, fileset.first_ioproc, 0);
    H5E_BEGIN_TRY {
      file.file = H5Fopen (file.filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    } H5E_END_TRY;
  }

  // if that still failed, try an unchunked file
  // (which is always written on processor 0)
  if (file.file < 0) {
    free (file.filename);
    file.filename = IOUtil_AssembleFilename (NULL, basefilename, "", ".h5",
                                             called_from, fileset.first_ioproc, 1);
    H5E_BEGIN_TRY {
      file.file = H5Fopen (file.filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    } H5E_END_TRY;
  }

  // return if no valid checkpoint could be found
  if (file.file < 0) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No valid HDF5 file with basename '%s' found", basefilename);
    free (file.filename);
    return (filesets.end());
  }

  if (CCTK_Equals (verbose, "full")) {
    CCTK_VInfo (CCTK_THORNSTRING, "opening %s file '%s'",
                called_from == CP_RECOVER_PARAMETERS ?
                "checkpoint" : "input", file.filename);
  }

  // read all the metadata information
  ReadMetadata (fileset, file.file);

  // browse through all datasets contained in this file
  HDF5_ERROR (H5Giterate (file.file, "/", NULL, BrowseDatasets, &file));
  assert (file.patches.size() > 0);

  // recover parameters
  if (called_from == CP_RECOVER_PARAMETERS) {
    if (not CCTK_Equals (verbose, "none")) {
      CCTK_VInfo (CCTK_THORNSTRING, "Recovering parameters from checkpoint "
                  "file '%s'", file.filename);
    }
    hid_t dataset;
    HDF5_ERROR (dataset = H5Dopen (file.file,
                                   METADATA_GROUP "/" ALL_PARAMETERS));
    char *parameters = new char[H5Dget_storage_size (dataset) + 1];
    HDF5_ERROR (H5Dread (dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, parameters));
    HDF5_ERROR (H5Dclose (dataset));
    IOUtil_SetAllParameters (parameters);
    delete[] parameters;

    // use refinement levels parameter from checkpointing file ?
    if (use_reflevels_from_checkpoint) {
      char buffer[32];

      snprintf (buffer, sizeof (buffer), "%d", fileset.num_reflevels);
      CCTK_ParameterSet ("refinement_levels", "CarpetRegrid", buffer);

      if (not CCTK_Equals (verbose, "none")) {
        CCTK_VInfo (CCTK_THORNSTRING, "Using %i reflevels read from "
                    "checkpoint file. Ignoring value in parameter file.",
                    fileset.num_reflevels);
      }
    }
  }

  // allocate and initialise the list of input files for this set
  fileset.files.resize (fileset.nioprocs);
  for (int i = 0; i < fileset.nioprocs; i++) {
    fileset.files[i].file = -1;
  }
  fileset.files[fileset.first_ioproc] = file;
  fileset.setname = setname;
  fileset.basefilename = basefilename;

  // add the new fileset to the list of sets and return an iterator to it
  filesets.push_front (fileset);
  return (filesets.begin());
}


//////////////////////////////////////////////////////////////////////////////
// Read the metadata for a set of input files
//////////////////////////////////////////////////////////////////////////////
static void ReadMetadata (fileset_t& fileset, hid_t file)
{
  DECLARE_CCTK_PARAMETERS;

  fileset.nioprocs = 1;
  hid_t metadata, attr, dataspace;
  H5E_BEGIN_TRY {
    metadata = H5Gopen (file, METADATA_GROUP);
  } H5E_END_TRY;
  if (metadata < 0) {
    // no metadata at all - this must be old-fashioned data output file
    return;
  }

  // in old times (before June 2005) the metadata attributes were attached
  // to the parameters dataset which was stored in checkpoint files only
  H5E_BEGIN_TRY {
    attr = H5Aopen_name (metadata, "nioprocs");
  } H5E_END_TRY;
  const bool is_old_fashioned_file = attr < 0;
  if (is_old_fashioned_file) {
    HDF5_ERROR (H5Gclose (metadata));
    HDF5_ERROR (metadata = H5Dopen (file,
                                    METADATA_GROUP "/" ALL_PARAMETERS));
  } else {
    HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &fileset.nioprocs));
    HDF5_ERROR (H5Aclose (attr));
  }
  HDF5_ERROR (attr = H5Aopen_name (metadata, "carpet_reflevels"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &fileset.num_reflevels));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (metadata, "numberofmgtimes"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &fileset.num_mglevels));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (metadata, "GH$iteration"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &fileset.cctk_iteration));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (metadata, "main loop index"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &fileset.main_loop_index));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (metadata, "carpet_global_time"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, &fileset.global_time));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (metadata, "carpet_delta_time"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, &fileset.delta_time));
  HDF5_ERROR (H5Aclose (attr));

  fileset.mgleveltimes.resize (fileset.num_mglevels * fileset.num_reflevels);
  for (int i = 0; i < fileset.num_mglevels; i++) {
    char buffer[32];

    snprintf (buffer, sizeof (buffer), "mgleveltimes %d", i);
    HDF5_ERROR (attr = H5Aopen_name (metadata, buffer));
    HDF5_ERROR (dataspace = H5Aget_space (attr));
    assert (H5Sget_simple_extent_npoints (dataspace) == fileset.num_reflevels);
    HDF5_ERROR (H5Sclose (dataspace));
    HDF5_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE,
                         &fileset.mgleveltimes[i * fileset.num_reflevels]));
    HDF5_ERROR (H5Aclose (attr));
  }

  if (is_old_fashioned_file) {
    HDF5_ERROR (H5Dclose (metadata));
  } else {
    HDF5_ERROR (H5Gclose (metadata));
  }
}


//////////////////////////////////////////////////////////////////////////////
// Browses through all the datasets in a checkpoint/data file
// and stores this information in a list of patch_t objects
//////////////////////////////////////////////////////////////////////////////
static herr_t BrowseDatasets (hid_t group, const char *objectname, void *arg)
{
  file_t *file = (file_t *) arg;
  patch_t patch;
  hid_t dataset, dataspace, attr, attrtype;
  H5G_stat_t object_info;


  // skip everything in the metadata group
  if (strcmp (objectname, METADATA_GROUP) == 0) {
    return (0);
  }
  // we are interested only in datasets
  HDF5_ERROR (H5Gget_objinfo (group, objectname, 0, &object_info));
  if (object_info.type != H5G_DATASET) {
    return (0);
  }

  HDF5_ERROR (dataset = H5Dopen (group, objectname));
  HDF5_ERROR (dataspace = H5Dget_space (dataset));
  patch.rank = H5Sget_simple_extent_ndims (dataspace);
  assert (patch.rank > 0 and patch.rank <= ::dim);
  HDF5_ERROR (H5Sclose (dataspace));
  HDF5_ERROR (attr = H5Aopen_name (dataset, "name"));
  HDF5_ERROR (attrtype = H5Aget_type (attr));
  size_t length = H5Tget_size (attrtype);
  assert (length > 0);
  vector<char> varname(length + 1);
  HDF5_ERROR (H5Aread (attr, attrtype, &varname[0]));
  HDF5_ERROR (H5Tclose (attrtype));
  HDF5_ERROR (H5Aclose (attr));
  patch.vindex = CCTK_VarIndex (&varname[0]);
  HDF5_ERROR (attr = H5Aopen_name (dataset, "carpet_mglevel"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &patch.mglevel));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (dataset, "level"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &patch.reflevel));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (dataset, "timestep"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &patch.timestep));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Aopen_name (dataset, "group_timelevel"));
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &patch.timelevel));
  HDF5_ERROR (H5Aclose (attr));
  // in old days (before February 2005)
  // timelevels were stored as negative numbers
  patch.timelevel = abs (patch.timelevel);
  HDF5_ERROR (attr = H5Aopen_name (dataset, "iorigin"));
  HDF5_ERROR (dataspace = H5Aget_space (attr));
  assert (H5Sget_simple_extent_npoints (dataspace) == patch.rank);
  HDF5_ERROR (H5Sclose (dataspace));
  patch.iorigin = 0;
  HDF5_ERROR (H5Aread (attr, H5T_NATIVE_INT, &patch.iorigin[0]));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (H5Dclose (dataset));

  // add this patch to our list
  if (patch.vindex >=0 and patch.vindex < CCTK_NumVars ()) {
    patch.patchname = objectname;
    file->patches.push_back (patch);
  } else {
    CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Ignoring dataset '%s' (invalid variable name)", objectname);
  }

  return (0);
}


//////////////////////////////////////////////////////////////////////////////
// Reads a single grid variable from a checkpoint/data file
//////////////////////////////////////////////////////////////////////////////
static int ReadVar (const cGH* const cctkGH,
                    hid_t file,
                    list<patch_t>::const_iterator patch,
                    vector<ibset> &bboxes_read,
                    bool in_recovery)
{
  DECLARE_CCTK_PARAMETERS;

  const int gindex = CCTK_GroupIndexFromVarI (patch->vindex);
  assert (gindex >= 0 and gindex < Carpet::arrdata.size());
  const int var = patch->vindex - CCTK_FirstVarIndexI(gindex);
  cGroup group;
  CCTK_GroupData (gindex, &group);


  // grid arrays and scalars are read once on the coarsest reflevel
  assert (group.grouptype == CCTK_GF or (mglevel == 0 and reflevel == 0));

  if (CCTK_Equals (verbose, "full")) {
    char *fullname = CCTK_FullName (patch->vindex);
    CCTK_VInfo (CCTK_THORNSTRING, "  reading '%s'", fullname);
    free (fullname);
  }

  // Check for storage
  if (not CCTK_QueryGroupStorageI(cctkGH, gindex)) {
    char *fullname = CCTK_FullName (patch->vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable '%s' (no storage)", fullname);
    free (fullname);
    return 0;
  }

  hid_t dataset;
  HDF5_ERROR (dataset = H5Dopen (file, patch->patchname.c_str()));
  assert (dataset >= 0);

  // filereader reads the current timelevel
  int timelevel = in_recovery ? patch->timelevel : 0;

  // get dataset dimensions
  hid_t filespace;
  HDF5_ERROR (filespace = H5Dget_space (dataset));
  assert (group.dim == (group.grouptype == CCTK_SCALAR ?
                        patch->rank-1 : patch->rank));
  vector<hsize_t> h5shape(patch->rank);
  HDF5_ERROR (H5Sget_simple_extent_dims (filespace, &h5shape[0], NULL));

  ivect shape(1);
  int elems = 1;
  for (int i = 0; i < patch->rank; i++) {
    shape[i] = h5shape[patch->rank-i-1];
    elems   *= shape[i];
  }
  const hid_t datatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype, 0);

  const ivect stride = group.grouptype == CCTK_GF ?
                       maxspacereflevelfact/spacereflevelfact : 1;
  ivect lower = patch->iorigin * stride;
  ivect upper = lower + (shape - 1) * stride;

  // Traverse all local components on all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype) {
    struct arrdesc& data = arrdata.at(gindex).at(Carpet::map);

    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, group.grouptype) {

      // reset bounds for DISTRIB=CONST variables
      if (group.disttype == CCTK_DISTRIB_CONSTANT) {
        assert (group.grouptype==CCTK_SCALAR or group.grouptype==CCTK_ARRAY);
        const int newlower = lower[patch->rank-1] +
          (upper[patch->rank-1]-lower[patch->rank-1]+1)*
          (data.hh->processors().at(reflevel).at(component));
        const int newupper = upper[patch->rank-1] +
          (upper[patch->rank-1]-lower[patch->rank-1]+1)*
          (data.hh->processors().at(reflevel).at(component));
        lower[patch->rank-1] = newlower;
        upper[patch->rank-1] = newupper;
      }
      const ibbox filebox(lower, upper, stride);

      ibbox& interior_membox =
        data.dd->boxes.at(mglevel).at(reflevel).at(component).interior;

      // skip this dataset if it doesn't overlap with this component's interior
      const ibbox interior_overlap = interior_membox & filebox;
      if (interior_overlap.empty()) {
        continue;
      }

      // get the overlap with this component's exterior
      ibbox& membox =
        data.dd->boxes.at(mglevel).at(reflevel).at(component).exterior;
      const ibbox overlap = membox & filebox;
      bboxes_read.at(Carpet::map) |= overlap;

      // calculate hyperslab selection parameters
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 6 && H5_VERS_RELEASE >= 4)
      hsize_t memorigin[dim], fileorigin[dim];
#else
      hssize_t memorigin[dim], fileorigin[dim];
#endif
      hsize_t memdims[dim], count[dim];
      for (int i = 0; i < patch->rank; i++) {
        memdims[patch->rank-i-1] =
          (membox.upper() - membox.lower())[i] / stride[i] + 1;
        count[patch->rank-i-1] =
          (overlap.upper() - overlap.lower())[i] / stride[i] + 1;
        memorigin[patch->rank-i-1] =
          (overlap.lower() - membox.lower())[i] / stride[i];
        fileorigin[patch->rank-i-1] =
          group.disttype == CCTK_DISTRIB_CONSTANT ?
          0 : (overlap.lower() - lower)[i] / stride[i];
      }
      // read the hyperslab
      hid_t memspace;
      HDF5_ERROR (memspace  = H5Screate_simple (patch->rank, memdims, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, fileorigin,
                                       NULL, count, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (memspace, H5S_SELECT_SET, memorigin,
                                       NULL, count, NULL));
      HDF5_ERROR (H5Dread (dataset, datatype, memspace, filespace, H5P_DEFAULT,
                           cctkGH->data[patch->vindex][timelevel]));
      HDF5_ERROR (H5Sclose (memspace));

    } END_LOCAL_COMPONENT_LOOP;

    if (in_recovery) {
      data.tt->set_time (reflevel, mglevel,
                         ((cctkGH->cctk_time - cctk_initial_time)
                          / (delta_time * mglevelfact)) );
    }

  } END_MAP_LOOP;

  HDF5_ERROR (H5Sclose (filespace));
  HDF5_ERROR (H5Dclose (dataset));

  return 1;
}

} // namespace CarpetIOHDF5
