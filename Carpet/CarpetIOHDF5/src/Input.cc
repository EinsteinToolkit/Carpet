#include <cassert>
#include <cstring>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Parameters.h"

#include "CarpetIOHDF5.hh"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "defs.hh"


namespace CarpetIOHDF5
{

using namespace std;
using namespace Carpet;

// structure describing a patch as a single dataset of an HDF5 file to read from
typedef struct {
  string patchname;
  int vindex;
  int map;
  int mglevel;
  int reflevel;
  int timestep;
  int timelevel;
  int component;
  int rank;
  ivect iorigin;
  vector<hsize_t> shape;
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
  double global_time;
  double delta_time;
  vector<double> mgleveltimes;  // [num_mglevels*num_reflevels]

  vector<vector<vector<region_t> > > grid_structure; // [map][reflevel][component]
  vector<vector<vector<CCTK_REAL> > > grid_times;    // [map][mglevel][reflevel]
  vector<vector<CCTK_REAL> > leveltimes;             // [mglevel][reflevel]

} fileset_t;

// list of checkpoint/filereader files
static list<fileset_t> filesets;

// number of reflevels in the checkpoint
static int num_reflevels = -1;


static list<fileset_t>::iterator OpenFileSet (const cGH* const cctkGH,
                                              const string setname,
                                              const char *basefilename,
                                              int called_from);
static void ReadMetadata (fileset_t& fileset, hid_t file);
static herr_t BrowseDatasets (hid_t group, const char *objectname, void *arg);
static int ReadVar (const cGH* const cctkGH,
                    hid_t file,
                    CCTK_REAL & io_bytes,
                    list<patch_t>::const_iterator patch,
                    vector<ibset> &bboxes_read,
                    bool in_recovery);


//////////////////////////////////////////////////////////////////////////////
// Register with the Cactus Recovery Interface
//////////////////////////////////////////////////////////////////////////////
int CarpetIOHDF5_RecoverParameters ()
{
  return IOUtil_RecoverParameters (Recover, ".h5", "HDF5");
}

//////////////////////////////////////////////////////////////////////////////
// Recover the grid structure
//////////////////////////////////////////////////////////////////////////////
void CarpetIOHDF5_RecoverGridStructure (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  fileset_t & fileset = * filesets.begin();
  
  // Abort with an error if there is no grid structure in the
  // checkpoint file, or if the number of maps is wrong
  if (fileset.grid_structure.empty()) {
    CCTK_WARN(0, "No grid structure information found in checkpoint file !");
  } else if (int(fileset.grid_structure.size()) != maps) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Number of maps in the checkpoint's grid structure information "
               "(%d) doesn't match current number of maps (%d) !",
               int(fileset.grid_structure.size()), maps);
  }
  
  vector<vector<vector<region_t> > > superregsss = fileset.grid_structure;
  vector<vector<vector<vector<region_t> > > > regssss (maps);
  
  int type;
  void const * const ptr =
    CCTK_ParameterGet ("regrid_in_level_mode", "Carpet", & type);
  assert (ptr != 0);
  assert (type == PARAMETER_BOOLEAN);
  CCTK_INT const regrid_in_level_mode = * static_cast<CCTK_INT const *> (ptr);
  
  if (not regrid_in_level_mode) {
    // Distribute each map independently
    
    for (int m = 0; m < maps; ++ m) {
      vector <vector <region_t> > & superregss = superregsss.at(m);
      
      // Make multiprocessor aware
      vector <vector <region_t> > regss(superregss.size());
      for (size_t rl = 0; rl < superregss.size(); ++ rl) {
        Carpet::SplitRegions (cctkGH, superregss.at(rl), regss.at(rl));
      } // for rl
      
      // Make multigrid aware
      Carpet::MakeMultigridBoxes (cctkGH, m, regss, regssss.at(m));
      
    } // for m
    
  } else {                     // if regrid_in_level_mode
    // Distribute all maps at the same time
    
    vector<vector<vector<region_t> > > regsss(maps);
    
    // Count levels
    vector <int> rls (maps);
    for (int m = 0; m < maps; ++ m) {
      rls.at(m) = superregsss.at(m).size();
    }
    int maxrl = 0;
    for (int m = 0; m < maps; ++ m) {
      maxrl = max (maxrl, rls.at(m));
    }
    // All maps must have the same number of levels
    for (int m = 0; m < maps; ++ m) {
      superregsss.at(m).resize (maxrl);
      regsss.at(m).resize (maxrl);
    }
    
    // Make multiprocessor aware
    for (int rl = 0; rl < maxrl; ++ rl) {
      vector <vector <region_t> > superregss (maps);
      for (int m = 0; m < maps; ++ m) {
        superregss.at(m) = superregsss.at(m).at(rl);
      }
      vector <vector <region_t> > regss (maps);
      Carpet::SplitRegionsMaps (cctkGH, superregss, regss);
      for (int m = 0; m < maps; ++ m) {
        superregsss.at(m).at(rl) = superregss.at(m);
        regsss.at(m).at(rl) = regss.at(m);
      }
    } // for rl
    
    // Make multigrid aware
    Carpet::MakeMultigridBoxesMaps (cctkGH, regsss, regssss);
    
  } // if regrid_in_level_mode
  
  for (int m = 0; m < maps; ++ m) {
    
    // Regrid
    RegridMap (cctkGH, m, superregsss.at(m), regssss.at(m));
    
    // Set time hierarchy correctly after RegridMap created it
    for (int ml = 0; ml < mglevels; ++ ml) {
      for (size_t rl = 0; rl < fileset.grid_times.at(m).at(mglevel).size(); ++ rl) {
        vtt.at(m)->set_time (rl, ml, fileset.grid_times.at(m).at(ml).at(rl));
      }
    }
    
  } // for m
  
  PostRegrid (cctkGH);
  
  // Set level times correctly after PostRegrid created them
  leveltimes = fileset.leveltimes;

  for (int rl = 0; rl < reflevels; ++ rl) {
    Recompose (cctkGH, rl, false);
  }
}

//////////////////////////////////////////////////////////////////////////////
// Overwrite the "CarpetRegrid::refinement_levels"
// with the number of levels given in the checkpoint file
//
// Note that this has to be done after parameter recovery in order to have
// any effect of steering "CarpetRegrid::refinement_levels".
//////////////////////////////////////////////////////////////////////////////
int CarpetIOHDF5_SetNumRefinementLevels ()
{
  DECLARE_CCTK_PARAMETERS;

  if (num_reflevels > 0) {
    if (not CCTK_Equals (verbose, "none")) {
      char *buffer = CCTK_ParameterValString ("refinement_levels",
                                              "CarpetRegrid");
      assert (buffer);
      CCTK_VInfo (CCTK_THORNSTRING, "Using %i reflevels from checkpoint file. "
                  "Ignoring value '%s' in parameter file.",
                  num_reflevels, buffer);
      free (buffer);
    }

    char buffer[32];
    snprintf (buffer, sizeof (buffer), "%d", num_reflevels);
    int const retval = CCTK_ParameterSet ("refinement_levels", "CarpetRegrid",
                                          buffer);
    assert (retval == 0);
  }
  
  return 0;
}


//////////////////////////////////////////////////////////////////////////////
// close all open checkpoint/filereader files after recovering grid variables
//////////////////////////////////////////////////////////////////////////////
void CarpetIOHDF5_CloseFiles (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int error_count = 0;


  for (list<fileset_t>::const_iterator set = filesets.begin();
       set != filesets.end(); set++) {
    for (unsigned int i = 0; i < set->files.size(); i++) {
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
}


//////////////////////////////////////////////////////////////////////////////
// Top-level recovery/filereader routine
// Opens the checkpoint/data file on the first time through,
// recovers parameters and variables
//////////////////////////////////////////////////////////////////////////////
int Recover (cGH* cctkGH, const char *basefilename, int called_from)
{
  int error_count = 0;
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
    int const idx = mglevel*fileset->num_reflevels + reflevel;
    cctkGH->cctk_time = fileset->mgleveltimes.at(idx);

    if (use_grid_structure_from_checkpoint) {
      // recover the grid structure only once
      static bool is_first = true;
      if (is_first) {
        is_first = false;
        CarpetIOHDF5_RecoverGridStructure (cctkGH);
      }
    }
  }

  if (not CCTK_Equals (verbose, "none")) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "reading grid variables on mglevel %d reflevel %d",
                mglevel, reflevel);
  }

  CCTK_REAL io_files = 0;
  CCTK_REAL io_bytes = 0;
  BeginTimingIO (cctkGH);

  // create a bbox set for each active timelevel of all variables
  // to mark how much has been read already
  const int numvars = CCTK_NumVars ();
  vector<vector<bool> > read_completely(numvars);
  vector<vector<vector<ibset> > > bboxes_read (numvars);
  for (unsigned int vindex = 0; vindex < bboxes_read.size(); vindex++) {
    const int timelevels = CCTK_ActiveTimeLevelsVI (cctkGH, vindex);
    read_completely[vindex].resize (timelevels, false);
    bboxes_read[vindex].resize (timelevels);
    for (int tl = 0; tl < timelevels; tl++) {
      bboxes_read[vindex][tl].resize (maps);
    }
  }

  // query for groups which have the 'CHECKPOINT = "no"' option set
  // Such groups are not checked against being read completely.
  vector<bool> not_checkpointed(numvars);
  if (in_recovery) {
    for (unsigned int vindex = 0; vindex < not_checkpointed.size(); vindex++) {
      int gindex = CCTK_GroupIndexFromVarI (vindex);
      int tagstable = CCTK_GroupTagsTableI (gindex);
      int const len = Util_TableGetString (tagstable, 0, NULL, "checkpoint");
      if (len > 0) {
        char* value = new char[len + 1];
        Util_TableGetString (tagstable, len + 1, value, "checkpoint");
        if (len == sizeof ("no") - 1 and CCTK_Equals (value, "no")) {
          not_checkpointed[vindex] = true;
        }
        delete[] value;
      }
    }
  }

  // create a bbox set for each group to list how much needs to be read
  const int numgroups = CCTK_NumGroups ();
  vector<vector<ibset> > group_bboxes (numgroups);
  for (unsigned int gindex = 0; gindex < group_bboxes.size(); gindex++) {
    group_bboxes[gindex].resize (maps);
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

  // mark variables in groups with no grid points (size 0) as already done
  for (int group = 0; group < (int)group_bboxes.size(); group++) {
    bool is_empty = true;
    for (int m = 0; m < (int)group_bboxes[group].size(); m++) {
      is_empty &= group_bboxes[group][m].empty();
    }
    if (is_empty) {
      int vindex = CCTK_FirstVarIndexI (group);
      const int last = vindex + CCTK_NumVarsInGroupI (group);
      while (vindex < last) {
        read_completely[vindex].assign (read_completely[vindex].size(), true);
        vindex++;
      }
    }
  }

  const ioGH* ioUtilGH = (const ioGH*) CCTK_GHExtension (cctkGH, "IO");
  CarpetIOHDF5GH* myGH =
                  (CarpetIOHDF5GH*) CCTK_GHExtension (cctkGH, CCTK_THORNSTRING);
  // allocate list of recovery filenames
  // if the recovery checkpoint should be removed eventually
  if (in_recovery and recover_and_remove and checkpoint_keep > 0) {
    if (not myGH->recovery_filename_list) {
      myGH->recovery_num_filenames = fileset->files.size();
      myGH->recovery_filename_list =
            (char **) calloc (fileset->files.size(), sizeof (char **));
      assert (myGH->cp_filename_index == 0);

      // add a dummy entry in the checkpoint filename ring buffer
      myGH->cp_filename_list[myGH->cp_filename_index] = "bla";
      myGH->cp_filename_index = (myGH->cp_filename_index+1) % checkpoint_keep;
    }
  }

  // loop over all input files of this set
  for (unsigned int i = 0; i < fileset->files.size(); i++) {

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
      io_files += 1;

      if (CCTK_Equals (verbose, "full")) {
        CCTK_VInfo (CCTK_THORNSTRING, "opening %s file '%s'",
                    in_recovery ? "checkpoint" : "input", file.filename);
      }

      // browse through all datasets contained in this file
      HDF5_ERROR (H5Giterate (file.file, "/", NULL, BrowseDatasets, &file));
    }
    assert (file.patches.size() > 0);
    if (myGH->recovery_filename_list and not myGH->recovery_filename_list[i]) {
      myGH->recovery_filename_list[i] = strdup (file.filename);
    }

    // some optimisation for the case when all processors recover
    // from a single chunked checkpoint file:
    // reorder the list so that processor-local components come first
    if (fileset->nioprocs == 1) {
      list<patch_t>::iterator patch = file.patches.begin();
      while (patch != file.patches.end()) {
        if (patch->component == dist::rank()) {
          file.patches.push_front (*patch);
          patch = file.patches.erase (patch);
        } else {
          patch++;
        }
      }
    }

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

      // check the timelevel
      if ((unsigned int) patch->timelevel >=
          read_completely.at(patch->vindex).size()) {
        CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Ignoring dataset '%s' (invalid timelevel %d)",
                    patch->patchname.c_str(), patch->timelevel);
        continue;
      }

      // actually read the patch
      if (not read_completely.at(patch->vindex).at(patch->timelevel)) {
        error_count += ReadVar (cctkGH, file.file, io_bytes, patch,
                                bboxes_read.at(patch->vindex).at(patch->timelevel),
                                in_recovery);

        // update the read_completely entry
        const int gindex = CCTK_GroupIndexFromVarI (patch->vindex);
        read_completely.at(patch->vindex).at(patch->timelevel) =
          bboxes_read.at(patch->vindex).at(patch->timelevel) ==
          group_bboxes.at(gindex);
      }
    }

    // check if all variables have been read completely already
    bool all_done = true;
    for (unsigned int vindex = 0; vindex < read_completely.size(); vindex++) {

      // skip all variables which aren't expected to be recovered
      if (not_checkpointed[vindex] or
          (CCTK_GroupTypeFromVarI (vindex) != CCTK_GF and reflevel > 0)) {
        continue;
      }

      for (unsigned int tl = 0; tl < read_completely[vindex].size(); tl++) {
        all_done &= read_completely[vindex][tl];
	if(not read_completely[vindex][tl]) {
	  char * const fullname = CCTK_FullName (vindex);
	  CCTK_VWarn(1,__LINE__,__FILE__, CCTK_THORNSTRING, "Variable %s on rl %d and tl %d not"
		      " read completely. Will have to look for it in other files.", 
		      fullname,reflevel,tl);
	  free(fullname);
	}
      }
    }
    if (all_done) {
      break;
    }

    // keep the file open if not requested otherwise by the user
    if (open_one_input_file_at_a_time) {
      HDF5_ERROR (H5Fclose (file.file));
      file.file = -1;
    }
  }

  // check that all variables have been read completely on this mglevel/reflevel
  int num_incomplete = 0;
  for (unsigned int vindex = 0; vindex < read_completely.size(); vindex++) {

    if (CCTK_GroupTypeFromVarI (vindex) != CCTK_GF and reflevel > 0) {
      continue;
    }

    for (unsigned int tl = 0; tl < read_completely[vindex].size(); tl++) {
      if (called_from == FILEREADER_DATA and not
          (ioUtilGH->do_inVars and ioUtilGH->do_inVars[vindex])) {
        continue;
      }

      if (not read_completely[vindex][tl]) {
        // check if the variable has been read partially
        size_t size = 0;
        for (int m = 0; m < maps; m++) {
          size += bboxes_read[vindex][tl][m].size();
        }
        char* fullname = CCTK_FullName (vindex);
        if (size == 0) {
          if (not_checkpointed[vindex]) {
            CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "variable '%s' timelevel %d has not been read "
                        "(variable has option tag \"CHECKPOINT = 'no'\")",
                        fullname, tl);
          } else {
            CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "variable '%s' timelevel %d has not been read",
                        fullname, tl);
          }
        } else {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "variable '%s' timelevel %d has been read only partially",
                      fullname, tl);
          num_incomplete++;
        }
        free (fullname);
      }
    }
  }
  if (num_incomplete) {
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%d variables on mglevel %d reflevel %d have been read "
                "only partially", num_incomplete, mglevel, reflevel);
  }
  
  {
    CCTK_REAL local[2], global[2];
    local[0] = io_files;
    local[1] = io_bytes;
    MPI_Allreduce (local, global, 2, dist::datatype (local[0]), MPI_SUM, dist::comm());
    io_files = global[0];
    io_bytes = global[1];
  }
  EndTimingIO (cctkGH, io_files, io_bytes, true);
  
  // Synchronise all variables which have been read
  {
    vector <bool> dosync (CCTK_NumGroups(), false);
    for (size_t vindex = 0; vindex < read_completely.size(); vindex++) {
      if (CCTK_GroupTypeFromVarI (vindex) == CCTK_GF or reflevel == 0) {
        for (size_t tl = 0; tl < read_completely[vindex].size(); tl++) {
          if (called_from != FILEREADER_DATA or
              (ioUtilGH->do_inVars and ioUtilGH->do_inVars[vindex]))
          {
            if (read_completely[vindex][tl]) {
              int const gindex = CCTK_GroupIndexFromVarI (vindex);
              dosync.at(gindex) = true;
            }
          }
        }
      }
    }
    for (comm_state state; not state.done(); state.step()) {
      for (size_t group = 0; group < dosync.size(); ++ group) {
        if (dosync.at(group)) {
          for (size_t m = 0; m < arrdata.at(group).size(); ++ m) {
            arrdesc & ad = arrdata.at(group).at(m);
            for (int ml = 0; ml < ad.hh->mglevels(); ++ ml) {
              for (int rl = 0; rl < ad.hh->reflevels(); ++ rl) {
                for (size_t v = 0; v < ad.data.size(); ++ v) {
                  ggf * const gf = ad.data.at(v);
                  for (int tl = 0; tl < gf->timelevels (ml, rl); ++ tl) {
                    gf->sync_all (state, tl, rl, ml);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (error_count and abort_on_io_errors) {
    CCTK_WARN(0, "Found errors while trying to restart from checkpoint, aborting.");
  }
  
  if (in_recovery and not CCTK_Equals (verbose, "none")) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "restarting simulation on mglevel %d reflevel %d "
                "at iteration %d (simulation time %g)", mglevel, reflevel,
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
  int error_count = 0;
  list<string> filenames;
  DECLARE_CCTK_PARAMETERS;

  // first try to open a chunked file written on this processor
  // (note that dist::rank() cannot be called yet during RECOVER_PARAMETERS)
  fileset.first_ioproc = CCTK_MyProc (cctkGH);
  file.filename = IOUtil_AssembleFilename (NULL, basefilename, "", ".h5",
                                           called_from, fileset.first_ioproc, 0);

  // try to open the file (prevent HDF5 error messages if it fails)
  H5E_BEGIN_TRY {
    filenames.push_back (string (file.filename));
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
      filenames.push_back (string (file.filename));
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
      filenames.push_back (string (file.filename));
      file.file = H5Fopen (file.filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    } H5E_END_TRY;
  }

  // return if no valid checkpoint could be found
  if (file.file < 0) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No valid HDF5 file with basename \"%s\" found", basefilename);
    free (file.filename);
    for (list<string>::const_iterator
           lsi = filenames.begin(); lsi != filenames.end(); ++ lsi)
    {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Tried filename \"%s\"", lsi->c_str());
    }
    return (filesets.end());
  }

  if (CCTK_Equals (verbose, "full")) {
    CCTK_VInfo (CCTK_THORNSTRING, "opening %s file '%s'",
                called_from == CP_RECOVER_PARAMETERS ?
                "checkpoint" : "input", file.filename);
  }

  // read all the metadata information
  ReadMetadata (fileset, file.file);

  // first try to open a chunked file written on this processor
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

    num_reflevels = fileset.num_reflevels;
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
  int error_count = 0;
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

  // Read grid structure if it is present
  hid_t dataset;
  H5E_BEGIN_TRY {
    dataset = H5Dopen (metadata, GRID_STRUCTURE);
  } H5E_END_TRY;
  if (dataset >= 0) {
    vector<char> gs_cstr (H5Dget_storage_size (dataset) + 1);
    HDF5_ERROR (H5Dread (dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, &gs_cstr.front()));
    HDF5_ERROR (H5Dclose (dataset));
    istringstream gs_buf (&gs_cstr.front());
    gs_buf >> fileset.grid_structure;
    gs_buf >> fileset.grid_times;
    gs_buf >> fileset.leveltimes;
  }

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
  int error_count = 0;
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
  patch.shape.resize (patch.rank);
  HDF5_ERROR (H5Sget_simple_extent_dims (dataspace, &patch.shape[0], NULL));
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

  // try to obtain the map and component number from the patch's name
  // (cleaner way would be to store attributes with the dataset)
  patch.map = 0;
  const char* map_string = strstr (objectname, " m=");
  if (map_string) {
    sscanf (map_string, " m=%d", &patch.map);
  }
  patch.component = -1;
  const char* component_string = strstr (objectname, " c=");
  if (component_string) {
    sscanf (component_string, " c=%d", &patch.component);
  }

  // add this patch to our list
  if (patch.vindex >= 0 and patch.vindex < CCTK_NumVars ()) {
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
                    CCTK_REAL & io_bytes,
                    list<patch_t>::const_iterator patch,
                    vector<ibset> &bboxes_read,
                    bool in_recovery)
{
  int error_count = 0;
  DECLARE_CCTK_PARAMETERS;

  const int gindex = CCTK_GroupIndexFromVarI (patch->vindex);
  assert (gindex >= 0 and (unsigned int) gindex < Carpet::arrdata.size());
  cGroup group;
  CCTK_GroupData (gindex, &group);


  // grid arrays and scalars are read once on the coarsest reflevel
  assert (group.grouptype == CCTK_GF or (mglevel == 0 and reflevel == 0));

  // Check for storage
  if (not CCTK_QueryGroupStorageI(cctkGH, gindex)) {
    char *fullname = CCTK_FullName (patch->vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable '%s' (no storage)", fullname);
    free (fullname);
    return 1;
  }

  // filereader reads the current timelevel
  int timelevel = in_recovery ? patch->timelevel : 0;

  // get dataset dimensions
  assert (group.dim == (group.grouptype == CCTK_SCALAR ?
                        patch->rank-1 : patch->rank));
  ivect shape(1);
  int elems = 1;
  for (int i = 0; i < patch->rank; i++) {
    shape[i] = patch->shape[patch->rank-i-1];
    elems   *= shape[i];
  }
  const hid_t datatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype, 0);

  const ivect stride = group.grouptype == CCTK_GF ?
                       maxspacereflevelfact/spacereflevelfact : 1;
  ivect lower = patch->iorigin * stride;
  ivect upper = lower + (shape - 1) * stride;

  // Traverse all local components on all maps
  hid_t filespace = -1, dataset = -1;
  hid_t xfer = H5P_DEFAULT;
  H5Z_EDC_t checksums = H5Z_NO_EDC;

  BEGIN_MAP_LOOP (cctkGH, group.grouptype) {

    // skip this dataset if it belongs to another map
    if (group.grouptype == CCTK_GF and patch->map != Carpet::map) {
      continue;
    }

    struct arrdesc& data = arrdata.at(gindex).at(Carpet::map);

    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, group.grouptype) {

      // reset bounds for DISTRIB=CONST variables
      if (group.disttype == CCTK_DISTRIB_CONSTANT) {
        assert (group.grouptype==CCTK_SCALAR or group.grouptype==CCTK_ARRAY);
        const int newlower = lower[patch->rank-1] +
          (upper[patch->rank-1]-lower[patch->rank-1]+1)*
          (data.hh->processor(reflevel,component));
        const int newupper = upper[patch->rank-1] +
          (upper[patch->rank-1]-lower[patch->rank-1]+1)*
          (data.hh->processor(reflevel,component));
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

      if (CCTK_Equals (verbose, "full")) {
        char *fullname = CCTK_FullName (patch->vindex);
        CCTK_VInfo (CCTK_THORNSTRING, "  reading '%s' from dataset '%s'",
                    fullname, patch->patchname.c_str());
        free (fullname);
      }

      // get the overlap with this component's exterior
      ibbox& membox =
        data.dd->boxes.at(mglevel).at(reflevel).at(component).exterior;
      const ibbox overlap = membox & filebox;
      bboxes_read.at(Carpet::map) |= overlap;

      // calculate hyperslab selection parameters

      // before HDF5-1.6.4 the H5Sselect_hyperslab() function expected
      // the 'start' argument to be of type 'hssize_t'
      slice_start_size_t memorigin[dim], fileorigin[dim];
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

      // open the dataset on the first time through
      if (dataset < 0) {
        HDF5_ERROR (dataset = H5Dopen (file, patch->patchname.c_str()));
        HDF5_ERROR (filespace = H5Dget_space (dataset));
        xfer = H5Pcreate (H5P_DATASET_XFER);
        checksums = H5Pget_edc_check (xfer);
        if (use_checksums && (checksums == H5Z_DISABLE_EDC))
          CCTK_WARN(1, "Checksums not enabled in HDF5 library, but "
                       "requested in parameter, reading data without "
                       "tests on checksums.");
        if (!use_checksums && (checksums == H5Z_ENABLE_EDC))
          H5Pset_edc_check(xfer, H5Z_DISABLE_EDC);
      }

      // read the hyperslab
      hid_t memspace;
      HDF5_ERROR (memspace  = H5Screate_simple (patch->rank, memdims, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, fileorigin,
                                       NULL, count, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (memspace, H5S_SELECT_SET, memorigin,
                                       NULL, count, NULL));
      HDF5_ERROR (H5Dread (dataset, datatype, memspace, filespace, xfer,
                           cctkGH->data[patch->vindex][timelevel]));
      hid_t datatype;
      HDF5_ERROR (datatype = H5Dget_type (dataset));
      io_bytes += H5Sget_select_npoints (filespace) * H5Tget_size (datatype);
      HDF5_ERROR (H5Tclose (datatype));
      HDF5_ERROR (H5Sclose (memspace));

    } END_LOCAL_COMPONENT_LOOP;

    if (in_recovery) {
      data.tt->set_time (reflevel, mglevel,
                         ((cctkGH->cctk_time - cctk_initial_time)
                          / (delta_time * mglevelfact)) );
    }

  } END_MAP_LOOP;

  if (dataset >= 0) {
    HDF5_ERROR (H5Dclose (dataset));
    HDF5_ERROR (H5Sclose (filespace));
  }

  return error_count;
}

} // namespace CarpetIOHDF5
