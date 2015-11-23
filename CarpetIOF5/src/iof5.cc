#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#include <hdf5.h>

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include <Timer.hh>

#include "iof5.hh"

namespace CarpetIOF5 {

using namespace std;
using namespace Carpet;

// Checkpointing state
static int last_checkpoint_iteration = -1;

// Scheduled startup routine
int CarpetIOF5_Startup() {
  CCTK_RegisterBanner("AMR F5 I/O provided by CarpetIOF5");

  int const GHExtension = CCTK_RegisterGHExtension(CCTK_THORNSTRING);
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  return 0;
}

// Registered GH extension setup routine
void *SetupGH(tFleshConfig *const fleshconfig, int const convLevel,
              cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // register I/O method
  int const ierr = IOUtil_RegisterRecover("CarpetIOF5 recovery", Input);
  assert(not ierr);

  int const IOMethod = CCTK_RegisterIOMethod("IOF5");
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);

  // there no actual extension data structure
  return NULL;
}

// Scheduled initialisation routine
void CarpetIOF5_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = -1;
  *next_output_iteration = 0;

  last_checkpoint_iteration = cctk_iteration;
}

// A mechanism to keep the HDF5 output file open across multiple
// write operations:
hid_t open_file = H5I_INVALID_HID;
int keep_file_open = 0;
void enter_keep_file_open() { ++keep_file_open; }
void leave_keep_file_open() {
  assert(keep_file_open > 0);
  --keep_file_open;
  if (keep_file_open == 0 and open_file >= 0) {
    herr_t const herr = H5Fclose(open_file);
    assert(not herr);
    open_file = H5I_INVALID_HID;
    H5garbage_collect();
  }
}

// Interpret Cactus parameters to decide which variables have been
// selected for output
class selection_t {
  int iteration;
  int times_set;
  vector<bool> selected;
  vector<int> last_output_iteration;
  static void do_output(int const vindex, char const *const optstring,
                        void *const callback_arg) {
    static_cast<selection_t *>(callback_arg)->selected.at(vindex) = true;
  }

public:
  selection_t() : iteration(-1), times_set(-1) {}
  bool is_selected(cGH const *const cctkGH, int const vindex) {
    DECLARE_CCTK_PARAMETERS;
    // Check whether the parameter out_vars changed only once per
    // iteration
    if (cctkGH->cctk_iteration != iteration) {
      iteration = cctkGH->cctk_iteration;
      int const current_times_set =
          CCTK_ParameterQueryTimesSet("out_vars", CCTK_THORNSTRING);
      // Re-scan out_vars (which is somewhat expensive) only if it
      // has changed
      if (current_times_set > times_set) {
        times_set = current_times_set;
        selected.resize(CCTK_NumVars());
        fill(selected.begin(), selected.end(), false);
        // for (int n=0; n<CCTK_NumVars(); ++n) selected.at(n) = false;
        CCTK_TraverseString(out_vars, do_output, this, CCTK_GROUP_OR_VAR);
      }
    }
    return selected.at(vindex);
  }
  bool should_output(cGH const *const cctkGH, int const vindex) {
    last_output_iteration.resize(CCTK_NumVars(), -1);
    return last_output_iteration.at(vindex) < cctkGH->cctk_iteration;
  }
  void did_output(cGH const *const cctkGH, int const vindex) {
    last_output_iteration.resize(CCTK_NumVars(), -1);
    last_output_iteration.at(vindex) = cctkGH->cctk_iteration;
  }
};
selection_t output_variables;

// Callbacks for CarpetIOHDF5's I/O method

int OutputGH(cGH const *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("F5::OutputGH");
  timer.start();

  enter_keep_file_open();
  for (int vindex = 0; vindex < CCTK_NumVars(); ++vindex) {
    if (TimeToOutput(cctkGH, vindex)) {
      TriggerOutput(cctkGH, vindex);
    }
  }
  leave_keep_file_open();

  timer.stop();

  return 0;
}

int TimeToOutput(cGH const *const cctkGH, int const vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const numvars = CCTK_NumVars();
  assert(vindex >= 0 and vindex < numvars);

  // Output only in global mode
  if (not do_global_mode)
    return 0;

  // Output only selected variables
  if (not output_variables.is_selected(cctkGH, vindex))
    return 0;

  if (not output_variables.should_output(cctkGH, vindex))
    return 0;

  // whether to output at this iteration
  bool output_this_iteration = false;

  int const my_out_every = out_every == -2 ? IO_out_every : out_every;
  int const my_out0D_every = out0D_every == -2 ? my_out_every : out0D_every;
  int const my_out1D_every = out1D_every == -2 ? my_out_every : out1D_every;
  int const my_out2D_every = out2D_every == -2 ? my_out_every : out2D_every;
  int const my_out3D_every = out3D_every == -2 ? my_out_every : out3D_every;
  if (my_out0D_every > 0 or my_out1D_every > 0 or my_out2D_every > 0 or
      my_out3D_every > 0) {
    if (*this_iteration == cctk_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else if (cctk_iteration >= *next_output_iteration) {
      // it is time for the next output
      output_this_iteration = true;
      *this_iteration = cctk_iteration;
      *next_output_iteration = cctk_iteration + my_out_every;
    }
  }

  return output_this_iteration ? 1 : 0;
}

int TriggerOutput(cGH const *const cctkGH, int const vindex) {
  DECLARE_CCTK_PARAMETERS;

  char *const fullname = CCTK_FullName(vindex);
  int const gindex = CCTK_GroupIndexFromVarI(vindex);
  char *const groupname = CCTK_GroupName(gindex);
  for (char *p = groupname; *p; ++p)
    *p = tolower(*p);
  int const retval = OutputVarAs(cctkGH, fullname, groupname);
  free(groupname);
  free(fullname);

  output_variables.did_output(cctkGH, vindex);

  return retval;
}

int OutputVarAs(cGH const *const cctkGH, const char *const varname,
                const char *const alias) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_level_mode());
  BEGIN_GLOBAL_MODE(cctkGH) {
    DECLARE_CCTK_ARGUMENTS;

    CCTK_VInfo(CCTK_THORNSTRING, "F5::OutputVarAs: iteration=%d, variable=%s",
               cctk_iteration, varname);

    // We don't know how to open multiple files yet
    assert(CCTK_EQUALS(file_content, "everything"));

    // Determine number of I/O processes
    int const myproc = CCTK_MyProc(cctkGH);
    int const nprocs = CCTK_nProcs(cctkGH);

    int const ioproc_every =
        max_nioprocs == 0 ? 1 : (nprocs + max_nioprocs - 1) / max_nioprocs;
    assert(ioproc_every > 0);
    int const nioprocs = nprocs / ioproc_every;
    assert(nioprocs > 0);
    assert(max_nioprocs == 0 or nioprocs <= max_nioprocs);
    assert(nioprocs <= nprocs);
    int const myioproc = myproc / ioproc_every * ioproc_every;
    // We split processes into "I/O groups" which range from
    // myioproc to myioproc+ioproc_every-1 (inclusive). Within each
    // group, at most one process can perform I/O.

    // If I am not the first process in my I/O group, wait for a
    // token from my predecessor
    if (myproc > myioproc) {
      MPI_Recv(NULL, 0, MPI_INT, myproc - 1, 0, dist::comm(),
               MPI_STATUS_IGNORE);
    }

    // Open file
    static bool first_time = true;

    // The file name doesn't matter since we currently write
    // everything into a single file
    int const vindex = CCTK_VarIndex(varname);
    assert(vindex >= 0);
    string const basename = generate_basename(cctkGH, vindex);
    int const ioproc = myioproc / ioproc_every;
    string const name =
        create_filename(cctkGH, basename, cctkGH->cctk_iteration, ioproc,
                        io_dir_output, first_time);

    indent_t indent;
    cout << indent << "I/O process=" << ioproc << "\n";

    enter_keep_file_open();
    bool const truncate_file =
        first_time and IO_TruncateOutputFiles(cctkGH) and myproc == myioproc;
    if (open_file < 0) {
      // Reuse file hid if file is already open
      hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
      open_file = truncate_file ? H5Fcreate(name.c_str(), H5F_ACC_TRUNC,
                                            H5P_DEFAULT, fapl)
                                : H5Fopen(name.c_str(), H5F_ACC_RDWR, fapl);
      assert(open_file >= 0);
      H5Pclose(fapl);
    }
    first_time = false;

    vector<bool> output_var(CCTK_NumVars());
    output_var.at(vindex) = true;
    // NOTE: We should output metadata at most once per iteration,
    // probably only once per restart (per file)
    output(cctkGH, open_file, output_var, false, true);

    // Close file
    leave_keep_file_open();

    // If I am not the last process in my I/O group, send a token to
    // my successor
    if (myproc < min(myioproc + ioproc_every, nprocs) - 1) {
      MPI_Send(NULL, 0, MPI_INT, myproc + 1, 0, dist::comm());
    }
  }
  END_GLOBAL_MODE;

  return 0; // no error
}

// Checkpointing

void Checkpoint(cGH const *const cctkGH, int const called_from) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_global_mode());

  CCTK_VInfo(CCTK_THORNSTRING, "F5::Checkpoint: iteration=%d",
             cctkGH->cctk_iteration);

#if 0
    // generate filenames for both the temporary and real checkpoint
    // files
    int const ioproc = CCTK_MyProc(cctkGH);
    int const parallel_io = 1;
    char *const filename =
      IOUtil_AssembleFilename(cctkGH, NULL, "", ".f5",
                              called_from, ioproc, not parallel_io);
    char *const tempname =
      IOUtil_AssembleFilename(cctkGH, NULL, ".tmp", ".f5",
                              called_from, ioproc, not parallel_io);
#endif
  int const myproc = CCTK_MyProc(cctkGH);
  int const proc = myproc;
  string const name =
      create_filename(cctkGH, "checkpoint", cctkGH->cctk_iteration, proc,
                      io_dir_checkpoint, true);
  string const tempname =
      create_filename(cctkGH, "checkpoint.tmp", cctkGH->cctk_iteration, proc,
                      io_dir_checkpoint, true);

  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  hid_t const file =
      H5Fcreate(tempname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  assert(file >= 0);
  H5Pclose(fapl);

  vector<bool> output_var(CCTK_NumVars(), false);
  for (int gindex = 0; gindex < CCTK_NumGroups(); ++gindex) {
    // We can't check for storage here, since this requires at least
    // level mode
    if (not CCTK_QueryGroupStorageI(cctkGH, gindex))
      continue;

    // do not checkpoint groups with a "checkpoint=no" tag
    cGroup gdata;
    CCTK_GroupData(gindex, &gdata);
    int const len = Util_TableGetString(gdata.tagstable, 0, NULL, "checkpoint");
    if (len > 0) {
      char buf[1000];
      Util_TableGetString(gdata.tagstable, sizeof buf, buf, "checkpoint");
      if (CCTK_EQUALS(buf, "no"))
        continue;
      assert(CCTK_EQUALS(buf, "yes"));
    }

    int const first_vindex = CCTK_FirstVarIndexI(gindex);
    int const num_vars = CCTK_NumVarsInGroupI(gindex);
    if (num_vars > 0) {
      for (int vindex = first_vindex; vindex < first_vindex + num_vars;
           ++vindex) {
        output_var.at(vindex) = true;
      }
    }
  }

  // NOTE: We could write the metadata into only one of the process
  // files (to save space), or write it into a separate metadata
  // file
  output(cctkGH, file, output_var, true, true);

  // Close file
  herr_t const herr = H5Fclose(file);
  assert(not herr);
  H5garbage_collect();

  // Wait until all files have been written, then rename the
  // checkpoint files
  // TODO: ensure there were no errors
  CCTK_Barrier(cctkGH);
  int const ierr = rename(tempname.c_str(), name.c_str());
  assert(not ierr);
}

void CarpetIOF5_InitialDataCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint and checkpoint_ID) {
    Checkpoint(cctkGH, CP_INITIAL_DATA);
  }
}

void CarpetIOF5_EvolutionCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  bool const checkpoint_by_iteration =
      checkpoint_every > 0 and
      cctk_iteration >= last_checkpoint_iteration + checkpoint_every;
  bool const do_checkpoint =
      checkpoint and (checkpoint_by_iteration or checkpoint_next);

  if (do_checkpoint) {
    Checkpoint(cctkGH, CP_EVOLUTION_DATA);
  }
}

void CarpetIOF5_TerminationCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint and checkpoint_on_terminate) {

    bool const did_checkpoint = cctk_iteration == last_checkpoint_iteration;
    bool const do_checkpoint = not did_checkpoint;

    if (do_checkpoint) {
      Checkpoint(cctkGH, CP_EVOLUTION_DATA);
    }
  }
}

// Recovery information
static int recovery_iteration = -1;

int Input(cGH *const cctkGH, char const *const basefilename,
          int const called_from) {
  DECLARE_CCTK_PARAMETERS;

  herr_t herr;

  assert(called_from == CP_RECOVER_PARAMETERS or
         called_from == CP_RECOVER_DATA or called_from == FILEREADER_DATA);
  bool const in_recovery =
      called_from == CP_RECOVER_PARAMETERS or called_from == CP_RECOVER_DATA;
  if (not in_recovery) {
    assert(is_level_mode());
    if (reflevel > 0)
      return 0; // no error
  }

  BEGIN_GLOBAL_MODE(cctkGH) {
    DECLARE_CCTK_ARGUMENTS;

    if (in_recovery) {
      CCTK_VInfo(CCTK_THORNSTRING, "F5::Input: recovering iteration %d",
                 recovery_iteration);
    } else {
      CCTK_VInfo(CCTK_THORNSTRING, "F5::Input: reading iteration %d",
                 cctk_iteration);
    }
    cout << "called_from=" << (called_from == CP_RECOVER_PARAMETERS
                                   ? "CP_RECOVER_PARAMETERS"
                                   : called_from == CP_RECOVER_DATA
                                         ? "CP_RECOVER_DATA"
                                         : called_from == FILEREADER_DATA
                                               ? "FILEREADER_DATA"
                                               : NULL)
         << "\n";

    // Determine which variables to read
    ioGH const *const ioUtilGH = (ioGH const *)CCTK_GHExtension(cctkGH, "IO");
    vector<bool> input_var(CCTK_NumVars(), true);
    if (ioUtilGH->do_inVars) {
      for (int n = 0; n < CCTK_NumVars(); ++n) {
        input_var.at(n) = ioUtilGH->do_inVars[n];
      }
    }

    // Keep track of which files could be read, and which could not
    int foundproc = -1, notfoundproc = -1;

    // string const basename =
    //   in_recovery
    //   ? "checkpoint"
    //   : generate_basename(cctkGH, CCTK_VarIndex("grid::r"));
    string const basename = basefilename;

    {
      scatter_t scatter(cctkGH);

      // Iterate over files

      // TODO: Store how many processes contributed to the output,
      // and expect exactly that many files
      int const myproc = CCTK_MyProc(cctkGH);
      int const nprocs = CCTK_nProcs(cctkGH);
      // Loop over all (possible) files
      for (int proc = myproc;; proc += nprocs) {
        string const name =
            create_filename(cctkGH, basename, cctkGH->cctk_iteration, proc,
                            in_recovery ? io_dir_recover : io_dir_input, false);

        bool file_exists;
        H5E_BEGIN_TRY { file_exists = H5Fis_hdf5(name.c_str()) > 0; }
        H5E_END_TRY;
        if (not file_exists) {
          notfoundproc = proc;
          break;
        }
        foundproc = proc;

        indent_t indent;
        cout << indent << "process=" << proc << "\n";

        hid_t const file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        assert(file >= 0);

        // Iterate over all time slices
        bool const input_past_timelevels = in_recovery;
        // Read metadata when recoverying parameters
        bool const input_metadata = called_from == CP_RECOVER_PARAMETERS;
        input(cctkGH, file, input_var, input_past_timelevels, input_metadata,
              scatter);

        // Close file
        herr = H5Fclose(file);
        assert(not herr);
        H5garbage_collect();
      }

      // Destroy scatter object
    }

    {
      int maxfoundproc;
      MPI_Allreduce(&foundproc, &maxfoundproc, 1, MPI_INT, MPI_MAX,
                    dist::comm());
      if (maxfoundproc == -1) {
        string const name = create_filename(
            cctkGH, basename, cctkGH->cctk_iteration, notfoundproc,
            in_recovery ? io_dir_recover : io_dir_input, false);
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not read input file \"%s\"", name.c_str());
        return 1;
      }
      if (notfoundproc > -1 and notfoundproc <= maxfoundproc) {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not read file of process %d (but could read file of "
                   "process %d)",
                   notfoundproc, maxfoundproc);
      }
    }

    // The scatter object is implicitly destroyed here
  }
  END_GLOBAL_MODE;

  return 0; // no error
}

int CarpetIOF5_RecoverParameters() {
#if 0
    return IOUtil_RecoverParameters(Input, ".f5", "F5");
#endif

  DECLARE_CCTK_PARAMETERS;

  char const *const IO_dir = IO_recover_dir;
  char const *const F5_dir = recover_dir;
  bool const use_IO_dir = strcmp(F5_dir, "") == 0;
  char const *const my_recover_dir = use_IO_dir ? IO_dir : F5_dir;

  DIR *const dir = opendir(my_recover_dir);
  if (not dir) {
    // The recovery directory does not exist
    if (CCTK_Equals(recover, "autoprobe")) {
      // This is harmless when "autoprobe" is used
      CCTK_VInfo(CCTK_THORNSTRING, "Recovery directory \"%s\" doesn't exist",
                 my_recover_dir);
      return 0;
    } else {
      // This is an error when "auto" is used
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Recovery directory \"%s\" doesn't exist", my_recover_dir);
      return -2;
    }
  }

  // Get the list of potential recovery files
  char const *const my_recover_file = "checkpoint";
  string const prefix = string(my_recover_file) + ".i";
  string const infix = ".p";
  string const suffix = ".f5";
  assert(recovery_iteration < 0);
  while (dirent *const file = readdir(dir)) {
    char *p = file->d_name;

    // First check the file prefix
    if (prefix.compare(0, prefix.length(), p, prefix.length()) != 0)
      continue;
    p += prefix.length();

    // Now check if there is an iteration number following the file
    // prefix
    int const iter = strtol(p, &p, 10);
    if (!*p)
      continue;

    // Read (and ignore) the process number
    if (infix.compare(0, infix.length(), p, infix.length()) != 0)
      continue;
    p += infix.length();
    /* int const proc = */ strtol(p, &p, 10);
    if (!*p)
      continue;

    // Finally check the file extension
    if (suffix.compare(0, suffix.length(), p, suffix.length()) != 0)
      continue;
    p += suffix.length();

    // Check whether we read the whole string
    if (*p)
      continue;

    // Found a recovery file by that basename
    recovery_iteration = max(recovery_iteration, iter);
  }
  closedir(dir);

  // There is no recovery file
  if (recovery_iteration < 0) {
    if (CCTK_Equals(recover, "autoprobe")) {
      // This is harmless when "autoprobe" is used
      CCTK_VInfo(CCTK_THORNSTRING,
                 "No F5 checkpoint files with basefilename \"%s\" found in "
                 "recovery directory \"%s\"",
                 my_recover_file, my_recover_dir);
      return 0;
    } else {
      // This is an error when "auto" is used
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "No F5 checkpoint files with basefilename \"%s\" found in "
                 "recovery directory \"%s\"",
                 my_recover_file, my_recover_dir);
      return -1;
    }
  }

  return Input(NULL, "checkpoint", CP_RECOVER_PARAMETERS);
}

} // end namespace CarpetIOF5
