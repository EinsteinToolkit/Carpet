#include "input.hpp"
#include "output.hpp"
#include "pthread_wrapper.hpp"
#include "util.hpp"

#include <CactusBase/IOUtil/src/ioGH.h>
#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>

#include <Timer.hh>
#include <carpet.hh>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <H5Cpp.h>
#include <H5public.h>
#include <SimulationIO/SimulationIO.hpp>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace CarpetSimulationIO {
using namespace std;

// Checkpointing state
int last_checkpoint_iteration = -1;
int last_checkpoint_walltime_seconds = -1;

// Active output threads
vector<unique_ptr<pthread_wrapper_t> > output_pthreads;

// Align upwards to the next multiple
template <typename T> T align_up(T i, T j) {
  assert(i >= 0);
  assert(j > 0);
  return (i + j - 1) / j * j;
}

void *SetupGH(tFleshConfig *fleshconfig, int convLevel, cGH *cctkGH);

// Scheduled startup routine
extern "C" int CarpetSimulationIO_Startup() {
  CCTK_RegisterBanner("SimulationIO-based I/O provided by CarpetSimulationIO");

  int GHExtension = CCTK_RegisterGHExtension(CCTK_THORNSTRING);
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  return 0;
}

int OutputGH(const cGH *cctkGH);
int TimeToOutput(const cGH *cctkGH, int vindex);
int TriggerOutput(const cGH *cctkGH, int vindex);
int OutputVarAs(const cGH *cctkGH, const char *varname, const char *alias);

int Input(cGH *cctkGH, const char *basefilename, int called_from);

// Registered GH extension setup routine
void *SetupGH(tFleshConfig *fleshconfig, int convLevel, cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // Register I/O method
  int IOMethod = CCTK_RegisterIOMethod("SimulationIO");
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);

  // Register file reader
  int ierr = IOUtil_RegisterRecover("CarpetIOF5 recovery", Input);
  assert(not ierr);

  return nullptr; // there is no actual extension data structure
}

// Scheduled initialisation routine
extern "C" void CarpetSimulationIO_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // *this_iteration = -1;
  // *next_output_iteration = 0;
  *last_out0D_iteration = -1;
  *last_out1D_iteration = -1;
  *last_out2D_iteration = -1;
  *last_out3D_iteration = -1;

  last_checkpoint_iteration = cctk_iteration;
  last_checkpoint_walltime_seconds = CCTK_RunTime();

  init_comm();
}

void finalise_input_file_ptrs();

extern "C" void CarpetSimulationIO_FinaliseInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  finalise_input_file_ptrs();
}

// Scheduled finalization routine
extern "C" void CarpetSimulationIO_Finalise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // Wait for all output threads before shutting down
  output_pthreads.clear();

  finalise_comm();
}

////////////////////////////////////////////////////////////////////////////////

// Output

// Interpret Cactus parameters to decide which variables have been
// selected for output
class selection_t {
  int iteration;
  int times_set;
  vector<bool> selected;
  vector<int> last_output_iteration;
  static void do_output(int vindex, const char *optstring, void *callback_arg) {
    static_cast<selection_t *>(callback_arg)->selected.at(vindex) = true;
  }

public:
  selection_t() : iteration(-1), times_set(-1) {}
  bool is_selected(const cGH *cctkGH, int vindex) {
    DECLARE_CCTK_PARAMETERS;
    // Check whether the parameter out_vars changed (only once per iteration)
    if (cctkGH->cctk_iteration != iteration) {
      iteration = cctkGH->cctk_iteration;
      int current_times_set =
          CCTK_ParameterQueryTimesSet("out_vars", CCTK_THORNSTRING);
      // Re-scan out_vars (which is somewhat expensive) only if it has changed
      if (current_times_set > times_set) {
        times_set = current_times_set;
        selected.resize(CCTK_NumVars());
        fill(selected.begin(), selected.end(), false);
        CCTK_TraverseString(out_vars, do_output, this, CCTK_GROUP_OR_VAR);
      }
    }
    return selected.at(vindex);
  }
  bool should_output(const cGH *cctkGH, int vindex) {
    last_output_iteration.resize(CCTK_NumVars(), -1);
    return last_output_iteration.at(vindex) < cctkGH->cctk_iteration;
  }
  void did_output(const cGH *cctkGH, int vindex) {
    last_output_iteration.resize(CCTK_NumVars(), -1);
    last_output_iteration.at(vindex) = cctkGH->cctk_iteration;
  }
};
selection_t output_variables;

struct output_state_t {
  unique_ptr<output_file_t> output_file_hdf5_ptr;
  unique_ptr<output_file_t> global_file_hdf5_ptr;
  unique_ptr<output_file_t> output_file_asdf_ptr;
  unique_ptr<output_file_t> global_file_asdf_ptr;
  bool did_write;
  output_state_t(const cGH *cctkGH, io_dir_t io_dir, const string &projectname)
      : did_write(false) {
    DECLARE_CCTK_PARAMETERS;

    // We split processes into "I/O groups" which range from myioproc to
    // myioproc + ioproc_every - 1 (inclusive). Within each group, only
    // one process performs I/O; the other processes then communicate
    // with that process.
    int myproc = CCTK_MyProc(cctkGH);
    int nprocs = CCTK_nProcs(cctkGH);
    int ioproc_every =
        max_nioprocs == 0 ? 1 : (nprocs + max_nioprocs - 1) / max_nioprocs;
    assert(ioproc_every > 0);
    int nioprocs = (nprocs + ioproc_every - 1) / ioproc_every;
    assert(nioprocs > 0);
    assert(max_nioprocs == 0 or nioprocs <= max_nioprocs);
    assert(nioprocs <= nprocs);
    int myioproc = myproc / ioproc_every * ioproc_every;

    if (not output_hdf5 and not output_asdf)
      CCTK_WARN(CCTK_WARN_ALERT, "No file format is active -- no SimulationIO "
                                 "output files be written");

    if (output_hdf5) {
      output_file_hdf5_ptr.reset(
          new output_file_t(cctkGH, io_dir, projectname, file_format::hdf5,
                            file_type::local, false, myioproc, ioproc_every));
      if (output_global_file and myproc == 0)
        global_file_hdf5_ptr.reset(new output_file_t(
            cctkGH, io_dir, projectname, file_format::hdf5, file_type::global,
            false, myioproc, ioproc_every));
    }
    if (output_asdf) {
      output_file_asdf_ptr.reset(
          new output_file_t(cctkGH, io_dir, projectname, file_format::asdf,
                            file_type::local, false, myioproc, ioproc_every));
      if (output_global_file and myproc == 0)
        global_file_asdf_ptr.reset(new output_file_t(
            cctkGH, io_dir, projectname, file_format::asdf, file_type::global,
            false, myioproc, ioproc_every));
    }
  }
  ~output_state_t() { assert(did_write); }
  void insert_vars(const vector<int> &varindices, int reflevel, int timelevel,
                   data_handling handle_data) {
    if (output_file_hdf5_ptr)
      output_file_hdf5_ptr->insert_vars(varindices, reflevel, timelevel,
                                        handle_data);
    if (global_file_hdf5_ptr)
      global_file_hdf5_ptr->insert_vars(varindices, reflevel, timelevel,
                                        handle_data);
    if (output_file_asdf_ptr)
      output_file_asdf_ptr->insert_vars(varindices, reflevel, timelevel,
                                        handle_data);
    if (global_file_asdf_ptr)
      global_file_asdf_ptr->insert_vars(varindices, reflevel, timelevel,
                                        handle_data);
  }
  void write() {
    assert(not did_write);
    if (output_file_hdf5_ptr)
      output_file_hdf5_ptr->write();
    if (global_file_hdf5_ptr)
      global_file_hdf5_ptr->write();
    if (output_file_asdf_ptr)
      output_file_asdf_ptr->write();
    if (global_file_asdf_ptr)
      global_file_asdf_ptr->write();
    did_write = true;
  }
};

unique_ptr<output_state_t> output_state_ptr;

int TriggerVar(const cGH *cctkGH, int vindex, int reflevel,
               data_handling handle_data);
int OutputVar(const cGH *cctkGH, output_state_t &output_state, int vindex,
              int reflevel, data_handling handle_data);

// Callbacks for CarpetSimulationIO's I/O method

int OutputGH(const cGH *cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("OutputGH");

  assert(Carpet::is_level_mode());

  // Output only in global mode
  if (not Carpet::do_global_mode)
    return 0;

  static Timers::Timer timer("SimulationIO::OutputGH");
  timer.start();

  const data_handling handle_data =
      async_output ? data_handling::attach : data_handling::write;

  int numvars = CCTK_NumVars();
  for (int vindex = 0; vindex < numvars; ++vindex)
    if (TimeToOutput(cctkGH, vindex))
      TriggerVar(cctkGH, vindex, -1, handle_data);

  if (output_state_ptr) {
    // output_state_ptr->write();
    // output_state_ptr.reset();
    shared_ptr<output_state_t> osp(std::move(output_state_ptr));
    auto task = [=]() { osp->write(); };
    assert(not output_state_ptr);
    if (async_output) {
      // Output in a new pthreads thread
      // Ensure that the HDF5 library is thread-safe
#ifndef H5_HAVE_THREADSAFE
      CCTK_ERROR("HDF5 library is not thread-safe");
#endif
      // Prevent overlapping I/O. (We could easily allow several
      // overlapping output threads if we wanted do.)
      output_pthreads.clear();
      // Create new thread
      string projectname = generate_projectname();
      output_pthreads.emplace_back(new pthread_wrapper_t(projectname, task));
    } else {
      // Output synchronously
      task();
    }
  }

  timer.stop();

  return 0;
}

int TimeToOutput(const cGH *cctkGH, int vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(Carpet::is_level_mode());

  if (vindex < 0)
    CCTK_ERROR("vindex is negative");

  // if (verbose)
  //   CCTK_VINFO("TimeToOutput variable=\"%s\"",
  //              charptr2string(CCTK_FullName(vindex)).c_str());

  // Output only selected variables
  if (not output_variables.is_selected(cctkGH, vindex))
    return 0;

  // TODO: check reflevel as well
  if (not output_variables.should_output(cctkGH, vindex))
    return 0;

  // whether to output at this iteration
  bool output_this_iteration = false;

  int my_out_every = out_every == -2 ? IO_out_every : out_every;
  int my_out0D_every = out0D_every == -2 ? my_out_every : out0D_every;
  int my_out1D_every = out1D_every == -2 ? my_out_every : out1D_every;
  int my_out2D_every = out2D_every == -2 ? my_out_every : out2D_every;
  int my_out3D_every = out3D_every == -2 ? my_out_every : out3D_every;
  if (my_out0D_every > 0) {
    if (cctk_iteration == *last_out0D_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else if (*last_out0D_iteration < 0 or
               cctk_iteration >=
                   align_up(*last_out0D_iteration + 1, my_out0D_every)) {
      // it is time for the next output
      output_this_iteration = true;
      *last_out0D_iteration = cctk_iteration;
    }
  }
  if (my_out1D_every > 0) {
    if (cctk_iteration == *last_out1D_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else if (*last_out1D_iteration < 0 or
               cctk_iteration >=
                   align_up(*last_out1D_iteration + 1, my_out1D_every)) {
      // it is time for the next output
      output_this_iteration = true;
      *last_out1D_iteration = cctk_iteration;
    }
  }
  if (my_out2D_every > 0) {
    if (cctk_iteration == *last_out2D_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else if (*last_out2D_iteration < 0 or
               cctk_iteration >=
                   align_up(*last_out2D_iteration + 1, my_out2D_every)) {
      // it is time for the next output
      output_this_iteration = true;
      *last_out2D_iteration = cctk_iteration;
    }
  }
  if (my_out3D_every > 0) {
    if (cctk_iteration == *last_out3D_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else if (*last_out3D_iteration < 0 or
               cctk_iteration >=
                   align_up(*last_out3D_iteration + 1, my_out3D_every)) {
      // it is time for the next output
      output_this_iteration = true;
      *last_out3D_iteration = cctk_iteration;
    }
  }

  return output_this_iteration;
}

int TriggerOutput(const cGH *cctkGH, int vindex) {
  DECLARE_CCTK_PARAMETERS;

  assert(Carpet::is_level_mode());

  if (vindex < 0)
    CCTK_ERROR("vindex is negative");

  int grouptype = CCTK_GroupTypeFromVarI(vindex);
  int reflevel = grouptype == CCTK_GF ? Carpet::reflevel : 0;

  if (verbose)
    CCTK_VINFO("TriggerOutput variable=\"%s\" reflevel=%d",
               charptr2string(CCTK_FullName(vindex)).c_str(), reflevel);

  int retval = TriggerVar(cctkGH, vindex, reflevel, data_handling::attach);

  return retval;
}

int OutputVarAs(const cGH *cctkGH, const char *varname, const char *alias) {
  DECLARE_CCTK_PARAMETERS;

  assert(Carpet::is_level_mode());

  if (verbose)
    CCTK_VINFO("OutputVar variable=\"%s\" alias=\"%s\"", varname, alias);

  int vindex = CCTK_VarIndex(varname);
  if (vindex < 0)
    CCTK_VERROR("Unkonwn variable \"%s\"", varname);

  int grouptype = CCTK_GroupTypeFromVarI(vindex);
  int reflevel = grouptype == CCTK_GF ? Carpet::reflevel : 0;

  output_state_t output_state(cctkGH, io_dir_t::output, alias);
  OutputVar(cctkGH, output_state, vindex, reflevel, data_handling::write);
  output_state.write();

  return 0; // no error
}

int TriggerVar(const cGH *cctkGH, int vindex, int reflevel,
               data_handling handle_data) {
  DECLARE_CCTK_PARAMETERS;

  int numvars = CCTK_NumVars();
  if (vindex < 0 or vindex >= numvars)
    CCTK_ERROR("vindex out of range");

  if (verbose)
    CCTK_VINFO("TriggerVar variable=\"%s\"",
               charptr2string(CCTK_FullName(vindex)).c_str());

  if (not output_state_ptr) {
    string projectname = generate_projectname();
    output_state_ptr.reset(
        new output_state_t(cctkGH, io_dir_t::output, projectname));
  }

  int retval =
      OutputVar(cctkGH, *output_state_ptr, vindex, reflevel, handle_data);
  output_variables.did_output(cctkGH, vindex);
  return retval;
}

int OutputVar(const cGH *cctkGH, output_state_t &output_state, int vindex,
              int reflevel, data_handling handle_data) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_VINFO("OutputVar variable=\"%s\" reflevel=%d",
               charptr2string(CCTK_FullName(vindex)).c_str(), reflevel);

  int timelevel = output_all_timelevels ? -1 : 0;
  vector<int> varindices{vindex};
  output_state.insert_vars(varindices, reflevel, timelevel, handle_data);

  return 0; // no error
}

////////////////////////////////////////////////////////////////////////////////

// Checkpointing

void Checkpoint(const cGH *cctkGH, int called_from);

extern "C" void CarpetSimulationIO_InitialDataCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint and checkpoint_ID)
    Checkpoint(cctkGH, CP_INITIAL_DATA);
}

extern "C" void CarpetSimulationIO_EvolutionCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  bool checkpoint_by_iteration =
      checkpoint_every > 0 and
      cctk_iteration >= last_checkpoint_iteration + checkpoint_every;
  bool checkpoint_by_walltime =
      checkpoint_every_walltime_hours > 0 and
      CCTK_RunTime() >= last_checkpoint_walltime_seconds +
                            checkpoint_every_walltime_hours * 3600;
  bool do_checkpoint =
      checkpoint and
      (checkpoint_by_iteration or checkpoint_by_walltime or checkpoint_next);

  if (do_checkpoint)
    Checkpoint(cctkGH, CP_EVOLUTION_DATA);
}

extern "C" void CarpetSimulationIO_TerminationCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint and checkpoint_on_terminate) {

    bool did_checkpoint = cctk_iteration == last_checkpoint_iteration;
    bool do_checkpoint = not did_checkpoint;

    if (do_checkpoint)
      Checkpoint(cctkGH, CP_EVOLUTION_DATA);
  }
}

void Checkpoint(const cGH *cctkGH, int called_from) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("Checkpoint");

  assert(Carpet::is_global_mode());

  int myproc = CCTK_MyProc(cctkGH);
  int nprocs = CCTK_nProcs(cctkGH);

  if (myproc == 0)
    CCTK_VINFO("Writing checkpoint for iteration %d...",
               cctkGH->cctk_iteration);

  vector<int> varindices;
  varindices.reserve(CCTK_NumVars());
  for (int gindex = 0; gindex < CCTK_NumGroups(); ++gindex) {
    // We can't check for storage here, since this requires at least
    // level mode
    if (not CCTK_QueryGroupStorageI(cctkGH, gindex))
      continue;

    cGroup gdata;
    CCTK_GroupData(gindex, &gdata);

    // Output distrib=constant groups only on process 0
    if (gdata.disttype == CCTK_DISTRIB_CONSTANT and myproc != 0)
      continue;

    // Do not checkpoint groups with a "checkpoint=no" tag
    int len = Util_TableGetString(gdata.tagstable, 0, nullptr, "checkpoint");
    if (len > 0) {
      array<char, 10> buf;
      Util_TableGetString(gdata.tagstable, buf.size(), buf.data(),
                          "checkpoint");
      if (CCTK_EQUALS(buf.data(), "no"))
        continue;
      assert(CCTK_EQUALS(buf.data(), "yes"));
    }

    int varindex0 = CCTK_FirstVarIndexI(gindex);
    int numvars = CCTK_NumVarsInGroupI(gindex);
    for (int vindex1 = 0; vindex1 < numvars; ++vindex1)
      varindices.push_back(varindex0 + vindex1);
  }

  string projectname =
      called_from == CP_INITIAL_DATA ? checkpoint_ID_file : checkpoint_file;

  // We split processes into "I/O groups" which range from myioproc to
  // myioproc + ioproc_every - 1 (inclusive). Within each group, only
  // one process performs I/O; the other processes then communicate
  // with that process.
  int ioproc_every =
      max_nioprocs == 0 ? 1 : (nprocs + max_nioprocs - 1) / max_nioprocs;
  assert(ioproc_every > 0);
  int nioprocs = (nprocs + ioproc_every - 1) / ioproc_every;
  assert(nioprocs > 0);
  assert(max_nioprocs == 0 or nioprocs <= max_nioprocs);
  assert(nioprocs <= nprocs);
  int myioproc = myproc / ioproc_every * ioproc_every;

  // TODO: Allow asynchronous checkpointing? We would start writing
  // checkpoint files, but only wait for them at a later time (e.g. a
  // few iterations or a certain time later, probably at least before
  // the next checkpoint).

  // Wait for all output threads to ensure all previous output has
  // been written when the checkpoint file is created
  output_pthreads.clear();

  vector<file_format> formats;
  if (output_hdf5)
    formats.push_back(file_format::hdf5);
  if (output_asdf)
    formats.push_back(file_format::asdf);
  if (formats.empty())
    CCTK_WARN(CCTK_WARN_ALERT, "No file format is active -- no SimulationIO "
                               "checkpoints will be written");

  for (const auto format : formats) {
    {
      output_file_t output_file(cctkGH, io_dir_t::checkpoint, projectname,
                                format, file_type::local, true, myioproc,
                                ioproc_every);
      output_file.insert_vars(varindices, -1, -1, data_handling::write);
      output_file.write();
    }

    if (output_global_file and myproc == 0) {
      output_file_t global_file(cctkGH, io_dir_t::checkpoint, projectname,
                                format, file_type::global, true, myioproc,
                                ioproc_every);
      global_file.insert_vars(varindices, -1, -1, data_handling::write);
      global_file.write();
    }
  }

  last_checkpoint_iteration = cctkGH->cctk_iteration;
  last_checkpoint_walltime_seconds = CCTK_RunTime();

  // Delete old checkpoint files
  if (called_from == CP_EVOLUTION_DATA) {
    static vector<int> checkpoint_iterations;
    // Safety check to prevent double-inserting iterations, which
    // would then keep too few checkpoints
    if (not checkpoint_iterations.empty())
      assert(cctkGH->cctk_iteration > checkpoint_iterations.back());
    checkpoint_iterations.push_back(cctkGH->cctk_iteration);
    if (checkpoint_keep > 0) {
      while (checkpoint_iterations.size() > size_t(checkpoint_keep)) {
        auto iterpos = checkpoint_iterations.begin();
        int iteration = *iterpos;
        if (myproc == 0)
          CCTK_VINFO("Deleting old checkpoint for iteration %d...", iteration);
        for (const auto format : formats) {
          auto filename = generate_filename(
              io_dir_t::checkpoint, projectname, "", iteration, format,
              file_type::local, myioproc, ioproc_every);
          remove(filename.c_str());
          if (output_global_file and myproc == 0) {
            auto filename = generate_filename(
                io_dir_t::checkpoint, projectname, "", iteration, format,
                file_type::global, myioproc, ioproc_every);
            remove(filename.c_str());
          }
        }
        checkpoint_iterations.erase(iterpos);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

// Input

map<string, unique_ptr<input_file_t> > input_file_hdf5_ptrs;
map<string, unique_ptr<input_file_t> > input_file_asdf_ptrs;
void finalise_input_file_ptrs() {
  input_file_hdf5_ptrs.clear();
  input_file_asdf_ptrs.clear();
}

int Input(cGH *cctkGH, const char *basefilename, int called_from) {
  DECLARE_CCTK_PARAMETERS;

  assert(called_from == CP_RECOVER_PARAMETERS or
         called_from == CP_RECOVER_DATA or called_from == FILEREADER_DATA);
  bool do_recover =
      called_from == CP_RECOVER_PARAMETERS or called_from == CP_RECOVER_DATA;
  bool read_parameters = called_from == CP_RECOVER_PARAMETERS;

  if (not do_recover) {
    assert(Carpet::is_level_mode());
    if (Carpet::reflevel > 0)
      return 0; // no error
  }

  static HighResTimer::HighResTimer timer("SimulationIO::input");
  auto timer_clock = timer.start();

  if (do_recover and not read_parameters) {
    // Initialise these in case they should be recovered, but are
    // missing in the checkpoint file. This can happen only if the
    // checkpoint file was written with an older version of this thorn
    // where these variables did not yet exist. Usually, these
    // variables would either be recovered or be initialised, and both
    // happen later, so that the settings here will be overwritten.
    CCTK_INT *restrict const last_out0D_iteration = static_cast<CCTK_INT *>(
        CCTK_VarDataPtr(cctkGH, 0, CCTK_THORNSTRING "::last_out0D_iteration"));
    CCTK_INT *restrict const last_out1D_iteration = static_cast<CCTK_INT *>(
        CCTK_VarDataPtr(cctkGH, 0, CCTK_THORNSTRING "::last_out1D_iteration"));
    CCTK_INT *restrict const last_out2D_iteration = static_cast<CCTK_INT *>(
        CCTK_VarDataPtr(cctkGH, 0, CCTK_THORNSTRING "::last_out2D_iteration"));
    CCTK_INT *restrict const last_out3D_iteration = static_cast<CCTK_INT *>(
        CCTK_VarDataPtr(cctkGH, 0, CCTK_THORNSTRING "::last_out3D_iteration"));
    assert(last_out0D_iteration);
    assert(last_out1D_iteration);
    assert(last_out2D_iteration);
    assert(last_out3D_iteration);
    *last_out0D_iteration = -1;
    *last_out1D_iteration = -1;
    *last_out2D_iteration = -1;
    *last_out3D_iteration = -1;
  }

  string projectname;
  int iteration;
  if (do_recover) {
    // basename is only passed for CP_RECOVER_PARAMETERS, and needs to
    // be remembered
    static string saved_basefilename;
    if (read_parameters)
      // Split basename into the actual basename and the iteration number
      saved_basefilename = basefilename;
    assert(not saved_basefilename.empty());
    auto pi = split_filename(saved_basefilename);
    projectname = get<0>(pi);
    iteration = get<1>(pi);
  } else {
    projectname = basefilename;
    iteration = filereader_ID_iteration;
  }

  if (do_recover)
    if (read_parameters)
      CCTK_VINFO("Recovering parameters for iteration %d from file \"%s\"",
                 iteration, projectname.c_str());
    else
      CCTK_VINFO("Recovering variables for iteration %d from file \"%s\"",
                 iteration, projectname.c_str());
  else
    CCTK_VINFO("Reading iteration %d from file \"%s\"", iteration,
               projectname.c_str());

  if (do_recover and not read_parameters) {
    // Set global Cactus variables
    // CCTK_SetMainLoopIndex(main_loop_index);
    cctkGH->cctk_iteration = iteration;
  }

  // Determine which variables to read
  auto ioUtilGH = static_cast<const ioGH *>(CCTK_GHExtension(cctkGH, "IO"));
  vector<int> input_vars;
  input_vars.reserve(CCTK_NumVars());
  if (do_recover) {
    for (int gindex = 0; gindex < CCTK_NumGroups(); ++gindex) {
      // We can't check for storage here, since this requires at least
      // level mode
      if (not CCTK_QueryGroupStorageI(cctkGH, gindex))
        continue;

      cGroup gdata;
      CCTK_GroupData(gindex, &gdata);

      // Do not recover groups with a "checkpoint=no" tag
      int len = Util_TableGetString(gdata.tagstable, 0, nullptr, "checkpoint");
      if (len > 0) {
        array<char, 10> buf;
        Util_TableGetString(gdata.tagstable, buf.size(), buf.data(),
                            "checkpoint");
        if (CCTK_EQUALS(buf.data(), "no"))
          continue;
        assert(CCTK_EQUALS(buf.data(), "yes"));
      }

      int varindex0 = CCTK_FirstVarIndexI(gindex);
      int numvars = CCTK_NumVarsInGroupI(gindex);
      for (int vindex1 = 0; vindex1 < numvars; ++vindex1)
        input_vars.push_back(varindex0 + vindex1);
    }
  } else {
    for (int vindex = 0; vindex < CCTK_NumVars(); ++vindex)
      if (not ioUtilGH->do_inVars or ioUtilGH->do_inVars[vindex])
        input_vars.push_back(vindex);
  }

  io_dir_t io_dir = do_recover ? io_dir_t::recover : io_dir_t::input;
  bool did_read_parameters = false;
  bool did_read_grid_structure = false;
  // TODO: Emit this warning only if both kinds of files were actually found, in
  // particular when recovering
  // if (bool(input_hdf5) + bool(input_asdf) > 1)
  //   CCTK_WARN(
  //       CCTK_WARN_ALERT,
  //       "Reading both HDF5 and ASDF files -- this usually does not make
  //       sense");
  if (input_hdf5) {
    if (not input_file_hdf5_ptrs.count(projectname)) {
      static HighResTimer::HighResTimer timer1("SimulationIO::read_file_hdf5");
      auto timer1_clock = timer1.start();
      input_file_hdf5_ptrs[projectname] = make_unique<input_file_t>(
          io_dir, projectname, file_format::hdf5, iteration, -1, -1);
      timer1_clock.stop(0);
    }
    const auto &input_file_ptr = input_file_hdf5_ptrs.at(projectname);
    if (read_parameters) {
      static HighResTimer::HighResTimer timer1(
          "SimulationIO::read_parameters_hdf5");
      auto timer1_clock = timer1.start();
      input_file_ptr->read_params();
      did_read_parameters = true;
      timer1_clock.stop(0);
    } else {
      if (not did_read_grid_structure) {
        static HighResTimer::HighResTimer timer1(
            "SimulationIO::read_grid_structure_hdf5");
        auto timer1_clock = timer1.start();
        input_file_ptr->read_grid_structure(cctkGH);
        did_read_grid_structure = true;
        timer1_clock.stop(0);
      }
      static HighResTimer::HighResTimer timer1(
          "SimulationIO::read_variables_hdf5");
      auto timer1_clock = timer1.start();
      input_file_ptr->read_vars(input_vars, -1, -1);
      timer1_clock.stop(0);
    }
  }
  if (input_asdf) {
    if (not input_file_asdf_ptrs.count(projectname)) {
      static HighResTimer::HighResTimer timer1("SimulationIO::read_file_asdf");
      auto timer1_clock = timer1.start();
      input_file_asdf_ptrs[projectname] = make_unique<input_file_t>(
          io_dir, projectname, file_format::asdf, iteration, -1, -1);
      timer1_clock.stop(0);
    }
    const auto &input_file_ptr = input_file_asdf_ptrs.at(projectname);
    if (read_parameters) {
      static HighResTimer::HighResTimer timer1(
          "SimulationIO::read_parameters_asdf");
      auto timer1_clock = timer1.start();
      input_file_ptr->read_params();
      did_read_parameters = true;
      timer1_clock.stop(0);
    } else {
      if (not did_read_grid_structure) {
        static HighResTimer::HighResTimer timer1(
            "SimulationIO::read_grid_structure_asdf");
        auto timer1_clock = timer1.start();
        input_file_ptr->read_grid_structure(cctkGH);
        did_read_grid_structure = true;
        timer1_clock.stop(0);
      }
      static HighResTimer::HighResTimer timer1(
          "SimulationIO::read_variables_asdf");
      auto timer1_clock = timer1.start();
      input_file_ptr->read_vars(input_vars, -1, -1);
      timer1_clock.stop(0);
    }
  }

  timer_clock.stop(0);

  if (read_parameters)
    return did_read_parameters ? 1 : 0;

  assert(did_read_grid_structure);
  return 0; // no error
}

////////////////////////////////////////////////////////////////////////////////

// Recovering

extern "C" int CarpetSimulationIO_RecoverParameters() {
  DECLARE_CCTK_PARAMETERS;
  int iret = 0;
  if (input_hdf5 and iret == 0)
    iret = IOUtil_RecoverParameters(Input, out_extension, "SimulationIO");
  if (input_asdf and iret == 0)
    iret = IOUtil_RecoverParameters(Input, out_extension_asdf, "SimulationIO");
  return iret;
}

} // end namespace CarpetSimulationIO
