#include "output.hpp"
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
#include <SimulationIO.hpp>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace CarpetSimulationIO {
using namespace std;

// Checkpointing state
static int last_checkpoint_iteration = -1;

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

  *this_iteration = -1;
  *next_output_iteration = 0;

  last_checkpoint_iteration = cctk_iteration;
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

int TriggerVar(const cGH *cctkGH, output_file_t &output_file, int vindex,
               int reflevel, bool global);
int OutputVar(const cGH *cctkGH, output_file_t &output_file, int vindex,
              int reflevel, bool global);

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

  int myproc = CCTK_MyProc(cctkGH); // TODO: use ioproc
  int nprocs = CCTK_nProcs(cctkGH); // TODO: use nioprocs
  string projectname = generate_projectname(cctkGH);
  unique_ptr<output_file_t> output_file_ptr;
  unique_ptr<output_file_t> global_file_ptr;

  int numvars = CCTK_NumVars();
  for (int vindex = 0; vindex < numvars; ++vindex) {
    if (TimeToOutput(cctkGH, vindex)) {
      if (not output_file_ptr)
        output_file_ptr.reset(new output_file_t(cctkGH, io_dir_output,
                                                projectname, myproc, nprocs));
      TriggerVar(cctkGH, *output_file_ptr, vindex, -1, false);
      if (myproc == 0) {
        if (not global_file_ptr)
          global_file_ptr.reset(
              new output_file_t(cctkGH, io_dir_output, projectname, -1, -1));
        TriggerVar(cctkGH, *global_file_ptr, vindex, -1, true);
      }
    }
  }

  if (output_file_ptr)
    output_file_ptr->write();
  if (global_file_ptr)
    global_file_ptr->write();

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

  return output_this_iteration;
}

int TriggerOutput(const cGH *cctkGH, int vindex) {
  DECLARE_CCTK_PARAMETERS;

  assert(Carpet::is_level_mode());

  if (vindex < 0)
    CCTK_ERROR("vindex is negative");

  int myproc = CCTK_MyProc(cctkGH); // TODO: use ioproc
  int nprocs = CCTK_nProcs(cctkGH); // TODO: use nioprocs
  int grouptype = CCTK_GroupTypeFromVarI(vindex);
  int reflevel = grouptype == CCTK_GF ? Carpet::reflevel : 0;

  if (verbose)
    CCTK_VINFO("TriggerOutput variable=\"%s\" refleve=%d",
               charptr2string(CCTK_FullName(vindex)).c_str(), reflevel);

  string projectname = generate_projectname(cctkGH, vindex);
  output_file_t output_file(cctkGH, io_dir_output, projectname, myproc, nprocs);
  int retval = TriggerVar(cctkGH, output_file, vindex, reflevel, false);
  output_file.write();
  if (myproc == 0) {
    output_file_t global_file(cctkGH, io_dir_output, projectname, -1, -1);
    TriggerVar(cctkGH, global_file, vindex, reflevel, true);
    global_file.write();
  }

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

  int myproc = CCTK_MyProc(cctkGH); // TODO: use ioproc
  int nprocs = CCTK_nProcs(cctkGH); // TODO: use nioprocs
  int grouptype = CCTK_GroupTypeFromVarI(vindex);
  int reflevel = grouptype == CCTK_GF ? Carpet::reflevel : 0;

  output_file_t output_file(cctkGH, io_dir_output, alias, myproc, nprocs);
  OutputVar(cctkGH, output_file, vindex, reflevel, false);
  output_file.write();
  if (myproc == 0) {
    output_file_t global_file(cctkGH, io_dir_output, alias, -1, -1);
    OutputVar(cctkGH, global_file, vindex, reflevel, true);
    global_file.write();
  }

  return 0; // no error
}

int TriggerVar(const cGH *cctkGH, output_file_t &output_file, int vindex,
               int reflevel, bool global) {
  DECLARE_CCTK_PARAMETERS;

  int numvars = CCTK_NumVars();
  if (vindex < 0 or vindex >= numvars)
    CCTK_ERROR("vindex out of range");

  if (verbose)
    CCTK_VINFO("TriggerVar variable=\"%s\"",
               charptr2string(CCTK_FullName(vindex)).c_str());

  int retval = OutputVar(cctkGH, output_file, vindex, reflevel, global);
  output_variables.did_output(cctkGH, vindex);
  return retval;
}

int OutputVar(const cGH *cctkGH, output_file_t &output_file, int vindex,
              int reflevel, bool global) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_VINFO("OutputVar variable=\"%s\" reflevel=%d",
               charptr2string(CCTK_FullName(vindex)).c_str(), reflevel);

  int gindex = CCTK_GroupIndexFromVarI(vindex);
  cGroup gdata;
  CCTK_GroupData(gindex, &gdata);

  // Determine number of I/O processes
  int myproc = CCTK_MyProc(cctkGH);
  // int nprocs = CCTK_nProcs(cctkGH);

  // Output distrib=constant groups only on process 0
  if (gdata.disttype == CCTK_DISTRIB_CONSTANT and myproc != 0)
    return 0;

#if 0
  // We split processes into "I/O groups" which range from myioproc
  // to myioproc + ioproc_every - 1 (inclusive). Within each group,
  // at most one process can perform I/O.
  int ioproc_every =
      max_nioprocs == 0 ? 1 : (nprocs + max_nioprocs - 1) / max_nioprocs;
  assert(ioproc_every > 0);
  int nioprocs = (nprocs + ioproc_every - 1) / ioproc_every;
  assert(nioprocs > 0);
  assert(max_nioprocs == 0 or nioprocs <= max_nioprocs);
  assert(nioprocs <= nprocs);
  int myioproc = myproc / ioproc_every * ioproc_every;

  // not yet implemented
  assert(ioproc_every == 1);
#endif

  int timelevel = output_all_timelevels ? -1 : 0;
  vector<int> varindices{vindex};
  output_file.insert_vars(varindices, reflevel, timelevel, global);

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
  bool do_checkpoint =
      checkpoint and (checkpoint_by_iteration or checkpoint_next);

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
      char buf[1000];
      Util_TableGetString(gdata.tagstable, sizeof buf, buf, "checkpoint");
      if (CCTK_EQUALS(buf, "no"))
        continue;
      assert(CCTK_EQUALS(buf, "yes"));
    }

    int varindex0 = CCTK_FirstVarIndexI(gindex);
    int numvars = CCTK_NumVarsInGroupI(gindex);
    for (int vindex1 = 0; vindex1 < numvars; ++vindex1)
      varindices.push_back(varindex0 + vindex1);
  }

  string projectname =
      called_from == CP_INITIAL_DATA ? checkpoint_ID_file : checkpoint_file;
  output_file_t output_file(cctkGH, io_dir_checkpoint, projectname,

                            myproc, nprocs);
  output_file.insert_vars(varindices, -1, -1, false);
  output_file.write();

  if (myproc == 0) {
    output_file_t global_file(cctkGH, io_dir_checkpoint, projectname, -1, -1);
    global_file.insert_vars(varindices, -1, -1, true);
    global_file.write();
  }

  last_checkpoint_iteration = cctkGH->cctk_iteration;

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
        auto filename = generate_filename(
            cctkGH, io_dir_checkpoint, projectname, iteration, myproc, nprocs);
        remove(filename.c_str());
        if (myproc == 0) {
          auto filename = generate_filename(cctkGH, io_dir_checkpoint,
                                            projectname, iteration, -1, -1);
          remove(filename.c_str());
        }
        checkpoint_iterations.erase(iterpos);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

// Input

// Recovery information
static int recovery_iteration = -1;

int Input(cGH *cctkGH, const char *basefilename, int called_from) {
  assert(called_from == CP_RECOVER_PARAMETERS or
         called_from == CP_RECOVER_DATA or called_from == FILEREADER_DATA);
  bool do_recover =
      called_from == CP_RECOVER_PARAMETERS or called_from == CP_RECOVER_DATA;

  if (not do_recover) {
    assert(Carpet::is_level_mode());
    if (Carpet::reflevel > 0)
      return 0; // no error
  }

  if (do_recover)
    CCTK_VINFO("Recovering iteration %d", recovery_iteration);
  else
    CCTK_VINFO("Reading iteration %d", cctkGH->cctk_iteration);

  // Determine which variables to read
  const ioGH *ioUtilGH = (const ioGH *)CCTK_GHExtension(cctkGH, "IO");
  vector<int> input_vars;
  int numvars = CCTK_NumVars();
  for (int vindex = 0; vindex < numvars; ++vindex)
    if (not ioUtilGH->do_inVars or ioUtilGH->do_inVars[vindex])
      input_vars.push_back(vindex);

  return 0; // no error
}

} // end namespace CarpetSimulationIO
