#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <hdf5.h>
#include <iof5.hh>

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"



namespace CarpetIOF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Checkpointing state
  static int last_checkpoint_iteration = -1;
  
  
  
  // Scheduled startup routine
  int CarpetIOF5_Startup ()
  {
    CCTK_RegisterBanner ("AMR F5 I/O provided by CarpetIOF5");
    
    int const GHExtension = CCTK_RegisterGHExtension (CCTK_THORNSTRING);
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
    
    return 0;
  }
  
  // Registered GH extension setup routine
  void* SetupGH (tFleshConfig* const fleshconfig,
                 int const convLevel, cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // register I/O method
    int const IOMethod = CCTK_RegisterIOMethod ("IOF5");
    CCTK_RegisterIOMethodOutputGH      (IOMethod, OutputGH     );
    CCTK_RegisterIOMethodTimeToOutput  (IOMethod, TimeToOutput );
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
    CCTK_RegisterIOMethodOutputVarAs   (IOMethod, OutputVarAs  );
    
    // there no actual extension data structure
    return NULL;
  }
  
  // Scheduled initialisation routine
  void CarpetIOF5_Init (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    *this_iteration        = -1;
    *next_output_iteration =  0;
    
    last_checkpoint_iteration = cctk_iteration;
  }
  
  
  
  // Callbacks for CarpetIOHDF5's I/O method

  struct callback_arg_t {
    cGH const* cctkGH;
  };
  void do_output (int const vindex,
                  char const* const optstring,
                  void* const callback_arg);
  
  int OutputGH (cGH const* const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    static Carpet::Timer timer ("F5::OutputGH");
    timer.start();
    
    callback_arg_t callback_arg;
    callback_arg.cctkGH = cctkGH;
    CCTK_TraverseString (out_vars, do_output, &callback_arg, CCTK_GROUP_OR_VAR);
    
    timer.stop(0);
    
    return 0;
  }
  
  void do_output (int const vindex,
                  char const* const optstring,
                  void* const callback_arg_)
  {
    callback_arg_t& callback_arg =
      * static_cast<callback_arg_t*>(callback_arg_);
    cGH const* const cctkGH = callback_arg.cctkGH;
    if (TimeToOutput (cctkGH, vindex)) {
      TriggerOutput (cctkGH, vindex);
    }
  }
  
  int TimeToOutput (cGH const* const cctkGH, int const vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    int const numvars = CCTK_NumVars();
    assert (vindex>=0 and vindex<numvars);
    
    // Output only in global mode
    if (not do_global_mode) return 0;
    
    // whether to output at this iteration
    bool output_this_iteration = false;
    
    int const my_out_every = out_every == -2 ? IO_out_every : out_every;
    if (my_out_every > 0) {
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
  
  int TriggerOutput (cGH const* const cctkGH, int const vindex)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char* const fullname = CCTK_FullName(vindex);
    int const gindex = CCTK_GroupIndexFromVarI(vindex);
    char* const groupname = CCTK_GroupName(gindex);
    for (char* p=groupname; *p; ++p) *p=tolower(*p);
    int const retval = OutputVarAs (cctkGH, fullname, groupname);
    free (groupname);
    free (fullname);
    
    return retval;
  }
  
  int OutputVarAs (cGH const* const cctkGH,
                   const char* const varname, const char* const alias)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_level_mode());
    BEGIN_GLOBAL_MODE(cctkGH) {
      DECLARE_CCTK_ARGUMENTS;
      
      CCTK_VInfo (CCTK_THORNSTRING,
                  "F5::OutputVarAs: iteration=%d", cctk_iteration);
      
      
      
      // We don't know how to open multiple files yet
      assert (CCTK_EQUALS (file_content, "everything"));
      
      // Open file
      static bool first_time = true;
      
      // The file name doesn't matter since we currently write
      // everything into a single file
      int const vindex = CCTK_VarIndex (varname);
      assert (vindex >= 0);
      string const basename = generate_basename (cctkGH, vindex);
      int const myproc = CCTK_MyProc(cctkGH);
      int const proc = myproc;
      string const name = create_filename (cctkGH, basename, proc, first_time);
      
      indent_t indent;
      cout << indent << "process=" << proc << "\n";
      
      bool const truncate_file = first_time and IO_TruncateOutputFiles(cctkGH);
      hid_t const file =
        truncate_file ?
        H5Fcreate (name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) :
        H5Fopen   (name.c_str(), H5F_ACC_RDWR , H5P_DEFAULT);
      assert (file >= 0);
      first_time = false;
      
      vector<bool> output_var(CCTK_NumVars());
      output_var.at(vindex) = true;
      output (cctkGH, file, output_var, false);
      
      // Close file
      herr_t const herr = H5Fclose (file);
      assert (not herr);
      
    } END_GLOBAL_MODE;
    
    return 0;                   // no error
  }
  
  
  
  // Checkpointing
  
  void Checkpoint (cGH const* const cctkGH, int const called_from)
  {
    assert (is_global_mode());
    
    // generate filenames for both the temporary and real checkpoint
    // files
    int const ioproc = CCTK_MyProc(cctkGH);
    int const parallel_io = 1;
    char* const filename =
      IOUtil_AssembleFilename (cctkGH, NULL, "", ".f5",
                               called_from, ioproc, not parallel_io);
    char* const tempname =
      IOUtil_AssembleFilename (cctkGH, NULL, ".tmp", ".f5",
                               called_from, ioproc, not parallel_io);
    
    hid_t const file =
      H5Fcreate (tempname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert (file >= 0);
    
    vector<bool> output_var(CCTK_NumVars());
    for (int gindex=0; gindex<CCTK_NumGroups(); ++gindex) {
      // only checkpoint groups with storage
      if (CCTK_QueryGroupStorageI(cctkGH, gindex) <= 0) continue;
      if (CCTK_NumVarsInGroupI(gindex) == 0) continue;
      
      // do not checkpoint groups with a "checkpoint=no" tag
      cGroup gdata;
      CCTK_GroupData (gindex, &gdata);
      int const len =
        Util_TableGetString (gdata.tagstable, 0, NULL, "checkpoint");
      if (len > 0) {
        char buf[1000];
        Util_TableGetString (gdata.tagstable, sizeof buf, buf, "checkpoint");
        if (CCTK_EQUALS(buf, "no")) continue;
        assert (CCTK_EQUALS(buf, "yes"));
      }
      
      int const first_vindex = CCTK_FirstVarIndexI (gindex);
      int const num_vars = CCTK_NumVarsInGroupI (gindex);
      for (int vindex=first_vindex; vindex<first_vindex+num_vars; ++vindex) {
        output_var.at(vindex) = true;
      }
    }
    
    output (cctkGH, file, output_var, true);
    
    // Close file
    herr_t const herr = H5Fclose (file);
    assert (not herr);
    
    // Wait until all files have been written, then rename the
    // checkpoint files
    // TODO: ensure there were no errors
    CCTK_Barrier (cctkGH);
    int const ierr = rename (tempname, filename);
    assert (not ierr);
  }
  
  void CarpetIOF5_InitialDataCheckpoint (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (checkpoint and checkpoint_ID) {
      Checkpoint (cctkGH, CP_INITIAL_DATA);
    }
  }
  
  void CarpetIOF5_EvolutionCheckpoint (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    bool const checkpoint_by_iteration =
      checkpoint_every > 0 and
      cctk_iteration >= last_checkpoint_iteration + checkpoint_every;
    bool const do_checkpoint =
      checkpoint and (checkpoint_by_iteration or checkpoint_next);
    
    if (do_checkpoint) {
      Checkpoint (cctkGH, CP_EVOLUTION_DATA);
    }
  }
  
  void CarpetIOF5_TerminationCheckpoint (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (checkpoint and checkpoint_on_terminate) {
      
      bool const did_checkpoint = cctk_iteration == last_checkpoint_iteration;
      bool const do_checkpoint = not did_checkpoint;
      
      if (do_checkpoint) {
        Checkpoint (cctkGH, CP_EVOLUTION_DATA);
      }
    }
  }
  
} // end namespace CarpetIOF5
