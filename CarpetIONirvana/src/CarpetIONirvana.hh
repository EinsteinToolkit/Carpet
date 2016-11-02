#ifndef CARPETIONIRVANA_HH
#define CARPETIONIRVANA_HH

#include "CarpetN5.hh"
#include "metadata.hh"

#include <vector>

#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "carpet.hh"
#include "cctk_Arguments.h"

// CarpetIONirvana GH extension structure
typedef struct {
  // default number of times to output
  int out_every_default;

  // list of variables to output
  char *out_vars;

  // stop on I/O parameter parsing errors ?
  int stop_on_parse_errors;

  // I/O request description list (for all variables)
  vector<ioRequest *> requests;

  // directory in which to output
  char *out_dir;

  // ring buffer for list of successfully created cp files
  int checkpoint_keep;
  int cp_filename_index;
  char **cp_filename_list;

  // list of recovery files to remove
  int recovery_num_filenames;
  char **recovery_filename_list;

  // iteration number of the last checkpoint
  int last_checkpoint_iteration;

  // hdf5 datatype for complex variables; to be set at run time
  hid_t HDF5_COMPLEX, HDF5_COMPLEX8, HDF5_COMPLEX16, HDF5_COMPLEX32;

} CarpetIONirvanaGH;

namespace CarpetIONirvana {

// worker routines to write a single variable
int WriteVar(const cGH *const cctkGH, const string &filename, const int filenum,
             CCTK_REAL &io_bytes, const ioRequest *const request);

// returns an HDF5 datatype corresponding to the given CCTK datatype
hid_t CCTKtoHDF5_Datatype(const cGH *const cctkGH, int cctk_type,
                          bool single_precision);

// Everything is a class template, so that it can easily be
// instantiated for all output dimensions

/*  template<int outdim>
  struct IONirvana {

    // name of the output directory
    static char* my_out_slice_dir;

    // list of variables to output
    static char* my_out_slice_vars;

    // I/O request description list (for all variables)
    static vector<ioRequest*> slice_requests;



    // Scheduled functions
    static int Startup();

    // Registered functions
    static void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cctkGH);

    static int OutputGH (const cGH* cctkGH);
    static int OutputVarAs (const cGH* cctkGH,
                            const char* varname, const char* alias);
    static int TimeToOutput (const cGH* cctkGH, int vindex);
    static int TriggerOutput (const cGH* cctkGH, int vindex);

    // Other functions
    static void CheckSteerableParameters (const cGH* cctkGH);

    static bool DidOutput (const cGH* cctkGH,
                           int vindex,
                           string basefilename,
                           bool& is_new_file, bool& truncate_file);

    static bool DirectionIsRequested (const vect<int,outdim>& dirs);

    static void OutputDirection (const cGH* cctkGH,
                                 int vindex,
                                 string alias,
                                 string basefilename,
                                 const vect<int,outdim>& dirs,
                                 bool is_new_file,
                                 bool truncate_file);


    static ivect GetOutputOffset (const cGH* cctkGH, int m,
                                  const vect<int,outdim>& dirs);

  };*/ // struct IONirvana

// scheduled routines (must be declared as C according to schedule.ccl)
extern "C" {

int CarpetIONirvana_RecoverParameters(void);
int CarpetIONirvana_SetNumRefinementLevels(void);
int CarpetIONirvana_Startup(void);
void CarpetIONirvana_Init(CCTK_ARGUMENTS);
void CarpetIONirvana_InitCheckpointingIntervals(CCTK_ARGUMENTS);
void CarpetIONirvana_RecoverGridStructure(CCTK_ARGUMENTS);
void CarpetIONirvana_CloseFiles(CCTK_ARGUMENTS);
void CarpetIONirvana_InitialDataCheckpoint(CCTK_ARGUMENTS);
void CarpetIONirvana_EvolutionCheckpoint(CCTK_ARGUMENTS);
void CarpetIONirvana_TerminationCheckpoint(CCTK_ARGUMENTS);

} // extern "C"

} // namespace CarpetIONirvana

#endif // !defined(CarpetIONirvana_HH)
