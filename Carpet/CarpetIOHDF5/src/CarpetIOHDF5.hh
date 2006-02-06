#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH


#include <hdf5.h>
#include "carpet.hh"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"


// some macros for HDF5 group names
#define METADATA_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"

// atomic HDF5 datatypes for the generic CCTK datatypes
// (the one for CCTK_COMPLEX is created at startup as a compound HDF5 datatype
#define HDF5_CHAR   H5T_NATIVE_CHAR

#ifdef  CCTK_REAL_PRECISION_16
#define HDF5_REAL   H5T_NATIVE_LDOUBLE
#elif   CCTK_REAL_PRECISION_8
#define HDF5_REAL   H5T_NATIVE_DOUBLE
#elif   CCTK_REAL_PRECISION_4
#define HDF5_REAL   H5T_NATIVE_FLOAT
#endif

#ifdef  CCTK_INTEGER_PRECISION_8
#define HDF5_INT    H5T_NATIVE_LLONG
#elif   CCTK_INTEGER_PRECISION_4
#define HDF5_INT    H5T_NATIVE_INT
#elif   CCTK_INTEGER_PRECISION_2
#define HDF5_INT    H5T_NATIVE_SHORT
#elif   CCTK_INTEGER_PRECISION_1
#define HDF5_INT    H5T_NATIVE_CHAR
#endif


// check return code of HDF5 call and print a warning in case of an error
#define HDF5_ERROR(fn_call)                                                   \
        {                                                                     \
          int _error_code = fn_call;                                          \
                                                                              \
                                                                              \
          if (_error_code < 0)                                                \
          {                                                                   \
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,              \
                        "HDF5 call '%s' returned error code %d",              \
                        #fn_call, _error_code);                               \
            error_count++;                                                    \
          }                                                                   \
        }


// CarpetIOHDF5 GH extension structure
typedef struct
{
  // default number of times to output
  int out_every_default;

  // list of variables to output
  char *out_vars;

  // stop on I/O parameter parsing errors ?
  int stop_on_parse_errors;

  // I/O request description list (for all variables)
  vector<ioRequest*> requests;

  // directory in which to output
  char *out_dir;

  // ring buffer for list of successfully created cp files
  int    checkpoint_keep;
  int    cp_filename_index;
  char **cp_filename_list;

  // iteration number of the last checkpoint
  int last_checkpoint_iteration;

  // hdf5 datatype for complex variables; to be set at run time
  hid_t HDF5_COMPLEX, HDF5_COMPLEX8, HDF5_COMPLEX16, HDF5_COMPLEX32;
  
} CarpetIOHDF5GH;


namespace CarpetIOHDF5
{
  // callback routine registered for recovery/filereader
  int Recover (cGH* cctkGH, const char *basefilename, int called_from);

  // worker routines to write a single variable
  int WriteVarUnchunked (const cGH* const cctkGH,
                         hid_t file,
                         const ioRequest* const request,
                         bool called_from_checkpoint);
  int WriteVarChunkedSequential (const cGH* const cctkGH,
                                 hid_t file,
                                 const ioRequest* const request,
                                 bool called_from_checkpoint);
  int WriteVarChunkedParallel (const cGH* const cctkGH,
                               hid_t file,
                               const ioRequest* const request,
                               bool called_from_checkpoint);

  // returns an HDF5 datatype corresponding to the given CCTK datatype
  hid_t CCTKtoHDF5_Datatype (const cGH* const cctkGH,
                             int cctk_type,
                             bool single_precision);

  // scheduled routines (must be declared as C according to schedule.ccl)
  extern "C" {

    int CarpetIOHDF5_Startup (void);
    int CarpetIOHDF5_Init (const cGH* const);
    int CarpetIOHDF5_SetNumRefinementLevels (void);
    int CarpetIOHDF5_CloseFiles (void);
    int CarpetIOHDF5_InitialDataCheckpoint (const cGH* const);
    int CarpetIOHDF5_EvolutionCheckpoint (const cGH* const);
    int CarpetIOHDF5_TerminationCheckpoint (const cGH* const);
    int CarpetIOHDF5_RecoverParameters (void);

  } // extern "C"

} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
