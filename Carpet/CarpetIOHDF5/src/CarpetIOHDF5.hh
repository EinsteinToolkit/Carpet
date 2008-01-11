#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH

// some macros to fix compatibility issues as long
// as 1.8.0 is in beta phase
#define H5_USE_16_API 1

#include <hdf5.h>

#if (H5_VERS_MAJOR == 1 && (H5_VERS_MINOR == 8) && (H5_VERS_RELEASE == 0))
#warning "Hacking HDF5 1.8.0 compatiblity with 1.6.x; fix once 1.8.0 stable"
#warning "Hacking HDF5 1.8.0 compatiblity with 1.6.x; fix once 1.8.0 stable"
#warning "Hacking HDF5 1.8.0 compatiblity with 1.6.x; fix once 1.8.0 stable"
#warning "Hacking HDF5 1.8.0 compatiblity with 1.6.x; fix once 1.8.0 stable"
#warning "Hacking HDF5 1.8.0 compatiblity with 1.6.x; fix once 1.8.0 stable"
#endif


#include "cctk_Arguments.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "carpet.hh"





// some macros for HDF5 group names
#define METADATA_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"
#define GRID_STRUCTURE "Grid Structure"

// atomic HDF5 datatypes for the generic CCTK datatypes
// (the one for CCTK_COMPLEX is created at startup as a compound HDF5 datatype)
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
        do {                                                                  \
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
        } while (0)


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

  // list of recovery files to remove
  int    recovery_num_filenames;
  char **recovery_filename_list;

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
                         long long & io_bytes,
                         const ioRequest* const request,
                         bool called_from_checkpoint);
  int WriteVarChunkedSequential (const cGH* const cctkGH,
                                 hid_t file,
                                 long long & io_bytes,
                                 const ioRequest* const request,
                                 bool called_from_checkpoint);
  int WriteVarChunkedParallel (const cGH* const cctkGH,
                               hid_t file,
                               long long & io_bytes,
                               const ioRequest* const request,
                               bool called_from_checkpoint);

  // returns an HDF5 datatype corresponding to the given CCTK datatype
  hid_t CCTKtoHDF5_Datatype (const cGH* const cctkGH,
                             int cctk_type,
                             bool single_precision);

  // scheduled routines (must be declared as C according to schedule.ccl)
  extern "C" {

    int CarpetIOHDF5_Startup (void);
    void CarpetIOHDF5_Init (CCTK_ARGUMENTS);
    int CarpetIOHDF5_SetNumRefinementLevels (void);
    void CarpetIOHDF5_CloseFiles (CCTK_ARGUMENTS);
    void CarpetIOHDF5_InitialDataCheckpoint (CCTK_ARGUMENTS);
    void CarpetIOHDF5_EvolutionCheckpoint (CCTK_ARGUMENTS);
    void CarpetIOHDF5_TerminationCheckpoint (CCTK_ARGUMENTS);
    int CarpetIOHDF5_RecoverParameters (void);
    void CarpetIOHDF5_RecoverGridStructure (CCTK_ARGUMENTS);

  } // extern "C"

} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
