#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "carpet.hh"

#include "CactusBase/IOUtil/src/ioutil_Utils.h"

/* some macros for HDF5 group names */
#define METADATA_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"

// Carpet version ID to tag data output and checkpoint files
#define CARPET_VERSION	1

// Some MPI Datatypes we need for Recovery
// Originally written by Thomas Radke.

#ifdef  CCTK_INT4
#define CARPET_MPI_INT4  (sizeof (CCTK_INT4) == sizeof (int) ? MPI_INT :        \
                        sizeof (CCTK_INT4) == sizeof (short) ? MPI_SHORT :    \
                        MPI_DATATYPE_NULL)
#endif

#define CARPET_MPI_CHAR      MPI_CHAR

/* floating point types are architecture-independent,
   ie. a float has always 4 bytes, and a double has 8 bytes

   PUGH_MPI_REAL  is used for communicating reals of the generic CCTK_REAL type
   PUGH_MPI_REALn is used to explicitely communicate n-byte reals */
#ifdef  CCTK_REAL4
#define CARPET_MPI_REAL4  MPI_FLOAT
#endif
#ifdef  CCTK_REAL8
#define CARPET_MPI_REAL8  MPI_DOUBLE
#endif
#ifdef  CCTK_REAL16
#define CARPET_MPI_REAL16  (sizeof (CCTK_REAL16) == sizeof (long double) ?      \
                          MPI_LONG_DOUBLE : MPI_DATATYPE_NULL)
#endif


#ifdef  CCTK_REAL_PRECISION_16
#define CARPET_MPI_REAL   CARPET_MPI_REAL16
#elif   CCTK_REAL_PRECISION_8
#define CARPET_MPI_REAL   CARPET_MPI_REAL8
#elif   CCTK_REAL_PRECISION_4
#define CARPET_MPI_REAL   CARPET_MPI_REAL4
#endif


/*** Define the different datatypes used for HDF5 I/O
     NOTE: the complex datatype SHOULD be [is] defined dynamically at runtime in Startup.c
     100% of the definitions below were taken from Thomas Radke's IOHDF5Util thorn for PUGH
 ***/
/* char type is easy */
#define HDF5_CHAR   H5T_NATIVE_CHAR

/* floating point types are architecture-independent,
   ie. a float has always 4 bytes, and a double has 8 bytes
     HDF5_REAL  is used for storing reals of the generic CCTK_REAL type
     HDF5_REALn is used to explicitely store n-byte reals */
#ifdef  CCTK_REAL4
#define HDF5_REAL4  H5T_NATIVE_FLOAT
#endif
#ifdef  CCTK_REAL8
#define HDF5_REAL8  H5T_NATIVE_DOUBLE
#endif
#ifdef  CCTK_REAL16
#define HDF5_REAL16 (sizeof (CCTK_REAL16) == sizeof (long double) ?           \
                     H5T_NATIVE_LDOUBLE : -1)
#endif


#ifdef  CCTK_REAL_PRECISION_16
#define HDF5_REAL   HDF5_REAL16
#elif   CCTK_REAL_PRECISION_8
#define HDF5_REAL   HDF5_REAL8
#elif   CCTK_REAL_PRECISION_4
#define HDF5_REAL   HDF5_REAL4
#endif


/* integer types are architecture-dependent:
     HDF5_INT  is used for communicating integers of the generic CCTK_INT type
     HDF5_INTn is used to explicitely communicate n-byte integers */
#ifdef  CCTK_INT8
#define HDF5_INT8   (sizeof (CCTK_INT8) == sizeof (int) ? H5T_NATIVE_INT :    \
                     sizeof (CCTK_INT8) == sizeof (long) ? H5T_NATIVE_LONG :  \
                     sizeof (CCTK_INT8) == sizeof (long long) ?               \
                     H5T_NATIVE_LLONG : -1)
#endif

#ifdef  CCTK_INT4
#define HDF5_INT4   (sizeof (CCTK_INT4) == sizeof (int) ? H5T_NATIVE_INT :    \
                     sizeof (CCTK_INT4) == sizeof (short) ?                   \
                     H5T_NATIVE_SHORT : -1)
#endif

#ifdef  CCTK_INT2
#define HDF5_INT2   (sizeof (CCTK_INT2) == sizeof (short) ?                   \
                     H5T_NATIVE_SHORT : -1)
#endif

#ifdef  CCTK_INT1
#define HDF5_INT1   H5T_NATIVE_CHAR
#endif

#ifdef  CCTK_INTEGER_PRECISION_8
#define HDF5_INT    HDF5_INT8
#elif   CCTK_INTEGER_PRECISION_4
#define HDF5_INT    HDF5_INT4
#elif   CCTK_INTEGER_PRECISION_2
#define HDF5_INT    HDF5_INT2
#elif   CCTK_INTEGER_PRECISION_1
#define HDF5_INT    HDF5_INT1
#endif

/* Nice error handling. Stolen from Thomas Radke */
/* check return code of HDF5 call and print a warning in case of an error */
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
          }                                                                   \
        }


/* CarpetIOHDF5 GH extension structure */
typedef struct
{
  /* default number of times to output */
  int out_every_default;

  /* number of times to output for each variable */
  CCTK_INT *out_every;

  /* the last iteration output for each variable */
  int *out_last;

  /* list of variables to output */
  char *out_vars;

  /* stop on I/O parameter parsing errors ? */
  int stop_on_parse_errors;

  /* I/O request description list (for all variables) */
  ioRequest **requests;

  /* directory in which to output */
  char *out_dir;

  /* timer array for checkpointing/recovery */
  // int timers[IOHDF5_NUM_TIMERS];

  /* flag to indicate request for timer output */
  // int print_timing_info;

  /* ring buffer for list of successfully created cp files */
  int    cp_filename_index;
  char **cp_filename_list;

  /* iteration number of the last checkpoint */
  int last_checkpoint_iteration;

  /* hdf5 datatype for stupid complex variables; to be set at run time */
  hid_t HDF5_COMPLEX, HDF5_COMPLEX8, HDF5_COMPLEX16, HDF5_COMPLEX32;
  
} CarpetIOHDF5GH;


namespace CarpetIOHDF5 {
  
  using namespace std;
  using namespace Carpet;
  
  // Variable definitions
  extern vector<bool> do_truncate; // [var]
  extern vector<vector<vector<int> > > last_output; // [ml][rl][var]
  
  int WriteVar (const cGH* const cctkGH, const hid_t writer, const ioRequest* request,
		   const int called_from_checkpoint);

  int InputGH (const cGH* const cctkGH);
  int ReadVar (const cGH* const cctkGH, const int vindex,
	       const hid_t currdataset, vector<ibset> &regions_read, 
	       const int called_from_recovery);

  int Recover (cGH* cgh, const char *basefilename, int called_from);

  // auxiliary functions defined in iohdf5utils.cc

  void WriteAttribute (const hid_t dataset, const char* name, int value);
  void WriteAttribute (const hid_t dataset, const char* name, const int* values, int nvalues);
  void WriteAttribute (const hid_t dataset, const char* name, double value);
  void WriteAttribute (const hid_t dataset, const char* name, const double* values, int nvalues);
  void WriteAttribute (const hid_t dataset, const char* name, char value);
  void WriteAttribute (const hid_t dataset, const char* name, const char* values);
  void WriteAttribute (const hid_t dataset, const char* name, const char* values, int nvalues);
  
  int ReadAttribute (const hid_t dataset, const char* name, int& value);
  int ReadAttribute (const hid_t dataset, const char* name, int* values, int nvalues);
  int ReadAttribute (const hid_t dataset, const char* name, double& value);
  int ReadAttribute (const hid_t dataset, const char* name, double* values, int nvalues);
  int ReadAttribute (const hid_t dataset, const char* name, char& value);
  int ReadAttribute (const hid_t dataset, const char* name, char*& values);
  int ReadAttribute (const hid_t dataset, const char* name, char* values, int nvalues);
  
  int GetnDatasets (const hid_t reader);
  void GetDatasetName (const hid_t reader, const int _index, char* name);

  hid_t h5DataType (const cGH* const cctkGH, int cctk_type,
                    int single_precision);

extern "C" {

// scheduled routines
int CarpetIOHDF5_Startup (void);
int CarpetIOHDF5_Init (const cGH* const);
int CarpetIOHDF5_ReadData (const cGH* const);
int CarpetIOHDF5_CloseFile (void);
int CarpetIOHDF5_InitialDataCheckpoint (const cGH* const);
int CarpetIOHDF5_EvolutionCheckpoint (const cGH* const);
int CarpetIOHDF5_TerminationCheckpoint (const cGH* const);

// routines registered for recovery
int CarpetIOHDF5_Recover (cGH* cgh, const char *basefilename, int called_from);
int CarpetIOHDF5_RecoverParameters (void);

} // extern "C"

} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
