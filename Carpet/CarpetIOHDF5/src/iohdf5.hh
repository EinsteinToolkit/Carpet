// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5.hh,v 1.9 2004/04/03 12:40:21 schnetter Exp $

#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "carpet.hh"

#include "iohdf5.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

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


namespace CarpetIOHDF5 {
  
  using namespace std;
  using namespace Carpet;
  
  // Variable definitions
  extern int GHExtension;
  extern int IOMethod;
  extern vector<bool> do_truncate; // [var]
  extern vector<vector<vector<int> > > last_output; // [ml][rl][var]
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cctkGH);
  
  int OutputGH (const cGH* const cctkGH);
  int WriteVar (const cGH* const cctkGH, const hid_t writer, const ioRequest* request,
		   const int called_from_checkpoint);
  int OutputVarAs (const cGH* const cctkGH, const char* const varname,
		   const char* const alias);

  int TimeToOutput (const cGH* const cctkGH, const int vindex);
  int TriggerOutput (const cGH* const cctkGH, const int vindex);
  
  int InputGH (const cGH* const cctkGH);
  int ReadVar (const cGH* const cctkGH, const hid_t reader, const char* const varname,
	       const hid_t currdataset, vector<ibset> &regions_read, 
	       const int called_from_recovery);

  int InputVarAs (const cGH* const cctkGH, const char* const varname,
		  const char* const alias);
  
  int Recover (cGH* const cctkGH, const char *basefilename,
               const int called_from);
  
  int CarpetIOHDF5_Recover (cGH* cgh, const char *basefilename, int called_from);

  // auxiliary functions defined in iohdf5utils.cc

  bool CheckForVariable (const cGH* const cctkGH,
                         const char* const varlist, const int vindex);
  void SetFlag (int index, const char* optstring, void* arg);
  
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

  hid_t h5DataType(int cctk_type);

} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
