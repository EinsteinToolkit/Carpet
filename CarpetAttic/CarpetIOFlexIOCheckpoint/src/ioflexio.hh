// $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/ioflexio.hh,v 1.8 2004/01/06 08:56:39 cott Exp $

#ifndef CARPETIOFLEXIO_HH
#define CARPETIOFLEXIO_HH


#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "ioflexio.h"
#include "ioflexioGH.h"

/* define the IOFlexIO datatypes according to CCTK_??? datatypes */
#define FLEXIO_CHAR    IObase::Char

#ifdef  CCTK_INT2
#define FLEXIO_INT2    IObase::Int16
#endif
#ifdef  CCTK_INT4
#define FLEXIO_INT4    IObase::Int32
#endif
#ifdef  CCTK_INT8
#define FLEXIO_INT8    IObase::Int64
#endif

#ifdef  CCTK_REAL4
#define FLEXIO_REAL4   IObase::Float32
#endif
#ifdef  CCTK_REAL8
#define FLEXIO_REAL8   IObase::Float64
#endif
#ifdef  CCTK_REAL16
#define FLEXIO_REAL16  -1
#endif

/* define the FlexIO types for the generic CCTK_INT and CCTK_REAL datatypes */
#ifdef  CCTK_INTEGER_PRECISION_8
#define FLEXIO_INT     IObase::Int64
#elif   CCTK_INTEGER_PRECISION_4
#define FLEXIO_INT     IObase::Int32
#elif   CCTK_INTEGER_PRECISION_2
#define FLEXIO_INT     IObase::Int16
#endif

#ifdef  CCTK_REAL_PRECISION_4
#define FLEXIO_REAL    FLEXIO_REAL4
#elif   CCTK_REAL_PRECISION_8
#define FLEXIO_REAL    FLEXIO_REAL8
#elif   CCTK_REAL_PRECISION_16
#define FLEXIO_REAL    FLEXIO_REAL16
#endif

/* some macros needed for recovery */
#ifdef CCTK_MPI

#define CACTUS_MPI_ERROR(fn_call)                                             \
          do {                                                                \
            int errcode;                                                      \
                                                                              \
            if ((errcode = fn_call) != MPI_SUCCESS)                           \
            {                                                                 \
              char mpi_error_string[MPI_MAX_ERROR_STRING+1];                  \
              int resultlen;                                                  \
                                                                              \
              MPI_Error_string (errcode, mpi_error_string, &resultlen);       \
              fprintf (stderr, "MPI call '%s' returned error code %d (%s)\n", \
                               #fn_call, errcode, mpi_error_string);          \
              fprintf(stderr, "At line %d of file %s\n", __LINE__, __FILE__); \
            }                                                                 \
          } while (0)


#ifdef  CCTK_INT4
#define CARPET_MPI_INT4  (sizeof (CCTK_INT4) == sizeof (int) ? MPI_INT :        \
                        sizeof (CCTK_INT4) == sizeof (short) ? MPI_SHORT :    \
                        MPI_DATATYPE_NULL)
#endif

#define CARPET_MPI_CHAR      MPI_CHAR

#endif

namespace CarpetIOFlexIO {
  
  // Variable definitions
  extern int GHExtension;
  extern int IOMethod;
  extern vector<bool> do_truncate;
  extern vector<vector<int> > last_output;
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cgh);
  
  int OutputGH (const cGH* const cgh);
  int OutputVarAs (const cGH* const cgh, const char* const varname,
		   const char* const alias);
  int TimeToOutput (const cGH* const cgh, const int vindex);
  int TriggerOutput (const cGH* const cgh, const int vindex);
  
  int InputGH (const cGH* const cgh);
  int InputVarAs (const cGH* const cgh, const char* const varname,
		  const char* const alias);
  
  static const char* GetStringParameter (const char* const parametername,
					 const char* const fallback);

  int WriteGF (const cGH* const cgh, IObase* writer, AMRwriter* amrwriter, ioRequest* request);
  int WriteGS (const cGH* const cgh, IObase* writer, ioRequest* request);

} // namespace CarpetIOFlexIO

namespace CarpetIOFlexIOUtil {

  IObase::DataType FlexIODataType (int cctk_type);
  void DumpCommonAttributes (const cGH *cgh, IObase* writer, ioRequest* request);

}
#endif // !defined(CARPETIOFLEXIO_HH)

/* structure holding necessary information about a recovery file */


