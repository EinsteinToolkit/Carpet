// $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIOCheckpoint/src/ioflexio.hh,v 1.1 2003/05/16 14:02:18 hawke Exp $

#ifndef CARPETIOFLEXIO_HH
#define CARPETIOFLEXIO_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "ioflexio.h"

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
  
  const char* GetStringParameter (const char* const parametername,
					 const char* const fallback);

  int WriteVarAs (const cGH* const cgh, IObase* writer,AMRwriter* amrwriter, int varindex, int group);

} // namespace CarpetIOFlexIO

#endif // !defined(CARPETIOFLEXIO_HH)

