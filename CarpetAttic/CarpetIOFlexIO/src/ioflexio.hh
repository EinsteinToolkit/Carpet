// $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIO/src/ioflexio.hh,v 1.7 2004/02/07 16:21:56 schnetter Exp $

#ifndef CARPETIOFLEXIO_HH
#define CARPETIOFLEXIO_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "ioflexio.h"

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
  
  int Recover (cGH* const cgh, const char *basefilename,
               const int called_from);
  
} // namespace CarpetIOFlexIO

#endif // !defined(CARPETIOFLEXIO_HH)
