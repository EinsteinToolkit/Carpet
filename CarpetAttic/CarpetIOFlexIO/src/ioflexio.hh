// $Header: /home/eschnett/C/carpet/Carpet/CarpetAttic/CarpetIOFlexIO/src/ioflexio.hh,v 1.3 2001/03/17 22:37:28 eschnett Exp $

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

namespace CarpetIOFlexIO {
  
  // Variable definitions
  extern int GHExtension;
  extern int IOMethod;
  extern vector<bool> do_truncate;
  extern vector<vector<int> > last_output;
  
  // Scheduled functions
  extern "C" {
    int CarpetIOFlexIOStartup ();
    int CarpetIOFlexIOReadData (CCTK_ARGUMENTS);
  }
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cgh);
  
  int OutputGH (cGH* const cgh);
  int OutputVarAs (cGH* const cgh, const char* const varname,
		   const char* const alias);
  int TimeToOutput (cGH* const cgh, const int vindex);
  int TriggerOutput (cGH* const cgh, const int vindex);
  
  int InputGH (cGH* const cgh);
  int InputVarAs (cGH* const cgh, const char* const varname,
		  const char* const alias);
  
} // namespace CarpetIOFlexIO
