// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5.hh,v 1.1 2004/03/03 09:44:26 schnetter Exp $

#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "iohdf5.h"

namespace CarpetIOHDF5 {
  
  // Variable definitions
  extern int GHExtension;
  extern int IOMethod;
  extern vector<bool> do_truncate; // [var]
  extern vector<vector<vector<int> > > last_output; // [ml][rl][var]
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cctkGH);
  
  int OutputGH (const cGH* const cctkGH);
  int OutputVarAs (const cGH* const cctkGH, const char* const varname,
		   const char* const alias);
  int TimeToOutput (const cGH* const cctkGH, const int vindex);
  int TriggerOutput (const cGH* const cctkGH, const int vindex);
  
  int InputGH (const cGH* const cctkGH);
  int InputVarAs (const cGH* const cctkGH, const char* const varname,
		  const char* const alias);
  
  int Recover (cGH* const cctkGH, const char *basefilename,
               const int called_from);
  
} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
