// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5.hh,v 1.4 2004/03/09 16:02:48 cott Exp $

#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "iohdf5.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

namespace CarpetIOHDF5 {
  
  // Variable definitions
  extern int GHExtension;
  extern int IOMethod;
  extern vector<bool> do_truncate; // [var]
  extern vector<vector<vector<int> > > last_output; // [ml][rl][var]
  
  void* SetupGH (tFleshConfig* const fc,
		 const int convLevel, cGH* const cctkGH);
  
  int OutputGH (const cGH* const cctkGH);
  int WriteVar (const cGH* const cctkGH, hid_t writer, const ioRequest* request,
		   const int called_from_checkpoint);
  int OutputVarAs (const cGH* const cctkGH, const char* const varname,
		   const char* const alias);

  int TimeToOutput (const cGH* const cctkGH, const int vindex);
  int TriggerOutput (const cGH* const cctkGH, const int vindex);
  
  int InputGH (const cGH* const cctkGH);
  int InputVarAs (const cGH* const cctkGH, const char* const varname,
		  const char* const alias);
  
  int Recover (cGH* const cctkGH, const char *basefilename,
               const int called_from);
  

  // auxiliary functions defined in iohdf5utils.cc

  const char* GetStringParameter (const char* const parametername,
                                  const char* const fallback);
  int GetIntParameter (const char* const parametername, int fallback);
  
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

} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
