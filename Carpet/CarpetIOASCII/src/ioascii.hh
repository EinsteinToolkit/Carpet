// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOASCII/src/ioascii.hh,v 1.11 2002/08/30 18:08:53 schnetter Exp $

#include <vector>

#include "cctk.h"



namespace CarpetIOASCII {
  
  using namespace std;
  
  
  
  // scheduled functions
  extern "C" {
    int CarpetIOASCIIStartup();
  }
  
  
  
  // Everything is a class template, so that it can easily be
  // instantiated for all output dimensions.
  
  template<int outdim>
  struct IOASCII {
    
    // handle from CCTK_RegisterGHExtension
    static int GHExtension;
    
    // handles from CCTK_RegisterIOMethed
    static int IOMethod;
    
    
    
    // Do truncate the output files for a variable
    static vector<bool> do_truncate;
    
    // Last iteration on which a refinement level of a variable was
    // output (INT_MIN for none)
    static vector<vector<int> > last_output; // [rl][var]
    
    
    
    // scheduled functions
    static int Startup();
    
    
    
    // registered functions
    
    static void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
    
    static int OutputGH (const cGH* cgh);
    static int OutputVarAs (const cGH* cgh,
			    const char* varname, const char* alias);
    static int TimeToOutput (const cGH* cgh, int vindex);
    static int TriggerOutput (const cGH* cgh, int vindex);
    
    static int GetGridOffset (const cGH* cgh, int dir,
			      const char* itempl, const char* iglobal,
			      const char* ctempl, const char* cglobal,
			      CCTK_REAL cfallback);
    static int CoordToOffset (const cGH* cgh, int dir, CCTK_REAL coord,
			      int ifallback);
    
    static const char* GetStringParameter
    (const char* parametertemplate, const char* fallback);
    static int GetIntParameter (const char* parametertemplate, int fallback);
    
  };				// struct IOASCII
  
} // namespace CarpetIOASCII
