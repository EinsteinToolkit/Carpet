// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOASCII/src/ioascii.hh,v 1.15 2004/02/18 15:12:29 schnetter Exp $

#ifndef CARPETIOASCII_HH
#define CARPETIOASCII_HH

#include <vector>

#include "cctk.h"

#include "ioascii.h"



namespace CarpetIOASCII {
  
  using namespace std;
  
  
  
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
    static vector<vector<vector<int> > > last_output; // [ml][rl][var]
    
    
    
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
    static CCTK_INT GetIntParameter
    (const char* parametertemplate, CCTK_INT fallback);
    static CCTK_REAL GetRealParameter
    (const char* parametertemplate, CCTK_REAL fallback);
    
  };                            // struct IOASCII
  
} // namespace CarpetIOASCII

#endif // !defined(CARPETIOASCII_HH)
