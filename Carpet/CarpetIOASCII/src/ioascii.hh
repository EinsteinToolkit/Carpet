// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOASCII/src/ioascii.hh,v 1.2 2001/03/05 14:30:27 eschnett Exp $

#include <vector>

#include "cctk.h"



namespace Carpet {
  
  // scheduled functions
  extern "C" {
    int CarpetIOASCIIStartup();
  }
  
  
  
  // Everything is a template class, so that it can easily be
  // instantiated for all output dimensions.
  
  template<int outdim>
  struct IOASCII {
    
    // handle from CCTK_RegisterGHExtension
    static int GHExtension;
    
    // handles from CCTK_RegisterIOMethed
    static int IOMethod;
    
    
    
    // Do truncate the output files for a variable
    static vector<bool> do_truncate;
    
    // Last iteration on which a variable was output (-1 for none)
    static vector<int> last_output;
    
    
    
    // scheduled functions
    static int Startup();
    
    
    
    // registered functions
    
    static void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
    
    static int OutputGH (cGH* cgh);
    static int OutputVarAs (cGH* cgh, const char* varname, const char* alias);
    static int TimeToOutput (cGH* cgh, int vindex);
    static int TriggerOutput (cGH* cgh, int vindex);
    
    static int GetGridOffset (cGH* cgh, int dir,
			      const char* itempl, const char* iglobal,
			      const char* ctempl, const char* cglobal,
			      double cfallback);
    static int CoordToOffset (cGH* cgh, int dir, double coord);
    
    static const char* GetStringParameter
    (const char* parametertemplate, const char* fallback);
    static int GetIntParameter (const char* parametertemplate, int fallback);
    
  };
  
  bool CheckForVariable (cGH* const cgh,
			 const char* const varlist, const int vindex);
  void SetFlag (int index, const char* optstring, void* arg);
  
} // namespace Carpet
