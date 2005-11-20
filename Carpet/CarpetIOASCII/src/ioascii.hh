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

    // name of the output directory
    static const char* my_out_dir;

    // list of variables to output
    static char* my_out_vars;

    // Last iteration on which a refinement level of a variable was
    // output (INT_MIN for none)
    static vector<vector<vector<int> > > last_output; // [ml][rl][var]

    // I/O request description list (for all variables)
    static vector<ioRequest*> requests;


    // scheduled functions
    static int Startup();



    // registered functions

    static void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cctkGH);

    static int OutputGH (const cGH* cctkGH);
    static int OutputVarAs (const cGH* cctkGH,
                            const char* varname, const char* alias);
    static int TimeToOutput (const cGH* cctkGH, int vindex);
    static int TriggerOutput (const cGH* cctkGH, int vindex);

    static int GetGridOffset (const cGH* cctkGH, int dir,
                              const char* itempl, const char* iglobal,
                              const char* ctempl, const char* cglobal,
                              CCTK_REAL cfallback);
    static int CoordToOffset (const cGH* cctkGH, int dir, CCTK_REAL coord,
                              int ifallback);

    static void CheckSteerableParameters (const cGH* cctkGH);
    static const char* GetStringParameter (const char* parametertemplate);
    static CCTK_INT GetIntParameter (const char* parametertemplate);
    static CCTK_REAL GetRealParameter (const char* parametertemplate);

  };                            // struct IOASCII

} // namespace CarpetIOASCII

#endif // !defined(CARPETIOASCII_HH)
