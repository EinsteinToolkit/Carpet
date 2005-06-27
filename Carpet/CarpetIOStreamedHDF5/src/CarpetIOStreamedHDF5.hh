#ifndef CARPETIOSTREAMEDHDF5_HH
#define CARPETIOSTREAMEDHDF5_HH

#include "CarpetIOHDF5.hh"
#include "SocketUtils.h"


// CarpetIOStreamed GH extension structure
typedef struct
{
  // port for clients to connect to
  unsigned int port;

  // socket to stream data over a TCP connection
  SOCKET socket;

  // default number of times to output
  int out_every_default;

  // the last iteration output for each variable
  vector<int> out_last;

  // list of variables to output
  char *out_vars;

  // stop on I/O parameter parsing errors ?
  int stop_on_parse_errors;

  // I/O request description list (for all variables)
  vector<ioRequest*> requests;

} CarpetIOStreamedHDF5GH;


namespace CarpetIOStreamedHDF5
{ 
  // scheduled routines (must be declared as C according to schedule.ccl)
  extern "C" {

    void CarpetIOStreamedHDF5_Startup (void);
    void CarpetIOStreamedHDF5_Init (const cGH* const);
    void CarpetIOStreamedHDF5_Terminate (const cGH* const);

  } // extern "C"

} // namespace CarpetIOStreamedHDF5

#endif // !defined(CARPETIOSTREAMEDHDF5_HH)
