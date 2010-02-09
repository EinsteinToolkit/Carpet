#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Network.h"

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "CactusBase/IOUtil/src/ioGH.h"

#include "CarpetIOStreamedHDF5.hh"


namespace CarpetIOStreamedHDF5
{

using namespace std;
using namespace Carpet;
using namespace CarpetIOHDF5;


// Variable definitions
static vector<vector<vector<int> > > last_output; // [ml][rl][var]

// registered GH extension setup routine
static void* SetupGH (tFleshConfig* const fleshconfig,
                      const int convLevel, cGH* const cctkGH);

// callbacks for CarpetIOStreamedHDF5's I/O method
static int OutputGH (const cGH* const cctkGH);
static int TimeToOutput (const cGH* const cctkGH, const int vindex);
static int OutputVar (const cGH* const cctkGH, const ioRequest* const request,
                      hid_t file);

static void CheckSteerableParameters (const cGH *const cctkGH,
                                      CarpetIOStreamedHDF5GH *myGH);

//////////////////////////////////////////////////////////////////////////////
// public routines
//////////////////////////////////////////////////////////////////////////////

void CarpetIOStreamedHDF5_Startup (void)
{
  CCTK_RegisterBanner ("AMR streamed HDF5 output "
                       "provided by CarpetIOStreamedHDF5");

  const int GHExtension = CCTK_RegisterGHExtension (CCTK_THORNSTRING);
  CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
}


void CarpetIOStreamedHDF5_Init (const cGH* const cctkGH)
{
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = -1;
  *next_output_iteration = 0;
  *next_output_time = cctk_time;
}


// close the socket used for HDF5 streaming
void CarpetIOStreamedHDF5_Terminate (const cGH* const cctkGH)
{
  CarpetIOStreamedHDF5GH* myGH =
    (CarpetIOStreamedHDF5GH*) CCTK_GHExtension (cctkGH, CCTK_THORNSTRING);

  Socket_CloseSocket (myGH->socket);
}


//////////////////////////////////////////////////////////////////////////////
// private routines
//////////////////////////////////////////////////////////////////////////////
static void* SetupGH (tFleshConfig* const fleshconfig,
                      const int convLevel, cGH* const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  // processor 0 opens a socket on the given port to listen for clients
  unsigned int real_port;
  SOCKET real_socket = INVALID_SOCKET;
  if (dist::rank() == 0) {

    real_socket = Socket_TCPOpenServerSocket (port, &real_port, 1);

    if (real_socket == INVALID_SOCKET) {

      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Couldn't open TCP server socket on output port %d. "
                  "No HDF5 streaming output will be available !", port);

    } else if (Socket_SetNonBlocking (real_socket) < 0) {

      CCTK_WARN (1, "Couldn't set output socket into non-blocking mode. "
                    "No HDF5 streaming output will be available !");
      Socket_CloseSocket (real_socket);
      real_socket = INVALID_SOCKET;

    }
  }

  // broadcast success of socket operation to all processors
  int have_socket = real_socket != INVALID_SOCKET;
  MPI_Bcast (&have_socket, 1, MPI_INT, 0, dist::comm());

  // do not continue here if no socket is available
  if (not have_socket) {
    return (NULL);
  }

  if (dist::rank() == 0) {
    char hostname[256];
    Util_GetHostName (hostname, sizeof (hostname));
    CCTK_VInfo (CCTK_THORNSTRING,
                "data streaming service started on '%s:%u'", hostname, port);
  }

  // register CarpetIOStreamedHDF5's routines as a new I/O method
  const int IOMethod = CCTK_RegisterIOMethod ("IOStreamedHDF5");
  CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
#if 0
  // for now only register the OutputGH() callback
  CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
#endif

  if (not CCTK_Equals (verbose, "none")) {
    CCTK_INFO ("I/O Method 'IOStreamedHDF5' registered: AMR streamed HDF5 "
               "output of grid variables");
  }

  const int numvars = CCTK_NumVars ();

  // allocate a new GH extension structure
  CarpetIOStreamedHDF5GH* myGH = new CarpetIOStreamedHDF5GH;

  myGH->out_last.resize(numvars);
  for (int i = 0; i < numvars; i++) {
    myGH->out_last[i] = -1;
  }
  myGH->requests.resize(numvars);
  myGH->out_vars = strdup ("");
  myGH->out_every_default = out_every - 1;

  // initial I/O parameter check
  myGH->stop_on_parse_errors = strict_io_parameter_check;
  CheckSteerableParameters (cctkGH, myGH);
  myGH->stop_on_parse_errors = 0;

  // no iterations have yet been output
  last_output.resize(mglevels);
  for (int i = 0; i < mglevels; i++) {
    last_output[i].resize (maxreflevels);
    for (int j = 0; j < maxreflevels; j++) {
      last_output[i][j].resize (numvars, INT_MIN);
    }
  }

  myGH->socket = real_socket;
  myGH->port = real_port;

  return (myGH);
}


static void CheckSteerableParameters (const cGH *const cctkGH,
                                      CarpetIOStreamedHDF5GH *myGH)
{
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOHDF5::out_vars' parameter if it has changed
  if (strcmp (out_vars, myGH->out_vars)) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                               "IOStreamedHDF5::out_vars",
                               myGH->stop_on_parse_errors, out_vars,
                               -1, -1.0, &myGH->requests[0]);
#else
    IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                               "IOStreamedHDF5::out_vars",
                               myGH->stop_on_parse_errors, out_vars,
                               -1, &myGH->requests[0]);
#endif

    // notify the user about the new setting
    if (not CCTK_Equals (verbose, "none")) {
      int count = 0;
      string msg ("Periodic streamed HDF5 output requested for '");
      for (int i = CCTK_NumVars () - 1; i >= 0; i--) {
        if (myGH->requests[i]) {
          if (count++) {
            msg += "', '";
          }
          char *fullname = CCTK_FullName (i);
          msg += fullname;
          free (fullname);
        }
      }
      if (count) {
        msg += "'";
        CCTK_INFO (msg.c_str());
      }
    }

    // save the last setting of 'IOHDF5::out_vars' parameter
    free (myGH->out_vars);
    myGH->out_vars = strdup (out_vars);
  }
}


static int OutputGH (const cGH* const cctkGH)
{
  int error_count = 0;
  DECLARE_CCTK_PARAMETERS;

  CarpetIOStreamedHDF5GH *myGH =
    (CarpetIOStreamedHDF5GH *) CCTK_GHExtension (cctkGH, CCTK_THORNSTRING);

  static hid_t file = -1;
  static int client_is_ready = 0;

  // loop over all variables
  for (int vindex = CCTK_NumVars () - 1; vindex >= 0; vindex--) {
    if (TimeToOutput (cctkGH, vindex)) {

      // check if any client is ready to receive streamed data
      //
      // Even though this requires a broadcast operation, it should be
      // by far less expensive than gathering all the data on processor 0
      // and find out afterwards that no client wants to have it.
      if (not client_is_ready) {
        if (file < 0 and dist::rank() == 0) {
          fd_set read_set;
          FD_ZERO (&read_set);
          FD_SET (myGH->socket, &read_set);

          struct timeval timeout = {0, 0};
          client_is_ready =
            select (FD_SETSIZE, &read_set, NULL, NULL, &timeout) == 1;
        }
        MPI_Bcast (&client_is_ready, 1, MPI_INT, 0, dist::comm());
      }

      // short cut if there is nothing to do
      if (not client_is_ready) {
        continue;
      }

      // when called the first time during the current iteration,
      // open the file on processor 0
      if (file < 0 and dist::rank() == 0) {

        if (CCTK_Equals (verbose, "full")) {
          CCTK_VInfo (CCTK_THORNSTRING, "Opening HDF5 output file on output "
                                        "port %u", myGH->port);
        }

        // set file access property list to use the Stream VFD and open the file
        H5FD_stream_fapl_t fapl;
        fapl.increment = 0;
        fapl.socket = myGH->socket;
        fapl.do_socket_io = 1;
        fapl.backlog = max_num_clients;
        fapl.broadcast_fn = NULL;
        fapl.broadcast_arg = NULL;

        hid_t plist;
        HDF5_ERROR (plist = H5Pcreate (H5P_FILE_ACCESS));
        HDF5_ERROR (H5Pset_fapl_stream (plist, &fapl));

        // a filename is not used by Stream VFD
        assert (file < 0);
        HDF5_ERROR (file = H5Fcreate ("unused", H5F_ACC_TRUNC, H5P_DEFAULT,
                                      plist));
        HDF5_ERROR (H5Pclose (plist));
      }
      assert (dist::rank() or file >= 0);

      OutputVar (cctkGH, myGH->requests[vindex], file);
    }
  }

  // close an open file if the finest level has been output
  if (reflevel == reflevels - 1) {
    if (dist::rank() == 0 and file >= 0) {
      if (CCTK_Equals (verbose, "full")) {
        CCTK_VInfo (CCTK_THORNSTRING, "Closing HDF5 output file on port %u",
                    myGH->port);
      }
      HDF5_ERROR (H5Fclose (file));
      file = -1;
    }
    client_is_ready = 0;
  }

  // return negative number of errors occured during this output
  return (-error_count);
}


static int TimeToOutput (const cGH* const cctkGH, const int vindex)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int numvars = CCTK_NumVars();
  assert (vindex>=0 and vindex<numvars);

  if (CCTK_GroupTypeFromVarI (vindex) != CCTK_GF and not do_global_mode) {
    return 0;
  }

  CarpetIOStreamedHDF5GH *myGH =
    (CarpetIOStreamedHDF5GH *) CCTK_GHExtension (cctkGH, CCTK_THORNSTRING);
  CheckSteerableParameters (cctkGH, myGH);

  // check if output for this variable was requested
  if (not myGH->requests[vindex]) {
    return (0);
  }

  // check whether this refinement level should be output
  if (not (myGH->requests[vindex]->refinement_levels & (1 << reflevel))) {
    return (0);
  }

  // check if output for this variable was requested individually
  // by a "<varname>{ out_every = <number> }" option string
  // this will overwrite the output criterion setting
  const char *myoutcriterion = CCTK_EQUALS (out_criterion, "default") ?
                               io_out_criterion : out_criterion;
  if (myGH->requests[vindex]->out_every >= 0) {
    myoutcriterion = "divisor";
  }

  if (CCTK_EQUALS (myoutcriterion, "never")) {
    return (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration = false;

  if (CCTK_EQUALS (myoutcriterion, "iteration")) {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myoutevery > 0) {
      if (*this_iteration == cctk_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >= *next_output_iteration) {
        // it is time for the next output
        output_this_iteration = true;
        *this_iteration = cctk_iteration;
        *next_output_iteration = cctk_iteration + myoutevery;
      }
    }
  } else if (CCTK_EQUALS (myoutcriterion, "divisor")) {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myGH->requests[vindex]->out_every >= 0) {
      myoutevery = myGH->requests[vindex]->out_every;
    }
    if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
  } else if (CCTK_EQUALS (myoutcriterion, "time")) {
    CCTK_REAL myoutdt = out_dt == -2 ? io_out_dt : out_dt;
    if (myoutdt == 0 or *this_iteration == cctk_iteration) {
      output_this_iteration = true;
    } else if (myoutdt > 0 and (cctk_time / cctk_delta_time
                             >= *next_output_time / cctk_delta_time - 1.0e-12)) {
      // it is time for the next output
      output_this_iteration = true;
      *this_iteration = cctk_iteration;
      *next_output_time = cctk_time + myoutdt;
    }
  }

  if (not output_this_iteration) {
    return 0;
  }

  if (last_output.at(mglevel).at(reflevel).at(vindex) == cctk_iteration) {
    // Has already been output during this iteration
    char* varname = CCTK_FullName(vindex);
    CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Skipping output for variable \"%s\", because this variable "
                "has already been output during the current iteration -- "
                "probably via a trigger during the analysis stage",
                varname);
    free (varname);
    return 0;
  }

  assert (last_output.at(mglevel).at(reflevel).at(vindex) < cctk_iteration);

  // Should be output during this iteration
  return 1;
}


static int OutputVar (const cGH* const cctkGH, const ioRequest* const request,
                      hid_t file)
{
  DECLARE_CCTK_PARAMETERS;

  const int group = CCTK_GroupIndexFromVarI (request->vindex);
  assert (group >= 0);
  cGroup groupdata;
  CCTK_GroupData (group, &groupdata);
  if (groupdata.grouptype == CCTK_SCALAR or groupdata.grouptype == CCTK_ARRAY) {
    assert (do_global_mode);
  }

  // check for storage
  if (not CCTK_QueryGroupStorageI (cctkGH, group)) {
    char* fullname = CCTK_FullName (request->vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot output variable '%s' because it has no storage",
                fullname);
    free (fullname);
    return (0);
  }

  if (out_unchunked) {
    WriteVarUnchunked (cctkGH, file, request, false);
  } else {
    WriteVarChunkedSequential (cctkGH, file, request, false);
  }

  last_output.at(mglevel).at(reflevel).at(request->vindex) =
    cctkGH->cctk_iteration;

  return (0);
}


} // namespace CarpetIOStreamedHDF5
