// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.3 2001/03/10 20:55:03 eschnett Exp $

#include <vector>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"

namespace Carpet {
  
  
  
  const int dim = 3;
  
  
  
  // handle from CCTK_RegisterGHExtension
  extern int GHExtension;
  
  // data for scalars
  extern vector<vector<vector<void*> > > scdata;// [group][var][tl]
  
  // data for arrays
  struct arrdesc {
    gh<dim>* hh;
    th<dim>* tt;
    dh<dim>* dd;
    vector<generic_gf<dim>* > data; // [var]
    int size[dim];
  };
  extern vector<arrdesc> arrdata; // [group]
  
  // data for grid functions
  
  // the grid hierarchy
  extern gh<dim>* hh;
  extern th<dim>* tt;
  extern dh<dim>* dd;
  extern int gfsize[dim];
  
  struct gfdesc {
    vector<generic_gf<dim>* > data; // [var]
  };
  extern vector<gfdesc> gfdata;	// [group]
  
  // active time level
  extern int activetimelevel;	// 0 for current, 1 for next
  
  // current position on the grid hierarchy
  extern int mglevel;
  extern int reflevel;
  extern int component;
  
  
  
  // scheduled functions
  extern "C" {
    int CarpetStartup();
  }
  
  // registered functions
  void* SetupGH (tFleshConfig *fc, int convLevel, cGH *cgh);
  
  int Initialise (tFleshConfig *config);
  int Evolve (tFleshConfig *config);
  int Shutdown (tFleshConfig *config);
  int CallFunction (void *function, cFunctionData *attribute, void *data);
  
  void reflevel_up (cGH* cgh);
  void reflevel_down (cGH* cgh);
  
  int SyncGroup (cGH *cgh, const char *groupname);
  int EnableGroupStorage (cGH *cgh, const char *groupname);
  int DisableGroupStorage (cGH *cgh, const char *groupname); 
  int EnableGroupComm (cGH *cgh, const char *groupname);
  int DisableGroupComm (cGH *cgh, const char *groupname);
  int Barrier (cGH *cgh);
  int Exit (cGH *cgh, int retval);
  int Abort (cGH *cgh, int retval);
  int myProc (cGH *cgh);
  int nProcs (cGH *cgh);
  const int* ArrayGroupSizeB (cGH *cgh, int dir, int group,
			      const char *groupname);
  int QueryGroupStorageB (cGH *cgh, int group, const char *groupname);
  
  
  
  // Helper functions
  extern "C" {
    MPI_Comm CarpetMPICommunicator();
  }
  
} // namespace Carpet
