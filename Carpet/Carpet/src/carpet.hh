// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet.hh,v 1.9 2001/06/12 14:56:56 schnetter Exp $

#include <vector>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "Carpet/CarpetLib/src/dh.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/th.hh"

namespace Carpet {
  
  
  
  const int dim = 3;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  extern int GHExtension;
  
  // Maximum number of refinement levels
  extern int maxreflevels;
  
  // Refinement factor on finest grid
  extern int maxreflevelfact;
  
  // Current iteration per refinement level
  extern vector<int> iteration;
  
  // Current position on the grid hierarchy
  extern int reflevel;
  extern int mglevel;
  extern int component;
  
  // Current refinement factor
  extern int reflevelfact;
  
  // Time step on base grid
  extern CCTK_REAL base_delta_time;
  
  
  
  // Data for scalars
  struct scdesc {
    vector<vector<vector<void*> > > data; // [var][rl][tl]
  };
  extern vector<scdesc> scdata;	// group
  
  // Data for arrays
  struct arrdesc {
    gh<dim>* hh;
    th* tt;
    dh<dim>* dd;
    vector<generic_gf<dim>* > data; // [var]
    int size[dim];
  };
  extern vector<arrdesc> arrdata; // [group]
  
  // Data for grid functions
  
  // The grid hierarchy
  extern gh<dim>* hh;
  extern th* tt;
  extern dh<dim>* dd;
  extern int gfsize[dim];
  
  struct gfdesc {
    vector<generic_gf<dim>* > data; // [var]
  };
  extern vector<gfdesc> gfdata;	// [group]
  
  // Checksums
  struct ckdesc {
    bool valid;
    int sum;
  };
  extern vector<vector<vector<vector<ckdesc> > > > checksums; // [n][rl][tl][c]
  
  
  
  // Scheduled functions
  extern "C" {
    int CarpetStartup();
  }
  
  // Registered functions
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
  
  int Initialise (tFleshConfig* config);
  int Evolve (tFleshConfig* config);
  int Shutdown (tFleshConfig* config);
  int CallFunction (void* function, cFunctionData* attribute, void* data);
  
  int SyncGroup (cGH* cgh, const char* groupname);
  int EnableGroupStorage (cGH* cgh, const char* groupname);
  int DisableGroupStorage (cGH* cgh, const char* groupname); 
  int EnableGroupComm (cGH* cgh, const char* groupname);
  int DisableGroupComm (cGH* cgh, const char* groupname);
  int Barrier (cGH* cgh);
  int Exit (cGH* cgh, int retval);
  int Abort (cGH* cgh, int retval);
  int MyProc (cGH* cgh);
  int nProcs (cGH* cgh);
  const int* ArrayGroupSizeB (cGH* cgh, int dir, int group,
			      const char* groupname);
  int QueryGroupStorageB (cGH* cgh, int group, const char* groupname);
  
  
  
  // Helper functions
  void set_reflevel (cGH* cgh, int rl);
  void set_mglevel (cGH* cgh, int ml);
  void set_component (cGH* cgh, int c);
  
  
  
  // Refinement level iterator
  
#define BEGIN_REFLEVEL_LOOP(cgh)		\
  do {						\
    int _rl;					\
    assert (reflevel==0);			\
    for (;;) {					\
      {
#define END_REFLEVEL_LOOP(cgh)			\
      }						\
      if (reflevel==maxreflevels-1) break;	\
      set_reflevel ((cgh), reflevel+1);		\
    }						\
    set_reflevel ((cgh), 0);			\
    assert (reflevel==0);			\
    _rl = 0;					\
  } while (0)
  
  
  
  // Reverse refinement level iterator
  
#define BEGIN_REVERSE_REFLEVEL_LOOP(cgh)	\
  do {						\
    int _rrl;					\
    assert (reflevel==0);			\
    set_reflevel ((cgh), maxreflevels-1);	\
    for (;;) {					\
      {
#define END_REVERSE_REFLEVEL_LOOP(cgh)		\
      }						\
      if (reflevel==0) break;			\
      set_reflevel ((cgh), reflevel-1);		\
    }						\
    assert (reflevel==0);			\
    _rrl = 0;					\
  } while (0)
  
  
  
  // Component iterator
  
#define BEGIN_COMPONENT_LOOP(cgh)		\
  do {						\
    int _cl;					\
    assert (component==-1);			\
    set_component ((cgh), 0);			\
    for (;;) {					\
      {
#define END_COMPONENT_LOOP(cgh)				\
      }							\
      if (component==hh->components(reflevel)-1) break;	\
      set_component ((cgh), component+1);		\
    }							\
    set_component ((cgh), -1);				\
    assert (component==-1);				\
    _cl = 0;						\
  } while (0)
  
  
  
  extern "C" {
    MPI_Comm CarpetMPICommunicator();
  }
  
} // namespace Carpet
