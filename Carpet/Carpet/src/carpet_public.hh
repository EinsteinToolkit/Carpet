// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/carpet_public.hh,v 1.36 2003/09/19 16:08:37 schnetter Exp $

// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

#ifndef CARPET_PUBLIC_HH
#define CARPET_PUBLIC_HH

#include <vector>

#include "cctk.h"
#include "cctk_Schedule.h"

#include "dh.hh"
#include "gf.hh"
#include "ggf.hh"
#include "gh.hh"
#include "th.hh"
  
#include "carpet_public.h"
  
  
  
// Stuff with C linkage
#include "carpet_public.h"

  
  
namespace Carpet {
  
  
  
  const int dim = 3;
  
  
  
  // Handle from CCTK_RegisterGHExtension
  extern int GHExtension;
  
  // Refinement factor
  extern int reffact;
  
  // Refinement factor on finest grid
  extern int maxreflevelfact;
  
  // Multigrid levels
  extern int mglevels;
  
  // Multigrid factor
  extern int mgfact;
  
  // Maximum number of refinement levels
  extern int maxreflevels;
  
  // Multigrid factor on coarsest grid
  extern int maxmglevelfact;
  
  
  
  // Current position on the grid hierarchy
  extern int reflevel;
  extern int mglevel;
  extern int component;
  
  // Current refinement factor
  extern int reflevelfact;
  
  // Current multigrid factor
  extern int mglevelfact;
  
  // Is this the time for a global mode call?
  extern bool do_global_mode;
  
  // Is prolongation enabled?
  extern bool do_prolongate;
  
  // Current times on the refinement levels
  extern vector<CCTK_REAL> refleveltimes;
  extern CCTK_REAL delta_time;
  
  
  
  // Data for grid functions
  
  // The grid hierarchy
  extern gh<dim>* hh;
  extern th<dim>* tt;
  extern dh<dim>* dd;
  
  // Data for everything
  struct arrdesc {
    // points to hh etc. for GF, and is unique for SCALAR and ARRAY
    cGroupDynamicData info;
    gh<dim>* hh;
    th<dim>* tt;
    dh<dim>* dd;
    vector<ggf<dim>*> data;     // [var]
    bool do_transfer;           // prolongate and restrict
    // VGF
  };
  extern vector<arrdesc> arrdata; // [group]
  
  
  
  // Checksums
  struct ckdesc {
    bool valid;
    unsigned int sum;
  };
  // [rl][c][group][var][tl]
  extern vector<vector<vector<vector<vector<ckdesc> > > > > checksums;
  
  
  
  // Registered functions
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh);
  
  int Initialise (tFleshConfig* config);
  int Evolve (tFleshConfig* config);
  int Shutdown (tFleshConfig* config);
  int CallFunction (void* function, cFunctionData* attribute, void* data);
  
  int SyncGroup (const cGH* cgh, const char* groupname);
  int EnableGroupStorage (const cGH* cgh, const char* groupname);
  int DisableGroupStorage (const cGH* cgh, const char* groupname); 
  int EnableGroupComm (const cGH* cgh, const char* groupname);
  int DisableGroupComm (const cGH* cgh, const char* groupname);
  int Barrier (const cGH* cgh);
  int Exit (cGH* cgh, int retval);
  int Abort (cGH* cgh, int retval);
  int MyProc (const cGH* cgh);
  int nProcs (const cGH* cgh);
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname);
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname);
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data);
  
  
  
  // Functions for recomposing the grid hierarchy
  void RegisterRegridRoutine (int (*routine)(const cGH * cckgGH,
					     gh<dim>::rexts& bbsss,
					     gh<dim>::rbnds& obss,
					     gh<dim>::rprocs& pss));
  
  void SplitRegions (const cGH* cgh, vector<bbox<int,dim> >& bbs,
		     vector<vect<vect<bool,2>,dim> >& obs);
  void SplitRegions_AlongZ (const cGH* cgh, vector<bbox<int,dim> >& bbs,
			    vector<vect<vect<bool,2>,dim> >& obs);
  
  void MakeProcessors (const cGH* cgh, const gh<dim>::rexts& bbsss,
		       gh<dim>::rprocs& pss);
  
  
  
  // Helper functions
  void set_reflevel (cGH* cgh, int rl);
  void set_mglevel (cGH* cgh, int ml);
  void set_component (cGH* cgh, int c);
  
  
  
  // Refinement level iterator
  
#define BEGIN_REFLEVEL_LOOP(cgh)                        \
  do {                                                  \
    int _rll;                                           \
    cGH * const _cgh = const_cast<cGH*>(cgh);           \
    assert (reflevel==-1);                              \
    for (int _rl=0; _rl<hh->reflevels(); ++_rl) {       \
      set_reflevel (_cgh, _rl);                         \
      {
#define END_REFLEVEL_LOOP                       \
      }                                         \
    }                                           \
    set_reflevel (_cgh, -1);                    \
    _rll = 0;                                   \
  } while (0)
  
  
  
  // Reverse refinement level iterator
  
#define BEGIN_REVERSE_REFLEVEL_LOOP(cgh)                \
  do {                                                  \
    int _rrll;                                          \
    cGH * const _cgh = const_cast<cGH*>(cgh);           \
    assert (reflevel==-1);                              \
    for (int _rl=hh->reflevels()-1; _rl>=0; --_rl) {    \
      set_reflevel (_cgh, _rl);                         \
      {
#define END_REVERSE_REFLEVEL_LOOP               \
      }                                         \
    }                                           \
    set_reflevel (_cgh, -1);                    \
    _rrll = 0;                                  \
  } while (0)
  
  
  
  // Multigrid level iterator
  
#define BEGIN_MGLEVEL_LOOP(cgh)                         \
  do {                                                  \
    int _mgl;                                           \
    cGH * const _cgh = const_cast<cGH*>(cgh);           \
    assert (reflevel>=0 && reflevel<hh->reflevels());   \
    assert (mglevel==-1);                               \
    for (int _ml=mglevels-1; _ml>=0; --_ml) {           \
      set_mglevel (_cgh, _ml);                          \
      {
#define END_MGLEVEL_LOOP                        \
      }                                         \
    }                                           \
    set_mglevel (_cgh, -1);                     \
    _mgl = 0;                                   \
  } while (0)
  
  
  
  // Component iterator
  
#define BEGIN_COMPONENT_LOOP(cgh, grouptype)            \
  do {                                                  \
    int _cl;                                            \
    cGH * const _cgh = const_cast<cGH*>(cgh);           \
    int const _grouptype = (grouptype);                 \
    int const _savec = (component);                     \
    int _minc, _maxc;                                   \
    if (_grouptype == CCTK_GF) {                        \
      assert (reflevel>=0 && reflevel<hh->reflevels()); \
      assert (mglevel>=0 && mglevel<mglevels);          \
      assert (component==-1);                           \
      _minc=0;                                          \
      _maxc=hh->components(reflevel);                   \
    } else {                                            \
      _minc=component;                                  \
      _maxc=component+1;                                \
    }                                                   \
    for (int _c=_minc; _c<_maxc; ++_c) {                \
      if (component!=_c) set_component (_cgh, _c);      \
      {
#define END_COMPONENT_LOOP                              \
      }                                                 \
    }                                                   \
    if (component!=_savec) set_component (_cgh, -1);    \
    _cl = 0;                                            \
  } while (0)


  
#define BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype)                      \
  do {                                                                  \
    int _lcl;                                                           \
    cGH * const _cgh = const_cast<cGH*>(cgh);                           \
    int const _grouptype = (grouptype);                                 \
    int const _savec = (component);                                     \
    int _minc, _maxc;                                                   \
    if (_grouptype == CCTK_GF) {                                        \
      assert (reflevel>=0 && reflevel<hh->reflevels());                 \
      assert (mglevel>=0 && mglevel<mglevels);                          \
      assert (hh->local_components(reflevel)==1 || component==-1);      \
      _minc=0;                                                          \
      _maxc=hh->components(reflevel);                                   \
    } else {                                                            \
      _minc=component;                                                  \
      _maxc=component+1;                                                \
    }                                                                   \
    for (int _c=_minc; _c<_maxc; ++_c) {                                \
      if (_grouptype!=CCTK_GF || hh->is_local(reflevel,_c)) {           \
        if (component!=_c) set_component (_cgh, _c);                    \
        {
#define END_LOCAL_COMPONENT_LOOP                                \
        }                                                       \
      }                                                         \
    }                                                           \
    if (component!=_savec) set_component (_cgh, _savec);        \
    _lcl = 0;                                                   \
  } while (0)
  
} // namespace Carpet

#endif // !defined(CARPET_PUBLIC_HH)
