#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include "cctk.h"
#include "cctk_FortranString.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"

#include "carpet.hh"



namespace Carpet {

  using namespace std;



  // Get Carpet's GH extension
  CarpetGH const * GetCarpetGH (const cGH * const cgh)
  {
    assert (cgh);
    return &carpetGH;
  }

  // Get current refinement level
  extern "C"
  CCTK_INT Carpet_GetRefinementLevel (CCTK_POINTER_TO_CONST const cctkGH)
  {
    return reflevel;
  }

  // Get number of refinement levels
  extern "C"
  CCTK_INT Carpet_GetRefinementLevels (CCTK_POINTER_TO_CONST const cctkGH)
  {
    return reflevels;
  }

  // Get number of local components
  extern "C"
  CCTK_INT Carpet_GetLocalComponents (CCTK_POINTER_TO_CONST const cctkGH)
  {
    assert (reflevel >= 0);
    assert (map >= 0);
    return vhh.AT(map)->local_components(reflevel);
  }



  // Enable or disable prolongating
  CCTK_INT CarpetEnableProlongating (const CCTK_INT flag)
  {
    assert (flag==0 or flag==1);
    do_prolongate = flag;
    if (do_prolongate) {
      Checkpoint ("Prolongating enabled");
    } else {
      Checkpoint ("Prolongating disabled");
    }
    return 0;
  }

  CCTK_INT CarpetQueryProlongating ()
  {
    return do_prolongate;
  }
  
  
  
  // Get pointer to grid variable for a specific map and refinement level
  CCTK_POINTER
  Carpet_VarDataPtrI (CCTK_POINTER_TO_CONST const cctkGH,
                      CCTK_INT const m,
                      CCTK_INT const rl,
                      CCTK_INT const c,
                      CCTK_INT const tl,
                      CCTK_INT const varindex)
  {
    assert (cctkGH);
    assert (varindex >= 0 and varindex < CCTK_NumVars());
    int const groupindex = CCTK_GroupIndexFromVarI (varindex);
    assert (groupindex >= 0);
    int const grouptype = CCTK_GroupTypeI (groupindex);
    assert (mglevel >= 0);
    if (grouptype == CCTK_GF) {
      assert (m >= 0 and m < maps);
      assert (rl >= 0 and rl < reflevels);
      assert (c >= 0 and c < arrdata.AT(groupindex).AT(m).hh->components(reflevel));
      assert (arrdata.AT(groupindex).AT(m).hh->is_local(reflevel, c));
    } else {
      assert (m == 0);
      assert (rl == 0);
      assert (c == CCTK_MyProc (NULL));
    }
    int const maxtimelevels = CCTK_MaxTimeLevelsGI (groupindex);
    assert (tl >= 0 and tl < maxtimelevels);
    
    int const activetimelevels =
      groupdata.AT(groupindex).activetimelevels.AT(mglevel).AT(rl);
    if (tl < activetimelevels) {
      int const var = varindex - CCTK_FirstVarIndexI (varindex);
      ggf * const ff = arrdata.at(groupindex).at(m).data.at(var);
      gdata * const data = (*ff) (tl, rl, c, mglevel);
      return data->storage();
    } else {
      return NULL;
    }
  }
  
  
  
  // Multi-Model
  CCTK_POINTER_TO_CONST
  Carpet_GetMPICommUniverse (CCTK_POINTER_TO_CONST const cctkGH)
  {
    assert (comm_universe != MPI_COMM_NULL);
    return & comm_universe;
  }

  CCTK_POINTER_TO_CONST
  Carpet_GetMPICommWorld (CCTK_POINTER_TO_CONST const cctkGH)
  {
    assert (comm_world != MPI_COMM_NULL);
    return & comm_world;
  }
  


  CCTK_INT
  Carpet_GetCoordRange (CCTK_POINTER_TO_CONST const cctkGH_,
                        CCTK_INT              const m,
                        CCTK_INT              const ml,
                        CCTK_INT              const size,
                        CCTK_INT            * const gsh,
                        CCTK_REAL           * const lower,
                        CCTK_REAL           * const upper,
                        CCTK_REAL           * const delta)
  {
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    assert (cctkGH);
    assert (m >= 0 and m < maps);
    assert (ml >= 0 and ml < mglevels);
    assert (size >= 0);
    assert (not size or gsh);
    assert (not size or lower);
    assert (not size or upper);
    assert (not size or delta);
    
    assert (size == dim);
    
    ibbox const & baseext = vhh.at(m)->baseextents.at(ml).at(0);
    ivect const igsh = baseext.shape() / baseext.stride();
    
    for (int d = 0; d < dim; ++ d) {
      gsh[d] = igsh[d];
      lower[d] = origin_space.at(m).at(ml)[d];
      delta[d] = delta_space.at(m)[d] * mglevelfact;
      upper[d] = lower[d] + delta[d] * (gsh[d] - 1);
    }
    
    return 0;
  }



  // Communication

  int Barrier (const cGH* cgh)
  {
    void *dummy = &dummy;
    dummy = &cgh;

    MPI_Barrier (dist::comm());
    return 0;
  }



  int Exit (cGH* cgh, int retval)
  {
    CCTK_Barrier (cgh);
    dist::finalize();
    exit (retval);
    return -1;
  }

  int Abort (cGH* cgh, int retval)
  {
    void *dummy = &dummy;
    dummy = &cgh;
    
    assert (0);
    
    MPI_Abort (MPI_COMM_WORLD, retval);
    exit (retval);
    return -1;
  }



  int MyProc (const cGH* cgh)
  {
    // This may be called very early, before dist:comm() is valid
    int rank;
    MPI_Comm_rank (dist::goodcomm(), & rank);
    return rank;
  }

  int nProcs (const cGH* cgh)
  {
    // This may be called very early, before dist:comm() is valid
    int size;
    MPI_Comm_size (dist::goodcomm(), & size);
    return size;
  }



  MPI_Comm CarpetMPIComm ()
  {
    return dist::comm();
  }



  // Datatypes

  MPI_Datatype CarpetMPIDatatype (const int vartype)
  {
    switch (vartype) {
#define TYPECASE(N,T)				\
    case N: {					\
      T dummy;					\
      return dist::datatype(dummy);		\
    }
#include "typecase"
#undef TYPECASE
    default:
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Carpet does not support the variable type %d.", vartype);
    }
    // notreached
    return MPI_CHAR;
  }

  MPI_Datatype CarpetSimpleMPIDatatype (const int vartype)
  {
    switch (vartype) {
#ifdef CARPET_COMPLEX
    case CCTK_VARIABLE_COMPLEX:
      return CarpetMPIDatatype (CCTK_VARIABLE_REAL);
#endif
#ifdef CARPET_COMPLEX8
#  ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
      return CarpetMPIDatatype (CCTK_VARIABLE_REAL4);
#  endif
#endif
#ifdef CARPET_COMPLEX16
#  ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
      return CarpetMPIDatatype (CCTK_VARIABLE_REAL8);
#  endif
#endif
#ifdef CARPET_COMPLEX32
#  ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
      return CarpetMPIDatatype (CCTK_VARIABLE_REAL16);
#  endif
#endif
    default:
      return CarpetMPIDatatype (vartype);
    }
    // notreached
    return MPI_CHAR;
  }

  int CarpetSimpleMPIDatatypeLength (const int vartype)
  {
    switch (vartype) {
#ifdef CARPET_COMPLEX
    case CCTK_VARIABLE_COMPLEX:
#endif
#ifdef CARPET_COMPLEX8
#  ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#  endif
#endif
#ifdef CARPET_COMPLEX16
#  ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#  endif
#endif
#ifdef CARPET_COMPLEX32
#  ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#  endif
#endif
      return 2;
    default:
      return 1;
    }
    // notreached
    return 0;
  }



  // Timelevels

  int min_timelevel (const checktimes where, const int num_tl)
  {
    assert (num_tl>0);
    switch (where) {
    case currenttime:
      return 0;
    case currenttimebutnotifonly:
      // don't include current time if there is only one time level
      return num_tl>1 ? 0 : 1;
    case previoustime:
      return 1;
    case allbutlasttime:
      // do include current time if there is only one time level
      return 0;
    case allbutcurrenttime:
      return 1;
    case alltimes:
      return 0;
    default:
      assert (0);
    }
    return -999;
  }

  int max_timelevel (const checktimes where, const int num_tl)
  {
    assert (num_tl>0);
    switch (where) {
    case currenttime:
      return 0;
    case currenttimebutnotifonly:
      return 0;
    case previoustime:
      return num_tl>1 ? 1 : 0;
    case allbutlasttime:
      return num_tl-2;
    case allbutcurrenttime:
      return num_tl-1;
    case alltimes:
      return num_tl-1;
    default:
      assert (0);
    }
    return -999;
  }



  // Diagnostic output

  static void prepend_id (char* const msg, size_t const msglen)
  {
    if (mglevel!=-1) {
      snprintf (msg+strlen(msg), msglen-strlen(msg), "[%d]", mglevel);
      if (reflevel!=-1) {
        snprintf (msg+strlen(msg), msglen-strlen(msg), "[%d]", reflevel);
        if (map!=-1) {
          snprintf (msg+strlen(msg), msglen-strlen(msg), "[%d]", map);
          if (component!=-1) {
            snprintf (msg+strlen(msg), msglen-strlen(msg), "[%d]", component);
          }
        }
      }
      snprintf (msg+strlen(msg), msglen-strlen(msg), " ");
    }
  }

  void Output (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    va_list args;
    char msg[1000];
    snprintf (msg, sizeof msg, "%s", "");
    prepend_id (msg + strlen(msg), sizeof msg - strlen(msg));
    va_start (args, fmt);
    vsnprintf (msg + strlen(msg), sizeof msg - strlen(msg), fmt, args);
    va_end (args);
    CCTK_INFO (msg);
    if (barriers) {
      MPI_Barrier (dist::comm());
    }
  }

  void Waypoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (verbose or veryverbose) {
      va_list args;
      char msg[1000];
      snprintf (msg, sizeof msg, "%s", "");
      prepend_id (msg + strlen(msg), sizeof msg - strlen(msg));
      va_start (args, fmt);
      vsnprintf (msg + strlen(msg), sizeof msg - strlen(msg), fmt, args);
      va_end (args);
      CCTK_INFO (msg);
    }
    if (barriers) {
      MPI_Barrier (dist::comm());
    }
  }

  void Checkpoint (const char* fmt, ...)
  {
    DECLARE_CCTK_PARAMETERS;
    if (veryverbose) {
      va_list args;
      char msg[1000];
      snprintf (msg, sizeof msg, "%s", "");
      prepend_id (msg + strlen(msg), sizeof msg - strlen(msg));
      va_start (args, fmt);
      vsnprintf (msg + strlen(msg), sizeof msg - strlen(msg), fmt, args);
      va_end (args);
      CCTK_INFO (msg);
    }
    if (barriers) {
      MPI_Barrier (dist::comm());
    }
  }



  void UnsupportedVarType (const int vindex)
  {
    assert (vindex>=0 and vindex<CCTK_NumVars());
    CCTK_VWarn
      (0, __LINE__, __FILE__, CCTK_THORNSTRING,
       "Carpet does not support the type of the variable \"%s\".\n"
       "Either enable support for this type, "
       "or change the type of this variable.", CCTK_FullName(vindex));
  }

} // namespace Carpet
