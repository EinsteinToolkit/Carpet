#include <cassert>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "dh.hh"
#include "gf.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  int EnableGroupStorage (const cGH* cgh, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("EnableGroupStorage \"%s\"", groupname);
    
    // TODO: Enabling storage for one refinement level has to enable
    // it for all other refinement levels as well.  Disabling must
    // wait until all refinement levels have been disabled.
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (CCTK_NumVarsInGroupI(group)==0) return 0;
    
    const int grouptype = CCTK_GroupTypeI(group);
    
    // No storage change in local mode
    if (grouptype == CCTK_GF) {
      assert ((map == -1 || maps == 1)
              && (component == -1
                  || vhh.at(0)->local_components(reflevel) == 1));
    }
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was enabled previously
      const int n0 = CCTK_FirstVarIndexI(group);
      assert (n0>=0);
      const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
      assert (num_tl>0);
      return num_tl;
    }
    
    // There is a difference between the Cactus time levels and the
    // Carpet time levels.  If there are n time levels, then the
    // Cactus time levels are numbered 0 ... n-1, with the current
    // time level being 0.  In Carpet, the time levels are numbered
    // -(n-1) ... 0, where the current time level is also 0.
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tmin = 1 - num_tl;
    const int tmax = 0;
    
    if (grouptype == CCTK_GF) {
      if (max_refinement_levels > 1) {
        if (groupdata.at(group).transport_operator != op_none) {
          if (num_tl <= prolongation_order_time) {
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "There are not enough time levels for the desired temporal prolongation order in the grid function group \"%s\".  With Carpet::prolongation_order_time=%d, you need at least %d time levels.",
                        CCTK_GroupName(group),
                        prolongation_order_time, prolongation_order_time+1);
          }
        }
      }
    }
    
    cGroup gp;
    const int ierr = CCTK_GroupData (group, &gp);
    assert (!ierr);
    const int vectorlength = gp.vectorgroup ? gp.vectorlength : 1;
    
    assert (arrdata.at(group).at(0).data.size()==0
	    || arrdata.at(group).at(0).data.at(0) == 0);
    for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
      for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
        const int vectorindex = gp.vectorgroup ? var % vectorlength : 0;
        ggf<dim>* vectorleader
          = (gp.vectorgroup && vectorindex>0
             ? arrdata.at(group).at(m).data.at(var - vectorindex)
             : NULL);
        const int n = n0 + var;
        switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)                                                   \
        case N:                                                         \
          arrdata.at(group).at(m).data.at(var) = new gf<T,dim>          \
            (n, groupdata.at(group).transport_operator,                 \
             *arrdata.at(group).at(m).tt, *arrdata.at(group).at(m).dd,  \
             tmin, tmax, prolongation_order_time,                       \
             vectorlength, vectorindex, (gf<T,dim>*)vectorleader);      \
          break;
#include "typecase"
#undef TYPECASE
        default:
          UnsupportedVarType(n);
        } // switch
        
        if (grouptype != CCTK_GF) {
          for (int tl=0; tl<num_tl; ++tl) {
            assert (m == 0);
            int const c = CCTK_MyProc(cgh);
            cgh->data[n][tl] = ((*arrdata.at(group).at(m).data.at(var))
                                (-tl, 0, c, 0)->storage());
          }
        }
        
      } // for
    } // for
    
//     PoisonGroup (cgh, group, alltimes);
    
    // storage was not enabled previously
    return 0;
  }
  
  
  
  int DisableGroupStorage (const cGH* cgh, const char* groupname)
  {
    Checkpoint ("DisableGroupStorage \"%s\"", groupname);
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    // No storage change in local mode
    assert (! (component!=-1 && CCTK_GroupTypeI(group)==CCTK_GF));
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      // storage was disabled previously
      return 0;
    }
    
    // TODO
    CCTK_WARN (2, "Cannot disable storage -- storage management is not yet consistent for FMR");
    return 1;
    
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    
    for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
      for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
        const int n = n0 + var;
        switch (CCTK_VarTypeI(n)) {
#define TYPECASE(N,T)                                                   \
        case N:                                                         \
          assert (arrdata.at(group).at(m).data.at(var));                \
          delete (gf<T,dim>*)arrdata.at(group).at(m).data.at(var);      \
          arrdata.at(group).at(m).data.at(var) = NULL;                  \
          break;
#include "typecase"
#undef TYPECASE
        default:
          UnsupportedVarType(n);
        } // switch
        arrdata.at(group).at(m).data.at(var) = NULL;
        
        if (CCTK_GroupTypeI(group) != CCTK_GF) {
          for (int tl=0; tl<num_tl; ++tl) {
            cgh->data[n][tl] = NULL;
          }
        }
        
      } // for
    } // for
    
    // storage was enabled previously
    return 1;
  }
  
  
  
  int QueryGroupStorageB (const cGH* cgh, int group, const char* groupname)
  {
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (CCTK_NumVarsInGroupI(group)==0) return 0;
    
    const int n = CCTK_FirstVarIndexI(group);
    assert (n>=0 && n<CCTK_NumVars());
    const int var = 0;
    
    return arrdata.at(group).at(0).data.at(var) != NULL;
  }
  
  
  
  const int* ArrayGroupSizeB (const cGH* cgh, int dir, int group,
			      const char* groupname)
  {
    static const int zero = 0;
    static const int error = 0;
    
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (mglevel == -1) {
      return &error;            // meta mode
    }
    
    const int gptype = CCTK_GroupTypeI (group);
    if (gptype == CCTK_GF && map == -1) {
      return &error;            // global or level mode for a GF
    }
    
    const int gpdim = groupdata.at(group).info.dim;
    assert (dir>=0 && dir<gpdim);
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      
      return &groupdata.at(group).info.lsh[dir];
      
    } else {
      
      // no storage
      return &zero;
      
    }
  }
  
  
  
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data)
  {
    assert (group>=0 && group<CCTK_NumGroups());
    *data = groupdata.at(group).info;
    return 0;
  }
  
} // namespace Carpet
