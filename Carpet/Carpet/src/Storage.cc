#include <cassert>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "dh.hh"
#include "gf.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  static int
  GroupStorageCrease (const cGH* cgh, int n_groups, const int* groups,
                      const int* timelevels, int* status,
                      const bool inc);
  
  int
  GroupStorageCrease (const cGH* cgh, int n_groups, const int* groups,
                      const int* timelevels, int* status,
                      const bool inc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cgh);
    assert (n_groups >= 0);
    assert (groups);
    assert (timelevels);
    for (int n=0; n<n_groups; ++n) {
      assert (groups[n] >= 0 and groups[n] < CCTK_NumGroups());
#if 0
      for (int nn=0; nn<n; ++nn) {
        assert (groups[nn] != groups[n]);
      }
#endif
      assert (timelevels[n] >= 0);
    }
    
    int total_num_timelevels = 0;
    
    for (int n=0; n<n_groups; ++n) {
      int const group = groups[n];
      
      cGroup gp;
      const int ierr = CCTK_GroupData (group, &gp);
      assert (! ierr);
      
      int const firstvarindex = CCTK_FirstVarIndexI (group);
      assert (gp.numvars == 0
              or (firstvarindex >= 0 and firstvarindex < CCTK_NumVars()));
      
      // Check an assumption
      if (! gp.vectorgroup) assert (gp.vectorlength == 1);
      
      // No storage change in local mode
      if (gp.grouptype == CCTK_GF) {
        assert ((map == -1 or maps == 1)
                and (component == -1
                     or vhh.at(0)->local_components(reflevel) == 1));
      }
      
      // Record previous number of allocated time levels
      if (status) {
        status[n] = groupdata.at(group).activetimelevels;
      }
      
      // Only do something if the number of time levels actually needs
      // to be changed -- do nothing otherwise
      if ((inc and timelevels[n] > groupdata.at(group).activetimelevels)
          or
          (! inc and timelevels[n] < groupdata.at(group).activetimelevels))
      {
        
        // Set the new number of active time levels
        groupdata.at(group).activetimelevels = timelevels[n];
        
        // There is a difference between the Cactus time levels and
        // the Carpet time levels.  If there are n time levels, then
        // the Cactus time levels are numbered 0 ... n-1, with the
        // current time level being 0.  In Carpet, the time levels are
        // numbered -(n-1) ... 0, where the current time level is also
        // 0.
        const int tmin = - timelevels[n] + 1;
        const int tmax = 0;
        
        // Allocate the time levels as well
        for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
          for (int var=0; var<gp.numvars; ++var) {
            const int vectorindex
              = (gp.contiguous
                 ? var
                 : gp.vectorgroup ? var % gp.vectorlength : 0);
            const int vectorlength
              = (gp.contiguous
                 ? gp.numvars
                 : gp.vectorgroup ? gp.vectorlength : 1);
            assert (vectorindex>=0 and vectorindex<gp.numvars);
            assert (vectorlength>0 and vectorlength<=gp.numvars);
            ggf* vectorleader
              = (vectorindex>0
                 ? arrdata.at(group).at(m).data.at(var - vectorindex)
                 : NULL);
            const int varindex = firstvarindex + var;
            switch (gp.vartype) {
#define TYPECASE(N,T)                                                   \
              case N:                                                   \
                arrdata.at(group).at(m).data.at(var) = new gf<T>        \
                (varindex, groupdata.at(group).transport_operator,      \
                 *arrdata.at(group).at(m).tt, *arrdata.at(group).at(m).dd, \
                 tmin, tmax, prolongation_order_time,                   \
                 vectorlength, vectorindex, (gf<T>*)vectorleader);      \
              break;
#include "typecase"
#undef TYPECASE
            default:
              UnsupportedVarType (varindex);
            } // switch gp.vartype
            
            // Set the data pointers for grid arrays
            if (gp.grouptype != CCTK_GF) {
              assert (m == 0);
              int const c = CCTK_MyProc(cgh);
              for (int tl=0; tl<gp.numtimelevels; ++tl) {
                cgh->data[varindex][tl]
                  = (tl < groupdata.at(group).activetimelevels
                     ? ((*arrdata.at(group).at(m).data.at(var))
                        (-tl, 0, c, 0)->storage())
                     : NULL);
              }
            }
            
          } // for var
        } // for m
        
      } // if really change the number of active time levels
      
#if 0
      // Complain if there are not enough active time levels
      if (gp.grouptype == CCTK_GF) {
        if (max_refinement_levels > 1) {
          if (groupdata.at(group).transport_operator != op_none) {
            if (groupdata.at(group).activetimelevels != 0
                and (groupdata.at(group).activetimelevels
                     <= prolongation_order_time))
            {
              char * const groupname = CCTK_GroupName (group);
              CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "There are not enough time levels for the desired temporal prolongation order in the grid function group \"%s\".  With Carpet::prolongation_order_time=%d, you need at least %d time levels.",
                          groupname,
                          prolongation_order_time, prolongation_order_time+1);
              free (groupname);
            }
          }
        }
      }
#endif
      
      // Record current number of time levels
      total_num_timelevels += groupdata.at(group).activetimelevels;
      
    } // for n
    
    return total_num_timelevels;
  }
  
  
  
  int
  GroupStorageIncrease (const cGH* cgh, int n_groups, const int* groups,
                        const int* timelevels, int* status)
  {
    Checkpoint ("GroupStorageIncrease");
    return
      GroupStorageCrease (cgh, n_groups, groups, timelevels, status, true);
  }
  
  
  
  int
  GroupStorageDecrease (const cGH* cgh, int n_groups, const int* groups,
                        const int* timelevels, int* status)
  {
    Checkpoint ("GroupStorageDecrease");
    return
      GroupStorageCrease (cgh, n_groups, groups, timelevels, status, false);
  }
  
  
  
  int
  EnableGroupStorage (const cGH* cgh, const char* groupname)
  {
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 and group<CCTK_NumGroups());
    const int timelevels = CCTK_MaxTimeLevelsGI(group);
    int status;
    GroupStorageIncrease (cgh, 1, &group, &timelevels, &status);
    // Return whether storage was allocated previously
    return status;
  }
  
  
  
  int
  DisableGroupStorage (const cGH* cgh, const char* groupname)
  {
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 and group<CCTK_NumGroups());
    const int timelevels = 0;
    int status;
    GroupStorageDecrease (cgh, 1, &group, &timelevels, &status);
    // Return whether storage was allocated previously
    return status;
  }
  
  
  
  int
  QueryGroupStorageB (const cGH* cgh, int group, const char* groupname)
  {
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 and group<CCTK_NumGroups());
    const int timelevels = 0;
    int status;
    GroupStorageIncrease (cgh, 1, &group, &timelevels, &status);
    // Return whether storage is allocated
    return status;
  }
  
  
  
  const int*
  ArrayGroupSizeB (const cGH* cgh, int dir, int group, const char* groupname)
  {
    static const int zero = 0;
    static const int error = 0;
    
    if (groupname) {
      group = CCTK_GroupIndex(groupname);
    }
    assert (group>=0 and group<CCTK_NumGroups());
    
    if (mglevel == -1) {
      return &error;            // meta mode
    }
    
    const int gptype = CCTK_GroupTypeI (group);
    if (gptype == CCTK_GF and map == -1) {
      return &error;            // global or level mode for a GF
    }
    
    const int gpdim = groupdata.at(group).info.dim;
    assert (dir>=0 and dir<gpdim);
    
    if (CCTK_QueryGroupStorageI(cgh, group)) {
      
      return &groupdata.at(group).info.lsh[dir];
      
    } else {
      
      // no storage
      return &zero;
      
    }
  }
  
  
  
  int GroupDynamicData (const cGH* cgh, int group, cGroupDynamicData* data)
  {
    assert (group>=0 and group<CCTK_NumGroups());
    *data = groupdata.at(group).info;
    return 0;
  }
  
} // namespace Carpet
