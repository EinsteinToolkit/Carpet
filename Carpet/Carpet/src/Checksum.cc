#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Checksum.cc,v 1.12 2003/08/10 21:59:51 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Checksum_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void CalculateChecksums (const cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("%*sCalculateChecksums", 2*reflevel, "");
    
    checksums.resize(maxreflevels);
    checksums[reflevel].resize(hh->components(reflevel));
    BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
      checksums[reflevel][component].resize(CCTK_NumGroups());
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_QueryGroupStorageI(cgh, group)) {
          
          const int nvar = CCTK_NumVarsInGroupI(group);
          
          checksums[reflevel][component][group].resize(nvar);
          
          if (reflevel<arrdata[group].hh->reflevels()
              && component<arrdata[group].hh->components(reflevel)
              && arrdata[group].hh->is_local(reflevel, component)) {
            
            const int n0 = CCTK_FirstVarIndexI(group);
            const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
            assert (sz>0);
            
            vect<int,dim> size(1);
            const int gpdim = arrdata[group].info.dim;
            for (int d=0; d<gpdim; ++d) {
              size[d] = arrdata[group].info.lsh[d];
            }
            const int np = prod(size);
            
            const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
            assert (num_tl>0);
            const int min_tl = mintl(where, num_tl);
            const int max_tl = maxtl(where, num_tl);
            
            for (int var=0; var<nvar; ++var) {
              checksums[reflevel][component][group][var].resize(num_tl);
              for (int tl=min_tl; tl<=max_tl; ++tl) {
                
                const int n = n0 + var;
                const void* data = cgh->data[n][tl];
                unsigned int chk = 0;
                for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
                  chk += ((const unsigned int*)data)[i];
                }
                
                checksums[reflevel][component][group][var][tl].sum = chk;
                checksums[reflevel][component][group][var][tl].valid = true;
                
              } // for tl
            } // for var
          } // if local
        } // if storage
      } // for group
    } END_LOCAL_COMPONENT_LOOP;
    
  }
  
  
  
  void CheckChecksums (const cGH* cgh, const checktimes where)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (! checksum_timelevels) return;
    
    Checkpoint ("%*sCheckChecksums", 2*reflevel, "");
    
    assert ((int)checksums.size()==maxreflevels);
    assert ((int)checksums[reflevel].size()==hh->components(reflevel));
    BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
      assert ((int)checksums[reflevel][component].size()==CCTK_NumGroups());
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        if (CCTK_QueryGroupStorageI(cgh, group)) {
          
          const int nvar = CCTK_NumVarsInGroupI(group);
          
          assert ((int)checksums[reflevel][component][group].size()==nvar);
          
          if (reflevel<arrdata[group].hh->reflevels()
              && component<arrdata[group].hh->components(reflevel)
              && arrdata[group].hh->is_local(reflevel, component)) {
            
            const int n0 = CCTK_FirstVarIndexI(group);
            const int sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
            assert (sz>0);
          
            vect<int,dim> size(1);
            const int gpdim = arrdata[group].info.dim;
            for (int d=0; d<gpdim; ++d) {
              size[d] = arrdata[group].info.lsh[d];
            }
            const int np = prod(size);
            
            const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
            assert (num_tl>0);
            const int min_tl = mintl(where, num_tl);
            const int max_tl = maxtl(where, num_tl);
            
            for (int var=0; var<nvar; ++var) {
              assert ((int)checksums[reflevel][component][group][var].size()==num_tl);
              for (int tl=min_tl; tl<=max_tl; ++tl) {
                
                assert (checksums[reflevel][component][group][var][tl].valid);
                if (checksums[reflevel][component][group][var][tl].valid) {
                  
                  const int n = n0 + var;
                  const void* data = cgh->data[n][tl];
                  unsigned int chk = 0;
                  for (int i=0; i<np*sz/(int)sizeof chk; ++i) {
                    chk += ((const unsigned int*)data)[i];
                  }
                  const bool unexpected_change =
                    chk != checksums[reflevel][component][group][var][tl].sum;
                  
                  if (unexpected_change) {
                    char* fullname = CCTK_FullName(n);
                    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                                "Timelevel %d, component %d, refinement level %d of the variable \"%s\" has changed unexpectedly",
                                tl, component, reflevel, fullname);
                    free (fullname);
                  }
                  
                } // if valid
              } // for tl
            } // for var
          } // if local
        } // if storage
      } // for group
    } END_LOCAL_COMPONENT_LOOP;
  }
  
} // namespace Carpet
