#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "CAR.hh"

extern "C" {
  static const char* rcsid = "$Header:$";
  CCTK_FILEVERSION(Carpet_CarpetAdaptiveregrid_regrid_cc);
}



namespace CarpetAdaptiveRegrid {
  
  using namespace std;
  using namespace Carpet;

  extern "C" {
    void CCTK_FCALL CCTK_FNAME(copy_mask)
      (const CCTK_INT& snx, const CCTK_INT& sny, const CCTK_INT& snz,
       const CCTK_INT* smask, const CCTK_INT sbbox[3][3],
       const CCTK_INT& dnx, const CCTK_INT& dny, const CCTK_INT& dnz,
       CCTK_INT* dmask, const CCTK_INT dbbox[3][3]);
    void CCTK_FCALL CCTK_FNAME(check_box)
      (const CCTK_INT& nx, const CCTK_INT& ny, const CCTK_INT& nz,
       const CCTK_INT* mask,
       CCTK_INT* sum_x, CCTK_INT* sum_y, CCTK_INT* sum_z,
       CCTK_INT* sig_x, CCTK_INT* sig_y, CCTK_INT* sig_z,
       const CCTK_INT bbox[3][3],
       CCTK_INT newbbox1[3][3], CCTK_INT newbbox2[3][3],
       const CCTK_INT& min_width, const CCTK_REAL& min_fraction,
       CCTK_INT& didit);
  }
  
  CCTK_INT CarpetAdaptiveRegrid_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                                CCTK_POINTER const bbsss_,
                                CCTK_POINTER const obss_,
                                CCTK_POINTER const pss_,
				CCTK_INT force)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const cGH * const cctkGH = (const cGH *) cctkGH_;
    
    gh::mexts  & bbsss = * (gh::mexts  *) bbsss_;
    gh::rbnds  & obss  = * (gh::rbnds  *) obss_;
    gh::rprocs & pss   = * (gh::rprocs *) pss_;
    
    gh const & hh = *vhh.at(Carpet::map);
    
    assert (is_singlemap_mode());
    
    // In force mode (force == true) we do not check the
    // CarpetAdaptiveregrid parameters

    if (!force) {

      assert (regrid_every == -1 || regrid_every == 0
	      || regrid_every % maxmglevelfact == 0);
    
      // Return if no regridding is desired
      if (regrid_every == -1) return 0;
      
      // Return if we want to regrid during initial data only, and this
      // is not the time for initial data
      if (regrid_every == 0 && cctkGH->cctk_iteration != 0) return 0;

      // Return if we want to regrid regularly, but not at this time
      if (regrid_every > 0 && cctkGH->cctk_iteration != 0
	  && (cctkGH->cctk_iteration-1) % regrid_every != 0)
      {
	return 0;
      }

      // Return if it's initial data as we can't handle that yet.
      if (cctkGH->cctk_iteration ==0) {
        return 0;
      }

    }    

    if (reflevel == maxreflevels - 1) return 0;

    int do_recompose;
    do_recompose = 1;

    // We assume that we are called in fine->coarse order
    // So that on finer levels regridding has already been done

    // So the full algorithm should look something like:

    // Find how big the first bounding box should be on this level
    //   Do this by finding min lower and max upper bounds of all bboxes
    // Allocate box
    // Fill errors from local arrays
    // If grandchildren exist use their bboxes (expanded) to add to errors
    // Reduce errors (MPI sum)
    // Set errors to 0/1
    // Define vectors of bboxes final (empty) and todo (contains bbox)
    // Define vector of masks (contains error mask)
    // Loop over all entries in todo:
    //   Setup appropriate 1d array memory
    //   Call fortran routine
    //   If return is:
    //     zero: add bbox to final
    //     one:  add new bbox to todo and assoc mask to masklist
    //     two:  add both new bboxs to todo and assoc masks to masklist
    //   

    vector<ibbox> bbs = bbsss.at(mglevel).at(reflevel);
    
    ivect low = 100000, upp = -1;

    for ( vector<ibbox>::const_iterator bbi = bbs.begin(); 
          bbi != bbs.end();
          ++bbi) 
    {
      low = min(low, bbi->lower());
      upp = max(upp, bbi->upper());
    }
    
    // low and upp now define the starting bbox.

    ibbox bb(low, upp, bbs.at(0).stride());
    
    vector<CCTK_INT> mask(prod(bb.shape()/bb.stride()), 0);
    
    // Setup the mask.
      
    const ibbox& baseext = 
      vdd.at(Carpet::map)->bases.at(mglevel).at(reflevel).exterior;
    ivect imin = (bb.lower() - baseext.lower())/bb.stride(), 
      imax = (bb.upper() - baseext.lower())/bb.stride();

    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) 
    {
      const CCTK_REAL *error_var_ptr = 
        static_cast<const CCTK_REAL*>(CCTK_VarDataPtr(cctkGH, 0, error_var));
     
      assert(all(imin >= 0));
      assert(all(imax >= 0));
      assert(all(imin < ivect::ref(cctkGH->cctk_lsh)));
      assert(all(imax < ivect::ref(cctkGH->cctk_lsh)));
      assert(all(imin <= imax));
      
      for (CCTK_INT k = 0; k < cctkGH->cctk_lsh[2]; ++k) {
        for (CCTK_INT j = 0; j < cctkGH->cctk_lsh[1]; ++j) {
          for (CCTK_INT i = 0; i < cctkGH->cctk_lsh[0]; ++i) {
            CCTK_INT index = CCTK_GFINDEX3D(cctkGH, i, j, k);
            if (abs(error_var_ptr[index]) > max_error) {
              CCTK_INT ii = i + cctkGH->cctk_lbnd[0] - imin[0];
              CCTK_INT jj = j + cctkGH->cctk_lbnd[1] - imin[1];
              CCTK_INT kk = k + cctkGH->cctk_lbnd[2] - imin[2];
              assert (ii >= 0);
              assert (jj >= 0);
              assert (kk >= 0);
              assert (ii < imax[0] - imin[0]);
              assert (jj < imax[1] - imin[1]);
              assert (kk < imax[2] - imin[2]);
              CCTK_INT mindex = ii + 
                (imax[0]-imin[0])*(jj + (imax[1] - imin[1]) * kk);
              mask[mindex] = 1;
            }
          }
        }
      }
    } END_LOCAL_COMPONENT_LOOP;

    // FIXME    
    // If grandchildren exist use their bboxes (expanded) to add to errors
    // Reduce errors (MPI sum)
    // Set errors to 0/1

    // Pad the errors: stage 1 - buffer points marked as 2.

    for (CCTK_INT k = 0; k < imax[2] - imin[2]; k++) {
      for (CCTK_INT j = 0; j < imax[1] - imin[1]; j++) {
        for (CCTK_INT i = 0; i < imax[0] - imin[0]; i++) {
          CCTK_INT index =  i + 
                (imax[0]-imin[0])*(j + (imax[1] - imin[1]) * k);
          if (mask[index] == 1) {
            for (CCTK_INT kk = max(k-pad, 0); 
                 kk < min(k+pad, imax[2] - imin[2]);
                 ++kk)
            {
              for (CCTK_INT jj = max(j-pad, 0); 
                   jj < min(j+pad, imax[1] - imin[1]);
                   ++jj)
              {
                for (CCTK_INT ii = max(i-pad, 0); 
                     ii < min(i+pad, imax[0] - imin[0]);
                     ++ii)
                {
                  CCTK_INT mindex = ii + 
                    (imax[0]-imin[0])*(jj + (imax[1] - imin[1]) * kk);
                  if (!mask[mindex]) mask[mindex] = 2;
                }
              }
            }
          }
        }
      }
    }
    // stage 2: all buffer points marked truly in error.
    // Also mark if there are any errors.
    bool should_regrid = false;
    for (CCTK_INT k = 0; k < imax[2] - imin[2]; k++) {
      for (CCTK_INT j = 0; j < imax[1] - imin[1]; j++) {
        for (CCTK_INT i = 0; i < imax[0] - imin[0]; i++) {
          CCTK_INT index =  i + 
                (imax[0]-imin[0])*(j + (imax[1] - imin[1]) * k);
          if (mask[index] > 1) mask[index] = 1;
          should_regrid |= (mask[index]);
        }
      }
    }    

    // Define vectors of bboxes final (empty) and todo (contains bbox)
      
    stack<ibbox> final;

    vector<vector<ibbox> > bbss = bbsss.at(0);
      

    if (should_regrid) {
    
      stack<ibbox> todo;
  
      todo.push(bb);
      
      // Define vector of masks (contains error mask)
      
      stack<vector<CCTK_INT> > masklist;
      
      masklist.push(mask);
      
      // Loop over all entries in todo:
      //   Setup appropriate 1d array memory
      //   Call fortran routine
      //   If return is:
      //     zero: add bbox to final
      //     one:  add new bbox to todo and assoc mask to masklist
      //     two:  add both new bboxs to todo and assoc masks to masklist
      
      while (!todo.empty())
      {
      
        ibbox bb = todo.top(); todo.pop();
        vector<CCTK_INT> mask = masklist.top(); masklist.pop();
        
        CCTK_INT nx = bb.shape()[0]/bb.stride()[0];
        CCTK_INT ny = bb.shape()[1]/bb.stride()[1];
        CCTK_INT nz = bb.shape()[2]/bb.stride()[2];
        
        vector<CCTK_INT> sum_x(nx, 0);
        vector<CCTK_INT> sig_x(nx, 0);
        vector<CCTK_INT> sum_y(ny, 0);
        vector<CCTK_INT> sig_y(ny, 0);
        vector<CCTK_INT> sum_z(nz, 0);
        vector<CCTK_INT> sig_z(nz, 0);
        
        int fbbox[3][3], fbbox1[3][3], fbbox2[3][3];
        
        for (CCTK_INT d = 0; d < 3; ++d) {
          fbbox[0][d] = bb.lower()[d];
          fbbox[1][d] = bb.upper()[d];
          fbbox[2][d] = bb.stride()[d];
        }
        
        CCTK_INT didit;
        
        CCTK_FNAME(check_box)(nx, ny, nz,
                              &mask.front(),
                              &sum_x.front(), &sum_y.front(), &sum_z.front(),
                              &sig_x.front(), &sig_y.front(), &sig_z.front(),
                              fbbox,
                              fbbox1, fbbox2,
                              min_width, min_fraction,
                              didit);
        
        if (didit == 0) {
          
          final.push(bb);
          
        }
        else if (didit == 1) {
          
          ibbox newbbox1(ivect::ref(&fbbox1[0][0]),
                         ivect::ref(&fbbox1[1][0]),
                         ivect::ref(&fbbox1[2][0]));
          todo.push(newbbox1);
          
          CCTK_INT dnx = newbbox1.shape()[0]/newbbox1.stride()[0];
          CCTK_INT dny = newbbox1.shape()[1]/newbbox1.stride()[1];
          CCTK_INT dnz = newbbox1.shape()[2]/newbbox1.stride()[2];
          
          vector<CCTK_INT >  
            newmask1(prod(newbbox1.shape()/newbbox1.stride()), 0);
          
          CCTK_FNAME(copy_mask)(nx, ny, nz,
                                &mask.front(), fbbox,
                                dnx, dny, dnz,
                                &newmask1.front(), fbbox1);
          masklist.push(newmask1);
        }
        else if (didit == 2) {
          
          ibbox newbbox1(ivect::ref(&fbbox1[0][0]),
                         ivect::ref(&fbbox1[1][0]),
                         ivect::ref(&fbbox1[2][0]));
          todo.push(newbbox1);
          ibbox newbbox2(ivect::ref(&fbbox2[0][0]),
                         ivect::ref(&fbbox2[1][0]),
                         ivect::ref(&fbbox2[2][0]));
          todo.push(newbbox2);
          
          CCTK_INT dnx = newbbox1.shape()[0]/newbbox1.stride()[0];
          CCTK_INT dny = newbbox1.shape()[1]/newbbox1.stride()[1];
          CCTK_INT dnz = newbbox1.shape()[2]/newbbox1.stride()[2];
          
          vector<CCTK_INT >  
            newmask1(prod(newbbox1.shape()/newbbox1.stride()), 0);
          
          CCTK_FNAME(copy_mask)(nx, ny, nz,
                                &mask.front(), fbbox,
                                dnx, dny, dnz,
                                &newmask1.front(), fbbox1);
          masklist.push(newmask1);
          
          dnx = newbbox2.shape()[0]/newbbox2.stride()[0];
          dny = newbbox2.shape()[1]/newbbox2.stride()[1];
          dnz = newbbox2.shape()[2]/newbbox2.stride()[2];
          
          vector<CCTK_INT >  
            newmask2(prod(newbbox2.shape()/newbbox2.stride()), 0);
          
          CCTK_FNAME(copy_mask)(nx, ny, nz,
                                &mask.front(), fbbox,
                                dnx, dny, dnz,
                                &newmask2.front(), fbbox2);
          masklist.push(newmask2);
        }
        else {
          CCTK_WARN(0, "The fortran routine must be confused.");
        }
        
      } // loop over todo vector (boxes needing to be done).
      
      // Fixup the stride
      vector<ibbox> newbbs;
      vector<bbvect> obs;
      while (! final.empty()) {
        ibbox bb = final.top(); final.pop();
        ivect lo = bb.lower();
        ivect hi = bb.upper();
        ivect str = bb.stride() / reffact;
        newbbs.push_back (ibbox(lo,hi,str));
        
        bbvect ob(false);
        obs.push_back(ob);
      }
      // FIXME: check if the newbbs is really different from the current bbs
      //        if not, set do_recompose = 0
      bbs = newbbs;

      cout << bbs << endl;
      
      // make multiprocessor aware
      gh::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
    

      if (bbss.size() < reflevel+2) {
        bbss.resize(reflevel+2);
        obss.resize(reflevel+2);
        pss.resize(reflevel+2);
      }
      
      bbss.at(reflevel+1) = bbs;
      obss.at(reflevel+1) = obs;
      pss.at(reflevel+1) = ps;
    
    } // should_regrid?
    else
    {
      if (bbss.size() > reflevel+1) {
        bbss.resize(reflevel+1);
        obss.resize(reflevel+1);
        pss.resize(reflevel+1);
      }
    }
    
    // make multigrid aware
    MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);


    // Make a cboxlist (or set) from the boxes on the level.

    // cboxlist is list of cbox's.

    // cbox is a bbox 
    // plus a dim-D array of booleans mask(i,j,k)
    // plus 3 sum arrays (sigma_i = sum_{j,k} mask(i,:,:) etc)
    // plus 3 signature arrays (1D laplacian of sigma arrays)

    // First thing to do; make a cbox containing all active boxes
    // on this level.
    // cbox list contains this cbox.
    // Create the mask for this cbox by looking at the error GF.

    // Then recursively for every cbox in the list until no changes:

    //   Prune the cbox (remove all zeros from the ends)

    //   Split the cbox (find zero in sigma_:, split along that line)

    //   if fraction too low then find zero crossings of signature and
    //   split there

    // if none of these steps are taken the cbox is finished.

    // Most of the actual work is done by the utility functions in
    // CAR_utils.F77.
    // What we have to do is

    // loop over every box in the list
    //   for each box call check_box
    //   if it returns zero the box is done
    //   if it returns one then the box needs shrinking
    //   if it returns two then the box needs splitting

    // so, the method will be:

    // Create a cboxlist with one member

    // loop over the list:
    //   call check_box
    //   if the return is two then create a new box from newbbox2 and add it
    //     to the end of the list
    //   if the return value is one or two shrink the current box to the
    //     values given in newbbox1
    //   if the return value is one or two go redo this loop again.
    
    return do_recompose;
  }
  
} // namespace CarpetAdaptiveRegrid
