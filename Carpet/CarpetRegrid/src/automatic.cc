#include <assert.h>
#include <string.h>

#include <algorithm>
#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gf.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Automatic (cGH const * const cctkGH,
                 gh const & hh,
                 gh::rregs & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    
    
    const int vi = CCTK_VarIndex (errorvar);
    assert (vi>=0 and vi<CCTK_NumVars());
    const int gi = CCTK_GroupIndexFromVarI (vi);
    assert (gi>=0 and gi<CCTK_NumGroups());
    const int v1 = CCTK_FirstVarIndexI(gi);
    assert (v1>=0 and v1<=vi and v1<CCTK_NumVars());
    
    assert (CCTK_GroupTypeI(gi) == CCTK_GF);
    assert (CCTK_VarTypeI(vi) == CCTK_VARIABLE_REAL);
    assert (CCTK_GroupDimI(gi) == dim);
    
    assert (arrdata.at(gi).at(Carpet::map).data.at(vi-v1));
    const gf<CCTK_REAL>& errorgf
      = (*dynamic_cast<const gf<CCTK_REAL>*>
         (arrdata.at(gi).at(Carpet::map).data.at(vi-v1)));
    
    assert (not smart_outer_boundaries);
    
    vector<region_t> regs;
    Automatic_OneLevel
      (cctkGH, hh, reflevel, min(reflevels+1, (int)refinement_levels),
       errorgf, regs);
    
    if (regs.size() == 0) {
      // remove all finer levels
      regss.resize(reflevel+1);
    } else {
      assert (reflevel < (int)regss.size());
      if (reflevel+1 == (int)regss.size()) {
	// add a finer level
	regss.push_back (regs);
      } else {
	// change a finer level
	regss.at(reflevel+1) = regs;
      }
    }
    
    return 1;
  }
  
  
  
  void Automatic_OneLevel (const cGH * const cctkGH,
                           const gh & hh,
                           const int rl,
                           const int numrl,
                           const gf<CCTK_REAL> & errorgf,
                           vector<region_t> & regs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (rl+1 >= numrl) return;
    
    // Arbitrary
    const int tl = 0;
    const int ml = 0;
    
    // Determine refinement factor
    ivect const reffact = spacereffacts.at(rl+1) / spacereffacts.at(rl);
    
    if (veryverbose) {
      cout << endl << "MRA: Choosing regions to refine on level " << rl << " in " << hh.components(rl) << " components" << endl;
    }
    
    list<region_t> regl;
    for (int c=0; c<hh.components(rl); ++c) {
      const region_t region = hh.regions.at(ml).at(rl).at(c);
      assert (! region.extent.empty());
      
      const data<CCTK_REAL>& errordata =
        *dynamic_cast<const data<CCTK_REAL>*> (errorgf(tl,rl,c,ml));
      
      Automatic_Recursive (cctkGH, hh, errordata, regl, region, reffact);
    }
    
    if (veryverbose) {
      int numpoints = 0;
      for (list<region_t>::const_iterator ireg = regl.begin(); ireg != regl.end(); ++ireg) {
        numpoints += ireg->extent.size();
      }
      cout << "MRA: Chose " << regl.size() << " regions with a total size of " << numpoints << " to refine." << endl << endl;
    }
    
    // Create regs from regl
    assert (regs.size() == 0);
    regs.reserve (regl.size());
    for (list<region_t>::const_iterator it = regl.begin(); it != regl.end(); ++it) {
      regs.push_back (*it);
    }
    
    // Remove grid points outside the outer boundary
    b2vect const obp (false);
    for (size_t c=0; c<regs.size(); ++c) {
      const ivect lb = either (obp[0],
                               regs.at(c).extent.lower(),
                               max (regs.at(c).extent.lower(),
                                    hh.baseextents.at(0).at(0).lower()));
      const ivect ub = either (obp[1],
                               regs.at(c).extent.upper(),
                               min (regs.at(c).extent.upper(),
                                    hh.baseextents.at(0).at(rl).upper()));
      regs.at(c).extent = ibbox(lb, ub, regs.at(c).extent.stride());
      regs.at(c).outer_boundaries =
        b2vect (regs.at(c).extent.lower() == hh.baseextents.at(0).at(0).lower(),
                regs.at(c).extent.upper() == hh.baseextents.at(0).at(0).upper());
      regs.at(c).map = Carpet::map;
    }
    
  }
  
  
  
  void Automatic_Recursive (const cGH * const cctkGH,
                            const gh & hh,
                            const data<CCTK_REAL> & errordata,
                            list<region_t> & regl,
                            const region_t & region,
                            const ivect & reffact)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Just to be sure
    assert (! region.extent.empty());
    
    // Count grid points that need to be refined
    // (this doesn't work yet on multiple processors)
    assert (CCTK_nProcs(cctkGH) == 1);
    int cnt = 0;
    {
      ivect const off
        = (region.extent.lower() - errordata.extent().lower()) / region.extent.stride();
      ivect const len
        = region.extent.shape() / region.extent.stride();
      ivect const lsh
        = errordata.extent().shape() / region.extent.stride();
      CCTK_REAL const * const errorptr
        = (CCTK_REAL const *) errordata.storage();
      for (int k=0; k<len[2]; ++k) {
        for (int j=0; j<len[1]; ++j) {
          for (int i=0; i<len[0]; ++i) {
            int const ind
              = (off[0] + i) + lsh[0] * ((off[1] + j) + lsh[1] * (off[2] + k));
            if (errorptr[ind] > maxerror) ++cnt;
          }
        }
      }
    }
    const int width = maxval(region.extent.shape() / region.extent.stride());
    const CCTK_REAL fraction = (CCTK_REAL) cnt / region.extent.size();
    
    if (cnt == 0) {
      // Don't refine
    } else if (any (reffact*width < (int)(2*minwidth))
               or fraction >= minfraction) {
      // Refine the whole region
      const ivect lo(region.extent.lower());
      const ivect up(region.extent.upper());
      const ivect str(region.extent.stride());
      const ibbox ext(lo,up+str-str/reffact,str/reffact);
      region_t region0 (region);
      region0.extent = ext;
      regl.push_back (region0);
      if (veryverbose) {
        cout << "MRA: Refining to " << region0.extent << " size " << region0.extent.size() << " fraction " << fraction << endl;
      }
    } else {
      // Split the region and check recursively
      const int dir = maxloc(region.extent.shape());
      const ivect lo(region.extent.lower());
      const ivect up(region.extent.upper());
      const ivect str(region.extent.stride());
      ivect lo1(lo), lo2(lo);
      ivect up1(up), up2(up);
      const int mgstr = ipow(hh.mgfact, mglevels); // honour multigrid factors
      const int step = str[dir]*mgstr;
      lo2[dir] = ((up[dir]+lo[dir])/2/step)*step;
      up1[dir] = lo2[dir]-str[dir];
      const ibbox ext1(lo1,up1,str);
      const ibbox ext2(lo2,up2,str);
      assert (ext1.is_contained_in(region.extent));
      assert (ext2.is_contained_in(region.extent));
      assert ((ext1 & ext2).empty());
      assert (ext1 + ext2 == region.extent);
      region_t region1(region);
      region_t region2(region);
      region1.extent = ext1;
      region2.extent = ext2;
      list<region_t> regl1, regl2;
      Automatic_Recursive (cctkGH, hh, errordata, regl1, region1, reffact);
      Automatic_Recursive (cctkGH, hh, errordata, regl2, region2, reffact);
      // Combine regions if possible
      up2 += str-str/reffact;
      up2[dir] = lo2[dir];
      const ibbox iface(lo2,up2,str/reffact);
      Automatic_Recombine (regl1, regl2, regl, iface, dir);
    }
    
  }
  
  
  
  void Automatic_Recombine (list<region_t> & regl1,
                            list<region_t> & regl2,
                            list<region_t> & regl,
                            const ibbox & iface,
                            const int dir)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (!iface.empty());
    assert (iface.lower()[dir] == iface.upper()[dir]);
    
    const int oldnumboxes = regl.size() + regl1.size() + regl2.size();
    int numcombinedboxes = 0;
    
    int oldnumpoints = 0;
    for (list<region_t>::const_iterator ireg = regl.begin();
         ireg != regl.end();
         ++ireg)
    {
      oldnumpoints += ireg->extent.size();
    }
    for (list<region_t>::const_iterator ireg1 = regl1.begin();
         ireg1 != regl1.end();
         ++ireg1)
    {
      oldnumpoints += ireg1->extent.size();
    }
    for (list<region_t>::const_iterator ireg2 = regl2.begin();
         ireg2 != regl2.end();
         ++ireg2)
    {
      oldnumpoints += ireg2->extent.size();
    }
    
#if 0
    // remember old bounding boxes
    bboxset<int,dim> oldboxes;
    for (list<region_t>::const_iterator ireg1 = regl1.begin();
         ireg1 != regl1.end();
         ++ireg1)
    {
      oldboxes += ireg1->extent;
    }
    for (list<regino_t>::const_iterator ireg2 = regl2.begin();
         ireg2 != regl2.end();
         ++ireg2)
    {
      oldboxes += ireg2;->extent
    }
#endif
#if 0
    cout << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list:" << endl;
    cout << regl << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list 1:" << endl;
    cout << regl1 << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list 2:" << endl;
    cout << regl2 << endl;
#endif
    
    const ivect lo = iface.lower();
    const ivect up = iface.upper();
    const ivect str = iface.stride();
    
    {
      // prune boxes on the left
      list<region_t>::iterator ireg1 = regl1.begin();
      while (ireg1 != regl1.end()) {
	// is this bbox just to the left of the interface?
// 	const ivect lo1 = ireg1->extent.lower();
	const ivect up1 = ireg1->extent.upper();
	const ivect str1 = ireg1->extent.stride();
	assert (up1[dir]+str1[dir] <= lo[dir]);
	assert (all(str1 == str));
	if (up1[dir]+str1[dir] < lo[dir]) {
	  // no: forget it
	  regl.push_back (*ireg1);
	  ireg1 = regl1.erase(ireg1);
	  continue;
	}
	++ireg1;
      }	// while
    }
    
    {
      // prune boxes on the right
      list<region_t>::iterator ireg2 = regl2.begin();
      while (ireg2 != regl2.end()) {
	// is this regox just to the right of the interface?
	const ivect lo2 = ireg2->extent.lower();
// 	const ivect up2 = ireg2->extent.upper();
	const ivect str2 = ireg2->extent.stride();
	assert (up[dir] <= lo2[dir]);
	assert (all(str2 == str));
	if (up[dir] < lo2[dir]) {
	  // no: forget it
	  regl.push_back (*ireg2);
	  ireg2 = regl2.erase(ireg2);
	  continue;
	}
	++ireg2;
      }	// while
    }
    
    {
      // walk all boxes on the left
      list<region_t>::iterator ireg1 = regl1.begin();
      while (ireg1 != regl1.end()) {
	ivect lo1 = ireg1->extent.lower();
	ivect up1 = ireg1->extent.upper();
	ivect str1 = ireg1->extent.stride();
	assert (up1[dir]+str1[dir] == lo[dir]);
	lo1[dir] = lo[dir];
	up1[dir] = up[dir];
	const ibbox iface1 (lo1,up1,str1);
	
	{
	  // walk all boxes on the right
	  list<region_t>::iterator ireg2 = regl2.begin();
	  while (ireg2 != regl2.end()) {
	    ivect lo2 = ireg2->extent.lower();
	    ivect up2 = ireg2->extent.upper();
	    ivect str2 = ireg2->extent.stride();
	    assert (lo2[dir] == up[dir]);
	    lo2[dir] = lo[dir];
	    up2[dir] = up[dir];
	    const ibbox iface2 (lo2,up2,str2);
	    
	    // check for a match
	    if (iface1 == iface2) {
              region_t combined (*ireg1);
              combined.extent = ibbox (ireg1->extent.lower(), ireg2->extent.upper(), str);
              combined.outer_boundaries[0] = ireg1->outer_boundaries[0];
              combined.outer_boundaries[1] = ireg2->outer_boundaries[0];
	      regl.push_back (combined);
	      ireg1 = regl1.erase(ireg1);
	      ireg2 = regl2.erase(ireg2);
	      ++numcombinedboxes;
              if (veryverbose) {
                assert (dir<3);
                cout << "MRA: Combining along " << "xyz"[dir] << " to " << regl.back().extent << " size " << regl.back().extent.size() << endl;
              }
	      goto continue_search;
	    }
	    
	    ++ireg2;
	  } // while
	}
	
	++ireg1;
      continue_search:;
      }	// while
    }
    
    regl.splice (regl.end(), regl1);
    regl.splice (regl.end(), regl2);
    
    assert (regl1.empty() and regl2.empty());
    
    const int newnumboxes = regl.size();
    assert (newnumboxes + numcombinedboxes == oldnumboxes);
    
    int newnumpoints = 0;
    for (list<region_t>::const_iterator ireg = regl.begin(); ireg != regl.end(); ++ireg) {
      newnumpoints += ireg->extent.size();
    }
    assert (newnumpoints == oldnumpoints);
    
#if 0
    // find new bounding boxes
    bboxset<int,dim> newboxes;
    for (list<region_t>::const_iterator ireg = regl.begin();
         ireg != regl.end();
         ++ireg)
    {
      newboxes += ireg->extent;
    }
    
    // Check that they are equal
    assert (newboxes.size() <= oldboxes.size());
    assert ((newboxes.size()==0) == (oldboxes.size()==0));
    assert (oldboxes == newboxes);
#endif
#if 0
    cout << "MakeRegions_Adaptively_Recombine: final list:" << endl;
    cout << regl << endl;
    cout << endl;
#endif
  }
  
} // namespace CarpetRegrid
