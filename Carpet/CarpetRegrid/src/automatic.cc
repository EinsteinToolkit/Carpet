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
                 gh::mexts  & bbsss,
                 gh::rbnds  & obss,
                 gh::rprocs & pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    assert (bbsss.size() >= 1);
    vector<vector<ibbox> > bbss = bbsss.at(0);
    
    
    
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
    
    assert (! smart_outer_boundaries);
    
    vector<ibbox> bbs;
    gh::cbnds obs;
    Automatic_OneLevel
      (cctkGH, hh, reflevel, min(reflevels+1, refinement_levels),
       errorgf, bbs, obs);
    
    // make multiprocessor aware
    gh::cprocs ps;
    SplitRegions (cctkGH, bbs, obs, ps);
    
    if (bbs.size() == 0) {
      // remove all finer levels
      bbss.resize(reflevel+1);
      obss.resize(reflevel+1);
      pss.resize(reflevel+1);
    } else {
      assert (reflevel < (int)bbss.size());
      if (reflevel+1 == (int)bbss.size()) {
	// add a finer level
	bbss.push_back (bbs);
	obss.push_back (obs);
	pss.push_back (ps);
      } else {
	// change a finer level
	bbss.at(reflevel+1) = bbs;
	obss.at(reflevel+1) = obs;
	pss.at(reflevel+1) = ps;
      }
    }
    
    // make multigrid aware
    MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);
    
    return 1;
  }
  
  
  
  void Automatic_OneLevel (const cGH * const cctkGH,
                           const gh & hh,
                           const int rl,
                           const int numrl,
                           const gf<CCTK_REAL> & errorgf,
                           vector<ibbox> & bbs,
                           vector<bbvect> & obs)
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
    
    list<ibbox> bbl;
    for (int c=0; c<hh.components(rl); ++c) {
      const ibbox region = hh.extents().at(ml).at(rl).at(c);
      assert (! region.empty());
      
      const data<CCTK_REAL>& errordata = *errorgf(tl,rl,c,ml);
      
      Automatic_Recursive (cctkGH, hh, errordata, bbl, region, reffact);
    }
    
    if (veryverbose) {
      int numpoints = 0;
      for (list<ibbox>::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
        numpoints += ibb->size();
      }
      cout << "MRA: Chose " << bbl.size() << " regions with a total size of " << numpoints << " to refine." << endl << endl;
    }
    
    // Create bbs from bbl
    assert (bbs.size() == 0);
    bbs.reserve (bbl.size());
    for (list<ibbox>::const_iterator it = bbl.begin(); it != bbl.end(); ++it) {
      bbs.push_back (*it);
    }
    
    // Remove grid points outside the outer boundary
    bbvect const obp (false);
    for (size_t c=0; c<bbs.size(); ++c) {
      const ivect lb = xpose(obp)[0].ifthen
        (bbs.at(c).lower(), max (bbs.at(c).lower(), hh.baseextent.lower()));
      const ivect ub = xpose(obp)[1].ifthen
        (bbs.at(c).upper(), min (bbs.at(c).upper(), hh.baseextent.upper()));
      bbs.at(c) = ibbox(lb, ub, bbs.at(c).stride());
    }
    
    // Create obs from bbs
    obs.resize (bbs.size());
    for (size_t c=0; c<bbs.size(); ++c) {
      assert (hh.bases().size()>0 and hh.bases().at(0).size()>0);
      obs.at(c) = zip ((vect<bool,2> (*) (bool, bool)) &vect<bool,2>::make,
                    bbs.at(c).lower() == hh.baseextent.lower(),
                    bbs.at(c).upper() == hh.baseextent.upper());
    }
    
  }
  
  
  
  void Automatic_Recursive (const cGH * const cctkGH,
                            const gh & hh,
                            const data<CCTK_REAL> & errordata,
                            list<ibbox> & bbl,
                            const ibbox & region,
                            const ivect & reffact)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Just to be sure
    assert (! region.empty());
    
    // Count grid points that need to be refined
    // (this doesn't work yet on multiple processors)
    assert (CCTK_nProcs(cctkGH) == 1);
    int cnt = 0;
    {
      ivect const off
        = (region.lower() - errordata.extent().lower()) / region.stride();
      ivect const len
        = region.shape() / region.stride();
      ivect const lsh
        = errordata.extent().shape() / region.stride();
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
    const int width = maxval(region.shape() / region.stride());
    const CCTK_REAL fraction = (CCTK_REAL) cnt / region.size();
    
    if (cnt == 0) {
      // Don't refine
    } else if (any (reffact*width < 2*minwidth) or fraction >= minfraction) {
      // Refine the whole region
      const ivect lo(region.lower());
      const ivect up(region.upper());
      const ivect str(region.stride());
      const ibbox region0(lo,up+str-str/reffact,str/reffact);
      bbl.push_back (region0);
      if (veryverbose) {
        cout << "MRA: Refining to " << region0 << " size " << region0.size() << " fraction " << fraction << endl;
      }
    } else {
      // Split the region and check recursively
      const int dir = maxloc(region.shape());
      const ivect lo(region.lower());
      const ivect up(region.upper());
      const ivect str(region.stride());
      ivect lo1(lo), lo2(lo);
      ivect up1(up), up2(up);
      const int mgstr = ipow(hh.mgfact, mglevels); // honour multigrid factors
      const int step = str[dir]*mgstr;
      lo2[dir] = ((up[dir]+lo[dir])/2/step)*step;
      up1[dir] = lo2[dir]-str[dir];
      const ibbox region1(lo1,up1,str);
      const ibbox region2(lo2,up2,str);
      assert (region1.is_contained_in(region));
      assert (region2.is_contained_in(region));
      assert ((region1 & region2).empty());
      assert (region1 + region2 == region);
      list<ibbox> bbl1, bbl2;
      Automatic_Recursive (cctkGH, hh, errordata, bbl1, region1, reffact);
      Automatic_Recursive (cctkGH, hh, errordata, bbl2, region2, reffact);
      // Combine regions if possible
      up2 += str-str/reffact;
      up2[dir] = lo2[dir];
      const ibbox iface(lo2,up2,str/reffact);
      Automatic_Recombine (bbl1, bbl2, bbl, iface, dir);
    }
    
  }
  
  
  
  void Automatic_Recombine (list<ibbox> & bbl1,
                            list<ibbox> & bbl2,
                            list<ibbox> & bbl,
                            const ibbox & iface,
                            const int dir)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (!iface.empty());
    assert (iface.lower()[dir] == iface.upper()[dir]);
    
    const int oldnumboxes = bbl.size() + bbl1.size() + bbl2.size();
    int numcombinedboxes = 0;
    
    int oldnumpoints = 0;
    for (list<ibbox>::const_iterator ibb = bbl.begin();
         ibb != bbl.end();
         ++ibb)
    {
      oldnumpoints += ibb->size();
    }
    for (list<ibbox>::const_iterator ibb1 = bbl1.begin();
         ibb1 != bbl1.end();
         ++ibb1)
    {
      oldnumpoints += ibb1->size();
    }
    for (list<ibbox>::const_iterator ibb2 = bbl2.begin();
         ibb2 != bbl2.end();
         ++ibb2)
    {
      oldnumpoints += ibb2->size();
    }
    
#if 0
    // remember old bounding boxes
    bboxset<int,dim> oldboxes;
    for (list<ibbox>::const_iterator ibb1 = bbl1.begin();
         ibb1 != bbl1.end();
         ++ibb1)
    {
      oldboxes += *ibb1;
    }
    for (list<ibbox>::const_iterator ibb2 = bbl2.begin();
         ibb2 != bbl2.end();
         ++ibb2)
    {
      oldboxes += *ibb2;
    }
#endif
#if 0
    cout << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list:" << endl;
    cout << bbl << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list 1:" << endl;
    cout << bbl1 << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list 2:" << endl;
    cout << bbl2 << endl;
#endif
    
    const ivect lo = iface.lower();
    const ivect up = iface.upper();
    const ivect str = iface.stride();
    
    {
      // prune boxes on the left
      list<ibbox>::iterator ibb1 = bbl1.begin();
      while (ibb1 != bbl1.end()) {
	// is this bbox just to the left of the interface?
// 	const ivect lo1 = ibb1->lower();
	const ivect up1 = ibb1->upper();
	const ivect str1 = ibb1->stride();
	assert (up1[dir]+str1[dir] <= lo[dir]);
	assert (all(str1 == str));
	if (up1[dir]+str1[dir] < lo[dir]) {
	  // no: forget it
	  bbl.push_back (*ibb1);
	  ibb1 = bbl1.erase(ibb1);
	  continue;
	}
	++ibb1;
      }	// while
    }
    
    {
      // prune boxes on the right
      list<ibbox>::iterator ibb2 = bbl2.begin();
      while (ibb2 != bbl2.end()) {
	// is this bbox just to the right of the interface?
	const ivect lo2 = ibb2->lower();
// 	const ivect up2 = ibb2->upper();
	const ivect str2 = ibb2->stride();
	assert (up[dir] <= lo2[dir]);
	assert (all(str2 == str));
	if (up[dir] < lo2[dir]) {
	  // no: forget it
	  bbl.push_back (*ibb2);
	  ibb2 = bbl2.erase(ibb2);
	  continue;
	}
	++ibb2;
      }	// while
    }
    
    {
      // walk all boxes on the left
      list<ibbox>::iterator ibb1 = bbl1.begin();
      while (ibb1 != bbl1.end()) {
	ivect lo1 = ibb1->lower();
	ivect up1 = ibb1->upper();
	ivect str1 = ibb1->stride();
	assert (up1[dir]+str1[dir] == lo[dir]);
	lo1[dir] = lo[dir];
	up1[dir] = up[dir];
	const ibbox iface1 (lo1,up1,str1);
	
	{
	  // walk all boxes on the right
	  list<ibbox>::iterator ibb2 = bbl2.begin();
	  while (ibb2 != bbl2.end()) {
	    ivect lo2 = ibb2->lower();
	    ivect up2 = ibb2->upper();
	    ivect str2 = ibb2->stride();
	    assert (lo2[dir] == up[dir]);
	    lo2[dir] = lo[dir];
	    up2[dir] = up[dir];
	    const ibbox iface2 (lo2,up2,str2);
	    
	    // check for a match
	    if (iface1 == iface2) {
	      const ibbox combined (ibb1->lower(), ibb2->upper(), str);
	      bbl.push_back (combined);
	      ibb1 = bbl1.erase(ibb1);
	      ibb2 = bbl2.erase(ibb2);
	      ++numcombinedboxes;
              if (veryverbose) {
                assert (dir<3);
                cout << "MRA: Combining along " << "xyz"[dir] << " to " << bbl.back() << " size " << bbl.back().size() << endl;
              }
	      goto continue_search;
	    }
	    
	    ++ibb2;
	  } // while
	}
	
	++ibb1;
      continue_search:;
      }	// while
    }
    
    bbl.splice (bbl.end(), bbl1);
    bbl.splice (bbl.end(), bbl2);
    
    assert (bbl1.empty() and bbl2.empty());
    
    const int newnumboxes = bbl.size();
    assert (newnumboxes + numcombinedboxes == oldnumboxes);
    
    int newnumpoints = 0;
    for (list<ibbox>::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      newnumpoints += ibb->size();
    }
    assert (newnumpoints == oldnumpoints);
    
#if 0
    // find new bounding boxes
    bboxset<int,dim> newboxes;
    for (list<ibbox>::const_iterator ibb = bbl.begin();
         ibb != bbl.end();
         ++ibb)
    {
      newboxes += *ibb;
    }
    
    // Check that they are equal
    assert (newboxes.size() <= oldboxes.size());
    assert ((newboxes.size()==0) == (oldboxes.size()==0));
    assert (oldboxes == newboxes);
#endif
#if 0
    cout << "MakeRegions_Adaptively_Recombine: final list:" << endl;
    cout << bbl << endl;
    cout << endl;
#endif
  }
  
} // namespace CarpetRegrid
