// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gh.cc,v 1.22 2003/09/19 16:06:41 schnetter Exp $

#include <assert.h>
#include <stdlib.h>
#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dh.hh"
#include "th.hh"

#include "gh.hh"

using namespace std;



  // Constructors
template<int D>
gh<D>::gh (const int reffact, const centering refcent,
	   const int mgfact, const centering mgcent,
	   const ibbox& baseextent)
  : reffact(reffact), refcent(refcent),
    mgfact(mgfact), mgcent(mgcent),
    baseextent(baseextent)
{
}

// Destructors
template<int D>
gh<D>::~gh () { }

// Modifiers
template<int D>
void gh<D>::recompose (const rexts& exts,
                       const rbnds& outer_bounds,
		       const rprocs& procs,
                       const int initialise_upto)
{
  DECLARE_CCTK_PARAMETERS;
  
  extents = exts;
  outer_boundaries = outer_bounds;
  processors = procs;
  
  // Consistency checks
  
  // nota bene: there might be 0 refinement levels.
  
  // Check processor number consistency
  for (int rl=0; rl<reflevels(); ++rl) {
    assert (processors.size() == extents.size());
    assert (outer_boundaries.size() == extents.size());
    for (int c=0; c<components(rl); ++c) {
      assert (processors[rl].size() == extents[rl].size());
      assert (outer_boundaries[rl].size() == extents[rl].size());
    }
  }
  
  // Check multigrid consistency
  for (int rl=0; rl<reflevels(); ++rl) {
    for (int c=0; c<components(rl); ++c) {
      assert (mglevels(rl,c)>0);
      for (int ml=1; ml<mglevels(rl,c); ++ml) {
	assert (all(extents[rl][c][ml].stride()
		    == ivect(mgfact) * extents[rl][c][ml-1].stride()));
	assert (extents[rl][c][ml]
		.contracted_for(extents[rl][c][ml-1])
		.is_contained_in(extents[rl][c][ml-1]));
      }
    }
  }
  
  // Check component consistency
  for (int rl=0; rl<reflevels(); ++rl) {
    assert (components(rl)>0);
    for (int c=0; c<components(rl); ++c) {
      for (int ml=0; ml<mglevels(rl,c); ++ml) {
	assert (all(extents[rl][c][ml].stride()
		    == extents[rl][0][ml].stride()));
	assert (extents[rl][c][ml].is_aligned_with(extents[rl][0][ml]));
      }
    }
  }
  
  // Check base grid extent
  if (reflevels()>0) {
    for (int c=0; c<components(0); ++c) {
      assert (extents[0][c][0].is_contained_in(baseextent));
    }
  }
  
  // Check refinement levels
  for (int rl=1; rl<reflevels(); ++rl) {
    assert (all(extents[rl-1][0][0].stride()
		== ivect(reffact) * extents[rl][0][0].stride()));
    // Check contained-ness:
    // first take all coarse grids ...
    bboxset<int,D> all;
    for (int c=0; c<components(rl-1); ++c) {
      all |= extents[rl-1][c][0];
    }
    // ... remember their size ...
    const int sz = all.size();
    // ... then add the coarsified fine grids ...
    for (int c=0; c<components(rl); ++c) {
      all |= extents[rl][c][0].contracted_for(extents[rl-1][0][0]);
    }
    // ... and then check the sizes:
    assert (all.size() == sz);
  }
  
  // Calculate base extents of all levels
  bases.resize(reflevels());
  for (int rl=0; rl<reflevels(); ++rl) {
    if (components(rl)==0) {
      bases[rl].resize(0);
    } else {
      bases[rl].resize(mglevels(rl,0));
      for (int ml=0; ml<mglevels(rl,0); ++ml) {
	bases[rl][ml] = ibbox();
	for (int c=0; c<components(rl); ++c) {
	  bases[rl][ml]
	    = bases[rl][ml].expanded_containing(extents[rl][c][ml]);
	}
      }
    }
  }
  
  if (output_bboxes) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
	for (int ml=0; ml<mglevels(rl,c); ++ml) {
	  cout << endl;
          cout << "gh bboxes:" << endl;
	  cout << "rl=" << rl << " c=" << c << " ml=" << ml << endl;
          cout << "extent=" << extents[rl][c][ml] << endl;
          cout << "outer_boundary=" << outer_boundaries[rl][c] << endl;
          cout << "processor=" << processors[rl][c] << endl;
        }
      }
    }
    for (int rl=0; rl<reflevels(); ++rl) {
      if (components(rl)>0) {
	for (int ml=0; ml<mglevels(rl,0); ++ml) {
	  cout << endl;
          cout << "gh bases:" << endl;
	  cout << "rl=" << rl << " ml=" << ml << endl;
          cout << "base=" << bases[rl][ml] << endl;
        }
      }
    }
  }
  
  // Recompose the other hierarchies
  
  for (typename list<th<D>*>::iterator t=ths.begin(); t!=ths.end(); ++t) {
    (*t)->recompose();
  }
  
  for (typename list<dh<D>*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
    (*d)->recompose(initialise_upto);
  }
}

// Helpers
template<int D>
typename gh<D>::cexts
gh<D>::make_reflevel_multigrid_boxes (const vector<ibbox>& exts,
				      const int mglevels)
  const
{
  assert (mglevels>0);
  
  cexts mexts (exts.size());
  for (int c=0; c<(int)exts.size(); ++c) {
    
    mexts[c].resize(mglevels);
    
    ibbox ext = exts[c];
    for (int ml=0; ml<mglevels; ++ml) {
      
      mexts[c][ml] = ext;
      
      if (ml == mglevels-1) break;
      
      // This level's characteristics
      ivect str = ext.stride();
      ivect lo = ext.lower();
      ivect up = ext.upper();
      
      // Transform to next (coarser) level
      switch (mgcent) {
      case vertex_centered:
	break;
      case cell_centered:
	for (int d=0; d<D; ++d) assert (str[d]%2 == 0);
	lo += str/2;
	break;
      default:
	assert (0);
      }
      str *= mgfact;
      up = up - (up - lo) % str;
      
      ext = ibbox(lo,up,str);
    } // for ml
  } // for c
  
  return mexts;
}

template<int D>
typename gh<D>::rexts
gh<D>::make_multigrid_boxes (const vector<vector<ibbox> >& exts,
			     const int mglevels)
  const
{
  assert (mglevels>0);
  
  rexts mexts (exts.size());
  for (int rl=0; rl<(int)exts.size(); ++rl) {
    mexts[rl] = make_reflevel_multigrid_boxes (exts[rl], mglevels);
  }
  
  return mexts;
}



// Accessors
template<int D>
int gh<D>::local_components (const int rl) const {
  int lc = 0;
  for (int c=0; c<components(rl); ++c) {
    if (is_local(rl,c)) ++lc;
  }
  return lc;
}



// Time hierarchy management
template<int D>
void gh<D>::add (th<D>* t) {
  ths.push_back(t);
}

template<int D>
void gh<D>::remove (th<D>* t) {
  ths.remove(t);
}



// Data hierarchy management
template<int D>
void gh<D>::add (dh<D>* d) {
  dhs.push_back(d);
}

template<int D>
void gh<D>::remove (dh<D>* d) {
  dhs.remove(d);
}



template<int D>
ostream& gh<D>::output (ostream& os) const {
  os << "gh<" << D << ">:"
     << "reffactor=" << reffact << ",refcentering=" << refcent << ","
     << "mgfactor=" << mgfact << ",mgcentering=" << mgcent << ","
     << "baseextent=" << baseextent << ","
     << "extents=" << extents << ","
     << "outer_boundaries=" << outer_boundaries << ","
     << "processors=" << processors << ","
     << "dhs={";
  int cnt=0;
  for (typename list<dh<D>*>::const_iterator d = dhs.begin();
       d != dhs.end(); ++d) {
    if (cnt++) os << ",";
    (*d)->output(os);
  }
  os << "}";
  return os;
}



template class gh<3>;
