/***************************************************************************
                          gh.cc  -  Grid Hierarchy
                          bounding boxes for each multigrid level of each
                          component of each refinement level
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gh.cc,v 1.14 2002/08/02 15:11:16 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <iostream>

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
  : dimgeneric_gh (reffact, refcent, mgfact, mgcent),
    baseextent(baseextent)
{
}

// Destructors
template<int D>
gh<D>::~gh () { }

// Modifiers
template<int D>
void gh<D>::recompose (const rexts& exts, const rbnds& outer_bounds,
		       const rprocs& procs) {
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
    const int s = all.size();
    // ... then add the coarsified fine grids ...
    for (int c=0; c<components(rl); ++c) {
      all |= extents[rl][c][0].contracted_for(extents[rl-1][0][0]);
    }
    // ... and then check the sizes:
    assert (all.size() == s);
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
  
  // Recompose the other hierarchies
  
  for (list<th*>::iterator t=ths.begin(); t!=ths.end(); ++t) {
    (*t)->recompose();
  }
  
  for (list<dh<D>*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
    (*d)->recompose();
  }
}

// Helpers
template<int D>
gh<D>::cexts gh<D>::make_reflevel_multigrid_boxes (const vector<ibbox>& exts,
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
gh<D>::rexts gh<D>::make_multigrid_boxes (const vector<vector<ibbox> >& exts,
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
     << "dhs={";
  int cnt=0;
  for (list<dh<D>*>::const_iterator d = dhs.begin();
       d != dhs.end(); ++d) {
    if (cnt++) os << ",";
    (*d)->output(os);
  }
  os << "}";
  return os;
}



template class gh<3>;
