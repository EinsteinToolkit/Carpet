/***************************************************************************
                          dh.cc  -  Data Hierarchy
                          A grid hierarchy plus ghost zones
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dh.cc,v 1.21 2002/10/14 20:40:38 schnetter Exp $

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

#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "vect.hh"

#include "dh.hh"

#undef DEBUG_OUTPUT

using namespace std;



// Constructors
template<int D>
dh<D>::dh (gh<D>& h, const ivect& lghosts, const ivect& ughosts,
	   int prolongation_order_space)
  : dimgeneric_dh(prolongation_order_space),
    h(h),
    lghosts(lghosts), ughosts(ughosts)
{
  assert (all(lghosts>=0 && ughosts>=0));
  h.add(this);
  CHECKPOINT;
  recompose();
}

// Destructors
template<int D>
dh<D>::~dh ()
{
  CHECKPOINT;
  h.remove(this);
}

// Modifiers
template<int D>
void dh<D>::recompose () {
  CHECKPOINT;
  
  boxes.clear();
  
  boxes.resize(h.reflevels());
  for (int rl=0; rl<h.reflevels(); ++rl) {
    boxes[rl].resize(h.components(rl));
    for (int c=0; c<h.components(rl); ++c) {
      boxes[rl][c].resize(h.mglevels(rl,c));
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
       	const ibbox intr = h.extents[rl][c][ml];
	
       	// Interior
       	// (the interior of the grid has the extent as specified by
       	// the user)
       	boxes[rl][c][ml].interior = intr;
	
       	// Exterior (add ghost zones)
       	// (the content of the exterior is completely determined by
       	// the interior of this or other components; the content of
       	// the exterior is redundant)
	ivect ldist(lghosts), udist(ughosts);
	for (int d=0; d<D; ++d) {
	  if (h.outer_boundaries[rl][c][d][0]) ldist[d] = 0;
	  if (h.outer_boundaries[rl][c][d][1]) udist[d] = 0;
	}
       	boxes[rl][c][ml].exterior = intr.expand(ldist, udist);
	
       	// Boundaries (ghost zones only)
       	// (interior + boundaries = exterior)
       	boxes[rl][c][ml].boundaries = boxes[rl][c][ml].exterior - intr;
	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	
      	// Sync boxes
      	const int cs = h.components(rl);
      	boxes[rl][c][ml].send_sync.resize(cs);
      	boxes[rl][c][ml].recv_sync.resize(cs);
      	
      	// Refinement boxes
      	if (rl>0) {
      	  const int csm1 = h.components(rl-1);
      	  boxes[rl][c][ml].send_ref_coarse.resize(csm1);
      	  boxes[rl][c][ml].recv_ref_coarse.resize(csm1);
      	  boxes[rl][c][ml].recv_ref_bnd_coarse.resize(csm1);
      	}
      	if (rl<h.reflevels()-1) {
      	  const int csp1 = h.components(rl+1);
      	  boxes[rl][c][ml].recv_ref_fine.resize(csp1);
      	  boxes[rl][c][ml].send_ref_fine.resize(csp1);
      	  boxes[rl][c][ml].send_ref_bnd_fine.resize(csp1);
      	}
	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	const ibbox intr = boxes[rl][c][ml].interior;
      	const ibbox extr = boxes[rl][c][ml].exterior;

      	const ibset bnds = boxes[rl][c][ml].boundaries;
	
      	// Sync boxes
      	for (int cc=0; cc<h.components(rl); ++cc) {
      	  assert (ml<h.mglevels(rl,cc));
      	  // intersect boundaries with interior of that component
      	  const ibset ovlp = bnds & boxes[rl][cc][ml].interior;
      	  for (typename ibset::const_iterator b=ovlp.begin();
	       b!=ovlp.end(); ++b) {
	    boxes[rl][c ][ml].recv_sync[cc].push_back(*b);
	    boxes[rl][cc][ml].send_sync[c ].push_back(*b);
	  }
      	}
      	
      	// Multigrid boxes
      	if (ml>0) {
      	  const ibbox intrf = boxes[rl][c][ml-1].interior;
      	  const ibbox extrf = boxes[rl][c][ml-1].exterior;
      	  // Restriction (interior)
      	  {
      	    // (the restriction must fill all of the interior of the
      	    // coarse grid, and may use the exterior of the fine grid)
      	    const ibbox recv = intr;
      	    const ibbox send = recv.expanded_for(extrf);
      	    assert (send.is_contained_in(extrf));
      	    boxes[rl][c][ml-1].send_mg_coarse.push_back(send);
      	    boxes[rl][c][ml  ].recv_mg_fine  .push_back(recv);
      	  }
      	  // Prolongation (interior)
      	  {
      	    // (the prolongation may use the exterior of the coarse
      	    // grid, and may fill only the interior of the fine grid,
      	    // and the bbox must be as large as possible)
      	    const ibbox recv = extr.contracted_for(intrf) & intrf;
      	    const ibbox send = recv.expanded_for(extr);
      	    boxes[rl][c][ml-1].recv_mg_coarse.push_back(recv);
      	    boxes[rl][c][ml  ].send_mg_fine  .push_back(send);
      	  }
      	}
      	
      	// Refinement boxes
      	if (rl<h.reflevels()-1) {
      	  for (int cc=0; cc<h.components(rl+1); ++cc) {
      	    const ibbox intrf = boxes[rl+1][cc][ml].interior;
// 	    const ibbox extrf = boxes[rl+1][cc][ml].exterior;
      	    // Restriction (interior)
      	    {
      	      // (the restriction may fill the interior of the of the
      	      // coarse grid, and may use the interior of the fine
      	      // grid, and the bbox must be as large as possible)
      	      const ibbox recv = intrf.contracted_for(intr) & intr;
      	      const ibbox send = recv.expanded_for(intrf);
      	      boxes[rl+1][cc][ml].send_ref_coarse[c ].push_back(send);
      	      boxes[rl  ][c ][ml].recv_ref_fine  [cc].push_back(recv);
      	    }
      	    // Prolongation (interior)
      	    {
      	      // (the prolongation may use the exterior of the coarse
      	      // grid, and must fill all of the interior of the fine
      	      // grid)
      	      const ibbox recv = extr.contracted_for(intrf) & intrf;
      	      const ibbox send = recv.expanded_for(extr);
      	      assert (send.is_contained_in(extr));
      	      boxes[rl+1][cc][ml].recv_ref_coarse[c ].push_back(recv);
      	      boxes[rl  ][c ][ml].send_ref_fine  [cc].push_back(send);
      	    }
      	    // Prolongation (boundaries)
      	    {
      	      ibset bndsf = boxes[rl+1][cc][ml].boundaries;
      	      // coarsify boundaries of fine component
      	      for (typename ibset::const_iterator bi=bndsf.begin();
		   bi!=bndsf.end(); ++bi) {
		const ibbox& bndf = *bi;
		// (the prolongation may use the exterior of the
		// coarse grid, and must fill all of the boundary of
		// the fine grid)
		const int pss = prolongation_stencil_size();
		const ibbox recv = (extr.expand(-pss,-pss).contracted_for(bndf)
				    & bndf);
		const ibbox send = recv.expanded_for(extr);
		assert (send.is_contained_in(extr));
		boxes[rl+1][cc][ml].recv_ref_bnd_coarse[c ].push_back(recv);
		boxes[rl  ][c ][ml].send_ref_bnd_fine  [cc].push_back(send);
      	      }
      	    }
      	
      	  }
      	}
	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	
	// Boundaries that are not synced, or are neither synced nor
	// prolonged to from coarser grids (outer boundaries)
	ibset& sync_not = boxes[rl][c][ml].sync_not;
	ibset& recv_not = boxes[rl][c][ml].recv_not;
	
	// The whole boundary
	sync_not = boxes[rl][c][ml].boundaries;
	recv_not = boxes[rl][c][ml].boundaries;
	
	// Subtract boxes received during synchronisation
	const iblistvect& recv_sync = boxes[rl][c][ml].recv_sync;
	for (typename iblistvect::const_iterator lvi=recv_sync.begin();
	     lvi!=recv_sync.end(); ++lvi) {
	  for (typename iblist::const_iterator li=lvi->begin();
	       li!=lvi->end(); ++li) {
	    sync_not -= *li;
	    recv_not -= *li;
	  }
	}
	
	// Subtract all boxes received
	const iblistvect& recv_ref_bnd_coarse
	  = boxes[rl][c][ml].recv_ref_bnd_coarse;
	for (typename iblistvect::const_iterator
	       lvi=recv_ref_bnd_coarse.begin();
	     lvi!=recv_ref_bnd_coarse.end(); ++lvi) {
	  for (typename iblist::const_iterator li=lvi->begin();
	       li!=lvi->end(); ++li) {
	    recv_not -= *li;
	  }
	}
	
	// Check that no boundaries are left over
	if (rl==0) assert (sync_not.empty());
	assert (recv_not.empty());
	
      } // for ml
    } // for c
  } // for rl
  
  bases.resize(h.reflevels());
  for (int rl=0; rl<h.reflevels(); ++rl) {
    if (h.components(rl)==0) {
      bases[rl].resize(0);
    } else {
      bases[rl].resize(h.mglevels(rl,0));
      for (int ml=0; ml<h.mglevels(rl,0); ++ml) {
	bases[rl][ml].exterior = ibbox();
	bases[rl][ml].interior = ibbox();
	for (int c=0; c<h.components(rl); ++c) {
	  bases[rl][ml].exterior
	    = (bases[rl][ml].exterior
	       .expanded_containing(boxes[rl][c][ml].exterior));
	  bases[rl][ml].interior
	    = (bases[rl][ml].interior
	       .expanded_containing(boxes[rl][c][ml].interior));
	}
	bases[rl][ml].boundaries
	  = bases[rl][ml].exterior - bases[rl][ml].interior;
      }
    }
  }
  
#ifdef DEBUG_OUTPUT
  cout << endl << h << endl;
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	cout << endl;
      	cout << "dh bboxes:" << endl;
      	cout << "rl=" << rl << " c=" << c << " ml=" << ml << endl;
      	cout << "exterior=" << boxes[rl][c][ml].exterior << endl;
      	cout << "interior=" << boxes[rl][c][ml].interior << endl;
      	cout << "send_mg_fine=" << boxes[rl][c][ml].send_mg_fine << endl;
      	cout << "send_mg_coarse=" << boxes[rl][c][ml].send_mg_coarse << endl;
      	cout << "recv_mg_fine=" << boxes[rl][c][ml].recv_mg_fine << endl;
      	cout << "recv_mg_coarse=" << boxes[rl][c][ml].recv_mg_coarse << endl;
      	cout << "send_ref_fine=" << boxes[rl][c][ml].send_ref_fine << endl;
      	cout << "send_ref_coarse=" << boxes[rl][c][ml].send_ref_coarse << endl;
      	cout << "recv_ref_fine=" << boxes[rl][c][ml].recv_ref_fine << endl;
      	cout << "recv_ref_coarse=" << boxes[rl][c][ml].recv_ref_coarse << endl;
      	cout << "send_sync=" << boxes[rl][c][ml].send_sync << endl;
      	cout << "send_ref_bnd_fine=" << boxes[rl][c][ml].send_ref_bnd_fine << endl;
      	cout << "boundaries=" << boxes[rl][c][ml].boundaries << endl;
      	cout << "recv_sync=" << boxes[rl][c][ml].recv_sync << endl;
      	cout << "recv_ref_bnd_coarse=" << boxes[rl][c][ml].recv_ref_bnd_coarse << endl;
      	cout << "sync_not=" << boxes[rl][c][ml].sync_not << endl;
      	cout << "recv_not=" << boxes[rl][c][ml].recv_not << endl;
      }
    }
  }
  for (int rl=0; rl<h.reflevels(); ++rl) {
    if (h.components(rl)>0) {
      for (int ml=0; ml<h.mglevels(rl,0); ++ml) {
      	cout << endl;
      	cout << "dh bases:" << endl;
      	cout << "rl=" << rl << " ml=" << ml << endl;
      	cout << "exterior=" << bases[rl][ml].exterior << endl;
      	cout << "interior=" << bases[rl][ml].interior << endl;
      	cout << "boundaries=" << bases[rl][ml].boundaries << endl;
      }
    }
  }
#endif
  
  for (typename list<generic_gf<D>*>::iterator f=gfs.begin();
       f!=gfs.end(); ++f) {
    (*f)->recompose();
  }
}



// Grid function management
template<int D>
void dh<D>::add (generic_gf<D>* f) {
  CHECKPOINT;
  gfs.push_back(f);
}

template<int D>
void dh<D>::remove (generic_gf<D>* f) {
  CHECKPOINT;
  gfs.remove(f);
}



// Output
template<int D>
void dh<D>::output (ostream& os) const {
  os << "dh<" << D << ">:"
     << "ghosts=[" << lghosts << "," << ughosts << "],"
     << "gfs={";
  int cnt=0;
  for (typename list<generic_gf<D>*>::const_iterator f = gfs.begin();
       f != gfs.end(); ++f) {
    if (cnt++) os << ",";
    (*f)->output(os);
  }
  os << "}";
}



template class dh<3>;
