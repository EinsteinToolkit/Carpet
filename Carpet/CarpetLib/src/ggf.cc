/***************************************************************************
                          ggf.cc  -  Generic Grid Function
                          grid function without type information
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/ggf.cc,v 1.11 2001/07/02 13:22:14 schnetter Exp $

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
#include <string>

#include "defs.hh"
#include "dh.hh"
#include "th.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GGF_HH)
#  include "ggf.hh"
#endif

using namespace std;



// Constructors
template<int D>
generic_gf<D>::generic_gf (const string name, th& t, dh<D>& d,
                           const int tmin, const int tmax)
  : dimgeneric_gf(name, t, tmin, tmax),
    h(d.h), d(d),
    storage(tmax-tmin+1)
{
  assert (t.h == &d.h);

  d.add(this);

//   recompose();
}

// Destructors
template<int D>
generic_gf<D>::~generic_gf () {
  d.remove(this);
}

// Comparison
template<int D>
bool generic_gf<D>::operator== (const generic_gf<D>& f) const {
  return this == &f;
}



// Modifiers
template<int D>
void generic_gf<D>::recompose () {
  // Retain storage that might be needed
  fdata oldstorage(tmax-tmin+1);
  for (int tl=tmin; tl<=tmax; ++tl) {
    oldstorage[tl-tmin].resize
      (min(h.reflevels(), (int)storage[tl-tmin].size()));
    for (int rl=0; rl<(int)storage[tl-tmin].size(); ++rl) {
      oldstorage[tl-tmin][rl].resize
        (min(h.components(rl), (int)storage[tl-tmin][rl].size()));
      for (int c=0; c<(int)storage[tl-tmin][rl].size(); ++c) {
        oldstorage[tl-tmin][rl][c].resize
	  (min(h.mglevels(rl,c), (int)storage[tl-tmin][rl][c].size()));
        for (int ml=0; ml<(int)storage[tl-tmin][rl][ml].size(); ++ml) {
          oldstorage[tl-tmin][rl][c][ml]->transfer_from
            (storage[tl-tmin][tl][c][ml]);
        } // for ml
      } // for c
    } // for rl
  } // for tl
  
  // Delete storage
  storage.clear();
  
  storage.resize(tmax-tmin+1);
  for (int tl=tmin; tl<=tmax; ++tl) {
    storage[tl-tmin].resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      storage[tl-tmin][rl].resize(h.components(rl));
      for (int c=0; c<h.components(rl); ++c) {
      	storage[tl-tmin][rl][c].resize(h.mglevels(rl,c));
      	for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	  
	  storage[tl-tmin][rl][c][ml] = typed_data();
	  
      	  // Allocate storage
      	  storage[tl-tmin][rl][c][ml]->allocate
      	    (d.boxes[rl][c][ml].exterior, h.proc(rl,c));
	  
      	  // Initialise from coarser level, if possible
      	  // TODO: init only un-copied regions
      	  if (rl>0) {
      	    ref_prolongate (tl,rl,c,ml);
      	  } // if rl
	  
      	  // Copy from old storage, if possible
      	  if (rl<(int)oldstorage[tl-tmin].size()) {
      	    for (int cc=0; cc<(int)oldstorage[tl-tmin][rl].size(); ++cc) {
      	      if (ml<(int)oldstorage[tl-tmin][rl][cc].size()) {
      		const ibbox ovlp =
      		  (d.boxes[rl][c][ml].exterior
      		   & oldstorage[tl-tmin][rl][cc][ml]->extent());
      		storage[tl-tmin][rl][c][ml]->copy_from
      		  (oldstorage[tl-tmin][rl][cc][ml], ovlp);
	      } // if ml
	    } // for cc
	  } // if rl
	  
	} // for ml
      } // for c
      
      // Delete old storage
      if (rl<(int)oldstorage[tl-tmin].size()) {
      	oldstorage[tl-tmin][rl].clear();
      }
      
    } // for rl
  } // for tl
  
  for (int tl=tmin; tl<=tmax; ++tl) {
    for (int rl=0; rl<h.reflevels(); ++rl) {
      
      // Set boundaries
      for (int c=0; c<h.components(rl); ++c) {
      	for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	  sync (tl,rl,c,ml);
      	  // TODO: assert that reflevel 0 boundaries are copied
      	  if (rl>0) {
      	    ref_bnd_prolongate (tl,rl,c,ml);
      	  } // if rl
      	} // for ml
      } // for c
      
    } // for rl
  } // for tl
}

// Cycle the time levels by rotating the data sets
template<int D>
void generic_gf<D>::cycle (int rl, int c, int ml) {
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  generic_data<D>* tmpdata = storage[tmin-tmin][rl][c][ml];
  for (int tl=tmin; tl<=tmax-1; ++tl) {
    storage[tl-tmin][rl][c][ml] = storage[tl+1-tmin][rl][c][ml];
  }
  storage[tmax-tmin][rl][c][ml] = tmpdata;
}



// Operations

// Copy a region
template<int D>
void generic_gf<D>::copycat (int tl1, int rl1, int c1, int ml1,
			     const ibbox dh<D>::dboxes::* recv_list,
			     int tl2, int rl2, int ml2,
			     const ibbox dh<D>::dboxes::* send_list)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (rl2>=0 && rl2<h.reflevels());
  const int c2=c1;
  assert (ml2<h.mglevels(rl2,c2));
  const ibbox recv = d.boxes[rl1][c1][ml1].*recv_list;
  const ibbox send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // copy the content
  assert (recv==send);
  storage[tl1-tmin][rl1][c1][ml1]->copy_from
    (storage[tl2-tmin][rl2][c2][ml2], recv);
}

// Copy regions
template<int D>
void generic_gf<D>::copycat (int tl1, int rl1, int c1, int ml1,
			     const iblist dh<D>::dboxes::* recv_list,
			     int tl2, int rl2, int ml2,
			     const iblist dh<D>::dboxes::* send_list)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (rl2>=0 && rl2<h.reflevels());
  const int c2=c1;
  assert (ml2<h.mglevels(rl2,c2));
  const iblist recv = d.boxes[rl1][c1][ml1].*recv_list;
  const iblist send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // walk all boxes
  for (iblist::const_iterator r=recv.begin(), s=send.begin();
       r!=recv.end(); ++r, ++s) {
    // (use the send boxes for communication)
    // copy the content
    storage[tl1-tmin][rl1][c1][ml1]->copy_from
      (storage[tl2-tmin][rl2][c2][ml2], *r);
  }
}

// Copy regions
template<int D>
void generic_gf<D>::copycat (int tl1, int rl1, int c1, int ml1,
			     const iblistvect dh<D>::dboxes::* recv_listvect,
			     int tl2, int rl2, int ml2,
			     const iblistvect dh<D>::dboxes::* send_listvect)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (rl2>=0 && rl2<h.reflevels());
  // walk all components
  for (int c2=0; c2<h.components(rl2); ++c2) {
    assert (ml2<h.mglevels(rl2,c2));
    const iblist recv = (d.boxes[rl1][c1][ml1].*recv_listvect)[c2];
    const iblist send = (d.boxes[rl2][c2][ml2].*send_listvect)[c1];
    assert (recv.size()==send.size());
    // walk all boxes
    for (iblist::const_iterator r=recv.begin(), s=send.begin();
      	 r!=recv.end(); ++r, ++s) {
      // (use the send boxes for communication)
      // copy the content
      storage[tl1-tmin][rl1][c1][ml1]->copy_from
      	(storage[tl2-tmin][rl2][c2][ml2], *r);
    }
  }
}

// Interpolate a region
template<int D>
void generic_gf<D>::intercat (int tl1, int rl1, int c1, int ml1,
			      const ibbox dh<D>::dboxes::* recv_list,
			      const vector<int> tl2s, int rl2, int ml2,
			      const ibbox dh<D>::dboxes::* send_list,
			      int order_space)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s[i]>=tmin && tl2s[i]<=tmax);
  }
  assert (rl2>=0 && rl2<h.reflevels());
  const int c2=c1;
  assert (ml2<h.mglevels(rl2,c2));
  
  vector<const generic_data<D>*> gsrcs(tl2s.size());
  vector<int> tls(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    gsrcs[i] = storage[tl2s[i]-tmin][rl2][c2][ml2];
    tls[i] = t.time(tl2s[i],rl2,ml2);
  }
  const int tl = t.time(tl1,rl1,ml1);
  
  const ibbox recv = d.boxes[rl1][c1][ml1].*recv_list;
  const ibbox send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // interpolate the content
  assert (recv==send);
  storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
    (gsrcs, tls, recv, tl, order_space);
}

// Interpolate regions
template<int D>
void generic_gf<D>::intercat (int tl1, int rl1, int c1, int ml1,
			      const iblist dh<D>::dboxes::* recv_list,
			      const vector<int> tl2s, int rl2, int ml2,
			      const iblist dh<D>::dboxes::* send_list,
			      int order_space)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s[i]>=tmin && tl2s[i]<=tmax);
  }
  assert (rl2>=0 && rl2<h.reflevels());
  const int c2=c1;
  assert (ml2<h.mglevels(rl2,c2));
  
  vector<const generic_data<D>*> gsrcs(tl2s.size());
  vector<int> tls(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    gsrcs[i] = storage[tl2s[i]-tmin][rl2][c2][ml2];
    tls[i] = t.time(tl2s[i],rl2,ml2);
  }
  const int tl = t.time(tl1,rl1,ml1);
  
  const iblist recv = d.boxes[rl1][c1][ml1].*recv_list;
  const iblist send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // walk all boxes
  for (iblist::const_iterator r=recv.begin(), s=send.begin();
       r!=recv.end(); ++r, ++s) {
    // (use the send boxes for communication)
    // interpolate the content
    storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      (gsrcs, tls, *r, tl, order_space);
  }
}

// Interpolate regions
template<int D>
void generic_gf<D>::intercat (int tl1, int rl1, int c1, int ml1,
			      const iblistvect dh<D>::dboxes::* recv_listvect,
			      const vector<int> tl2s, int rl2, int ml2,
			      const iblistvect dh<D>::dboxes::* send_listvect,
			      int order_space)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s[i]>=tmin && tl2s[i]<=tmax);
  }
  assert (rl2>=0 && rl2<h.reflevels());
  // walk all components
  for (int c2=0; c2<h.components(rl2); ++c2) {
    assert (ml2<h.mglevels(rl2,c2));
    
    vector<const generic_data<D>*> gsrcs(tl2s.size());
    vector<int> tls(tl2s.size());
    for (int i=0; i<(int)gsrcs.size(); ++i) {
      gsrcs[i] = storage[tl2s[i]-tmin][rl2][c2][ml2];
      tls[i] = t.time(tl2s[i],rl2,ml2);
    }
    const int tl = t.time(tl1,rl1,ml1);
    
    const iblist recv = (d.boxes[rl1][c1][ml1].*recv_listvect)[c2];
    const iblist send = (d.boxes[rl2][c2][ml2].*send_listvect)[c1];
    assert (recv.size()==send.size());
    // walk all boxes
    for (iblist::const_iterator r=recv.begin(), s=send.begin();
      	 r!=recv.end(); ++r, ++s) {
      // (use the send boxes for communication)
      // interpolate the content
      storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      	(gsrcs, tls, *r, tl, order_space);
    }
  }
}



// Copy a component from the next time level
template<int D>
void generic_gf<D>::copy (int tl, int rl, int c, int ml) {
  // Copy
  copycat (tl  ,rl,c,ml, &dh<D>::dboxes::exterior,
      	   tl+1,rl,  ml, &dh<D>::dboxes::exterior);
}

// Synchronise the boundaries a component
template<int D>
void generic_gf<D>::sync (int tl, int rl, int c, int ml) {
  // Copy
  copycat (tl,rl,c,ml, &dh<D>::dboxes::recv_sync,
      	   tl,rl,  ml, &dh<D>::dboxes::send_sync);
}

// Prolongate the boundaries of a component
template<int D>
void generic_gf<D>::ref_bnd_prolongate (int tl, int rl, int c, int ml,
					int order_space, int order_time) {
  // Interpolate
  assert (rl>=1);
  const int tmod
    = ((t.time(tl,rl,ml) - t.get_time(rl-1,ml)) % t.get_delta(rl-1, ml)
       + t.get_delta(rl-1, ml)) % t.get_delta(rl-1, ml);
  vector<int> tl2s;
  if (tmod == 0) {
    // No interpolation in time
    tl2s.resize(1);
    tl2s[0] = tl;
  } else {
    // Interpolation in time
    assert (tmax-tmin+1 >= order_time+1);
    tl2s.resize(order_time+1);
    for (int i=0; i<=order_time; ++i) tl2s[i] = tmax - i;
  }
  intercat (tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_bnd_coarse,
	    tl2s,rl-1,  ml, &dh<D>::dboxes::send_ref_bnd_fine,
	    order_space);
}

// Restrict a multigrid level
template<int D>
void generic_gf<D>::mg_restrict (int tl, int rl, int c, int ml,
				 int order_space) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl,ml-1));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml-1, &dh<D>::dboxes::send_mg_fine,
	    order_space);
}

// Prolongate a multigrid level
template<int D>
void generic_gf<D>::mg_prolongate (int tl, int rl, int c, int ml,
				   int order_space) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl,ml+1));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml+1, &dh<D>::dboxes::send_mg_fine,
	    order_space);
}

// Restrict a refinement level
template<int D>
void generic_gf<D>::ref_restrict (int tl, int rl, int c, int ml,
				  int order_space) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl+1,ml));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_fine,
	    tl2s,rl+1,  ml, &dh<D>::dboxes::send_ref_coarse,
	    order_space);
}

// Prolongate a refinement level
template<int D>
void generic_gf<D>::ref_prolongate (int tl, int rl, int c, int ml,
				    int order_space) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl-1,ml));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_coarse,
	    tl2s,rl-1,  ml, &dh<D>::dboxes::send_ref_fine,
	    order_space);
}



#if defined(TMPL_EXPLICIT)
template class generic_gf<1>;
template class generic_gf<2>;
template class generic_gf<3>;
#endif
