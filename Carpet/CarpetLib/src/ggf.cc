/***************************************************************************
                          ggf.cc  -  Generic Grid Function
                          grid function without type information
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/ggf.cc,v 1.5 2001/03/19 21:30:19 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <cassert>
#include <iostream>
#include <string>

#include "defs.hh"
#include "dh.hh"
#include "th.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GGF_HH)
#  include "ggf.hh"
#endif



// Constructors
template<int D>
generic_gf<D>::generic_gf (const string name, th<D>& t, dh<D>& d,
                           const int tmin, const int tmax)
  : name(name), h(d.h), t(t), d(d), tmin(tmin), tmax(tmax),
    storage(tmax-tmin+1)
{
  assert (&t.h == &d.h);
  assert (tmin<=tmax+1);

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



// Operations

// Copy region for a component (between time levels)
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

// Copy regions for a component (between multigrid levels)
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
    storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      (storage[tl2-tmin][rl2][c2][ml2], *r);
  }
}

// Copy regions for a level (between refinement levels)
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
      storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      	(storage[tl2-tmin][rl2][c2][ml2], *r);
    }
  }
}

// Interpolate a component (between time levels)
template<int D>
void generic_gf<D>::intercat (int tl1, int rl1, int c1, int ml1,
			      const ibbox dh<D>::dboxes::* recv_list,
			      int tl2, const double fact2,
			      int tl3, const double fact3,
			      int rl2, int ml2,
			      const ibbox dh<D>::dboxes::* send_list)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (tl3>=tmin && tl3<=tmax);
  assert (rl2>=0 && rl2<h.reflevels());
  const int c2=c1;
  assert (ml2<h.mglevels(rl2,c2));
  const ibbox recv = d.boxes[rl1][c1][ml1].*recv_list;
  const ibbox send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // interpolate the content
  assert (recv==send);
  storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
    (storage[tl2-tmin][rl2][c2][ml2], fact2,
     storage[tl3-tmin][rl2][c2][ml2], fact3, recv);
}

// Interpolate a component (between multigrid levels)
template<int D>
void generic_gf<D>::intercat (int tl1, int rl1, int c1, int ml1,
			      const iblist dh<D>::dboxes::* recv_list,
			      int tl2, const double fact2,
			      int tl3, const double fact3,
			      int rl2, int ml2,
			      const iblist dh<D>::dboxes::* send_list)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (tl3>=tmin && tl3<=tmax);
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
    // interpolate the content
    storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      (storage[tl2-tmin][rl2][c2][ml2], fact2,
       storage[tl3-tmin][rl2][c2][ml2], fact3, *r);
  }
}

// Interpolate a level (between refinement levels)
template<int D>
void generic_gf<D>::intercat (int tl1, int rl1, int c1, int ml1,
			      const iblistvect dh<D>::dboxes::* recv_listvect,
			      int tl2, const double fact2,
			      int tl3, const double fact3,
			      int rl2, int ml2,
			      const iblistvect dh<D>::dboxes::* send_listvect)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (rl2>=0 && rl2<h.reflevels());
  assert (tl3>=tmin && tl3<=tmax);
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
      // interpolate the content
      storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      	(storage[tl2-tmin][rl2][c2][ml2], fact2,
      	 storage[tl3-tmin][rl2][c2][ml2], fact3, *r);
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
void generic_gf<D>::ref_bnd_prolongate (int tl, int rl, int c, int ml) {
  // Interpolate
  assert (rl>=1);
  double time =
    (t.time(tl,rl,ml) - t.get_time(rl-1,ml)) / (double)t.get_delta(rl-1, ml);
  const int tl2 = (int)floor(time + 1e-6);
  assert (tl2>=tmin && tl2<=tmax);
  if (time==tl2) {
    copycat (tl ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_bnd_coarse,
      	     tl2,rl-1,  ml, &dh<D>::dboxes::send_ref_bnd_fine);
  } else {
    int tl3=tl2+1;
    assert (tl3>=tmin && tl3<=tmax);
    const double fact2 = 1 - (time - tl2);
    const double fact3 = 1 - fact2;
    intercat (tl,rl,c,ml, &dh<D>::dboxes::recv_ref_bnd_coarse,
      	      tl2,fact2, tl3,fact3,
      	      rl-1,ml, &dh<D>::dboxes::send_ref_bnd_fine);
  }
}

// Restrict a multigrid level
template<int D>
void generic_gf<D>::mg_restrict (int tl, int rl, int c, int ml) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl,ml-1));
  copycat (tl,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
      	   tl,rl,  ml-1, &dh<D>::dboxes::send_mg_fine);
}

// Prolongate a multigrid level
template<int D>
void generic_gf<D>::mg_prolongate (int tl, int rl, int c, int ml) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl,ml+1));
  copycat (tl,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
      	   tl,rl,  ml+1, &dh<D>::dboxes::send_mg_fine);
}

// Restrict a refinement level
template<int D>
void generic_gf<D>::ref_restrict (int tl, int rl, int c, int ml) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl+1,ml));
  copycat (tl,rl  ,c,ml, &dh<D>::dboxes::recv_ref_fine,
      	   tl,rl+1,  ml, &dh<D>::dboxes::send_ref_coarse);
}

// Prolongate a refinement level
template<int D>
void generic_gf<D>::ref_prolongate (int tl, int rl, int c, int ml) {
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl-1,ml));
  copycat (tl,rl  ,c,ml, &dh<D>::dboxes::recv_ref_coarse,
      	   tl,rl-1,  ml, &dh<D>::dboxes::send_ref_fine);
}



// Output
template<int D>
ostream& operator<< (ostream& os, const generic_gf<D>& f) {
  return f.out(os);
}



#if defined(TMPL_EXPLICIT)
template class generic_gf<3>;
template ostream& operator<< (ostream& os, const generic_gf<3>& f);
#endif
