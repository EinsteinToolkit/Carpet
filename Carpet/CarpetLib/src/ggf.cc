// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/ggf.cc,v 1.22 2003/03/26 17:34:43 schnetter Exp $

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>

#include "cctk.h"

#include "defs.hh"
#include "dh.hh"
#include "th.hh"

#include "ggf.hh"

using namespace std;



// Constructors
template<int D>
ggf<D>::ggf (const string name, th<D>& t, dh<D>& d,
             const int tmin, const int tmax,
             const int prolongation_order_time)
  : name(name), t(t),
    tmin(tmin), tmax(tmax),
    prolongation_order_time(prolongation_order_time),
    h(d.h), d(d),
    storage(tmax-tmin+1)
{
  assert (t.h == &d.h);

  d.add(this);
}

// Destructors
template<int D>
ggf<D>::~ggf () {
  d.remove(this);
}

// Comparison
template<int D>
bool ggf<D>::operator== (const ggf<D>& f) const {
  return this == &f;
}



// Modifiers
template<int D>
void ggf<D>::recompose () {
  
  // Retain storage that might be needed
  fdata oldstorage = storage;
  
  // Resize structure
  storage.resize(tmax-tmin+1);
  for (int tl=tmin; tl<=tmax; ++tl) {
    storage[tl-tmin].resize(h.reflevels());
    for (int rl=0; rl<h.reflevels(); ++rl) {
      storage[tl-tmin][rl].resize(h.components(rl));
      for (int c=0; c<h.components(rl); ++c) {
      	storage[tl-tmin][rl][c].resize(h.mglevels(rl,c));
	for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	  storage[tl-tmin][rl][c][ml] = 0;
	} // for ml
      } // for c
    } // for rl
  } // for tl
  
  // Initialise the new storage
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int tl=tmin; tl<=tmax; ++tl) {
      for (int c=0; c<h.components(rl); ++c) {
      	for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	  
	  storage[tl-tmin][rl][c][ml] = typed_data();
	  
      	  // Allocate storage
      	  storage[tl-tmin][rl][c][ml]->allocate
      	    (d.boxes[rl][c][ml].exterior, h.proc(rl,c));
	  
      	  // Initialise from coarser level, if possible
#warning "TODO: init only un-copied regions"
      	  if (rl>0) {
	    const CCTK_REAL time = t.time(tl,rl,ml);
      	    ref_prolongate (tl,rl,c,ml,time);
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
      
    } // for tl
  } // for rl
  
  // Delete old storage
  for (int tl=tmin; tl<=tmax; ++tl) {
    for (int rl=0; rl<(int)oldstorage[tl-tmin].size(); ++rl) {
      for (int c=0; c<(int)oldstorage[tl-tmin][rl].size(); ++c) {
        for (int ml=0; ml<(int)oldstorage[tl-tmin][rl][c].size(); ++ml) {
	  delete oldstorage[tl-tmin][rl][c][ml];
        } // for ml
      } // for c
    } // for rl
  } // for tl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int tl=tmin; tl<=tmax; ++tl) {
      
      // Set boundaries
      for (int c=0; c<h.components(rl); ++c) {
      	for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	  sync (tl,rl,c,ml);
      	  // TODO: assert that reflevel 0 boundaries are copied
      	  if (rl>0) {
	    const CCTK_REAL time = t.time(tl,rl,ml);
      	    ref_bnd_prolongate (tl,rl,c,ml,time);
      	  } // if rl
      	} // for ml
      } // for c
      
    } // for tl
  } // for rl
}

// Cycle the time levels by rotating the data sets
template<int D>
void ggf<D>::cycle (int rl, int c, int ml) {
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  gdata<D>* tmpdata = storage[tmin-tmin][rl][c][ml];
  for (int tl=tmin; tl<=tmax-1; ++tl) {
    storage[tl-tmin][rl][c][ml] = storage[tl+1-tmin][rl][c][ml];
  }
  storage[tmax-tmin][rl][c][ml] = tmpdata;
}

// Flip the time levels by exchanging the data sets
template<int D>
void ggf<D>::flip (int rl, int c, int ml) {
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  for (int t=0; t<(tmax-tmin)/2; ++t) {
    const int tl1 = tmin + t;
    const int tl2 = tmax - t;
    assert (tl1 < tl2);
    gdata<D>* tmpdata = storage[tl1-tmin][rl][c][ml];
    storage[tl1-tmin][rl][c][ml] = storage[tl2-tmin][rl][c][ml];
    storage[tl2-tmin][rl][c][ml] = tmpdata;
  }
}

// Copy data from current time level to all previous levels
template<int D>
void ggf<D>::copytoprevs (int rl, int c, int ml) {
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  for (int tl=tmin; tl<=tmax-1; ++tl) {
    storage[tl-tmin][rl][c][ml]->copy_from 
      (storage[tmax-tmin][rl][c][ml], storage[tmax-tmin][rl][c][ml]->extent());
  }
}





// Operations

// Copy a region
template<int D>
void ggf<D>::copycat (int tl1, int rl1, int c1, int ml1,
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
void ggf<D>::copycat (int tl1, int rl1, int c1, int ml1,
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
  for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
       r!=recv.end(); ++r, ++s) {
    // (use the send boxes for communication)
    // copy the content
    storage[tl1-tmin][rl1][c1][ml1]->copy_from
      (storage[tl2-tmin][rl2][c2][ml2], *r);
  }
}

// Copy regions
template<int D>
void ggf<D>::copycat (int tl1, int rl1, int c1, int ml1,
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
    for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
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
void ggf<D>::intercat (int tl1, int rl1, int c1, int ml1,
                       const ibbox dh<D>::dboxes::* recv_list,
                       const vector<int> tl2s, int rl2, int ml2,
                       const ibbox dh<D>::dboxes::* send_list,
                       CCTK_REAL time)
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
  assert (ml2>=0 && ml2<h.mglevels(rl2,c2));
  
  vector<const gdata<D>*> gsrcs(tl2s.size());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    assert (rl2<(int)storage[tl2s[i]-tmin].size());
    assert (c2<(int)storage[tl2s[i]-tmin][rl2].size());
    assert (ml2<(int)storage[tl2s[i]-tmin][rl2][ml2].size());
    gsrcs[i] = storage[tl2s[i]-tmin][rl2][c2][ml2];
    times[i] = t.time(tl2s[i],rl2,ml2);
  }
//   const CCTK_REAL time = t.time(tl1,rl1,ml1);
  
  const ibbox recv = d.boxes[rl1][c1][ml1].*recv_list;
  const ibbox send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // interpolate the content
  assert (recv==send);
  storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
    (gsrcs, times, recv, time,
     d.prolongation_order_space, prolongation_order_time);
}

// Interpolate regions
template<int D>
void ggf<D>::intercat (int tl1, int rl1, int c1, int ml1,
                       const iblist dh<D>::dboxes::* recv_list,
                       const vector<int> tl2s, int rl2, int ml2,
                       const iblist dh<D>::dboxes::* send_list,
                       const CCTK_REAL time)
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
  assert (ml2>=0 && ml2<h.mglevels(rl2,c2));
  
  vector<const gdata<D>*> gsrcs(tl2s.size());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    assert (rl2<(int)storage[tl2s[i]-tmin].size());
    assert (c2<(int)storage[tl2s[i]-tmin][rl2].size());
    assert (ml2<(int)storage[tl2s[i]-tmin][rl2][ml2].size());
    gsrcs[i] = storage[tl2s[i]-tmin][rl2][c2][ml2];
    times[i] = t.time(tl2s[i],rl2,ml2);
  }
//   const CCTK_REAL time = t.time(tl1,rl1,ml1);
  
  const iblist recv = d.boxes[rl1][c1][ml1].*recv_list;
  const iblist send = d.boxes[rl2][c2][ml2].*send_list;
  assert (recv.size()==send.size());
  // walk all boxes
  for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
       r!=recv.end(); ++r, ++s) {
    // (use the send boxes for communication)
    // interpolate the content
    storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      (gsrcs, times, *r, time,
       d.prolongation_order_space, prolongation_order_time);
  }
}

// Interpolate regions
template<int D>
void ggf<D>::intercat (int tl1, int rl1, int c1, int ml1,
                       const iblistvect dh<D>::dboxes::* recv_listvect,
                       const vector<int> tl2s, int rl2, int ml2,
                       const iblistvect dh<D>::dboxes::* send_listvect,
                       const CCTK_REAL time)
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
    assert (ml2>=0 && ml2<h.mglevels(rl2,c2));
    
    vector<const gdata<D>*> gsrcs(tl2s.size());
    vector<CCTK_REAL> times(tl2s.size());
    for (int i=0; i<(int)gsrcs.size(); ++i) {
      assert (rl2<(int)storage[tl2s[i]-tmin].size());
      assert (c2<(int)storage[tl2s[i]-tmin][rl2].size());
      assert (ml2<(int)storage[tl2s[i]-tmin][rl2][ml2].size());
      gsrcs[i] = storage[tl2s[i]-tmin][rl2][c2][ml2];
      times[i] = t.time(tl2s[i],rl2,ml2);
    }
//     const CCTK_REAL time = t.time(tl1,rl1,ml1);
    
    const iblist recv = (d.boxes[rl1][c1][ml1].*recv_listvect)[c2];
    const iblist send = (d.boxes[rl2][c2][ml2].*send_listvect)[c1];
    assert (recv.size()==send.size());
    // walk all boxes
    for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
      	 r!=recv.end(); ++r, ++s) {
      // (use the send boxes for communication)
      // interpolate the content
      storage[tl1-tmin][rl1][c1][ml1]->interpolate_from
      	(gsrcs, times, *r, time,
	 d.prolongation_order_space, prolongation_order_time);
    }
  }
}



// Copy a component from the next time level
template<int D>
void ggf<D>::copy (int tl, int rl, int c, int ml)
{
  // Copy
  copycat (tl  ,rl,c,ml, &dh<D>::dboxes::exterior,
      	   tl+1,rl,  ml, &dh<D>::dboxes::exterior);
}

// Synchronise the boundaries a component
template<int D>
void ggf<D>::sync (int tl, int rl, int c, int ml)
{
  // Copy
  copycat (tl,rl,c,ml, &dh<D>::dboxes::recv_sync,
      	   tl,rl,  ml, &dh<D>::dboxes::send_sync);
}

// Prolongate the boundaries of a component
template<int D>
void ggf<D>::ref_bnd_prolongate (int tl, int rl, int c, int ml,
                                 CCTK_REAL time)
{
  // Interpolate
  assert (rl>=1);
  vector<int> tl2s;
  // Interpolation in time
  assert (tmax-tmin+1 >= prolongation_order_time+1);
  tl2s.resize(prolongation_order_time+1);
  for (int i=0; i<=prolongation_order_time; ++i) tl2s[i] = tmax - i;
  intercat (tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_bnd_coarse,
	    tl2s,rl-1,  ml, &dh<D>::dboxes::send_ref_bnd_fine,
	    time);
}

// Restrict a multigrid level
template<int D>
void ggf<D>::mg_restrict (int tl, int rl, int c, int ml,
                          CCTK_REAL time)
{
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl,ml-1));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml-1, &dh<D>::dboxes::send_mg_fine,
	    time);
}

// Prolongate a multigrid level
template<int D>
void ggf<D>::mg_prolongate (int tl, int rl, int c, int ml,
                            CCTK_REAL time)
{
  // Require same times
  assert (t.get_time(rl,ml) == t.get_time(rl,ml+1));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml+1, &dh<D>::dboxes::send_mg_fine,
	    time);
}

// Restrict a refinement level
template<int D>
void ggf<D>::ref_restrict (int tl, int rl, int c, int ml,
                           CCTK_REAL time)
{
  // Require same times
  // SHH: removed assert and added warning
//   if (t.get_time(rl,ml) != t.get_time(rl+1,ml)) {
//     printf("WARNING: Time on rl %d is %d, time on rl %d is %d\n",rl,t.get_time(rl,ml),rl+1,t.get_time(rl+1,ml));
//   }
  //assert (t.get_time(rl,ml) == t.get_time(rl+1,ml));
  const vector<int> tl2s(1,tl);
  intercat (tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_fine,
	    tl2s,rl+1,  ml, &dh<D>::dboxes::send_ref_coarse,
	    time);
}

// Prolongate a refinement level
template<int D>
void ggf<D>::ref_prolongate (int tl, int rl, int c, int ml,
                             CCTK_REAL time)
{
  assert (rl>=1);
  vector<int> tl2s;
  // Interpolation in time
  assert (tmax-tmin+1 >= prolongation_order_time+1);
  tl2s.resize(prolongation_order_time+1);
  for (int i=0; i<=prolongation_order_time; ++i) tl2s[i] = tmax - i;
  intercat (tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_coarse,
	    tl2s,rl-1,  ml, &dh<D>::dboxes::send_ref_fine,
	    time);
}



template class ggf<3>;
