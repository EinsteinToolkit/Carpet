#include <cassert>
#include <cmath>
#include <cstdlib>
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
ggf<D>::ggf (const int varindex, const operator_type transport_operator,
             th<D>& t, dh<D>& d,
             const int tmin, const int tmax,
             const int prolongation_order_time,
             const int vectorlength, const int vectorindex,
             ggf* const vectorleader)
  : varindex(varindex), transport_operator(transport_operator), t(t),
    tmin(tmin), tmax(tmax),
    prolongation_order_time(prolongation_order_time),
    h(d.h), d(d),
    storage(tmax-tmin+1),
    vectorlength(vectorlength), vectorindex(vectorindex),
    vectorleader(vectorleader)
{
  assert (&t.h == &d.h);
  
  assert (vectorlength >= 1);
  assert (vectorindex >= 0 && vectorindex < vectorlength);
  assert ((vectorindex==0 && !vectorleader)
          || (vectorindex!=0 && vectorleader));
  
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
void ggf<D>::recompose_crop ()
{
  // Free storage that will not be needed
  storage.resize(tmax-tmin+1);
  for (int tl=tmin; tl<=tmax; ++tl) {
    for (int rl=h.reflevels(); rl<(int)storage.at(tl-tmin).size(); ++rl) {
      for (int c=0; c<(int)storage.at(tl-tmin).at(rl).size(); ++c) {
        for (int ml=0; ml<(int)storage.at(tl-tmin).at(rl).at(c).size(); ++ml) {
          delete storage.at(tl-tmin).at(rl).at(c).at(ml);
        } // for ml
      } // for c
    } // for rl
    storage.at(tl-tmin).resize(h.reflevels());
  } // for tl
}

template<int D>
void ggf<D>::recompose_allocate (const int rl)
{
  // TODO: restructure storage only when needed
  
  // Retain storage that might be needed
  oldstorage.resize(tmax-tmin+1);
  for (int tl=tmin; tl<=tmax; ++tl) {
    oldstorage.at(tl-tmin).resize(h.reflevels());
    assert (oldstorage.at(tl-tmin).at(rl).size() == 0);
    oldstorage.at(tl-tmin).at(rl) = storage.at(tl-tmin).at(rl);
    storage.at(tl-tmin).at(rl).resize(0);
  }
  
  // Resize structure and allocate storage
  storage.resize(tmax-tmin+1);
  for (int tl=tmin; tl<=tmax; ++tl) {
    storage.at(tl-tmin).resize(h.reflevels());
    storage.at(tl-tmin).at(rl).resize(h.components(rl));
    for (int c=0; c<h.components(rl); ++c) {
      storage.at(tl-tmin).at(rl).at(c).resize(h.mglevels(rl,c));
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
        storage.at(tl-tmin).at(rl).at(c).at(ml) = typed_data(tl,rl,c,ml);
        storage.at(tl-tmin).at(rl).at(c).at(ml)->allocate
          (d.boxes.at(rl).at(c).at(ml).exterior, h.proc(rl,c));
      } // for ml
    } // for c
  } // for tl
}

template<int D>
void ggf<D>::recompose_fill (comm_state<D>& state, const int rl,
                             const bool do_prolongate)
{
  // Initialise the new storage
  for (int c=0; c<h.components(rl); ++c) {
    for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      for (int tl=tmin; tl<=tmax; ++tl) {
        
        // Find out which regions need to be prolongated
        // (Copy the exterior because some variables are not prolongated)
        // TODO: do this once in the dh instead of for each variable here
        ibset work (d.boxes.at(rl).at(c).at(ml).exterior);
        
        // Copy from old storage, if possible
        // TODO: copy only from interior regions?
        if (rl<(int)oldstorage.at(tl-tmin).size()) {
          for (int cc=0; cc<(int)oldstorage.at(tl-tmin).at(rl).size(); ++cc) {
            if (ml<(int)oldstorage.at(tl-tmin).at(rl).at(cc).size()) {
              // TODO: prefer same processor, etc., see dh.cc
              ibset ovlp
                = work & oldstorage.at(tl-tmin).at(rl).at(cc).at(ml)->extent();
              ovlp.normalize();
              work -= ovlp;
              for (typename ibset::const_iterator r=ovlp.begin(); r!=ovlp.end(); ++r) {
                storage.at(tl-tmin).at(rl).at(c).at(ml)->copy_from
                  (state, oldstorage.at(tl-tmin).at(rl).at(cc).at(ml), *r);
              }
            } // if ml
          } // for cc
        } // if rl
        
        if (do_prolongate) {
          // Initialise from coarser level, if possible
          if (rl>0) {
            if (transport_operator != op_none) {
              const int numtl = prolongation_order_time+1;
              assert (tmax-tmin+1 >= numtl);
              vector<int> tls(numtl);
              vector<CCTK_REAL> times(numtl);
              for (int i=0; i<numtl; ++i) {
                tls.at(i) = tmax - i;
                times.at(i) = t.time(tls.at(i),rl-1,ml);
              }
              for (int cc=0; cc<(int)storage.at(tl-tmin).at(rl-1).size(); ++cc) {
                vector<const gdata<D>*> gsrcs(numtl);
                for (int i=0; i<numtl; ++i) {
                  gsrcs.at(i)
                    = storage.at(tls.at(i)-tmin).at(rl-1).at(cc).at(ml);
                  assert (gsrcs.at(i)->extent() == gsrcs.at(0)->extent());
                }
                const CCTK_REAL time = t.time(tl,rl,ml);
                
                // TODO: choose larger regions first
                // TODO: prefer regions from the same processor
                const iblist& list
                  = d.boxes.at(rl).at(c).at(ml).recv_ref_coarse.at(cc);
                for (typename iblist::const_iterator iter=list.begin(); iter!=list.end(); ++iter) {
                  ibset ovlp = work & *iter;
                  ovlp.normalize();
                  work -= ovlp;
                  for (typename ibset::const_iterator r=ovlp.begin(); r!=ovlp.end(); ++r) {
                    storage.at(tl-tmin).at(rl).at(c).at(ml)->interpolate_from
                      (state, gsrcs, times, *r, time,
                       d.prolongation_order_space, prolongation_order_time);
                  } // for r
                } // for iter
              } // for cc
            } // if transport_operator
          } // if rl
        } // if do_prolongate
        
        // Note that work need not be empty here; in this case, not
        // everything could be initialised.  This is okay on outer
        // boundaries.
        // TODO: check this.
        
      } // for tl
    } // for ml
  } // for c
}

template<int D>
void ggf<D>::recompose_free (const int rl)
{
  // Delete old storage
  for (int tl=tmin; tl<=tmax; ++tl) {
    for (int c=0; c<(int)oldstorage.at(tl-tmin).at(rl).size(); ++c) {
      for (int ml=0; ml<(int)oldstorage.at(tl-tmin).at(rl).at(c).size(); ++ml) {
        delete oldstorage.at(tl-tmin).at(rl).at(c).at(ml);
      } // for ml
    } // for c
    oldstorage.at(tl-tmin).at(rl).resize(0);
  } // for tl
}

template<int D>
void ggf<D>::recompose_bnd_prolongate (comm_state<D>& state, const int rl,
                                       const bool do_prolongate)
{
  if (do_prolongate) {
    // Set boundaries
    if (rl>0) {
      for (int c=0; c<h.components(rl); ++c) {
        for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
          for (int tl=tmin; tl<=tmax; ++tl) {
            
            // TODO: assert that reflevel 0 boundaries are copied
            const CCTK_REAL time = t.time(tl,rl,ml);
            ref_bnd_prolongate (state,tl,rl,c,ml,time);
            
          } // for tl
        } // for ml
      } // for c
    } // if rl
  } // if do_prolongate
}

template<int D>
void ggf<D>::recompose_sync (comm_state<D>& state, const int rl,
                             const bool do_prolongate)
{
  if (do_prolongate) {
    // Set boundaries
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
        for (int tl=tmin; tl<=tmax; ++tl) {
          
          sync (state,tl,rl,c,ml);
          
        } // for tl
      } // for ml
    } // for c
  } // if do_prolongate
}



// Cycle the time levels by rotating the data sets
template<int D>
void ggf<D>::cycle (int rl, int c, int ml) {
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  gdata<D>* tmpdata = storage.at(tmin-tmin).at(rl).at(c).at(ml);
  for (int tl=tmin; tl<=tmax-1; ++tl) {
    storage.at(tl-tmin).at(rl).at(c).at(ml) = storage.at(tl+1-tmin).at(rl).at(c).at(ml);
  }
  storage.at(tmax-tmin).at(rl).at(c).at(ml) = tmpdata;
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
    gdata<D>* tmpdata = storage.at(tl1-tmin).at(rl).at(c).at(ml);
    storage.at(tl1-tmin).at(rl).at(c).at(ml) = storage.at(tl2-tmin).at(rl).at(c).at(ml);
    storage.at(tl2-tmin).at(rl).at(c).at(ml) = tmpdata;
  }
}



// Operations

// Copy a region
template<int D>
void ggf<D>::copycat (comm_state<D>& state,
                      int tl1, int rl1, int c1, int ml1,
                      const ibbox dh<D>::dboxes::* recv_box,
                      int tl2, int rl2, int ml2,
                      const ibbox dh<D>::dboxes::* send_box)
{
  assert (tl1>=tmin && tl1<=tmax);
  assert (rl1>=0 && rl1<h.reflevels());
  assert (c1>=0 && c1<h.components(rl1));
  assert (ml1>=0 && ml1<h.mglevels(rl1,c1));
  assert (tl2>=tmin && tl2<=tmax);
  assert (rl2>=0 && rl2<h.reflevels());
  const int c2=c1;
  assert (ml2<h.mglevels(rl2,c2));
  const ibbox recv = d.boxes.at(rl1).at(c1).at(ml1).*recv_box;
  const ibbox send = d.boxes.at(rl2).at(c2).at(ml2).*send_box;
  assert (all(recv.shape()==send.shape()));
  // copy the content
  assert (recv==send);
  storage.at(tl1-tmin).at(rl1).at(c1).at(ml1)->copy_from
    (state, storage.at(tl2-tmin).at(rl2).at(c2).at(ml2), recv);
}

// Copy regions
template<int D>
void ggf<D>::copycat (comm_state<D>& state,
                      int tl1, int rl1, int c1, int ml1,
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
  const iblist recv = d.boxes.at(rl1).at(c1).at(ml1).*recv_list;
  const iblist send = d.boxes.at(rl2).at(c2).at(ml2).*send_list;
  assert (recv.size()==send.size());
  // walk all boxes
  for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
       r!=recv.end(); ++r, ++s) {
    // (use the send boxes for communication)
    // copy the content
    storage.at(tl1-tmin).at(rl1).at(c1).at(ml1)->copy_from
      (state, storage.at(tl2-tmin).at(rl2).at(c2).at(ml2), *r);
  }
}

// Copy regions
template<int D>
void ggf<D>::copycat (comm_state<D>& state,
                      int tl1, int rl1, int c1, int ml1,
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
    const iblist recv = (d.boxes.at(rl1).at(c1).at(ml1).*recv_listvect).at(c2);
    const iblist send = (d.boxes.at(rl2).at(c2).at(ml2).*send_listvect).at(c1);
    assert (recv.size()==send.size());
    // walk all boxes
    for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
      	 r!=recv.end(); ++r, ++s) {
      // (use the send boxes for communication)
      // copy the content
      storage.at(tl1-tmin).at(rl1).at(c1).at(ml1)->copy_from
      	(state, storage.at(tl2-tmin).at(rl2).at(c2).at(ml2), *r);
    }
  }
}

// Interpolate a region
template<int D>
void ggf<D>::intercat (comm_state<D>& state,
                       int tl1, int rl1, int c1, int ml1,
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
    assert (tl2s.at(i)>=tmin && tl2s.at(i)<=tmax);
  }
  assert (rl2>=0 and rl2<h.reflevels());
  const int c2=c1;
  assert (ml2>=0 && ml2<h.mglevels(rl2,c2));
  
  vector<const gdata<D>*> gsrcs(tl2s.size());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    assert (rl2<(int)storage.at(tl2s.at(i)-tmin).size());
    assert (c2<(int)storage.at(tl2s.at(i)-tmin).at(rl2).size());
    assert (ml2<(int)storage.at(tl2s.at(i)-tmin).at(rl2).at(c2).size());
    gsrcs.at(i) = storage.at(tl2s.at(i)-tmin).at(rl2).at(c2).at(ml2);
    times.at(i) = t.time(tl2s.at(i),rl2,ml2);
  }
  
  const ibbox recv = d.boxes.at(rl1).at(c1).at(ml1).*recv_list;
  const ibbox send = d.boxes.at(rl2).at(c2).at(ml2).*send_list;
  assert (all(recv.shape()==send.shape()));
  // interpolate the content
  assert (recv==send);
  storage.at(tl1-tmin).at(rl1).at(c1).at(ml1)->interpolate_from
    (state, gsrcs, times, recv, time,
     d.prolongation_order_space, prolongation_order_time);
}

// Interpolate regions
template<int D>
void ggf<D>::intercat (comm_state<D>& state,
                       int tl1, int rl1, int c1, int ml1,
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
    assert (tl2s.at(i)>=tmin && tl2s.at(i)<=tmax);
  }
  assert (rl2>=0 and rl2<h.reflevels());
  const int c2=c1;
  assert (ml2>=0 && ml2<h.mglevels(rl2,c2));
  
  vector<const gdata<D>*> gsrcs(tl2s.size());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    assert (rl2<(int)storage.at(tl2s.at(i)-tmin).size());
    assert (c2<(int)storage.at(tl2s.at(i)-tmin).at(rl2).size());
    assert (ml2<(int)storage.at(tl2s.at(i)-tmin).at(rl2).at(c2).size());
    gsrcs.at(i) = storage.at(tl2s.at(i)-tmin).at(rl2).at(c2).at(ml2);
    times.at(i) = t.time(tl2s.at(i),rl2,ml2);
  }
  
  const iblist recv = d.boxes.at(rl1).at(c1).at(ml1).*recv_list;
  const iblist send = d.boxes.at(rl2).at(c2).at(ml2).*send_list;
  assert (recv.size()==send.size());
  // walk all boxes
  for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
       r!=recv.end(); ++r, ++s) {
    // (use the send boxes for communication)
    // interpolate the content
    storage.at(tl1-tmin).at(rl1).at(c1).at(ml1)->interpolate_from
      (state, gsrcs, times, *r, time,
       d.prolongation_order_space, prolongation_order_time);
  }
}

// Interpolate regions
template<int D>
void ggf<D>::intercat (comm_state<D>& state,
                       int tl1, int rl1, int c1, int ml1,
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
    assert (tl2s.at(i)>=tmin && tl2s.at(i)<=tmax);
  }
  assert (rl2>=0 and rl2<h.reflevels());
  // walk all components
  for (int c2=0; c2<h.components(rl2); ++c2) {
    assert (ml2>=0 && ml2<h.mglevels(rl2,c2));
    
    vector<const gdata<D>*> gsrcs(tl2s.size());
    vector<CCTK_REAL> times(tl2s.size());
    for (int i=0; i<(int)gsrcs.size(); ++i) {
      assert (rl2<(int)storage.at(tl2s.at(i)-tmin).size());
      assert (c2<(int)storage.at(tl2s.at(i)-tmin).at(rl2).size());
      assert (ml2<(int)storage.at(tl2s.at(i)-tmin).at(rl2).at(c2).size());
      gsrcs.at(i) = storage.at(tl2s.at(i)-tmin).at(rl2).at(c2).at(ml2);
      times.at(i) = t.time(tl2s.at(i),rl2,ml2);
    }
    
    const iblist recv = (d.boxes.at(rl1).at(c1).at(ml1).*recv_listvect).at(c2);
    const iblist send = (d.boxes.at(rl2).at(c2).at(ml2).*send_listvect).at(c1);
    assert (recv.size()==send.size());
    // walk all boxes
    for (typename iblist::const_iterator r=recv.begin(), s=send.begin();
      	 r!=recv.end(); ++r, ++s) {
      // (use the send boxes for communication)
      // interpolate the content
      storage.at(tl1-tmin).at(rl1).at(c1).at(ml1)->interpolate_from
      	(state, gsrcs, times, *r, time,
	 d.prolongation_order_space, prolongation_order_time);
    }
  }
}



// Copy a component from the next time level
template<int D>
void ggf<D>::copy (comm_state<D>& state, int tl, int rl, int c, int ml)
{
  // Copy
  copycat (state,
           tl  ,rl,c,ml, &dh<D>::dboxes::exterior,
      	   tl+1,rl,  ml, &dh<D>::dboxes::exterior);
}

// Synchronise the boundaries a component
template<int D>
void ggf<D>::sync (comm_state<D>& state, int tl, int rl, int c, int ml)
{
  // Copy
  copycat (state,
           tl,rl,c,ml, &dh<D>::dboxes::recv_sync,
      	   tl,rl,  ml, &dh<D>::dboxes::send_sync);
}

// Prolongate the boundaries of a component
template<int D>
void ggf<D>::ref_bnd_prolongate (comm_state<D>& state, 
                                 int tl, int rl, int c, int ml,
                                 CCTK_REAL time)
{
  // Interpolate
  assert (rl>=1);
  if (transport_operator == op_none) return;
  vector<int> tl2s;
  // Interpolation in time
  assert (tmax-tmin+1 >= prolongation_order_time+1);
  tl2s.resize(prolongation_order_time+1);
  for (int i=0; i<=prolongation_order_time; ++i) tl2s.at(i) = tmax - i;
  intercat (state,
            tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_bnd_coarse,
	    tl2s,rl-1,  ml, &dh<D>::dboxes::send_ref_bnd_fine,
	    time);
}

// Restrict a multigrid level
template<int D>
void ggf<D>::mg_restrict (comm_state<D>& state,
                          int tl, int rl, int c, int ml,
                          CCTK_REAL time)
{
  // Require same times
  assert (abs(t.get_time(rl,ml) - t.get_time(rl,ml-1))
	  <= 1.0e-8 * abs(t.get_time(rl,ml)));
  const vector<int> tl2s(1,tl);
  intercat (state,
            tl  ,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml-1, &dh<D>::dboxes::send_mg_fine,
	    time);
}

// Prolongate a multigrid level
template<int D>
void ggf<D>::mg_prolongate (comm_state<D>& state,
                            int tl, int rl, int c, int ml,
                            CCTK_REAL time)
{
  // Require same times
  assert (abs(t.get_time(rl,ml) - t.get_time(rl,ml+1))
	  <= 1.0e-8 * abs(t.get_time(rl,ml)));
  const vector<int> tl2s(1,tl);
  intercat (state,
            tl  ,rl,c,ml,   &dh<D>::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml+1, &dh<D>::dboxes::send_mg_fine,
	    time);
}

// Restrict a refinement level
template<int D>
void ggf<D>::ref_restrict (comm_state<D>& state,
                           int tl, int rl, int c, int ml,
                           CCTK_REAL time)
{
  // Require same times
  assert (abs(t.get_time(rl,ml) - t.get_time(rl+1,ml))
	  <= 1.0e-8 * abs(t.get_time(rl,ml)));
  if (transport_operator == op_none) return;
  const vector<int> tl2s(1,tl);
  intercat (state,
            tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_fine,
	    tl2s,rl+1,  ml, &dh<D>::dboxes::send_ref_coarse,
	    time);
}

// Prolongate a refinement level
template<int D>
void ggf<D>::ref_prolongate (comm_state<D>& state,
                             int tl, int rl, int c, int ml,
                             CCTK_REAL time)
{
  assert (rl>=1);
  if (transport_operator == op_none) return;
  vector<int> tl2s;
  // Interpolation in time
  assert (tmax-tmin+1 >= prolongation_order_time+1);
  tl2s.resize(prolongation_order_time+1);
  for (int i=0; i<=prolongation_order_time; ++i) tl2s.at(i) = tmax - i;
  intercat (state,
            tl  ,rl  ,c,ml, &dh<D>::dboxes::recv_ref_coarse,
	    tl2s,rl-1,  ml, &dh<D>::dboxes::send_ref_fine,
	    time);
}



template class ggf<3>;
