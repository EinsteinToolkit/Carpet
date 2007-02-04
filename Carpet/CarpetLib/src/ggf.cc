#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "cctk.h"

#include "defs.hh"
#include "dh.hh"
#include "th.hh"
#include "timestat.hh"

#include "ggf.hh"

using namespace std;
using namespace CarpetLib;



// Constructors
ggf::ggf (const int varindex_, const operator_type transport_operator_,
          th& t_, dh& d_,
          const int prolongation_order_time_,
          const int vectorlength_, const int vectorindex_,
          ggf* const vectorleader_)
  : varindex(varindex_), transport_operator(transport_operator_), t(t_),
    prolongation_order_time(prolongation_order_time_),
    h(d_.h), d(d_),
    storage(h.mglevels()),
    vectorlength(vectorlength_), vectorindex(vectorindex_),
    vectorleader(vectorleader_)
{
  assert (&t.h == &d.h);
  
  assert (vectorlength >= 1);
  assert (vectorindex >= 0 and vectorindex < vectorlength);
  assert ((vectorindex==0 and !vectorleader)
          or (vectorindex!=0 and vectorleader));
  
  timelevels_.resize(d.h.mglevels());
  for (int ml=0; ml<d.h.mglevels(); ++ml) {
    timelevels_.AT(ml).resize(d.h.reflevels(), 0);
  }
  
  d.add(this);
}

// Destructors
ggf::~ggf () {
  d.remove(this);
}

// Comparison
bool ggf::operator== (const ggf& f) const {
  return this == &f;
}



// Modifiers
void ggf::set_timelevels (const int ml, const int rl, const int new_timelevels)
{
  assert (ml>=0 and ml<(int)storage.size());
  assert (rl>=0 and rl<(int)storage.AT(ml).size());
  
  assert (new_timelevels >= 0);
  
  if (new_timelevels < timelevels(ml,rl)) {
    
    for (int c=0; c<(int)storage.AT(ml).AT(rl).size(); ++c) {
      for (int tl=new_timelevels; tl<timelevels(ml,rl); ++tl) {
        delete storage.AT(ml).AT(rl).AT(c).AT(tl);
      }
      storage.AT(ml).AT(rl).AT(c).resize (new_timelevels);
    } // for c
    
  } else if (new_timelevels > timelevels(ml,rl)) {
    
    for (int c=0; c<(int)storage.AT(ml).AT(rl).size(); ++c) {
      storage.AT(ml).AT(rl).AT(c).resize (new_timelevels);
      for (int tl=timelevels(ml,rl); tl<new_timelevels; ++tl) {
        storage.AT(ml).AT(rl).AT(c).AT(tl) = typed_data(tl,rl,c,ml);
        storage.AT(ml).AT(rl).AT(c).AT(tl)->allocate
          (d.boxes.AT(ml).AT(rl).AT(c).exterior, h.processor(rl,c));
      } // for tl
    } // for c
    
  }
  
  timelevels_.AT(ml).AT(rl) = new_timelevels;
}



void ggf::recompose_crop ()
{
  // Free storage that will not be needed
  for (int ml=0; ml<h.mglevels(); ++ml) {
    for (int rl=h.reflevels(); rl<(int)storage.AT(ml).size(); ++rl) {
      for (int c=0; c<(int)storage.AT(ml).AT(rl).size(); ++c) {
        for (int tl=0; tl<(int)storage.AT(ml).AT(rl).AT(c).size(); ++tl) {
          delete storage.AT(ml).AT(rl).AT(c).AT(tl);
        } // for tl
      } // for c
    } // for rl
    storage.AT(ml).resize(h.reflevels());
  } // for ml
}

void ggf::recompose_allocate (const int rl)
{
  // Retain storage that might be needed
  oldstorage.resize(storage.size());
  for (int ml=0; ml<(int)storage.size(); ++ml) {
    oldstorage.AT(ml).resize(storage.AT(ml).size());
    oldstorage.AT(ml).AT(rl) = storage.AT(ml).AT(rl);
    storage.AT(ml).AT(rl).resize(0);
  }
  
  for (int ml=0; ml<d.h.mglevels(); ++ml) {
    timelevels_.AT(ml).resize(d.h.reflevels(), timelevels_.AT(ml).AT(0));
  }
  
  // Resize structure and allocate storage
  storage.resize(h.mglevels());
  for (int ml=0; ml<h.mglevels(); ++ml) {
    storage.AT(ml).resize(h.reflevels());
    storage.AT(ml).AT(rl).resize(h.components(rl));
    for (int c=0; c<h.components(rl); ++c) {
      storage.AT(ml).AT(rl).AT(c).resize(timelevels(ml,rl));
      for (int tl=0; tl<timelevels(ml,rl); ++tl) {
        storage.AT(ml).AT(rl).AT(c).AT(tl) = typed_data(tl,rl,c,ml);
        storage.AT(ml).AT(rl).AT(c).AT(tl)->allocate
          (d.boxes.AT(ml).AT(rl).AT(c).exterior, h.processor(rl,c));
      } // for tl
    } // for c
  } // for ml
}

void ggf::recompose_fill (comm_state& state, const int rl,
                          const bool do_prolongate)
{
  // Initialise the new storage
  for (int ml=0; ml<h.mglevels(); ++ml) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int tl=0; tl<timelevels(ml,rl); ++tl) {
        
        // Find out which regions need to be prolongated
        // (Copy the exterior because some variables are not prolongated)
        // TODO: do this once in the dh instead of for each variable here
        ibset work (d.boxes.AT(ml).AT(rl).AT(c).exterior);
        
        // Copy from old storage, if possible
        // TODO: copy only from interior regions?
        if (rl<(int)oldstorage.AT(ml).size()) {
          for (int cc=0; cc<(int)oldstorage.AT(ml).AT(rl).size(); ++cc) {
            // TODO: prefer same processor, etc., see dh.cc
            ibset ovlp
              = work & oldstorage.AT(ml).AT(rl).AT(cc).AT(tl)->extent();
            ovlp.normalize();
            work -= ovlp;
            for (ibset::const_iterator r=ovlp.begin(); r!=ovlp.end(); ++r) {
              storage.AT(ml).AT(rl).AT(c).AT(tl)->copy_from
                (state, oldstorage.AT(ml).AT(rl).AT(cc).AT(tl), *r);
            }
          } // for cc
        } // if rl
        
        if (do_prolongate) {
          // Initialise from coarser level, if possible
          if (rl>0) {
            if (transport_operator != op_none) {
              const int pos = d.prolongation_order_space;
              const int pot = (transport_operator != op_copy
                               ? prolongation_order_time
                               : 0);
              const int numtl = pot+1;
              assert (timelevels(ml,rl) >= numtl);
              vector<int> tls(numtl);
              vector<CCTK_REAL> times(numtl);
              for (int i=0; i<numtl; ++i) {
                tls.AT(i) = i;
                times.AT(i) = t.time(tls.AT(i),rl-1,ml);
              }
              for (int cc=0; cc<(int)storage.AT(ml).AT(rl-1).size(); ++cc) {
                vector<const gdata*> gsrcs(numtl);
                for (int i=0; i<numtl; ++i) {
                  gsrcs.AT(i) = storage.AT(ml).AT(rl-1).AT(cc).AT(tls.AT(i));
                  assert (gsrcs.AT(i)->extent() == gsrcs.AT(0)->extent());
                }
                const CCTK_REAL time = t.time(tl,rl,ml);
                
                // TODO: choose larger regions first
                // TODO: prefer regions from the same processor
                const iblist& list
                  = d.boxes.AT(ml).AT(rl).AT(c).recv_ref_coarse.AT(cc);
                for (iblist::const_iterator iter=list.begin();
                     iter!=list.end(); ++iter)
                {
                  ibset ovlp = work & *iter;
                  ovlp.normalize();
                  work -= ovlp;
                  for (ibset::const_iterator r=ovlp.begin();
                       r!=ovlp.end(); ++r)
                  {
                    storage.AT(ml).AT(rl).AT(c).AT(tl)->interpolate_from
                      (state, gsrcs, times, *r, time, pos, pot);
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
    } // for c
  } // for ml
}

void ggf::recompose_free_old (const int rl)
{
  // Delete old storage
  for (int ml=0; ml<(int)oldstorage.size(); ++ml) {
    for (int c=0; c<(int)oldstorage.AT(ml).AT(rl).size(); ++c) {
      for (int tl=0; tl<timelevels(ml,rl); ++tl) {
        delete oldstorage.AT(ml).AT(rl).AT(c).AT(tl);
      } // for tl
    } // for c
    oldstorage.AT(ml).AT(rl).resize(0);
  } // for ml
}

void ggf::recompose_bnd_prolongate (comm_state& state, const int rl,
                                    const bool do_prolongate)
{
  if (do_prolongate) {
    // Set boundaries
    if (rl>0) {
      for (int ml=0; ml<h.mglevels(); ++ml) {
        for (int c=0; c<h.components(rl); ++c) {
          for (int tl=0; tl<timelevels(ml,rl); ++tl) {
            
            // TODO: assert that reflevel 0 boundaries are copied
            const CCTK_REAL time = t.time(tl,rl,ml);
            ref_bnd_prolongate (state,tl,rl,c,ml,time);
            
          } // for tl
        } // for c
      } // for ml
    } // if rl
  } // if do_prolongate
}

void ggf::recompose_sync (comm_state& state, const int rl,
                          const bool do_prolongate)
{
  if (do_prolongate) {
    // Set boundaries
    for (int ml=0; ml<h.mglevels(); ++ml) {
      for (int c=0; c<h.components(rl); ++c) {
        for (int tl=0; tl<timelevels(ml,rl); ++tl) {
          
          sync (state,tl,rl,c,ml);
          
        } // for tl
      } // for c
    } // for ml
  } // if do_prolongate
}



// Cycle the time levels by rotating the data sets
void ggf::cycle (int const rl, int const c, int const ml) {
  assert (rl>=0 and rl<h.reflevels());
  assert (c>=0 and c<h.components(rl));
  assert (ml>=0 and ml<h.mglevels());
  int const ntl = timelevels(ml,rl);
  assert (ntl > 0);
  fdata & fdatas = storage.AT(ml).AT(rl).AT(c);
  gdata * const tmpdata = fdatas.AT(ntl-1);
  for (int tl=ntl-1; tl>0; --tl) {
    fdatas.AT(tl) = fdatas.AT(tl-1);
  }
  fdatas.AT(0) = tmpdata;
}

// Flip the time levels by exchanging the data sets
void ggf::flip (int rl, int c, int ml) {
  assert (rl>=0 and rl<h.reflevels());
  assert (c>=0 and c<h.components(rl));
  assert (ml>=0 and ml<h.mglevels());
  for (int tl=0; tl<(timelevels(ml,rl)-1)/2; ++tl) {
    const int tl1 =                  tl;
    const int tl2 = timelevels(ml,rl)-1 - tl;
    assert (tl1 < tl2);
    gdata* tmpdata = storage.AT(ml).AT(rl).AT(c).AT(tl1);
    storage.AT(ml).AT(rl).AT(c).AT(tl1) = storage.AT(ml).AT(rl).AT(c).AT(tl2);
    storage.AT(ml).AT(rl).AT(c).AT(tl2) = tmpdata;
  }
}



// Operations

// Copy a region
void ggf::copycat (comm_state& state,
                   int tl1, int rl1, int c1, int ml1,
                   const ibbox dh::dboxes::* recv_box,
                   int tl2, int rl2, int ml2)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (rl2>=0 and rl2<h.reflevels());
  int const c2= c1;
  assert (ml2<h.mglevels());
  assert (tl2>=0 and tl2<timelevels(ml2,rl2));
  static Timer copycat1 ("copycat_1");
  copycat1.start ();
  ibbox const & recv = d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_box;
  // copy the content
  gdata * const dst = storage.AT(ml1).AT(rl1).AT(c1).AT(tl1);
  gdata const * const src = storage.AT(ml2).AT(rl2).AT(c2).AT(tl2);
  dst->copy_from(state, src, recv);
  copycat1.stop (0);
}

// Copy regions
void ggf::copycat (comm_state& state,
                   int tl1, int rl1, int c1, int ml1,
                   const iblist dh::dboxes::* recv_list,
                   int tl2, int rl2, int ml2)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (           ml2<h.mglevels());
  assert (rl2>=0 and rl2<h.reflevels());
  int const c2 = c1;
  assert (tl2>=0 and tl2<timelevels(ml2,rl2));
  static Timer copycat1 ("copycat_list_1");
  copycat1.start ();
  iblist const & recv = d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_list;
  // walk all boxes
  for (iblist::const_iterator r=recv.begin(); r!=recv.end(); ++r) {
    // (use the send boxes for communication)
    // copy the content
    gdata * const dst = storage.AT(ml1).AT(rl1).AT(c1).AT(tl1);
    gdata const * const src = storage.AT(ml2).AT(rl2).AT(c2).AT(tl2);
    dst->copy_from(state, src, *r);
  }
  copycat1.stop (0);
}

// Copy regions
void ggf::copycat (comm_state& state,
                   int tl1, int rl1, int c1, int ml1,
                   const iblistvect dh::dboxes::* recv_listvect,
                   int tl2, int rl2, int ml2)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (           ml2<h.mglevels());
  assert (rl2>=0 and rl2<h.reflevels());
  assert (tl2>=0 and tl2<timelevels(ml2,rl2));
  // walk all components
  static Timer copycat1 ("copycat_listvect_1");
  copycat1.start ();
  iblistvect::const_iterator recvi =
    (d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_listvect).begin();
  gdata * const dst = storage.AT(ml1).AT(rl1).AT(c1).AT(tl1);
  cdata::const_iterator srci = storage.AT(ml2).AT(rl2).begin();
  int const nc = h.components(rl2);
  for (int c2=0; c2<nc; ++c2) {
    iblist const & recv = (*recvi);
    gdata const * const src = (*srci).AT(tl2);
    // walk all boxes
    for (iblist::const_iterator r=recv.begin(); r!=recv.end(); ++r) {
      // (use the send boxes for communication)
      // copy the content
      dst->copy_from(state, src, *r);
    }
    ++ recvi;
    ++ srci;
  }
  copycat1.stop (0);
}

// Copy regions
void ggf::copycat (comm_state& state,
                   int tl1, int rl1, int c1, int ml1,
                   const pvect dh::dboxes::* recv_pvect,
                   int tl2, int rl2, int ml2)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (           ml2<h.mglevels());
  assert (rl2>=0 and rl2<h.reflevels());
  assert (tl2>=0 and tl2<timelevels(ml2,rl2));
  // walk all components
  static Timer copycat1 ("copycat_pvect_1");
  copycat1.start ();
  gdata * const dst = storage.AT(ml1).AT(rl1).AT(c1).AT(tl1);
  cdata const & srcs = storage.AT(ml2).AT(rl2);
  pvect const & prs = d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_pvect;
  for (pvect::const_iterator ipr=prs.begin(); ipr!=prs.end(); ++ipr) {
    pseudoregion const & pr = * ipr;
    ibbox const & r = pr.extent;
    int const c2 = pr.processor;
    gdata const * const src = srcs.AT(c2).AT(tl2);
    dst->copy_from(state, src, r);
  }
  copycat1.stop (0);
}

// Interpolate a region
void ggf::intercat (comm_state& state,
                    int tl1, int rl1, int c1, int ml1,
                    const ibbox dh::dboxes::* recv_list,
                    const vector<int> tl2s, int rl2, int ml2,
                    CCTK_REAL time)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (rl2>=0 and rl2<h.reflevels());
  int const c2 = c1;
  assert (ml2>=0 and ml2<h.mglevels());
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s.AT(i)>=0 and tl2s.AT(i)<timelevels(ml2,rl2));
  }
  
  vector<const gdata*> gsrcs(tl2s.size());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    gsrcs.AT(i) = storage.AT(ml2).AT(rl2).AT(c2).AT(tl2s.AT(i));
    times.AT(i) = t.time(tl2s.AT(i),rl2,ml2);
  }
  
  static Timer intercat1 ("intercat_1");
  intercat1.start ();
  ibbox const & recv = d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_list;
  // interpolate the content
  storage.AT(ml1).AT(rl1).AT(c1).AT(tl1)->interpolate_from
    (state, gsrcs, times, recv, time,
     d.prolongation_order_space, prolongation_order_time);
  intercat1.stop (0);
}

// Interpolate regions
void ggf::intercat (comm_state& state,
                    int tl1, int rl1, int c1, int ml1,
                    const iblist dh::dboxes::* recv_list,
                    const vector<int> tl2s, int rl2, int ml2,
                    const CCTK_REAL time)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (rl2>=0 and rl2<h.reflevels());
  int const c2 = c1;
  assert (ml2>=0 and ml2<h.mglevels());
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s.AT(i)>=0 and tl2s.AT(i)<timelevels(ml2,rl2));
  }
  
  vector<const gdata*> gsrcs(tl2s.size());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)gsrcs.size(); ++i) {
    gsrcs.AT(i) = storage.AT(ml2).AT(rl2).AT(c2).AT(tl2s.AT(i));
    times.AT(i) = t.time(tl2s.AT(i),rl2,ml2);
  }
  
  static Timer intercat1 ("intercat_list_1");
  intercat1.start ();
  iblist const & recv = d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_list;
  // walk all boxes
  for (iblist::const_iterator r=recv.begin(); r!=recv.end(); ++r)
  {
    // (use the send boxes for communication)
    // interpolate the content
    storage.AT(ml1).AT(rl1).AT(c1).AT(tl1)->interpolate_from
      (state, gsrcs, times, *r, time,
       d.prolongation_order_space, prolongation_order_time);
  }
  intercat1.stop (0);
}

// Interpolate regions
void ggf::intercat (comm_state& state,
                    int tl1, int rl1, int c1, int ml1,
                    const iblistvect dh::dboxes::* recv_listvect,
                    const vector<int> tl2s, int rl2, int ml2,
                    const CCTK_REAL time)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (rl2>=0 and rl2<h.reflevels());
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s.AT(i)>=0 and tl2s.AT(i)<timelevels(ml2,rl2));
  }
  // walk all components
  static Timer intercat1 ("intercat_listvect_1");
  intercat1.start ();
  for (int c2=0; c2<h.components(rl2); ++c2) {
    assert (ml2>=0 and ml2<h.mglevels());
    
    vector<const gdata*> gsrcs(tl2s.size());
    vector<CCTK_REAL> times(tl2s.size());
    for (int i=0; i<(int)gsrcs.size(); ++i) {
      gsrcs.AT(i) = storage.AT(ml2).AT(rl2).AT(c2).AT(tl2s.AT(i));
      times.AT(i) = t.time(tl2s.AT(i),rl2,ml2);
    }
    
    iblist const &recv =
      (d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_listvect).AT(c2);
    // walk all boxes
    for (iblist::const_iterator r=recv.begin(); r!=recv.end(); ++r)
    {
      // (use the send boxes for communication)
      // interpolate the content
      int const pos = d.prolongation_order_space;
      int const pot =
        transport_operator == op_copy ? 0 : prolongation_order_time;
      storage.AT(ml1).AT(rl1).AT(c1).AT(tl1)->interpolate_from
      	(state, gsrcs, times, *r, time, pos, pot);
    }
  }
  intercat1.stop (0);
}

// Interpolate regions
void ggf::intercat (comm_state& state,
                    int tl1, int rl1, int c1, int ml1,
                    const pvect dh::dboxes::* recv_pvect,
                    const vector<int> tl2s, int rl2, int ml2,
                    const CCTK_REAL time)
{
  assert (rl1>=0 and rl1<h.reflevels());
  assert (c1>=0 and c1<h.components(rl1));
  assert (ml1>=0 and ml1<h.mglevels());
  assert (tl1>=0 and tl1<timelevels(ml1,rl1));
  assert (rl2>=0 and rl2<h.reflevels());
  assert (ml2>=0 and ml2<h.mglevels());
  vector<CCTK_REAL> times(tl2s.size());
  for (int i=0; i<(int)tl2s.size(); ++i) {
    assert (tl2s.AT(i)>=0 and tl2s.AT(i)<timelevels(ml2,rl2));
    times.AT(i) = t.time(tl2s.AT(i),rl2,ml2);
  }
  // walk all components
  static Timer intercat1 ("intercat_pvect_1");
  intercat1.start ();
  int const pos = d.prolongation_order_space;
  int const pot = transport_operator == op_copy ? 0 : prolongation_order_time;
  gdata * const dst = storage.AT(ml1).AT(rl1).AT(c1).AT(tl1);
  cdata const & srcs = storage.AT(ml2).AT(rl2);
  vector<const gdata*> gsrcs(tl2s.size());
  pvect const & prs = d.boxes.AT(ml1).AT(rl1).AT(c1).*recv_pvect;
  for (pvect::const_iterator ipr=prs.begin(); ipr!=prs.end(); ++ipr) {
    pseudoregion const & pr = * ipr;
    ibbox const & r = pr.extent;
    int const c2 = pr.processor;
    for (int i=0; i<(int)gsrcs.size(); ++i) {
      gsrcs.AT(i) = srcs.AT(c2).AT(tl2s.AT(i));
    }
    dst->interpolate_from (state, gsrcs, times, r, time, pos, pot);
  }
  intercat1.stop (0);
}



// Copy a component from the next time level
void ggf::copy (comm_state& state, int tl, int rl, int c, int ml)
{
  // Copy
  static Timer timer ("copy");
  timer.start ();
  copycat (state,
           tl  ,rl,c,ml, &dh::dboxes::exterior,
      	   tl+1,rl,  ml);
  timer.stop (0);
}

// Synchronise the boundaries a component
void ggf::sync (comm_state& state, int tl, int rl, int c, int ml)
{
  // Copy
  static Timer timer ("sync");
  timer.start ();
  copycat (state,
           tl,rl,c,ml, &dh::dboxes::recv_sync_fast,
      	   tl,rl,  ml);
  timer.stop (0);
}

// Prolongate the boundaries of a component
void ggf::ref_bnd_prolongate (comm_state& state, 
                              int tl, int rl, int c, int ml,
                              CCTK_REAL time)
{
  // Interpolate
  assert (rl>=1);
  if (transport_operator == op_none) return;
  vector<int> tl2s;
  static Timer timer ("ref_bnd_prolongate");
  timer.start ();
  if (transport_operator != op_copy) {
    // Interpolation in time
    if (not (timelevels(ml,rl) >= prolongation_order_time+1)) {
      char * const fullname = CCTK_FullName (varindex);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The variable \"%s\" has only %d active time levels, which is not enough for boundary prolongation of order %d",
                  fullname ? fullname : "<unknown variable>",
                  timelevels(ml,rl), prolongation_order_time);
      free (fullname);
    }
    assert (timelevels(ml,rl) >= prolongation_order_time+1);
    tl2s.resize(prolongation_order_time+1);
    for (int i=0; i<=prolongation_order_time; ++i) tl2s.AT(i) = i;
  } else {
    assert (timelevels(ml,rl) >= 1);
    tl2s.resize(1);
    tl2s.AT(0) = 0;
  }
  intercat (state,
            tl  ,rl  ,c,ml, &dh::dboxes::recv_ref_bnd_coarse_fast,
            tl2s,rl-1,  ml, time);
  timer.stop (0);
}

// Restrict a multigrid level
void ggf::mg_restrict (comm_state& state,
                       int tl, int rl, int c, int ml,
                       CCTK_REAL time)
{
  // Require same times
  assert (abs(t.get_time(rl,ml) - t.get_time(rl,ml-1))
	  <= 1.0e-8 * abs(t.get_time(rl,ml)));
  const vector<int> tl2s(1,tl);
  static Timer timer ("mg_restrict");
  timer.start ();
  intercat (state,
            tl  ,rl,c,ml,   &dh::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml-1, time);
  timer.stop (0);
}

// Prolongate a multigrid level
void ggf::mg_prolongate (comm_state& state,
                         int tl, int rl, int c, int ml,
                         CCTK_REAL time)
{
  // Require same times
  assert (abs(t.get_time(rl,ml) - t.get_time(rl,ml+1))
	  <= 1.0e-8 * abs(t.get_time(rl,ml)));
  static Timer timer ("mg_prolongate");
  timer.start ();
  const vector<int> tl2s(1,tl);
  intercat (state,
            tl  ,rl,c,ml,   &dh::dboxes::recv_mg_coarse,
	    tl2s,rl,  ml+1, time);
  timer.stop (0);
}

// Restrict a refinement level
void ggf::ref_restrict (comm_state& state,
                        int tl, int rl, int c, int ml,
                        CCTK_REAL time)
{
  // Require same times
  assert (abs(t.get_time(rl,ml) - t.get_time(rl+1,ml))
	  <= 1.0e-8 * abs(t.get_time(rl,ml)));
  if (transport_operator == op_none) return;
  static Timer timer ("ref_restrict");
  timer.start ();
  vector<int> const tl2s(1,tl);
  intercat (state,
            tl  ,rl  ,c,ml, &dh::dboxes::recv_ref_fine_fast,
	    tl2s,rl+1,  ml, time);
  timer.stop (0);
}

// Prolongate a refinement level
void ggf::ref_prolongate (comm_state& state,
                          int tl, int rl, int c, int ml,
                          CCTK_REAL time)
{
  assert (rl>=1);
  if (transport_operator == op_none) return;
  static Timer timer ("ref_prolongate");
  timer.start ();
  vector<int> tl2s;
  // Interpolation in time
  assert (timelevels(ml,rl) >= prolongation_order_time+1);
  tl2s.resize(prolongation_order_time+1);
  for (int i=0; i<=prolongation_order_time; ++i) tl2s.AT(i) = i;
  intercat (state,
            tl  ,rl  ,c,ml, &dh::dboxes::recv_ref_coarse_fast,
	    tl2s,rl-1,  ml, time);
  timer.stop (0);
}
