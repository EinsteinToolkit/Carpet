// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.cc,v 1.21 2003/01/03 15:49:36 schnetter Exp $

#include <assert.h>

#include <iostream>

#include "cctk.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "gdata.hh"

using namespace std;



// Constructors
template<int D>
gdata<D>::gdata ()
  : _has_storage(false)
{ }

// Destructors
template<int D>
gdata<D>::~gdata () { }



// Data manipulators
template<int D>
void gdata<D>::copy_from (const gdata* src, const ibbox& box)
{
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      && box.stride()==src->extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0));
  
  if (proc() == src->proc()) {
    // copy on same processor
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == proc()) {
      copy_from_innerloop (src, box);
    }
    
  } else {
    
    // copy to different processor
    gdata* const tmp = make_typed();
    tmp->allocate (box, src->proc());
    tmp->copy_from (src, box);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



template<int D>
void gdata<D>
::interpolate_from (const vector<const gdata*> srcs,
		    const vector<CCTK_REAL> times,
		    const ibbox& box, const CCTK_REAL time,
		    const int order_space,
		    const int order_time)
{
  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs[t]->has_storage());
    assert (all(box.lower()>=srcs[t]->extent().lower()));
    assert (all(box.upper()<=srcs[t]->extent().upper()));
  }
  
  if (proc() == srcs[0]->proc()) {
    // interpolate on same processor
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == proc()) {
      interpolate_from_innerloop
	(srcs, times, box, time, order_space, order_time);
    }
    
  } else {
    // interpolate from other processor
    
    gdata* const tmp = make_typed();
    tmp->allocate (box, srcs[0]->proc());
    tmp->interpolate_from (srcs, times, box, time, order_space, order_time);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



template class gdata<3>;
