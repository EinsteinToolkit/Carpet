/***************************************************************************
                          gdata.cc  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.cc,v 1.18 2002/01/09 23:42:42 schnetter Exp $

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

#include <iostream>

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GDATA_HH)
#  include "gdata.hh"
#endif

using namespace std;



// Constructors
template<int D>
generic_data<D>::generic_data ()
{ }

// Destructors
template<int D>
generic_data<D>::~generic_data () { }



// Data manipulators
template<int D>
void generic_data<D>::copy_from (const generic_data* src, const ibbox& box)
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
    generic_data* const tmp = make_typed();
    tmp->allocate (box, src->proc());
    tmp->copy_from (src, box);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



template<int D>
void generic_data<D>
::interpolate_from (const vector<const generic_data*> srcs,
		    const vector<int> tls,
		    const ibbox& box, const int tl,
		    const int order_space,
		    const int order_time)
{
  assert (has_storage());
  assert (all(box.lower()>=extent().lower()));
  assert (all(box.upper()<=extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  assert (srcs.size() == tls.size() && srcs.size()>0);
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
      interpolate_from_innerloop (srcs, tls, box, tl, order_space, order_time);
    }
    
  } else {
    // interpolate from other processor
    
    generic_data* const tmp = make_typed();
    tmp->allocate (box, srcs[0]->proc());
    tmp->interpolate_from (srcs, tls, box, tl, order_space, order_time);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



#if defined(TMPL_EXPLICIT)

template class generic_data<3>;

#endif
