/***************************************************************************
                          gdata.cc  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.cc,v 1.9 2001/03/22 18:42:05 eschnett Exp $

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

#include <fstream>
#include <iomanip>

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
  : _has_storage(false)
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
::interpolate_from (const generic_data* src, const ibbox& box)
{
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.upper()<=extent().upper()));
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  
  if (proc() == src->proc()) {
    // interpolate on same processor
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == proc()) {
      interpolate_from_innerloop (src, box);
    }
    
  } else {
    // interpolate from other processor
    
    generic_data* const tmp = make_typed();
    tmp->allocate (box, src->proc());
    tmp->interpolate_from (src, box);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



template<int D>
void generic_data<D>
::interpolate_from (const generic_data* src1, const int t1,
		    const generic_data* src2, const int t2,
		    const ibbox& box, const int t)
{
  assert (has_storage() && src1->has_storage() && src2->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.upper()<=extent().upper()));
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src1->extent().lower()
	      && box.lower()>=src2->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src1->extent().upper()
	      && box.upper()<=src2->extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src1->extent().lower())%box.stride() == 0
	      && (box.lower()-src2->extent().lower())%box.stride() == 0));
  assert (src1->proc() == src2->proc());
  
  if (proc() == src1->proc()) {
    // interpolate on same processor
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == proc()) {
      interpolate_from_innerloop (src1, t1, src2, t2, box, t);
    }
    
  } else {
    // interpolate from other processor
    
    generic_data* const tmp = make_typed();
    tmp->allocate (box, src1->proc());
    tmp->interpolate_from (src1, t1, src2, t2, box, t);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



// Output
template<int D>
template<int DD>
void generic_data<D>::write_ascii (const string name, const int time,
				   const vect<int,D>& org,
				   const vect<int,DD>& dirs,
				   const int tl, const int rl,
				   const int c, const int ml)
  const
{
  assert (_has_storage);
  CHECKPOINT;
  
  if (proc()==0) {
    // output on processor 0
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == 0) {
      
      ofstream file(name.c_str(), ios::app);
      assert (file.good());
      
      file << setprecision(15);
      
      file << "# iteration " << time << endl
	   << "# time level " << tl << "   refinement level " << rl
	   << "   component " << c << "   multigrid level " << ml << endl
	   << "# column format: it tl rl c ml";
      assert (D>=1 && D<=3);
      const char* const coords = "xyz";
      for (int d=0; d<D; ++d) file << " " << coords[d];
      file << " data" << endl;
      
      const vect<int,DD> lo = extent().lower()[dirs];
      const vect<int,DD> up = extent().upper()[dirs];
      const vect<int,DD> str = extent().stride()[dirs];
      const bbox<int,DD> ext(lo,up,str);
      
      // Check whether the output origin is contained in the extent of
      // the data
      ivect org1(org);
      for (int d=0; d<DD; ++d) org1[dirs[d]] = ext.lower()[d];
      if (extent().contains(org1)) {
	
	for (bbox<int,DD>::iterator it=ext.begin(); it!=ext.end(); ++it) {
	  ivect index(org);
	  for (int d=0; d<DD; ++d) index[dirs[d]] = (*it)[d];
	  file << time << " " << tl << " " << rl << " " << c << " " << ml
	       << " ";
	  for (int d=0; d<D; ++d) file << index[d] << " ";
	  write_ascii_output_element (file, index);
	  file << endl;
	  for (int d=0; d<DD; ++d) {
	    if (index[dirs[d]]!=extent().upper()[dirs[d]]) break;
	    file << endl;
	  }
	}
	
      }	else {
	
	file << "#" << endl;
	
      }	// if ! ext contains org
      
      file.close();
      assert (file.good());
      
    }
    
  } else {
    // copy to processor 0 and output there
    
    generic_data* const tmp = make_typed();
    tmp->allocate(extent(), 0);
    tmp->copy_from (this, extent());
    tmp->write_ascii (name, time, org, dirs, tl, rl, c, ml);
    delete tmp;
    
  }
}



template<int D>
ostream& operator<< (ostream& os, const generic_data<D>& f) {
  return f.out(os);
}



#if defined(TMPL_EXPLICIT)
template class generic_data<3>;
template ostream& operator<< (ostream& os, const generic_data<3>& d);

template void generic_data<3>
::write_ascii (const string name, const int time,
	       const vect<int,3>& org, const vect<int,1>& dirs,
	       const int tl, const int rl, const int c, const int ml) const;
template void generic_data<3>
::write_ascii (const string name, const int time,
	       const vect<int,3>& org, const vect<int,2>& dirs,
	       const int tl, const int rl, const int c, const int ml) const;
template void generic_data<3>
::write_ascii (const string name, const int time,
	       const vect<int,3>& org, const vect<int,3>& dirs,
	       const int tl, const int rl, const int c, const int ml) const;
#endif
