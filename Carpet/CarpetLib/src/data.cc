/***************************************************************************
                          data.cc  -  Data storage
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/data.cc,v 1.6 2001/03/17 00:35:52 eschnett Exp $

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
#include <fstream>
#include <string>

#include <mpi.h>

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#if !defined(TMPL_IMPLICIT) || !defined(DATA_HH)
#  include "data.hh"
#endif



// Constructors
template<class T, int D>
data<T,D>::data ()
  : _storage(0)
{ }

template<class T, int D>
data<T,D>::data (const ibbox& extent, const int proc)
  : _storage(0)
{
  allocate(extent, proc);
}

// Destructors
template<class T, int D>
data<T,D>::~data () {
  free();
}
  
// Pseudo constructors
template<class T, int D>
data<T,D>* data<T,D>::make_typed () const {
  return new data();
}



// Storage management
template<class T, int D>
void data<T,D>::allocate (const ibbox& extent, const int proc,
			  void* const mem=0) {
  assert (!_has_storage);
  _has_storage = true;
  // data
  _extent = extent;
  _shape = max(ivect(0), _extent.shape() / _extent.stride());
  _size = 1;
  for (int d=0; d<D; ++d) {
    _stride[d] = _size;
    _size *= _shape[d];
  }
  _proc = proc;
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  if (rank==_proc) {
    _owns_storage = !mem;
    if (_owns_storage) {
      _storage = new T[_size];
    } else {
      _storage = (T*)mem;
    }
  } else {
    assert (!mem);
  }
}

template<class T, int D>
void data<T,D>::free () {
  if (_storage && _owns_storage) delete [] _storage;
  _storage = 0;
  _has_storage = false;
}

template<class T, int D>
void data<T,D>::transfer_from (generic_data<D>* gsrc) {
  data* src = (data*)gsrc;
  assert (!_storage);
  *this = *src;
  *src = data();
}



// Processor management
template<class T, int D>
void data<T,D>::change_processor (const int newproc, void* const mem=0) {
  if (newproc == _proc) {
    assert (!mem);
    return;
  }
  
  if (_has_storage) {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == newproc) {
      // copy from other processor
      
      assert (!_storage);
      _owns_storage = !mem;
      if (_owns_storage) {
	_storage = new T[_size];
      } else {
	_storage = (T*)mem;
      }

      T dummy;
      MPI_Status status;
      MPI_Recv (_storage, _size, dist::datatype(dummy), _proc,
		dist::tag, dist::comm, &status);

    } else if (rank == _proc) {
      // copy to other processor
      
      assert (!mem);
      assert (_storage);
      T dummy;
      MPI_Send (_storage, _size, dist::datatype(dummy), newproc,
		dist::tag, dist::comm);

      if (_owns_storage) {
	delete [] _storage;
      }
      _storage = 0;
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }

  _proc = newproc;
}



// Data manipulators
template<class T, int D>
void data<T,D>::copy_from (const generic_data<D>* gsrc, const ibbox& box) {
  const data* src = (const data*)gsrc;
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
      
      for (ibbox::iterator it=box.begin(); it!=box.end(); ++it) {
        const ivect index = *it;
        (*this)[index] = (*src)[index];
      }
      
    }
    
  } else {
    
    // copy to different processor
    data* const tmp = new data(box, src->proc());
    tmp->copy_from (src, box);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}

template<class T, int D>
void data<T,D>::interpolate_from (const generic_data<D>* gsrc,
                                  const ibbox& box)
{
  const data* src = (const data*)gsrc;
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.upper()<=extent().upper()));
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      /* && box.stride()<=src->extent().stride() */ ));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      /* && (box.lower()-src->extent().lower())%box.stride() == 0 */ ));
  
  if (proc() == src->proc()) {
    // interpolate on same processor
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == proc()) {
      
      for (ibbox::iterator posi=box.begin(); posi!=box.end(); ++posi) {
        const ivect& pos = *posi;
	
        // get box around current position
        const ibbox frombox
	  = ibbox(pos,pos,extent().stride()).expanded_for(src->extent());
	
        // interpolate from box to position
        T sum = 0;
        for (ibbox::iterator fromposi=frombox.begin();
	     fromposi!=frombox.end(); ++fromposi)
     	  {
	    const ivect& frompos = *fromposi;
	    
	    // interpolation weight
	    const ivect str = src->extent().stride();
	    const T f = prod(vect<T,D>(str - abs(pos - frompos))
			     / vect<T,D>(str));
	    sum += f * (*src)[frompos];
	  }
        (*this)[pos] = sum;
	
      } // for pos
      
    }
    
  } else {
    // interpolate from other processor
    
    data* const tmp = new data(box, src->proc());
    tmp->interpolate_from (src, box);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}

template<class T, int D>
void data<T,D>::interpolate_from (const generic_data<D>* gsrc,
                                  const double sfact,
				  const generic_data<D>* gtrc,
				  const double tfact,
				  const ibbox& box)
{
  const data* src = (const data*)gsrc;
  const data* trc = (const data*)gtrc;
  assert (has_storage() && src->has_storage() && trc->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.upper()<=extent().upper()));
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()
	      && box.lower()>=trc->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()
	      && box.upper()<=trc->extent().upper()));
  assert (all(box.stride()==extent().stride()
	      /* && box.stride()<=src->extent().stride() */
	      /* && trc->extent().stride()==src->extent().stride() */ ));
  assert (all((box.lower()-extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0
	      && (box.lower()-trc->extent().lower())%box.stride() == 0));
  assert (src->proc() == trc->proc());
  
  if (proc() == src->proc()) {
    // interpolate on same processor
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == proc()) {
      
      for (ibbox::iterator posi=box.begin(); posi!=box.end(); ++posi) {
        const ivect& pos = *posi;
	
        // get box around current position
        const ibbox frombox
	  = ibbox(pos,pos,extent().stride()).expanded_for(src->extent());
	
        // interpolate from box to position
        T src_sum = 0;
        T trc_sum = 0;
        for (ibbox::iterator fromposi=frombox.begin();
	     fromposi!=frombox.end(); ++fromposi)
  	    {
	      const ivect& frompos = *fromposi;
	      
	      // interpolation weight
	      const ivect str = src->extent().stride();
	      const T f = prod(vect<T,D>(str - abs(pos - frompos))
			       / vect<T,D>(str));
	      src_sum += f * (*src)[frompos];
	      trc_sum += f * (*trc)[frompos];
	    }
        (*this)[pos] = (T)sfact * src_sum + (T)tfact * trc_sum;
	
      } // for pos
      
    }
    
  } else {
    // interpolate from other processor
    
    data* const tmp = new data(box, src->proc());
    tmp->interpolate_from (src, sfact, trc, tfact, box);
    tmp->change_processor (proc());
    copy_from (tmp, box);
    delete tmp;
    
  }
}



// Output
template<class T, int D>
void data<T,D>::write_ascii_output_element (ofstream& file, const ivect& index)
  const
{
  file << (*this)[index];
}



// Output
template<class T,int D>
ostream& data<T,D>::out (ostream& os) const {
  os << "data<" STR(T) "," << D << ">:"
     << "extent=" << extent() << ","
     << "stride=" << stride() << ",size=" << size();
  return os;
}



#if defined(TMPL_EXPLICIT)

#define INSTANTIATE(T)				\
template data<T,1>;				\
template data<T,2>;				\
template data<T,3>;

#include "instantiate"

#undef INSTANTIATE

#endif
