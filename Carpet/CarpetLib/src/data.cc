/***************************************************************************
                          data.cc  -  Data storage
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/data.cc,v 1.8 2001/03/22 18:42:05 eschnett Exp $

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
#include <string>

#include <mpi.h>

#include "cctk.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#if !defined(TMPL_IMPLICIT) || !defined(DATA_HH)
#  include "data.hh"
#endif

using namespace std;



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
			  void* const mem) {
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
void data<T,D>::change_processor (const int newproc, void* const mem) {
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
void data<T,D>
::copy_from_innerloop (const generic_data<D>* gsrc, const ibbox& box)
{
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
  
  assert (proc() == src->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
  for (ibbox::iterator it=box.begin(); it!=box.end(); ++it) {
    const ivect index = *it;
    (*this)[index] = (*src)[index];
  }
  
}



template<class T, int D>
void data<T,D>
::interpolate_from_innerloop (const generic_data<D>* gsrc, const ibbox& box)
{
  const data* src = (const data*)gsrc;
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.upper()<=extent().upper()));
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  
  assert (proc() == src->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
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



template<class T, int D>
void data<T,D>
::interpolate_from_innerloop (const generic_data<D>* gsrc1, const int t1,
			      const generic_data<D>* gsrc2, const int t2,
			      const ibbox& box, const int t)
{
  const data* src1 = (const data*)gsrc1;
  const data* src2 = (const data*)gsrc2;
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
  
  assert (proc() == src1->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
  const T one = 1;
  const T src1_fac = (T)((t - t2) * one / (t1 - t2));
  const T src2_fac = (T)((t - t1) * one / (t2 - t1));
  
  for (ibbox::iterator posi=box.begin(); posi!=box.end(); ++posi) {
    const ivect& pos = *posi;
    
    // get box around current position
    const ibbox frombox
      = ibbox(pos,pos,extent().stride()).expanded_for(src1->extent());
    
    // interpolate from box to position
    T src1_sum = 0;
    T src2_sum = 0;
    for (ibbox::iterator fromposi=frombox.begin();
	 fromposi!=frombox.end(); ++fromposi)
      {
	const ivect& frompos = *fromposi;
	
	// interpolation weight
	const ivect str = src1->extent().stride();
	const T f = prod(vect<T,D>(str - abs(pos - frompos))
			 / vect<T,D>(str));
	src1_sum += f * (*src1)[frompos];
	src2_sum += f * (*src2)[frompos];
      }
    (*this)[pos] = src1_fac * src1_sum + src2_fac * src2_sum;
    
  } // for pos
  
}



extern "C" {
  void CCTK_FCALL CCTK_FNAME(copy_3d_real8)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
}

template<>
void data<CCTK_REAL8,3>
::copy_from_innerloop (const generic_data<3>* gsrc, const ibbox& box)
{
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
  
  assert (proc() == src->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
  const ibbox& sext = src->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  for (int d=0; d<3; ++d) {
    srcshp[d] = (sext.shape() / sext.stride())[d];
    dstshp[d] = (dext.shape() / dext.stride())[d];
    
    srcbbox[0][d] = sext.lower()[d];
    srcbbox[1][d] = sext.upper()[d];
    srcbbox[2][d] = sext.stride()[d];
    
    dstbbox[0][d] = dext.lower()[d];
    dstbbox[1][d] = dext.upper()[d];
    dstbbox[2][d] = dext.stride()[d];
    
    regbbox[0][d] = box.lower()[d];
    regbbox[1][d] = box.upper()[d];
    regbbox[2][d] = box.stride()[d];
  }
  
  assert (all(dext.stride() == box.stride()));
  if (all(sext.stride() == dext.stride())) {
    CCTK_FNAME(copy_3d_real8) ((const CCTK_REAL8*)src->storage(),
			       srcshp[0], srcshp[1], srcshp[2],
			       (CCTK_REAL8*)storage(),
			       dstshp[0], dstshp[1], dstshp[2],
			       srcbbox,
			       dstbbox,
			       regbbox);
    
  } else {
    abort();
  }
}



extern "C" {
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(restrict_3d_real8)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
}

template<>
void data<CCTK_REAL8,3>
::interpolate_from_innerloop (const generic_data<3>* gsrc, const ibbox& box)
{
  const data* src = (const data*)gsrc;
  assert (has_storage() && src->has_storage());
  assert (all(box.lower()>=extent().lower()
	      && box.upper()<=extent().upper()));
  assert (all(box.lower()>=extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==extent().stride()));
  assert (all((box.lower()-extent().lower())%box.stride() == 0));
  
  assert (proc() == src->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
  const ibbox& sext = src->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  for (int d=0; d<3; ++d) {
    srcshp[d] = (sext.shape() / sext.stride())[d];
    dstshp[d] = (dext.shape() / dext.stride())[d];
    
    srcbbox[0][d] = sext.lower()[d];
    srcbbox[1][d] = sext.upper()[d];
    srcbbox[2][d] = sext.stride()[d];
    
    dstbbox[0][d] = dext.lower()[d];
    dstbbox[1][d] = dext.upper()[d];
    dstbbox[2][d] = dext.stride()[d];
    
    regbbox[0][d] = box.lower()[d];
    regbbox[1][d] = box.upper()[d];
    regbbox[2][d] = box.stride()[d];
  }
  
  assert (all(dext.stride() == box.stride()));
  if (all(sext.stride() > dext.stride())) {
    CCTK_FNAME(prolongate_3d_real8) ((const CCTK_REAL8*)src->storage(),
				     srcshp[0], srcshp[1], srcshp[2],
				     (CCTK_REAL8*)storage(),
				     dstshp[0], dstshp[1], dstshp[2],
				     srcbbox,
				     dstbbox,
				     regbbox);
  } else if (all(sext.stride() < dext.stride())) {
    CCTK_FNAME(restrict_3d_real8) ((const CCTK_REAL8*)src->storage(),
				   srcshp[0], srcshp[1], srcshp[2],
				   (CCTK_REAL8*)storage(),
				   dstshp[0], dstshp[1], dstshp[2],
				   srcbbox,
				   dstbbox,
				   regbbox);
  } else {
    abort();
  }
}



extern "C" {
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl)
    (const CCTK_REAL8* src1, const int& t1,
     const CCTK_REAL8* src2, const int& t2,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const int& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
}

template<>
void data<CCTK_REAL8,3>
::interpolate_from_innerloop (const generic_data<3>* gsrc1, const int t1,
			      const generic_data<3>* gsrc2, const int t2,
			      const ibbox& box, const int t)
{
  const data* src1 = (const data*)gsrc1;
  const data* src2 = (const data*)gsrc2;
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
  
  assert (proc() == src1->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
  assert (src1->extent() == src2->extent());
  
  const ibbox& sext = src1->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  for (int d=0; d<3; ++d) {
    srcshp[d] = (sext.shape() / sext.stride())[d];
    dstshp[d] = (dext.shape() / dext.stride())[d];
    
    srcbbox[0][d] = sext.lower()[d];
    srcbbox[1][d] = sext.upper()[d];
    srcbbox[2][d] = sext.stride()[d];
    
    dstbbox[0][d] = dext.lower()[d];
    dstbbox[1][d] = dext.upper()[d];
    dstbbox[2][d] = dext.stride()[d];
    
    regbbox[0][d] = box.lower()[d];
    regbbox[1][d] = box.upper()[d];
    regbbox[2][d] = box.stride()[d];
  }
  
  assert (all(dext.stride() == box.stride()));
  if (all(sext.stride() > dext.stride())) {
    CCTK_FNAME(prolongate_3d_real8_2tl)
      ((const CCTK_REAL8*)src1->storage(), t1,
       (const CCTK_REAL8*)src2->storage(), t2,
       srcshp[0], srcshp[1], srcshp[2],
       (CCTK_REAL8*)storage(), t,
       dstshp[0], dstshp[1], dstshp[2],
       srcbbox,
       dstbbox,
       regbbox);
  } else {
    abort();
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
template class data<T,3>;

#include "instantiate"

#undef INSTANTIATE

#endif
