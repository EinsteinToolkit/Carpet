// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/data.cc,v 1.33 2003/11/05 16:18:39 schnetter Exp $

#include <assert.h>
#include <limits.h>

#include <iostream>
#include <string>

#include <mpi.h>

#include "cctk.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "data.hh"

using namespace std;



// Hand out the next MPI tag
static int nexttag ()
{
  static int last = 100;
  ++last;
  if (last > 30000) last = 100;
  return last;
}



// Constructors
template<class T, int D>
data<T,D>::data (const int varindex_)
  : gdata<D>(varindex_),
    _storage(0),
    comm_active(false),
    tag(nexttag())
{ }

template<class T, int D>
data<T,D>::data (const int varindex_, const ibbox& extent_, const int proc_)
  : gdata<D>(varindex_),
    _storage(0),
    comm_active(false),
    tag(nexttag())
{
  allocate(extent_, proc_);
}

// Destructors
template<class T, int D>
data<T,D>::~data () {
  free();
}
  
// Pseudo constructors
template<class T, int D>
data<T,D>* data<T,D>::make_typed (const int varindex_) const {
  return new data(varindex_);
}



// Storage management
template<class T, int D>
void data<T,D>::allocate (const ibbox& extent_, const int proc_,
			  void* const mem) {
  assert (!this->_has_storage);
  this->_has_storage = true;
  // data
  this->_extent = extent_;
  this->_shape = max(ivect(0), this->_extent.shape() / this->_extent.stride());
  this->_size = 1;
  for (int d=0; d<D; ++d) {
    this->_stride[d] = this->_size;
    assert (this->_shape[d]==0 || this->_size <= INT_MAX / this->_shape[d]);
    this->_size *= this->_shape[d];
  }
  this->_proc = proc_;
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  if (rank==this->_proc) {
    this->_owns_storage = !mem;
    if (this->_owns_storage) {
      this->_storage = new T[this->_size];
    } else {
      this->_storage = (T*)mem;
    }
  } else {
    assert (!mem);
  }
}

template<class T, int D>
void data<T,D>::free () {
  if (this->_storage && this->_owns_storage) delete [] _storage;
  _storage = 0;
  this->_has_storage = false;
}

template<class T, int D>
void data<T,D>::transfer_from (gdata<D>* gsrc) {
  data* src = (data*)gsrc;
  assert (!_storage);
  *this = *src;
  *src = data(varindex);
}



// Processor management
template<class T, int D>
void data<T,D>::change_processor (comm_state<D>& state,
                                  const int newproc, void* const mem)
{
  switch (state.thestate) {
  case state_recv:
    change_processor_recv (newproc, mem);
    break;
  case state_send:
    change_processor_send (newproc, mem);
    break;
  case state_wait:
    change_processor_wait (newproc, mem);
    break;
  default:
    assert(0);
  }
}



template<class T, int D>
void data<T,D>::change_processor_recv (const int newproc, void* const mem)
{
  assert (!comm_active);
  comm_active = true;
  
  if (newproc == this->_proc) {
    assert (!mem);
    return;
  }
  
  if (this->_has_storage) {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == newproc) {
      // copy from other processor
      
      assert (!_storage);
      this->_owns_storage = !mem;
      if (this->_owns_storage) {
	_storage = new T[this->_size];
      } else {
	_storage = (T*)mem;
      }
      
      T dummy;
      MPI_Irecv (_storage, this->_size, dist::datatype(dummy), this->_proc,
                 this->tag, dist::comm, &request);
      
    } else if (rank == this->_proc) {
      // copy to other processor
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }
}



template<class T, int D>
void data<T,D>::change_processor_send (const int newproc, void* const mem)
{
  assert (comm_active);
  
  if (newproc == this->_proc) {
    assert (!mem);
    return;
  }
  
  if (this->_has_storage) {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == newproc) {
      // copy from other processor
      
    } else if (rank == this->_proc) {
      // copy to other processor
      
      assert (!mem);
      assert (_storage);
      
      T dummy;
      MPI_Isend (_storage, this->_size, dist::datatype(dummy), newproc,
                 this->tag, dist::comm, &request);
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }
}



template<class T, int D>
void data<T,D>::change_processor_wait (const int newproc, void* const mem)
{
  assert (comm_active);
  comm_active = false;
  
  if (newproc == this->_proc) {
    assert (!mem);
    return;
  }
  
  if (this->_has_storage) {
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == newproc) {
      // copy from other processor
      
      MPI_Status status;
      MPI_Wait (&request, &status);
      
    } else if (rank == this->_proc) {
      // copy to other processor
      
      assert (!mem);
      assert (_storage);
      
      MPI_Status status;
      MPI_Wait (&request, &status);
      
      if (this->_owns_storage) {
	delete [] _storage;
      }
      _storage = 0;
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }
  
  this->_proc = newproc;
}



// Data manipulators
template<class T, int D>
void data<T,D>
::copy_from_innerloop (const gdata<D>* gsrc, const ibbox& box)
{
  const data* src = (const data*)gsrc;
  assert (this->has_storage() && src->has_storage());
  assert (all(box.lower()>=this->extent().lower()
	      && box.lower()>=src->extent().lower()));
  assert (all(box.upper()<=this->extent().upper()
	      && box.upper()<=src->extent().upper()));
  assert (all(box.stride()==this->extent().stride()
	      && box.stride()==src->extent().stride()));
  assert (all((box.lower()-this->extent().lower())%box.stride() == 0
	      && (box.lower()-src->extent().lower())%box.stride() == 0));
  
  assert (this->proc() == src->proc());

  const int groupindex = CCTK_GroupIndexFromVarI(varindex);
  const int group_tags_table = CCTK_GroupTagsTableI(groupindex);
  assert (group_tags_table >= 0);
  
  // Disallow this.
  assert (0);
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == this->proc());
  
  for (typename ibbox::iterator it=box.begin(); it!=box.end(); ++it) {
    const ivect index = *it;
    (*this)[index] = (*src)[index];
  }
  
}



template<class T, int D>
void data<T,D>
::interpolate_from_innerloop (const vector<const gdata<D>*> gsrcs,
			      const vector<CCTK_REAL> times,
			      const ibbox& box, const CCTK_REAL time,
			      const int order_space,
			      const int order_time)
{
  assert (this->has_storage());
  assert (all(box.lower()>=this->extent().lower()));
  assert (all(box.upper()<=this->extent().upper()));
  assert (all(box.stride()==this->extent().stride()));
  assert (all((box.lower()-this->extent().lower())%box.stride() == 0));
  vector<const data*> srcs(gsrcs.size());
  for (int t=0; t<(int)srcs.size(); ++t) srcs[t] = (const data*)gsrcs[t];
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs[t]->has_storage());
    assert (all(box.lower()>=srcs[t]->extent().lower()));
    assert (all(box.upper()<=srcs[t]->extent().upper()));
    assert (this->proc() == srcs[t]->proc());
  }
  assert (order_space >= 0);
  assert (order_time >= 0);
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == this->proc());
  
  assert (varindex >= 0);
  const int groupindex = CCTK_GroupIndexFromVarI (varindex);
  assert (groupindex >= 0);
  char* groupname = CCTK_GroupName(groupindex);
  CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
	      "There is no interpolator available for the group \"%s\" with variable type %s, dimension %d, spatial interpolation order %d, temporal interpolation order %d.",
	      groupname, D, order_space, order_time);
  ::free (groupname);
}



extern "C" {
  void CCTK_FCALL CCTK_FNAME(copy_3d_int4)
    (const CCTK_INT4* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_INT4* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
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
void data<CCTK_INT4,3>
::copy_from_innerloop (const gdata<3>* gsrc, const ibbox& box)
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
    CCTK_FNAME(copy_3d_int4) ((const CCTK_INT4*)src->storage(),
                              srcshp[0], srcshp[1], srcshp[2],
                              (CCTK_INT4*)storage(),
                              dstshp[0], dstshp[1], dstshp[2],
                              srcbbox,
                              dstbbox,
                              regbbox);
    
  } else {
    assert (0);
  }
}

template<>
void data<CCTK_REAL8,3>
::copy_from_innerloop (const gdata<3>* gsrc, const ibbox& box)
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
    assert (0);
  }
}



extern "C" {
  
  void CCTK_FCALL CCTK_FNAME(restrict_3d_real8)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(restrict_3d_real8_rf2)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  
  
  
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_rf2)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_o3)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_o3_rf2)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_minmod)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_o5)
    (const CCTK_REAL8* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl_o3)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl_minmod)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl_o5)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const CCTK_REAL8* src3, const CCTK_REAL8& t3,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl_rf2)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const CCTK_REAL8* src3, const CCTK_REAL8& t3,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl_o3)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const CCTK_REAL8* src3, const CCTK_REAL8& t3,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl_o3_rf2)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const CCTK_REAL8* src3, const CCTK_REAL8& t3,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl_minmod)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const CCTK_REAL8* src3, const CCTK_REAL8& t3,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl_o5)
    (const CCTK_REAL8* src1, const CCTK_REAL8& t1,
     const CCTK_REAL8* src2, const CCTK_REAL8& t2,
     const CCTK_REAL8* src3, const CCTK_REAL8& t3,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_REAL8* dst, const CCTK_REAL8& t,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
}

template<>
void data<CCTK_REAL8,3>
::interpolate_from_innerloop (const vector<const gdata<3>*> gsrcs,
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
  vector<const data*> srcs(gsrcs.size());
  for (int t=0; t<(int)srcs.size(); ++t) srcs[t] = (const data*)gsrcs[t];
  assert (srcs.size() == times.size() && srcs.size()>0);
  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs[t]->has_storage());
    assert (all(box.lower()>=srcs[t]->extent().lower()));
    assert (all(box.upper()<=srcs[t]->extent().upper()));
  }
  
  assert (proc() == srcs[0]->proc());
  
  int rank;
  MPI_Comm_rank (dist::comm, &rank);
  assert (rank == proc());
  
  const ibbox& sext = srcs[0]->extent();
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
  if (all(sext.stride() < dext.stride())) {
    
    switch (transport_operator) {
      
    case op_Lagrange:
    case op_TVD:
      assert (srcs.size() == 1);
      if (all (dext.stride() == sext.stride() * 2)) {
        CCTK_FNAME(restrict_3d_real8_rf2)
          ((const CCTK_REAL8*)srcs[0]->storage(),
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(),
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
      } else {
        CCTK_FNAME(restrict_3d_real8)
          ((const CCTK_REAL8*)srcs[0]->storage(),
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(),
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
      }
      break;
      
    default:
      assert (0);
      
    } // switch (prolong_method)
    
  } else if (all(sext.stride() > dext.stride())) {
    
    switch (transport_operator) {
      
    case op_Lagrange:
      switch (order_time) {
        
      case 0:
        assert (srcs.size()>=1);
        switch (order_space) {
        case 0:
        case 1:
          if (all (sext.stride() == dext.stride() * 2)) {
            CCTK_FNAME(prolongate_3d_real8_rf2)
              ((const CCTK_REAL8*)srcs[0]->storage(),
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(),
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          } else {
            CCTK_FNAME(prolongate_3d_real8)
              ((const CCTK_REAL8*)srcs[0]->storage(),
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(),
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          }
          break;
        case 2:
        case 3:
          if (all (sext.stride() == dext.stride() * 2)) {
            CCTK_FNAME(prolongate_3d_real8_o3_rf2)
              ((const CCTK_REAL8*)srcs[0]->storage(),
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(),
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          } else {
            CCTK_FNAME(prolongate_3d_real8_o3)
              ((const CCTK_REAL8*)srcs[0]->storage(),
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(),
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          }
          break;
        case 4:
        case 5:
          CCTK_FNAME(prolongate_3d_real8_o5)
            ((const CCTK_REAL8*)srcs[0]->storage(),
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(),
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
          break;
        default:
          assert (0);
        }
        break;
        
      case 1:
        assert (srcs.size()>=2);
        switch (order_space) {
        case 0:
        case 1:
          CCTK_FNAME(prolongate_3d_real8_2tl)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
          break;
        case 2:
        case 3:
          CCTK_FNAME(prolongate_3d_real8_2tl_o3)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
          break;
        case 4:
        case 5:
          CCTK_FNAME(prolongate_3d_real8_2tl_o5)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
          break;
        default:
          assert (0);
        }
        break;
        
      case 2:
        assert (srcs.size()>=3);
        switch (order_space) {
        case 0:
        case 1:
          if (all (sext.stride() == dext.stride() * 2)) {
            CCTK_FNAME(prolongate_3d_real8_3tl_rf2)
              ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
               (const CCTK_REAL8*)srcs[1]->storage(), times[1],
               (const CCTK_REAL8*)srcs[2]->storage(), times[2],
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(), time,
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          } else {
            CCTK_FNAME(prolongate_3d_real8_3tl)
              ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
               (const CCTK_REAL8*)srcs[1]->storage(), times[1],
               (const CCTK_REAL8*)srcs[2]->storage(), times[2],
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(), time,
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          }
          break;
        case 2:
        case 3:
          if (all (sext.stride() == dext.stride() * 2)) {
            CCTK_FNAME(prolongate_3d_real8_3tl_o3_rf2)
              ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
               (const CCTK_REAL8*)srcs[1]->storage(), times[1],
               (const CCTK_REAL8*)srcs[2]->storage(), times[2],
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(), time,
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          } else {
            CCTK_FNAME(prolongate_3d_real8_3tl_o3)
              ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
               (const CCTK_REAL8*)srcs[1]->storage(), times[1],
               (const CCTK_REAL8*)srcs[2]->storage(), times[2],
               srcshp[0], srcshp[1], srcshp[2],
               (CCTK_REAL8*)storage(), time,
               dstshp[0], dstshp[1], dstshp[2],
               srcbbox, dstbbox, regbbox);
          }
          break;
        case 4:
        case 5:
          CCTK_FNAME(prolongate_3d_real8_3tl_o5)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             (const CCTK_REAL8*)srcs[2]->storage(), times[2],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
          break;
        default:
          assert (0);
        }
        break;
        
      default:
        assert (0);
      } // switch (order_time)
      break;
      
    case op_TVD:
      switch (order_time) {
      case 0: 
        CCTK_FNAME(prolongate_3d_real8_minmod)
          ((const CCTK_REAL8*)srcs[0]->storage(),
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(),
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
        break;
      case 1:
        CCTK_FNAME(prolongate_3d_real8_2tl_minmod)
          ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
           (const CCTK_REAL8*)srcs[1]->storage(), times[1],
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(), time,
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
      case 2: 
        CCTK_FNAME(prolongate_3d_real8_3tl_minmod)
          ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
           (const CCTK_REAL8*)srcs[1]->storage(), times[1],
           (const CCTK_REAL8*)srcs[2]->storage(), times[2],
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(), time,
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
        break;
      default:
        assert (0);
      }
      break;
      
    default:
      assert(0);
    } // switch (prolong_method)
    
  } else {
    assert (0);
  }
}

  



// Output
template<class T,int D>
ostream& data<T,D>::output (ostream& os) const {
  T Tdummy;
  os << "data<" << typestring(Tdummy) << "," << D << ">:"
     << "extent=" << this->extent() << ","
     << "stride=" << this->stride() << ",size=" << this->size();
  return os;
}



#define INSTANTIATE(T)				\
template class data<T,3>;

#include "instantiate"

#undef INSTANTIATE
