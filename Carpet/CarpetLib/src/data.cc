#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"
#include "vect.hh"

#include "data.hh"

using namespace std;



// Total number of currently allocated bytes and objects
static size_t total_allocated_bytes = 0;
static size_t total_allocated_objects = 0;



static const CCTK_REAL eps = 1.0e-10;

// Constructors
template<typename T>
data<T>::data (const int varindex_, const operator_type transport_operator_,
               const int vectorlength_, const int vectorindex_,
               data* const vectorleader_,
               const int tag_)
  : gdata(varindex_, transport_operator_, tag_),
    _storage(NULL), _allocated_bytes(0),
    vectorlength(vectorlength_), vectorindex(vectorindex_),
    vectorleader(vectorleader_)
{
  assert (vectorlength>=1);
  assert (vectorindex>=0 && vectorindex<vectorlength);
  assert ((vectorindex==0 && !vectorleader)
          || (vectorindex!=0 && vectorleader));
  if (vectorindex==0) vectorclients.resize (vectorlength);
  if (vectorindex!=0) vectorleader->register_client (vectorindex);
}

template<typename T>
data<T>::data (const int varindex_, const operator_type transport_operator_,
               const int vectorlength_, const int vectorindex_,
               data* const vectorleader_,
               const ibbox& extent_, const int proc_)
  : gdata(varindex_, transport_operator_),
    _storage(NULL), _allocated_bytes(0),
    vectorlength(vectorlength_), vectorindex(vectorindex_),
    vectorleader(vectorleader_)
{
  assert (vectorlength>=1);
  assert (vectorindex>=0 && vectorindex<vectorlength);
  assert ((vectorindex==0 && !vectorleader)
          || (vectorindex!=0 && vectorleader));
  if (vectorindex==0) vectorclients.resize (vectorlength);
  if (vectorindex!=0) vectorleader->register_client (vectorindex);
  allocate(extent_, proc_);
}

// Destructors
template<typename T>
data<T>::~data ()
{
  if (vectorleader) vectorleader->unregister_client (vectorindex);
  if (vectorindex==0) assert (! has_clients());
  free();
}
  
// Pseudo constructors
template<typename T>
data<T>* data<T>::make_typed (const int varindex_,
                              const operator_type transport_operator_,
                              const int tag_)
  const
{
  return new data(varindex_, transport_operator_, 1, 0, NULL, tag_);
}



// Vector mamagement
template<typename T>
void data<T>::register_client (const int index)
{
  assert (! vectorclients.at(index));
  vectorclients.at(index) = true;
}

template<typename T>
void data<T>::unregister_client (const int index)
{
  assert (vectorclients.at(index));
  vectorclients.at(index) = false;
}

template<typename T>
bool data<T>::has_clients () const
{
  return (find (vectorclients.begin(), vectorclients.end(), true)
          != vectorclients.end());
}



// Storage management
template<typename T>
void data<T>::getmem (const size_t nelems)
{
  const size_t nbytes = nelems * sizeof(T);
  try {
    assert (_allocated_bytes == 0);
    _storage = new T[nelems];
    _allocated_bytes = nbytes;
  } catch (...) {
    T Tdummy;
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to allocate %.0f bytes (%.3f MB) of memory for type %s.  %.0f bytes (%.3f MB) are currently allocated in %d objects",
                double(nbytes), double(nbytes/1.0e6),
                typestring(Tdummy),
                double(total_allocated_bytes),
                double(total_allocated_bytes/1.0e6),
                int(total_allocated_objects));
  }
  total_allocated_bytes += nbytes;
  ++ total_allocated_objects;
}



template<typename T>
void data<T>::freemem ()
{
  delete [] _storage;
  assert (total_allocated_bytes > _allocated_bytes);
  assert (total_allocated_objects > 0);
  total_allocated_bytes -= _allocated_bytes;
  -- total_allocated_objects;
  _allocated_bytes = 0;
}



template<typename T>
void data<T>::allocate (const ibbox& extent_,
                        const int proc_,
                        void* const mem)
{
  assert (!_has_storage);
  _has_storage = true;
  // prevent accidental wrap-around
  assert (all(extent_.lower() < numeric_limits<int>::max() / 2));
  assert (all(extent_.lower() > numeric_limits<int>::min() / 2));
  assert (all(extent_.upper() < numeric_limits<int>::max() / 2));
  assert (all(extent_.upper() > numeric_limits<int>::min() / 2));
  // data
  _extent = extent_;
  _shape = max(ivect(0), _extent.shape() / _extent.stride());
  _size = 1;
  for (int d=0; d<dim; ++d) {
    _stride[d] = _size;
    assert (_shape[d]==0 || _size <= INT_MAX / _shape[d]);
    _size *= _shape[d];
  }
  _proc = proc_;
  if (dist::rank() == _proc) {
    _owns_storage = !mem;
    if (_owns_storage) {
      if (vectorindex == 0) {
        assert (! vectorleader);
        getmem (vectorlength * _size);
      } else {
        assert (vectorleader);
        _storage = vectorleader->vectordata (vectorindex);
      }
    } else {
      _storage = (T*)mem;
    }
  } else {
    assert (!mem);
  }
}

template<typename T>
void data<T>::free ()
{
  if (_storage && _owns_storage && vectorindex==0) {
    freemem ();
  }
  _storage = 0;
  _has_storage = false;
}

template<typename T>
void data<T>::transfer_from (gdata* gsrc)
{
  assert (vectorlength==1);
  data* src = (data*)gsrc;
  assert (src->vectorlength==1);
  assert (!_storage);
  *this = *src;
  *src = data(varindex, transport_operator);
}

template<typename T>
T* data<T>::vectordata (const int vectorindex_) const
{
  assert (vectorindex==0);
  assert (! vectorleader);
  assert (vectorindex_>=0 && vectorindex_<vectorlength);
  assert (_storage && _owns_storage);
  return _storage + vectorindex_ * _size;
}



// Processor management
template<typename T>
void data<T>::change_processor_recv (comm_state& state,
                                     const int newproc,
                                     void* const mem)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (!comm_active);
  comm_active = true;
  
  if (newproc == _proc) {
    assert (!mem);
    return;
  }
  
  wtime_changeproc_recv.start();

  if (_has_storage) {
    if (dist::rank() == newproc) {
      // copy from other processor
      
      assert (!_storage);
      _owns_storage = !mem;
      if (_owns_storage) {
        getmem (_size);
      } else {
	_storage = (T*)mem;
      }
      
      wtime_irecv.start();
      T dummy;
      MPI_Irecv (_storage, _size, dist::datatype(dummy), proc(),
                 tag, dist::comm, &request);
      wtime_irecv.stop();
      if (use_waitall) {
        state.requests.push_back (request);
      }
      
    } else if (dist::rank() == _proc) {
      // copy to other processor
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }
  
  wtime_changeproc_recv.stop();
}



template<typename T>
void data<T>::change_processor_send (comm_state& state,
                                     const int newproc,
                                     void* const mem)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (comm_active);
  
  if (newproc == _proc) {
    assert (!mem);
    return;
  }
  
  wtime_changeproc_send.start();
  
  if (_has_storage) {
    if (dist::rank() == newproc) {
      // copy from other processor
      
    } else if (dist::rank() == _proc) {
      // copy to other processor
      
      assert (!mem);
      assert (_storage);
      
      wtime_isend.start();
      T dummy;
      MPI_Isend (_storage, _size, dist::datatype(dummy), newproc,
                 tag, dist::comm, &request);
      wtime_isend.stop();
      if (use_waitall) {
        state.requests.push_back (request);
      }
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }
  
  wtime_changeproc_send.stop();
}



template<typename T>
void data<T>::change_processor_wait (comm_state& state,
                                     const int newproc,
                                     void* const mem)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (comm_active);
  comm_active = false;
  
  if (newproc == _proc) {
    assert (!mem);
    return;
  }
  
  wtime_changeproc_wait.start();

  if (use_waitall) {
    if (! state.requests.empty()) {
      // wait for all requests at once
      wtime_irecvwait.start();
      MPI_Waitall
        (state.requests.size(), &state.requests.front(), MPI_STATUSES_IGNORE);
      wtime_irecvwait.stop();
      state.requests.clear();
    }
  }
  
  if (_has_storage) {
    if (dist::rank() == newproc) {
      // copy from other processor
      
      if (! use_waitall) {
        wtime_irecvwait.start();
        MPI_Wait (&request, MPI_STATUS_IGNORE);
        wtime_irecvwait.stop();
      }
      
    } else if (dist::rank() == _proc) {
      // copy to other processor
      
      assert (!mem);
      assert (_storage);
      
      if (! use_waitall) {
        wtime_isendwait.start();
        MPI_Wait (&request, MPI_STATUS_IGNORE);
        wtime_isendwait.stop();
      }
      
      if (_owns_storage) {
	freemem ();
      }
      _storage = 0;
      
    } else {
      assert (!mem);
      assert (!_storage);
    }
  }
  
  _proc = newproc;
  
  wtime_changeproc_wait.stop();
}



#if 0
template<typename T>
void
data<T>::copy_from_recv_inner (comm_state& state,
                               const gdata* gsrc, const ibbox& box)
{
  DECLARE_CCTK_PARAMETERS;
  
  wtime_copyfrom_recvinner_allocate.start();
  comm_state::commbuf<T> * b = new comm_state::commbuf<T>;
  b->am_receiver = true;
  b->am_sender = false;
  b->data.resize (prod (box.shape() / box.stride()));
  wtime_copyfrom_recvinner_allocate.stop();
  
  wtime_copyfrom_recvinner_recv.start();
  assert (dist::rank() == proc());
  T dummy;
  MPI_Irecv (&b->data.front(), b->data.size(),
             dist::datatype(dummy), gsrc->proc(),
             tag, dist::comm, &b->request);
  wtime_copyfrom_recvinner_recv.stop();
  if (use_waitall) {
    state.requests.push_back (b->request);
  }
  state.recvbufs.push (b);
}
#endif


#if 0
template<typename T>
void
data<T>::copy_from_send_inner (comm_state& state,
                               const gdata* gsrc, const ibbox& box)
{
  DECLARE_CCTK_PARAMETERS;
  
  wtime_copyfrom_sendinner_allocate.start();
  comm_state::gcommbuf * b = gsrc->make_typed_commbuf (box);
  b->am_receiver = false;
  b->am_sender = true;
  wtime_copyfrom_sendinner_allocate.stop();
  
  wtime_copyfrom_sendinner_copy.start();
  const data<T> * src = dynamic_cast<const data<T> *> (gsrc);
  assert (src->_has_storage);
  assert (src->_owns_storage);
  // copy src to b
#if 0
  {
    T * restrict p = & b->data.front();
    T const * restrict const q = src->_storage;
    ivect const imin = box.lower() / box.stride();
    ivect const imax = (box.upper() + box.stride()) / box.stride();
    ivect const lbnd = src->extent().lower() / src->extent().stride();
    ivect const lsh = src->extent().shape() / src->extent().stride();
    for (int k=imin[2]; k<imax[2]; ++k) {
      for (int j=imin[1]; j<imax[1]; ++j) {
        for (int i=imin[0]; i<imax[0]; ++i) {
          * p ++ = q [i - lbnd[0] + lsh[0] * (j - lbnd[1] + lsh[1] * (k - lbnd[2]))];
        }
      }
    }
  }
#endif
  {
    data<T> * tmp = src->make_typed (varindex, transport_operator, tag);
    tmp->allocate (box, src->proc(), &b->data.front());
    tmp->copy_from_innerloop (src, box);
    delete tmp;
  }
  wtime_copyfrom_sendinner_copy.stop();
  
  wtime_copyfrom_sendinner_send.start();
  assert (dist::rank() == src->proc());
  T dummy;
  MPI_Isend (b->pointer(), b->size(), b->datatype(), proc(),
             tag, dist::comm, &b->request);
  wtime_copyfrom_sendinner_send.stop();
  if (use_waitall) {
    state.requests.push_back (b->request);
  }
  state.sendbufs.push (b);
}
#endif



#if 0
template<typename T>
void
data<T>::copy_from_recv_wait_inner (comm_state& state,
                                    const gdata* gsrc, const ibbox& box)
{
  DECLARE_CCTK_PARAMETERS;
  
  comm_state::commbuf<T> * b
    = (comm_state::commbuf<T> *) state.recvbufs.front();
  state.recvbufs.pop();
  assert (b->am_receiver);
  assert (! b->am_sender);
  
  wtime_copyfrom_recvwaitinner_wait.start();
  if (use_waitall) {
    if (! state.requests.empty()) {
      // wait for all requests at once
      MPI_Waitall
        (state.requests.size(), &state.requests.front(), MPI_STATUSES_IGNORE);
      state.requests.clear();
    }
  }
  
  if (! use_waitall) {
    MPI_Wait (&b->request, MPI_STATUS_IGNORE);
  }
  wtime_copyfrom_recvwaitinner_wait.stop();
  
  wtime_copyfrom_recvwaitinner_copy.start();
  assert (_has_storage);
  assert (_owns_storage);
  // copy b to this
  {
    T * restrict const p = _storage;
    T const * restrict q = & b->data.front();
    ivect const imin = box.lower() / box.stride();
    ivect const imax = (box.upper() + box.stride()) / box.stride();
    ivect const lbnd = extent().lower() / extent().stride();
    ivect const lsh = extent().shape() / extent().stride();
    for (int k=imin[2]; k<imax[2]; ++k) {
      for (int j=imin[1]; j<imax[1]; ++j) {
        for (int i=imin[0]; i<imax[0]; ++i) {
          p [i - lbnd[0] + lsh[0] * (j - lbnd[1] + lsh[1] * (k - lbnd[2]))] = * q ++;
        }
      }
    }
  }
  wtime_copyfrom_recvwaitinner_copy.stop();
  
  wtime_copyfrom_recvwaitinner_delete.start();
  delete b;
  wtime_copyfrom_recvwaitinner_delete.stop();
}
#endif



#if 0
template<typename T>
void
data<T>::copy_from_send_wait_inner (comm_state& state,
                                    const gdata* gsrc, const ibbox& box)
{
  DECLARE_CCTK_PARAMETERS;
  
  comm_state::commbuf<T> * b
    = (comm_state::commbuf<T> *) state.sendbufs.front();
  state.sendbufs.pop();
  assert (! b->am_receiver);
  assert (b->am_sender);
  
  wtime_copyfrom_sendwaitinner_wait.start();
  if (use_waitall) {
    if (! state.requests.empty()) {
      // wait for all requests at once
      MPI_Waitall
        (state.requests.size(), &state.requests.front(), MPI_STATUSES_IGNORE);
      state.requests.clear();
    }
  }
  
  if (! use_waitall) {
    MPI_Wait (&b->request, MPI_STATUS_IGNORE);
  }
  wtime_copyfrom_sendwaitinner_wait.stop();
  
  wtime_copyfrom_sendwaitinner_delete.start();
  delete b;
  wtime_copyfrom_sendwaitinner_delete.stop();
}
#endif



// Data manipulators
template<typename T>
comm_state::gcommbuf *
data<T>::
make_typed_commbuf (const ibbox & box)
  const
{
  return new comm_state::commbuf<T> (box);
}



template<typename T>
void data<T>
::copy_from_innerloop (const gdata* gsrc, const ibbox& box)
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

  const int groupindex = CCTK_GroupIndexFromVarI(varindex);
  const int group_tags_table = CCTK_GroupTagsTableI(groupindex);
  assert (group_tags_table >= 0);
  
  // Disallow this.
  T Tdummy;
  CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
              "There is no copy operator available for the variable type %s",
              typestring(Tdummy));
  
  assert (dist::rank() == proc());
  
  for (typename ibbox::iterator it=box.begin(); it!=box.end(); ++it) {
    const ivect index = *it;
    (*this)[index] = (*src)[index];
  }
  
}



template<class T>
void data<T>
::fill_bbox_arrays (int srcshp[dim],
                    int dstshp[dim],
                    int srcbbox[dim][dim],
                    int dstbbox[dim][dim],
                    int regbbox[dim][dim],
                    const ibbox & box,
                    const ibbox & sext,
                    const ibbox & dext)
{
  for (int d=0; d<dim; ++d) {
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
}

template<typename T>
void data<T>
::interpolate_from_innerloop (const vector<const gdata*> gsrcs,
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
    assert (proc() == srcs[t]->proc());
  }
  assert (order_space >= 0);
  assert (order_time >= 0);
  
  assert (dist::rank() == proc());
  
  assert (varindex >= 0);
  const int groupindex = CCTK_GroupIndexFromVarI (varindex);
  assert (groupindex >= 0);
  char* groupname = CCTK_GroupName(groupindex);
  T Tdummy;
  CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
	      "There is no interpolator available for the group \"%s\" with variable type %s, spatial interpolation order %d, temporal interpolation order %d",
	      groupname, typestring(Tdummy), order_space, order_time);
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
  void CCTK_FCALL CCTK_FNAME(copy_3d_complex16)
    (const CCTK_COMPLEX16* src,
     const int& srciext, const int& srcjext, const int& srckext,
     CCTK_COMPLEX16* dst,
     const int& dstiext, const int& dstjext, const int& dstkext,
     const int srcbbox[3][3],
     const int dstbbox[3][3],
     const int regbbox[3][3]);
}

template<>
void data<CCTK_INT4>
::copy_from_innerloop (const gdata* gsrc, const ibbox& box)
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
  
  assert (dist::rank() == proc());
  
  const ibbox& sext = src->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  fill_bbox_arrays( srcshp, dstshp, srcbbox, dstbbox, regbbox,
	box, sext, dext );
  
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
void data<CCTK_REAL8>
::copy_from_innerloop (const gdata* gsrc, const ibbox& box)
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
  
  assert (dist::rank() == proc());
  
  const ibbox& sext = src->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  fill_bbox_arrays( srcshp, dstshp, srcbbox, dstbbox, regbbox,
	box, sext, dext );
  
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

template<>
void data<CCTK_COMPLEX16>
::copy_from_innerloop (const gdata* gsrc, const ibbox& box)
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
  
  assert (dist::rank() == proc());
  
  const ibbox& sext = src->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  fill_bbox_arrays( srcshp, dstshp, srcbbox, dstbbox, regbbox,
	box, sext, dext );
  
  assert (all(dext.stride() == box.stride()));
  if (all(sext.stride() == dext.stride())) {
    CCTK_FNAME(copy_3d_complex16) ((const CCTK_COMPLEX16*)src->storage(),
                                   srcshp[0], srcshp[1], srcshp[2],
                                   (CCTK_COMPLEX16*)storage(),
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
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_eno)
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
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl_rf2)
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
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl_o3_rf2)
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
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_2tl_eno)
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
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_3tl_eno)
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

template<typename T>
void data<T>
::interpolate_restrict (const vector<const data<T>*> & srcs,
                        const vector<CCTK_REAL> & times,
                        const ibbox& box)
{
  const ibbox& sext = srcs[0]->extent();
  const ibbox& dext = extent();
  
  int srcshp[3], dstshp[3];
  int srcbbox[3][3], dstbbox[3][3], regbbox[3][3];
  
  fill_bbox_arrays (srcshp, dstshp, srcbbox, dstbbox, regbbox,
                    box, sext, dext );
  
  switch (transport_operator) {
      
  case op_Lagrange:
  case op_TVD:
  case op_ENO:
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
  }
}

template<typename T>
void data<T>
::interpolate_prolongate (const vector<const data<T>*> & srcs,
                          const vector<CCTK_REAL> & times,
                          const ibbox& box, const CCTK_REAL time,
                          const int order_space,
                          const int order_time)
{
  const ibbox& sext = srcs[0]->extent();
  const ibbox& dext = extent();
  
  int srcshp[dim], dstshp[dim];
  int srcbbox[dim][dim], dstbbox[dim][dim], regbbox[dim][dim];
  
  fill_bbox_arrays (srcshp, dstshp, srcbbox, dstbbox, regbbox,
                    box, sext, dext);
  switch (transport_operator) {
    
  case op_Lagrange:
    switch (order_time) {
      
    case 0:
      assert (times.size() == 1);
      assert (abs(times[0] - time) < eps);
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
        if (all (sext.stride() == dext.stride() * 2)) {
          CCTK_FNAME(prolongate_3d_real8_2tl_rf2)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
        } else {
          CCTK_FNAME(prolongate_3d_real8_2tl)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
        }
        break;
      case 2:
      case 3:
        if (all (sext.stride() == dext.stride() * 2)) {
          CCTK_FNAME(prolongate_3d_real8_2tl_o3_rf2)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
        } else {
          CCTK_FNAME(prolongate_3d_real8_2tl_o3)
            ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
             (const CCTK_REAL8*)srcs[1]->storage(), times[1],
             srcshp[0], srcshp[1], srcshp[2],
             (CCTK_REAL8*)storage(), time,
             dstshp[0], dstshp[1], dstshp[2],
             srcbbox, dstbbox, regbbox);
        }
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
      assert (times.size() == 1);
      assert (abs(times[0] - time) < eps);
      switch (order_space) {
      case 0:
      case 1:
        CCTK_WARN (0, "There is no stencil for op=\"TVD\" with order_space=1");
        break;
      case 2:
      case 3:
        CCTK_FNAME(prolongate_3d_real8_rf2)
          ((const CCTK_REAL8*)srcs[0]->storage(),
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(),
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
//            CCTK_FNAME(prolongate_3d_real8_minmod)
//              ((const CCTK_REAL8*)srcs[0]->storage(),
//               srcshp[0], srcshp[1], srcshp[2],
//               (CCTK_REAL8*)storage(),
//               dstshp[0], dstshp[1], dstshp[2],
//               srcbbox, dstbbox, regbbox);
        break;
      default:
        assert (0);
      }
      break;
    case 1:
      switch (order_space) {
      case 0:
      case 1:
        CCTK_WARN (0, "There is no stencil for op=\"TVD\" with order_space=1");
        break;
      case 2:
      case 3:
        CCTK_FNAME(prolongate_3d_real8_2tl_rf2)
          ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
           (const CCTK_REAL8*)srcs[1]->storage(), times[1],
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(), time,
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
//            CCTK_FNAME(prolongate_3d_real8_2tl_minmod)
//              ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
//               (const CCTK_REAL8*)srcs[1]->storage(), times[1],
//               srcshp[0], srcshp[1], srcshp[2],
//               (CCTK_REAL8*)storage(), time,
//               dstshp[0], dstshp[1], dstshp[2],
//               srcbbox, dstbbox, regbbox);
        break;
      default:
        assert (0);
      }
      break;
    case 2: 
      switch (order_space) {
      case 0:
      case 1:
        CCTK_WARN (0, "There is no stencil for op=\"TVD\" with order_space=1");
        break;
      case 2:
      case 3:
        CCTK_FNAME(prolongate_3d_real8_3tl_rf2)
          ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
           (const CCTK_REAL8*)srcs[1]->storage(), times[1],
           (const CCTK_REAL8*)srcs[2]->storage(), times[2],
           srcshp[0], srcshp[1], srcshp[2],
           (CCTK_REAL8*)storage(), time,
           dstshp[0], dstshp[1], dstshp[2],
           srcbbox, dstbbox, regbbox);
//            CCTK_FNAME(prolongate_3d_real8_3tl_minmod)
//              ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
//               (const CCTK_REAL8*)srcs[1]->storage(), times[1],
//               (const CCTK_REAL8*)srcs[2]->storage(), times[2],
//               srcshp[0], srcshp[1], srcshp[2],
//               (CCTK_REAL8*)storage(), time,
//               dstshp[0], dstshp[1], dstshp[2],
//               srcbbox, dstbbox, regbbox);
        break;
      default:
        assert (0);
      }
      break;
    default:
      assert (0);
    } // switch (order_time)
    break;
    
  case op_ENO:
    switch (order_time) {
    case 0: 
      assert (times.size() == 1);
      assert (abs(times[0] - time) < eps);
      switch (order_space) {
      case 0:
      case 1:
        CCTK_WARN (0, "There is no stencil for op=\"ENO\" with order_space=1");
        break;
      case 2:
      case 3:
        CCTK_FNAME(prolongate_3d_real8_eno)
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
      switch (order_space) {
      case 0:
      case 1:
        CCTK_WARN (0, "There is no stencil for op=\"ENO\" with order_space=1");
        break;
      case 2:
      case 3:
        CCTK_FNAME(prolongate_3d_real8_2tl_eno)
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
      switch (order_space) {
      case 0:
      case 1:
        CCTK_WARN (0, "There is no stencil for op=\"ENO\" with order_space=1");
        break;
      case 2:
      case 3:
        CCTK_FNAME(prolongate_3d_real8_3tl_eno)
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
    
  default:
    assert(0);
  } // switch (transport_operator)
}

template<>
void data<CCTK_REAL8>
::Check_that_the_times_are_consistent (const vector<CCTK_REAL> & times,
                                       const CCTK_REAL time)
{
  assert (times.size() > 0);
  CCTK_REAL min_time = times[0];
  CCTK_REAL max_time = times[0];
  for (size_t tl=1; tl<times.size(); ++tl) {
    min_time = min(min_time, times[tl]);
    max_time = max(max_time, times[tl]);
  }
  if (time < min_time - eps || time > max_time + eps) {
    ostringstream buf;
    buf << "Internal error: extrapolation in time."
        << "  time=" << time
        << "  times=" << times;
    CCTK_WARN (0, buf.str().c_str());
  }
}

template<>
bool data<CCTK_REAL8>
::try_without_time_interpolation (const vector<const gdata*> & gsrcs,
                                  const vector<CCTK_REAL> & times,
                                  const ibbox& box, const CCTK_REAL time,
                                  const int order_space,
                                  const int order_time)
{
  for (size_t tl=0; tl<times.size(); ++tl) {
    if (abs(times[tl] - time) < eps) {
      vector<const gdata*> my_gsrcs(1);
      vector<CCTK_REAL> my_times(1);
      my_gsrcs[0] = gsrcs[tl];
      my_times[0] = times[tl];
      const int my_order_time = 0;
      interpolate_from_innerloop
        (my_gsrcs, my_times, box, time, order_space, my_order_time);
      return true;
    }
  }
  return false;
}

template<>
void data<CCTK_REAL8>
::interpolate_from_innerloop (const vector<const gdata*> gsrcs,
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

  for (int t=0; t<(int)srcs.size(); ++t)
    srcs[t] = (const data*)gsrcs[t];

  assert (srcs.size() == times.size() && srcs.size()>0);

  for (int t=0; t<(int)srcs.size(); ++t) {
    assert (srcs[t]->has_storage());
    assert (all(box.lower()>=srcs[t]->extent().lower()));
    assert (all(box.upper()<=srcs[t]->extent().upper()));
  }
  
  assert (proc() == srcs[0]->proc());
  
  assert (dist::rank() == proc());
  
  Check_that_the_times_are_consistent (times, time);

  bool did_time_interpolation = false;

  if (times.size() > 1) {
    // try to avoid time interpolation if possible
    did_time_interpolation = 
      try_without_time_interpolation
      (gsrcs, times, box, time, order_space, order_time);
  }
  
  if (!did_time_interpolation) {
    const ibbox& sext = srcs[0]->extent();
    const ibbox& dext = extent();
  
    assert (all(dext.stride() == box.stride()));

    if (all(sext.stride() < dext.stride())) {
    
      assert (times.size() == 1);
      assert (abs(times[0] - time) < eps);

      interpolate_restrict (srcs, times, box);
    
    } else if (all(sext.stride() > dext.stride())) {
    
      interpolate_prolongate (srcs, times, box, time, order_space, order_time);
    
    } else {
      assert (0);
    }
  }
}

// Output
template<typename T>
ostream& data<T>::output (ostream& os) const
{
  T Tdummy;
  os << "data<" << typestring(Tdummy) << ">:"
     << "extent=" << extent() << ","
     << "stride=" << stride() << ",size=" << size();
  return os;
}



#define INSTANTIATE(T)				\
template class data<T>;

#include "instantiate"

#undef INSTANTIATE
