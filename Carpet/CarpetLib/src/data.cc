// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/data.cc,v 1.31 2003/10/15 07:14:01 hawke Exp $

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

#include "util_ErrorCodes.h"
#include "util_Table.h"

using namespace std;



// Constructors
template<class T, int D>
data<T,D>::data (const int varindex_)
  : gdata<D>(varindex_),
    _storage(0)
{ }

template<class T, int D>
data<T,D>::data (const int varindex_, const ibbox& extent_, const int proc_)
  : gdata<D>(varindex_),
    _storage(0)
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
void data<T,D>::change_processor (const int newproc, void* const mem) {
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
      MPI_Status status;
      MPI_Recv (_storage, this->_size, dist::datatype(dummy), this->_proc,
		dist::tag, dist::comm, &status);

    } else if (rank == this->_proc) {
      // copy to other processor
      
      assert (!mem);
      assert (_storage);
      T dummy;
      MPI_Send (_storage, this->_size, dist::datatype(dummy), newproc,
		dist::tag, dist::comm);

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
  
  T Tdummy;
  CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
	      "There is no interpolator available for variable type %s, dimension %d, spatial interpolation order %d, temporal interpolation order %d.  The interpolation will not be done.",
	      typestring(Tdummy), D, order_space, order_time);
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
  
  
  
  void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8)
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

  const int varindex = (gsrcs[0])->var_index();
  const int groupindex = CCTK_GroupIndexFromVarI(varindex);
  const int group_tags_table = CCTK_GroupTagsTableI(groupindex);
  assert(group_tags_table >= 0);
  
  int typecode = -1;
  int keysize = -1;
  if (! Util_TableQueryValueInfo(group_tags_table,
                                 &typecode,
                                 &keysize,
                                 "Prolongation"))
  {
    CCTK_VWarn(4, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Tags table for group %s does not contain 'Prolongation'"
               " tag", CCTK_GroupName(groupindex));
  }
  else
  {  
    assert(typecode == CCTK_VARIABLE_CHAR);
  }

#define PROLONG_NONE 0
#define PROLONG_LAGRANGE 1
#define PROLONG_TVD 2

  int prolong_method;
  if (keysize < 0) { // No key in table - default to Lagrange.
    prolong_method = 1;
  }
  else {
    char prolong_string[keysize+10];
    const int error = Util_TableGetString(group_tags_table,
                                          keysize+10,
                                          prolong_string,
                                          "Prolongation");
    if (error < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error code %d getting prolongation method"
                 " from tags table for group %s.",
                 error, CCTK_GroupName(groupindex));
    }
    if (CCTK_Equals(prolong_string, "None")){
      prolong_method = PROLONG_NONE;
    }
    else if (CCTK_Equals(prolong_string, "Lagrange")){
      prolong_method = PROLONG_LAGRANGE;
    }
    else if (CCTK_Equals(prolong_string, "TVD")){
      prolong_method = PROLONG_TVD;
    }
    else {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Prolongation method %s for group %s not recognized.\n"
                 "Is this is a typo? Known values are:\n"
                 "None\nLagrange\nTVD",
                 prolong_method, CCTK_GroupName(groupindex));
    }
      
  }

//   CCTK_VInfo(CCTK_THORNSTRING,
//              "Variable %d is %s with method %d", 
//              varindex, CCTK_VarName(varindex), prolong_method);
  
  assert (all(dext.stride() == box.stride()));
  if (all(sext.stride() < dext.stride())) {
    
    assert (srcs.size() == 1);
    CCTK_FNAME(restrict_3d_real8)
      ((const CCTK_REAL8*)srcs[0]->storage(),
       srcshp[0], srcshp[1], srcshp[2],
       (CCTK_REAL8*)storage(),
       dstshp[0], dstshp[1], dstshp[2],
       srcbbox, dstbbox, regbbox);
    
  } else if (all(sext.stride() > dext.stride())) {

    switch (prolong_method) {

      case PROLONG_NONE:  // PROLONG_NONE: NOTHING
        break;
      case PROLONG_LAGRANGE:  // PROLONG_LAGRANGE: DEFAULT
        switch (order_time) {
          
          case 0:
            assert (srcs.size()>=1);
            switch (order_space) {
              case 0:
              case 1:
                CCTK_FNAME(prolongate_3d_real8)
                  ((const CCTK_REAL8*)srcs[0]->storage(),
                   srcshp[0], srcshp[1], srcshp[2],
                   (CCTK_REAL8*)storage(),
                   dstshp[0], dstshp[1], dstshp[2],
                   srcbbox, dstbbox, regbbox);
                break;
              case 2:
              case 3:
                CCTK_FNAME(prolongate_3d_real8_o3)
                  ((const CCTK_REAL8*)srcs[0]->storage(),
                   srcshp[0], srcshp[1], srcshp[2],
                   (CCTK_REAL8*)storage(),
                   dstshp[0], dstshp[1], dstshp[2],
                   srcbbox, dstbbox, regbbox);
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
                CCTK_FNAME(prolongate_3d_real8_3tl)
                  ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
                   (const CCTK_REAL8*)srcs[1]->storage(), times[1],
                   (const CCTK_REAL8*)srcs[2]->storage(), times[2],
                   srcshp[0], srcshp[1], srcshp[2],
                   (CCTK_REAL8*)storage(), time,
                   dstshp[0], dstshp[1], dstshp[2],
                   srcbbox, dstbbox, regbbox);
                break;
              case 2:
              case 3:
                CCTK_FNAME(prolongate_3d_real8_3tl_o3)
                  ((const CCTK_REAL8*)srcs[0]->storage(), times[0],
                   (const CCTK_REAL8*)srcs[1]->storage(), times[1],
                   (const CCTK_REAL8*)srcs[2]->storage(), times[2],
                   srcshp[0], srcshp[1], srcshp[2],
                   (CCTK_REAL8*)storage(), time,
                   dstshp[0], dstshp[1], dstshp[2],
                   srcbbox, dstbbox, regbbox);
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
      case PROLONG_TVD: // PROLONG_TVD
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
