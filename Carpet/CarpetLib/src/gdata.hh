// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.hh,v 1.19 2003/10/15 07:14:01 hawke Exp $

#ifndef GDATA_HH
#define GDATA_HH

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>

#include "cctk.h"

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "vect.hh"

using namespace std;



// A generic data storage without type information
template<int D>
class gdata {

  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;

protected:                      // should be readonly

  // Fields
  int varindex;                 // Cactus variable index, or -1
  
  bool _has_storage;		// has storage associated (on some processor)
  bool _owns_storage;		// owns the storage
  // (only valid if there is storage on this processor; it means that
  // the memory is allocated and freed by this class)
  int _size;			// size

  int _proc;			// stored on processor
  
  ivect _shape, _stride;      	// shape and index order
  
  ibbox _extent;		// bbox for all data
  
public:

  // Constructors
  gdata (const int varindex);

  // Destructors
  virtual ~gdata ();

  // Pseudo constructors
  virtual gdata<D>* make_typed (const int varindex) const = 0;
  
  // Processor management
  virtual void change_processor (const int newproc, void* const mem=0) = 0;
  
  // Storage management
  virtual void transfer_from (gdata<D>* src) = 0;
  
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0) = 0;
  virtual void free () = 0;
  
  // Accessors

  int var_index () const {
    return varindex;
  }
  
  bool has_storage () const {
    return _has_storage;
  }
  bool owns_storage () const {
    assert (_has_storage);
    return _owns_storage;
  }
  
  virtual const void* storage () const = 0;
  
  virtual void* storage () = 0;
  
  int size () const {
    assert (_has_storage);
    return _size;
  }

  int proc () const {
    assert (_has_storage);
    return _proc;
  }
  
  const ivect& shape () const {
    assert (_has_storage);
    return _shape;
  }

  const ivect& stride () const {
    assert (_has_storage);
    return _stride;
  }
  
  const ibbox& extent () const {
    assert (_has_storage);
    return _extent;
  }

  // Data accessors
  int offset (const ivect& index) const {
    assert (_has_storage);
    assert (all((index-extent().lower()) % extent().stride() == 0));
    ivect ind = (index-extent().lower()) / extent().stride();
    assert (all(ind>=0 && ind<=shape()));
    return dot(ind, stride());
  }

  // Data manipulators
  void copy_from (const gdata* src, const ibbox& box);
  void interpolate_from (const vector<const gdata*> srcs,
			 const vector<CCTK_REAL> times,
			 const ibbox& box, const CCTK_REAL time,
			 const int order_space,
			 const int order_time);
protected:
  virtual void
  copy_from_innerloop (const gdata* src, const ibbox& box) = 0;
  virtual void
  interpolate_from_innerloop (const vector<const gdata*> srcs,
			      const vector<CCTK_REAL> times,
			      const ibbox& box, const CCTK_REAL time,
			      const int order_space,
			      const int order_time) = 0;
  
};



#endif // GDATA_HH
