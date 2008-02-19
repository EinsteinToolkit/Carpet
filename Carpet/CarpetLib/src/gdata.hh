#ifndef GDATA_HH
#define GDATA_HH

#include <cassert>
#include <cstdlib>
#include <queue>
#include <iostream>
#include <string>
#include <vector>

#include "cctk.h"

#include "bbox.hh"
#include "commstate.hh"
#include "defs.hh"
#include "dist.hh"
#include "operators.hh"
#include "timestat.hh"
#include "vect.hh"

using namespace std;



// A generic data storage without type information
class gdata {

protected:                      // should be readonly

  // Fields
  void * _storage;              // A copy of the storage pointer
  
  const int varindex;           // Cactus variable index, or -1

  centering cent;
  operator_type transport_operator;
  
  bool _has_storage;		// has storage associated (on some processor)
  int _size;			// size

  int _proc;			// stored on processor
  
  ivect _shape, _stride;      	// shape and index order
  
  ibbox _extent;		// bbox for all data
  
  bool comm_active;             // a communication is going on
  MPI_Request request;          // outstanding MPI request
  
  int tag;                      // MPI tag for this object
  
private:
  // Forbid copying and passing by value
  gdata (gdata const &);
  gdata & operator= (gdata const &);
  
public:

  // Constructors
  gdata (const int varindex,
         const centering cent = error_centered,
         const operator_type transport_operator = op_error,
         const int tag = -1);

  // Destructors
  virtual ~gdata ();

  // Pseudo constructors
  virtual gdata*
  make_typed (const int varindex,
              const centering cent = error_centered,
              const operator_type transport_operator = op_error,
              const int tag = -1) const = 0;
  
  // Storage management
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const memptr = NULL, size_t const memsize = 0) = 0;
  virtual void free () = 0;
  virtual size_t allocsize (const ibbox& extent, const int proc) const = 0;
  
  // Accessors
  bool has_storage () const {
    return _has_storage;
  }
  
  void const *
  storage ()
    const
  {
    assert (_has_storage);
    return _storage;
  }
  
  void *
  storage ()
  {
    assert (_has_storage);
    return _storage;
  }
  
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

  int elementsize () const {
    return c_datatype_size();
  }

  // Data accessors
  int offset (const ivect& index) const {
    assert (_has_storage);
    assert (all((index-extent().lower()) % extent().stride() == 0));
    ivect ind = (index-extent().lower()) / extent().stride();
    assert (all(ind>=0 && ind<=shape()));
    return dot(ind, stride());
  }
  
private:
  // Datatype accessors
  // maps the C datatype of a data class object to a 0-based index
  virtual unsigned int c_datatype () const = 0;
  virtual size_t c_datatype_size () const = 0;
  
  // Data manipulators
  
public:
  void
  copy_from (comm_state & state,
             gdata const * src,
             ibbox const & box);
  
  void
  transfer_from (comm_state & state,
                 vector<gdata const *> const & srcs,
                 vector<CCTK_REAL>     const & times,
                 ibbox const & dstbox,
                 ibbox const & srcbox,
                 CCTK_REAL time,
                 int order_space,
                 int order_time);
  
protected:
  void
  find_source_timelevel (vector <CCTK_REAL> const & times,
                         CCTK_REAL time,
                         int order_time,
                         int & timelevel0,
                         int & ntimelevels)
    const;
  
private:
  virtual
  void
  copy_from_innerloop (gdata const * gsrc,
                       ibbox const & box)
    = 0;
  
  virtual
  void
  transfer_from_innerloop (vector <gdata const *> const & gsrcs,
                           vector <CCTK_REAL> const & times,
                           ibbox const & box,
                           CCTK_REAL time,
                           int order_space,
                           int order_time)
    = 0;
  
};



#endif // GDATA_HH
