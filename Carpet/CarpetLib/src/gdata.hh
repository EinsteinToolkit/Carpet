#ifndef GDATA_HH
#define GDATA_HH

#include <cassert>
#include <cstdlib>
#include <queue>
#include <iostream>
#include <string>
#include <vector>

#include "cctk.h"

#include "commstate.hh"
#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "operators.hh"
#include "timestat.hh"
#include "vect.hh"

using namespace std;



// A generic data storage without type information
class gdata {

protected:                      // should be readonly

  // Fields
  const int varindex;           // Cactus variable index, or -1
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
         const operator_type transport_operator = op_error,
         const int tag = -1);

  // Destructors
  virtual ~gdata ();

  // Pseudo constructors
  virtual gdata*
  make_typed (const int varindex,
              const operator_type transport_operator = op_error,
              const int tag = -1) const = 0;
  
  // Processor management
  void change_processor (comm_state& state,
                         const int newproc,
                         void* const mem=0);
 protected:
  virtual void change_processor_recv (comm_state& state,
                                      const int newproc,
                                      void* const mem=0)
    = 0;
  virtual void change_processor_send (comm_state& state,
                                      const int newproc,
                                      void* const mem=0)
    = 0;
  virtual void change_processor_wait (comm_state& state,
                                      const int newproc,
                                      void* const mem=0)
    = 0;
 public:
  
  // Storage management
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0) = 0;
  virtual void free () = 0;
  
  // Accessors

  bool has_storage () const {
    return _has_storage;
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
private:
  virtual comm_state::gcommbuf *
  make_typed_commbuf (const ibbox & box)
    const = 0;
  
 public:
    void copy_from (comm_state& state,
                    const gdata* src, const ibbox& box);
 private:
  void copy_from_nocomm (const gdata* src, const ibbox& box);
  void copy_from_recv (comm_state& state,
                       const gdata* src, const ibbox& box);
  void copy_from_send (comm_state& state,
                       const gdata* src, const ibbox& box);
  void copy_from_wait (comm_state& state,
                       const gdata* src, const ibbox& box);

  // Copy processor-local source data into communication send buffer 
  // of the corresponding destination processor
  // The case when both source and destination are local is also handled here.
  void copy_into_sendbuffer (comm_state& state,
                             const gdata* src, const ibbox& box);
  // Copy processor-local destination data from communication recv buffer
  // of the corresponding source processor
  void copy_from_recvbuffer (comm_state& state,
                             const gdata* src, const ibbox& box);

#if 0
 protected:
  virtual void
  copy_from_recv_inner (comm_state& state,
                        const gdata* src,
                        const ibbox& box)
    = 0;
  virtual void
  copy_from_send_inner (comm_state& state,
                        const gdata* src,
                        const ibbox& box)
    = 0;
  virtual void
  copy_from_recv_wait_inner (comm_state& state,
                             const gdata* src,
                             const ibbox& box)
    = 0;
  virtual void
  copy_from_send_wait_inner (comm_state& state,
                             const gdata* src,
                             const ibbox& box)
    = 0;
#endif
  
 public:
  void interpolate_from (comm_state& state,
                         const vector<const gdata*> srcs,
                         const vector<CCTK_REAL> times,
                         const ibbox& box, const CCTK_REAL time,
                         const int order_space,
                         const int order_time);
 private:
  void interpolate_from_nocomm (const vector<const gdata*> srcs,
                                const vector<CCTK_REAL> times,
                                const ibbox& box, const CCTK_REAL time,
                                const int order_space,
                                const int order_time);
  void interpolate_from_recv (comm_state& state,
                              const vector<const gdata*> srcs,
                              const vector<CCTK_REAL> times,
                              const ibbox& box, const CCTK_REAL time,
                              const int order_space,
                              const int order_time);
  void interpolate_from_send (comm_state& state,
                              const vector<const gdata*> srcs,
                              const vector<CCTK_REAL> times,
                              const ibbox& box, const CCTK_REAL time,
                              const int order_space,
                              const int order_time);
  void interpolate_from_wait (comm_state& state,
                              const vector<const gdata*> srcs,
                              const vector<CCTK_REAL> times,
                              const ibbox& box, const CCTK_REAL time,
                              const int order_space,
                              const int order_time);

  // Interpolate processor-local source data into communication send buffer
  // of the corresponding destination processor
  // The case when both source and destination are local is also handled here.
  void interpolate_into_sendbuffer (comm_state& state,
                                    const vector<const gdata*> srcs,
                                    const vector<CCTK_REAL> times,
                                    const ibbox& box,
                                    const CCTK_REAL time,
                                    const int order_space,
                                    const int order_time);
 public:
  
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
