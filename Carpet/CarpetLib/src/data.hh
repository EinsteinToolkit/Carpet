// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/data.hh,v 1.16 2004/01/25 14:57:29 schnetter Exp $

#ifndef DATA_HH
#define DATA_HH

#include <assert.h>

#include <iostream>
#include <string>

#include "cctk.h"

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "gdata.hh"
#include "vect.hh"

using namespace std;



// A real data storage
template<class T,int D>
class data: public gdata<D> {
  
  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;

  // Fields
  T* _storage;			// the data (if located on this processor)
  
  bool comm_active;
  MPI_Request request;
  
  int tag;                      // MPI tag for this object
  
public:
  
  // Constructors
  data (const int varindex = -1,
        const operator_type transport_operator = op_error);
  data (const int varindex, const operator_type transport_operator,
        const ibbox& extent, const int proc);

  // Destructors
  virtual ~data ();

  // Pseudo constructors
  virtual data* make_typed (const int varindex,
                            const operator_type transport_operator) const;

  // Storage management
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0);
  virtual void free ();
  virtual void transfer_from (gdata<D>* gsrc);

  // Processor management
  virtual void change_processor (comm_state<D>& state,
                                 const int newproc, void* const mem=0);
 private:
  virtual void change_processor_recv (const int newproc, void* const mem=0);
  virtual void change_processor_send (const int newproc, void* const mem=0);
  virtual void change_processor_wait (const int newproc, void* const mem=0);
 public:

  // Accessors
  virtual const void* storage () const {
    assert (this->_has_storage);
    return _storage;
  }

  virtual void* storage () {
    assert (this->_has_storage);
    return _storage;
  }
  
  // Data accessors
  const T& operator[] (const ivect& index) const {
    assert (_storage);
    return _storage[offset(index)];
  }
  
  T& operator[] (const ivect& index) {
    assert (_storage);
    return _storage[offset(index)];
  }
  
  // Data manipulators
  void copy_from_innerloop (const gdata<D>* gsrc,
			    const ibbox& box);
  void interpolate_from_innerloop (const vector<const gdata<D>*> gsrcs,
				   const vector<CCTK_REAL> times,
				   const ibbox& box, const CCTK_REAL time,
				   const int order_space,
				   const int order_time);
  
public:

  // Output
  ostream& output (ostream& os) const;
};



#endif // DATA_HH
