#ifndef DATA_HH
#define DATA_HH

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "cctk.h"

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "gdata.hh"
#include "mem.hh"
#include "vect.hh"

using namespace std;

template<typename T>
class data;

template<typename T>
ostream & operator << ( ostream & os, const data<T> & d );

// A distributed multi-dimensional array
template<typename T>
class data: public gdata
{
  
  // Fields
  mem<T> * _memory;             // the data (if located on this processor)
  
  // For vector groups with contiguous storage
  int vectorlength;             // number of vector elements
  int vectorindex;              // index of this vector element
  data* vectorleader;           // if index!=0: first vector element
   
private:
  // Forbid copying and passing by value
  data (data const &);
  data & operator= (data const &);
  
public:
  
  // Constructors
  data (const int varindex = -1,
        const operator_type transport_operator = op_error,
        const int vectorlength = 1, const int vectorindex = 0,
        data* const vectorleader = NULL,
        const int tag = -1);
  data (const int varindex, const operator_type transport_operator,
        const int vectorlength, const int vectorindex,
        data* const vectorleader,
        const ibbox& extent, const int proc);

  // Destructors
  virtual ~data ();

  // Pseudo constructors
  virtual data* make_typed (const int varindex,
                            const operator_type transport_operator,
                            const int tag) const;

  // Storage management
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const memptr = NULL);
  virtual void free ();

  // Processor management
private:
  virtual void change_processor_recv (comm_state& state,
                                      const int newproc,
                                      void* const memptr = NULL);
  virtual void change_processor_send (comm_state& state,
                                      const int newproc,
                                      void* const memptr = NULL);
  virtual void change_processor_wait (comm_state& state,
                                      const int newproc,
                                      void* const memptr = NULL);
public:

  // Accessors
  virtual const void* storage () const
  {
    assert (_has_storage);
    if (! _memory) return 0;
    return _memory->storage(vectorindex);
  }

  virtual void* storage () {
    assert (_has_storage);
    if (! _memory) return 0;
    return _memory->storage(vectorindex);
  }
  
  // Data accessors
  const T& operator[] (const ivect& index) const
  {
    assert (_has_storage);
    assert (_memory);
    return _memory->storage(vectorindex)[offset(index)];
  }
  
  T& operator[] (const ivect& index)
  {
    assert (_has_storage);
    assert (_memory);
    return _memory->storage(vectorindex)[offset(index)];
  }

#if 0
protected:
  virtual void
  copy_from_recv_inner (comm_state& state,
                        const gdata* src,
                        const ibbox& box);
  virtual void
  copy_from_send_inner (comm_state& state,
                        const gdata* src,
                        const ibbox& box);
  virtual void
  copy_from_recv_wait_inner (comm_state& state,
                             const gdata* src,
                             const ibbox& box);
  virtual void
  copy_from_send_wait_inner (comm_state& state,
                             const gdata* src,
                             const ibbox& box);
#endif
  
  // Datatype accessors
private:
  // maps the C datatype of a data class object to a 0-based index
  virtual unsigned int c_datatype () const
  {
    T dummy;
    return dist::c_datatype (dummy);
  }

  // Data manipulators
private:
  virtual comm_state::gcommbuf *
  make_typed_commbuf (const ibbox & box)
    const;
  
public:
  void copy_from_innerloop (const gdata* gsrc,
			    const ibbox& box);
  void interpolate_from_innerloop (const vector<const gdata*> gsrcs,
				   const vector<CCTK_REAL> times,
				   const ibbox& box, const CCTK_REAL time,
				   const int order_space,
				   const int order_time);
  
  
public:

  // Output
  ostream& output (ostream& os) const;
private:
  bool try_without_time_interpolation (const vector<const gdata*> & gsrcs,
                                       const vector<CCTK_REAL> & times,
                                       const ibbox& box, const CCTK_REAL time,
                                       const int order_space,
                                       const int order_time);
  void interpolate_restrict (const vector<const data<T>*> & gsrcs,
                              const vector<CCTK_REAL> & times,
                              const ibbox& box);
  void interpolate_prolongate (const vector<const data<T>*> & gsrcs,
                              const vector<CCTK_REAL> & times,
                              const ibbox& box, const CCTK_REAL time,
                              const int order_space,
                              const int order_time);
  void Check_that_the_times_are_consistent ( const vector<CCTK_REAL> & times,
                              const CCTK_REAL time );

  friend ostream & operator << <T> ( ostream & os, const data<T> & d );
};



// Declare a specialisation
template<>
void data<CCTK_REAL8>
::interpolate_from_innerloop (const vector<const gdata*> gsrcs,
                              const vector<CCTK_REAL> times,
                              const ibbox& box, const CCTK_REAL time,
                              const int order_space,
                              const int order_time);


#endif // DATA_HH
