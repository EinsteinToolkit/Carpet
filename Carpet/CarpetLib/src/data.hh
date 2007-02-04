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
        const centering cent = error_centered,
        const operator_type transport_operator = op_error,
        const int vectorlength = 1, const int vectorindex = 0,
        data* const vectorleader = NULL,
        const int tag = -1);
  data (const int varindex,
        const centering cent, const operator_type transport_operator,
        const int vectorlength, const int vectorindex,
        data* const vectorleader,
        const ibbox& extent, const int proc);

  // Destructors
  virtual ~data ();

  // Pseudo constructors
  virtual data* make_typed (const int varindex,
                            const centering cent,
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
  void interpolate_from_innerloop (const vector<const gdata*>& gsrcs,
				   const vector<CCTK_REAL>& times,
				   const ibbox& box, const CCTK_REAL time,
				   const int order_space,
				   const int order_time);
  
private:
  void interpolate_time (vector <data const *> const & srcs,
                         vector <CCTK_REAL> const & times,
                         ibbox const & box,
                         CCTK_REAL const time,
                         int const order_space,
                         int const order_time);
  void interpolate_p_r (data const * const src,
                        ibbox const & box,
                        int const order_space);
  void interpolate_p_vc_cc (data const * const src,
                            ibbox const & box,
                            int const order_space);
  void interpolate_prolongate (data const * src,
                               ibbox const & box,
                               int order_space);
  void interpolate_restrict (data const * src,
                             ibbox const & box,
                             int order_space);
  void time_interpolate (vector <data *> const & srcs,
                         ibbox const & box,
                         vector <CCTK_REAL> const & times,
                         CCTK_REAL time,
                         int order_time);
  
public:

  // Output
  ostream& output (ostream& os) const;

  friend ostream & operator << <T> ( ostream & os, const data<T> & d );
};

#endif // DATA_HH
