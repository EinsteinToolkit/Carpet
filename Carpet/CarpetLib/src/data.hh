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
#include "vect.hh"

using namespace std;



// A distributed multi-dimensional array
template<typename T>
class data: public gdata
{
  
  // Fields
  T* _storage;			// the data (if located on this processor)
  size_t _allocated_bytes;      // number of allocated bytes
  
  // For vector groups with contiguous storage
  int vectorlength;             // number of vector elements
  int vectorindex;              // index of this vector element
  data* vectorleader;           // if index!=0: first vector element
  vector<bool> vectorclients;   // if index==0: registered elements
  
  void register_client (int index);
  void unregister_client (int index);
  bool has_clients ();
  
public:
  
  // Constructors
  data (const int varindex = -1,
        const operator_type transport_operator = op_error,
        const int vectorlength = 1, const int vectorindex = 0,
        data* const vectorleader = NULL);
  data (const int varindex, const operator_type transport_operator,
        const int vectorlength, const int vectorindex,
        data* const vectorleader,
        const ibbox& extent, const int proc);

  // Destructors
  virtual ~data ();

  // Pseudo constructors
  virtual data* make_typed (const int varindex,
                            const operator_type transport_operator) const;

  // Storage management
private:
  void getmem (const size_t nelems);
  void freemem ();
public:
  virtual void allocate (const ibbox& extent, const int proc,
			 void* const mem=0);
  virtual void free ();
  virtual void transfer_from (gdata* gsrc);

private:
  T* vectordata (const int vectorindex) const;
public:

  // Processor management
private:
  virtual void change_processor_recv (comm_state& state,
                                      const int newproc,
                                      void* const mem=0);
  virtual void change_processor_send (comm_state& state,
                                      const int newproc,
                                      void* const mem=0);
  virtual void change_processor_wait (comm_state& state,
                                      const int newproc,
                                      void* const mem=0);
public:

  // Accessors
  virtual const void* storage () const
  {
    assert (_has_storage);
    return _storage;
  }

  virtual void* storage () {
    assert (_has_storage);
    return _storage;
  }
  
  // Data accessors
  const T& operator[] (const ivect& index) const
  {
    assert (_storage);
    return _storage[offset(index)];
  }
  
  T& operator[] (const ivect& index)
  {
    assert (_storage);
    return _storage[offset(index)];
  }
  
  // Data manipulators
private:
  static void
  fill_bbox_arrays (int srcshp[dim],
                    int dstshp[dim],
                    int srcbbox[dim][dim],
                    int dstbbox[dim][dim],
                    int regbbox[dim][dim],
                    const ibbox & box, const ibbox & sext, const ibbox & dext);
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
  bool interpolate_in_time (const vector<const gdata*> & gsrcs,
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
