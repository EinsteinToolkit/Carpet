// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/reduce.cc,v 1.36 2004/03/31 16:38:21 schnetter Exp $

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex>
#include <limits>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "dist.hh"
#include "vect.hh"

#include "carpet.hh"

#include "reduce.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/reduce.cc,v 1.36 2004/03/31 16:38:21 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetReduce_reduce_cc);
}



namespace CarpetReduce {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Helper functions and types
  
  // The minimum of two values
  template<typename T> inline T
  mymin (const T x, const T y)
  {
    return min(x, y);
  }
  
  // The maximum of two values
  template<typename T> inline T
  mymax (const T x, const T y)
  {
    return max(x, y);
  }
  
  // Square root
  template<typename T> inline T
  mysqrt (const T x)
  {
    return sqrt(x);
  }
  
  
  
  // Properties of numeric types
  template<typename T>
  struct my_numeric_limits {
    
    // The smallest possible value
    static T min ()
    {
      return mymin (numeric_limits<T>::min(), -numeric_limits<T>::max());
    }
    
    // The largest possible value
    static T max ()
    {
      return mymax (numeric_limits<T>::max(), -numeric_limits<T>::min());
    }
    
  };
  
  
  
  // Provide for each Cactus type a "good" type, translating between
  // CCTK_COMPLEX* and complex<CCTK_REAL*>
  template<typename T>
  struct typeconv {
    typedef T goodtype;
    typedef T badtype;
  };
  
  
  
  // Overload the above helper functions and types for complex values
  
#ifdef CCTK_REAL4
  
  template<> inline complex<CCTK_REAL4>
  mymin (const complex<CCTK_REAL4> x, const complex<CCTK_REAL4> y)
  {
    return complex<CCTK_REAL4> (mymin(x.real(), y.real()),
                                mymin(x.imag(), y.imag()));
  }
  
  template<> inline complex<CCTK_REAL4>
  mymax (const complex<CCTK_REAL4> x, const complex<CCTK_REAL4> y)
  {
    return complex<CCTK_REAL4> (mymax(x.real(), y.real()),
                                mymax(x.imag(), y.imag()));
  }
  
  template<>
  struct my_numeric_limits<complex<CCTK_REAL4> > {
    static complex<CCTK_REAL4> min ()
    {
      return complex<CCTK_REAL4> (my_numeric_limits<CCTK_REAL4>::min(),
                                  my_numeric_limits<CCTK_REAL4>::min());
    }
    static complex<CCTK_REAL4> max ()
    {
      return complex<CCTK_REAL4> (my_numeric_limits<CCTK_REAL4>::max(),
                                  my_numeric_limits<CCTK_REAL4>::max());
    }
  };

  template<>
  struct typeconv<CCTK_COMPLEX8> {
    typedef complex<CCTK_REAL4> goodtype;
    typedef CCTK_COMPLEX8 badtype;
  };
  
#endif
  
#ifdef CCTK_REAL8
  
  template<> inline complex<CCTK_REAL8>
  mymin (const complex<CCTK_REAL8> x, const complex<CCTK_REAL8> y)
  {
    return complex<CCTK_REAL8> (mymin(x.real(), y.real()),
                                mymin(x.imag(), y.imag()));
  }
  
  template<> inline complex<CCTK_REAL8>
  mymax (const complex<CCTK_REAL8> x, const complex<CCTK_REAL8> y)
  {
    return complex<CCTK_REAL8> (mymax(x.real(), y.real()),
                                mymax(x.imag(), y.imag()));
  }
  
  template<>
  struct my_numeric_limits<complex<CCTK_REAL8> > {
    static complex<CCTK_REAL8> min ()
    {
      return complex<CCTK_REAL8> (my_numeric_limits<CCTK_REAL8>::min(),
                                  my_numeric_limits<CCTK_REAL8>::min());
    }
    static complex<CCTK_REAL8> max ()
    {
      return complex<CCTK_REAL8> (my_numeric_limits<CCTK_REAL8>::max(),
                                  my_numeric_limits<CCTK_REAL8>::max());
    }
  };

  template<>
  struct typeconv<CCTK_COMPLEX16> {
    typedef complex<CCTK_REAL8> goodtype;
    typedef CCTK_COMPLEX16 badtype;
  };
  
#endif
  
#ifdef CCTK_REAL16
  
  template<> inline complex<CCTK_REAL16>
  mymin (const complex<CCTK_REAL16> x, const complex<CCTK_REAL16> y)
  {
    return complex<CCTK_REAL16> (mymin(x.real(), y.real()),
                                 mymin(x.imag(), y.imag()));
  }
  
  template<> inline complex<CCTK_REAL16>
  mymax (const complex<CCTK_REAL16> x, const complex<CCTK_REAL16> y)
  {
    return complex<CCTK_REAL16> (mymax(x.real(), y.real()),
                                 mymax(x.imag(), y.imag()));
  }
  
  template<>
  struct my_numeric_limits<complex<CCTK_REAL16> > {
    static complex<CCTK_REAL16> min ()
    {
      return complex<CCTK_REAL16> (my_numeric_limits<CCTK_REAL16>::min(),
                                   my_numeric_limits<CCTK_REAL16>::min());
    }
    static complex<CCTK_REAL16> max ()
    {
      return complex<CCTK_REAL16> (my_numeric_limits<CCTK_REAL16>::max(),
                                   my_numeric_limits<CCTK_REAL16>::max());
    }
  };

  template<>
  struct typeconv<CCTK_COMPLEX32> {
    typedef complex<CCTK_REAL16> goodtype;
    typedef CCTK_COMPLEX32 badtype;
  };
  
#endif
  
  
  
  // Provide a square root function for integer values
  
#ifdef CCTK_INT1
  template<> inline CCTK_INT1
  mysqrt (const CCTK_INT1 x)
  {
    return sqrt((CCTK_REAL)x);
  }
#endif
  
#ifdef CCTK_INT2
  template<> inline CCTK_INT2
  mysqrt (const CCTK_INT2 x)
  {
    return sqrt((CCTK_REAL)x);
  }
#endif
  
#ifdef CCTK_INT4
  template<> inline CCTK_INT4
  mysqrt (const CCTK_INT4 x)
  {
    return sqrt((CCTK_REAL)x);
  }
#endif
  
#ifdef CCTK_INT8
  template<> inline CCTK_INT8
  mysqrt (const CCTK_INT8 x)
  {
    return sqrt((CCTK_REAL)x);
  }
#endif
  
  
  
  // Poor man's RTTI
  enum ared { do_count, do_origin, do_minimum, do_maximum, do_product, do_sum,
              do_sum_abs, do_sum_squared, do_average, do_norm1, do_norm2,
              do_norm_inf };
  
  
  
  struct reduction {
    virtual ared thered () const = 0;
    virtual bool uses_cnt () const = 0;
    virtual MPI_Op mpi_op () const = 0;
  };
  
  
  
  // count: count the number of grid points
  struct count : reduction {
    count () { }
    ared thered () const { return do_count; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += T(1); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct origin : reduction {
    origin () { }
    ared thered () const { return do_origin; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { assert(0); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct minimum : reduction {
    minimum () { }
    ared thered () const { return do_minimum; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = my_numeric_limits<T>::max(); }
      static inline void reduce (T& accum, const T& val) { accum = mymin(accum, val); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_MIN; }
  };
  
  struct maximum : reduction {
    maximum () { }
    ared thered () const { return do_maximum; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = my_numeric_limits<T>::min(); }
      static inline void reduce (T& accum, const T& val) { accum = mymax(accum, val); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_MAX; }
  };
  
  struct product : reduction {
    product () { }
    ared thered () const { return do_product; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(1); }
      static inline void reduce (T& accum, const T& val) { accum *= val; }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_PROD; }
  };
  
  struct sum : reduction {
    sum () { }
    ared thered () const { return do_sum; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += val; }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct sum_abs : reduction {
    sum_abs () { }
    ared thered () const { return do_sum_abs; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += abs(val); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct sum_squared : reduction {
    sum_squared () { }
    ared thered () const { return do_sum_squared; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += val*val; }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct average : reduction {
    average () { }
    ared thered () const { return do_average; }
    bool uses_cnt () const { return true; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += val; }
      static inline void finalise (T& accum, const T& cnt) { accum /= cnt; }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct norm1 : reduction {
    norm1 () { }
    ared thered () const { return do_norm1; }
    bool uses_cnt () const { return true; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += abs(val); }
      static inline void finalise (T& accum, const T& cnt) { accum /= cnt; }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct norm2 : reduction {
    norm2 () { }
    ared thered () const { return do_norm2; }
    bool uses_cnt () const { return true; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum += abs(val)*abs(val); }
      static inline void finalise (T& accum, const T& cnt) { accum = mysqrt(accum / cnt); }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct norm_inf : reduction {
    norm_inf () { }
    ared thered () const { return do_norm_inf; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = T(0); }
      static inline void reduce (T& accum, const T& val) { accum = mymax(accum, (T)abs(val)); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  
  
  template<class T,class OP>
  void initialise (void* const outval, void* const cnt)
  {
    OP::initialise (*(T*)outval);
    *(T*)cnt = T(0);
  }
  
  template<class T,class OP>
  void reduce (const int* const lsh, const int* const bbox,
               const int* const nghostzones,
	       const void* const inarray, void* const outval, void* const cnt)
  {
    const T *myinarray = (const T*)inarray;
    T myoutval = *(T*)outval;
    T mycnt = *(T*)cnt;
    vect<int,dim> imin, imax;
    for (int d=0; d<dim; ++d) {
      imin[d] = bbox[2*d  ] ? 0      : nghostzones[d];
      imax[d] = bbox[2*d+1] ? lsh[d] : lsh[d] - nghostzones[d];
    }
    assert (dim==3);
    for (int k=imin[2]; k<imax[2]; ++k) {
      for (int j=imin[1]; j<imax[1]; ++j) {
        for (int i=imin[0]; i<imax[0]; ++i) {
	  const int index = i + lsh[0] * (j + lsh[1] * k);
	  OP::reduce (myoutval, myinarray[index]);
	  mycnt += 1;
	}
      }
    }
    *(T*)outval = myoutval;
    *(T*)cnt    = mycnt;
  }
  
  template<class T,class OP>
  void finalise (void* const outval, const void* const cnt)
  {
    OP::finalise (*(T*)outval, *(const T*)cnt);
  }
  
  
  
  void Initialise (const cGH* const cgh, const int proc,
		   const int num_outvals,
		   void* const myoutvals, const int outtype,
		   void* const mycounts,
		   const reduction* const red)
  {
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (num_outvals>=0);
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    assert (myoutvals);
    assert (mycounts);
    
    assert (red);
    
    for (int n=0; n<num_outvals; ++n) {
      
      switch (outtype) {
#define INITIALISE(OP,S)                                                \
      case do_##OP: {                                                   \
        typedef typeconv<S>::goodtype T;                                \
	initialise<T,OP::op<T> > (&((char*)myoutvals)[vartypesize*n],   \
				  &((char*)mycounts )[vartypesize*n]);  \
	break;                                                          \
      }
#define TYPECASE(N,T)				\
      case N: {					\
	switch (red->thered()) {		\
	  INITIALISE(count,T);			\
	  INITIALISE(origin,T);			\
	  INITIALISE(minimum,T);		\
	  INITIALISE(maximum,T);		\
	  INITIALISE(product,T);		\
	  INITIALISE(sum,T);			\
	  INITIALISE(sum_abs,T);		\
	  INITIALISE(sum_squared,T);		\
	  INITIALISE(average,T);		\
	  INITIALISE(norm1,T);			\
	  INITIALISE(norm2,T);			\
	  INITIALISE(norm_inf,T);		\
	default:				\
	  assert (0);				\
	}					\
	break;					\
      }
#include "Carpet/Carpet/src/typecase"
#undef TYPECASE
#undef INITIALISE
      default:
	assert (0);
      }
      
    } // for n
  }
  
  
  
  void Copy (const cGH* const cgh, const int proc,
             const int lsize,
             const int num_inarrays,
             const void* const* const inarrays, const int intype,
             const int num_outvals,
             void* const myoutvals, const int outtype,
             void* const mycounts)
  {
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (lsize >= 0);
    assert (num_outvals>=0);
    
    assert (num_inarrays>=0);
    assert (num_inarrays * lsize == num_outvals);
    assert (inarrays);
    for (int n=0; n<num_inarrays; ++n) {
      assert (inarrays[n]);
    }
    
    assert (myoutvals);
    assert (mycounts);
    
    assert (outtype == intype);
    
    for (int m=0; m<num_inarrays; ++m) {
      for (int n=0; n<lsize; ++n) {
        
        switch (outtype) {
#define COPY(S)                                                         \
        {                                                               \
          typedef typeconv<S>::goodtype T;                              \
          ((T*)myoutvals)[n+lsize*m] = ((const T*)inarrays[m])[n];      \
          ((T*)mycounts )[n+lsize*m] = T(1);                            \
        }
#define TYPECASE(N,T)                           \
        case N: {                               \
          COPY(T);                              \
          break;                                \
        }
#include "Carpet/Carpet/src/typecase"
#undef TYPECASE
#undef COPY
        default:
          assert (0);
        }
        
      } // for
    } // for
  }
  
  
  
  void Reduce (const cGH* const cgh, const int proc,
	       const int* const mylsh, const int* const mybbox,
               const int* const mynghostzones,
	       const int num_inarrays,
	       const void* const* const inarrays, const int intype,
	       const int num_outvals,
	       void* const myoutvals, const int outtype,
	       void* const mycounts,
	       const reduction* const red)
  {
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (num_outvals>=0);
    
    assert (num_inarrays>=0);
    assert (num_inarrays == num_outvals);
    assert (inarrays);
    for (int n=0; n<num_inarrays; ++n) {
      assert (inarrays[n]);
    }
    
    for (int d=0; d<dim; ++d) {
      assert (mylsh[d]>=0);
      assert (mynghostzones[d]>=0 && 2*mynghostzones[d]<=mylsh[d]);
    }
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    assert (myoutvals);
    assert (mycounts);
    
    assert (outtype == intype);
    
    for (int n=0; n<num_outvals; ++n) {
      
      switch (outtype) {
#define REDUCE(OP,S)                                                    \
      case do_##OP: {                                                   \
        typedef typeconv<S>::goodtype T;                                \
        reduce<T,OP::op<T> > (mylsh, mybbox, mynghostzones, inarrays[n], \
                              &((char*)myoutvals)[vartypesize*n],       \
                              &((char*)mycounts )[vartypesize*n]);      \
        break;                                                          \
      }
#define TYPECASE(N,T)				\
      case N: {					\
	switch (red->thered()) {		\
	  REDUCE(count,T);			\
	  REDUCE(origin,T);			\
	  REDUCE(minimum,T);			\
	  REDUCE(maximum,T);			\
	  REDUCE(product,T);			\
	  REDUCE(sum,T);			\
	  REDUCE(sum_abs,T);			\
	  REDUCE(sum_squared,T);		\
	  REDUCE(average,T);			\
	  REDUCE(norm1,T);			\
	  REDUCE(norm2,T);			\
	  REDUCE(norm_inf,T);			\
	default:				\
	  assert (0);				\
	}					\
	break;					\
      }
#include "Carpet/Carpet/src/typecase"
#undef TYPECASE
#undef REDUCE
      default:
	assert (0);
      }
      
    } // for n
  }
  
  
  
  void Finalise (const cGH* const cgh, const int proc,
		 const int num_outvals,
		 void* const outvals, const int outtype,
		 const void* const myoutvals,
		 const void* const mycounts,
		 const reduction* const red)
  {
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (num_outvals>=0);
    assert (outvals || (proc!=-1 && proc!=CCTK_MyProc(cgh)));
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    assert (myoutvals);
    assert (mycounts);
    
    vector<char> counts;
    if (proc==-1 || proc==CCTK_MyProc(cgh)) {
      counts.resize(vartypesize * num_outvals);
    }
    
    if (proc == -1) {
      MPI_Allreduce ((void*)myoutvals, outvals, num_outvals,
		     CarpetMPIDatatype(outtype), red->mpi_op(),
		     CarpetMPIComm());
      if (red->uses_cnt()) {
	MPI_Allreduce ((void*)mycounts, &counts[0], num_outvals,
		       CarpetMPIDatatype(outtype), MPI_SUM,
		       CarpetMPIComm());
      }
    } else {
      MPI_Reduce ((void*)myoutvals, outvals, num_outvals,
		  CarpetMPIDatatype(outtype), red->mpi_op(),
		  proc, CarpetMPIComm());
      if (red->uses_cnt()) {
	MPI_Reduce ((void*)mycounts, &counts[0], num_outvals,
		    CarpetMPIDatatype(outtype), MPI_SUM,
		    proc, CarpetMPIComm());
      }
    }
    
    if (proc==-1 || proc==CCTK_MyProc(cgh)) {
      
      for (int n=0; n<num_outvals; ++n) {
	
	assert (outvals);
	assert ((int)counts.size() == vartypesize * num_outvals);
	
	switch (outtype) {
#define FINALISE(OP,S)                                                  \
	case do_##OP: {                                                 \
          typedef typeconv<S>::goodtype T;                              \
	  finalise<T,OP::op<T> > (&((char*)outvals)[vartypesize*n],     \
                                  &        counts  [vartypesize*n]);    \
	  break;                                                        \
        }
#define TYPECASE(N,T)				\
	case N: {				\
	  switch (red->thered()) {		\
	    FINALISE(count,T);			\
	    FINALISE(origin,T);			\
	    FINALISE(minimum,T);		\
	    FINALISE(maximum,T);		\
	    FINALISE(product,T);		\
	    FINALISE(sum,T);			\
	    FINALISE(sum_abs,T);		\
	    FINALISE(sum_squared,T);		\
	    FINALISE(average,T);		\
	    FINALISE(norm1,T);			\
	    FINALISE(norm2,T);			\
	    FINALISE(norm_inf,T);		\
	  default:				\
	    assert (0);				\
	  }					\
	  break;				\
	}
#include "Carpet/Carpet/src/typecase"
#undef TYPECASE
#undef FINALISE
	default:
	  assert (0);
	}
	
      } // for n
      
    } // if
  }
  
  
  
  int ReduceArrays (const cGH* const cgh, const int proc,
		    const int num_dims, const int* const dims,
		    const int num_inarrays,
		    const void* const* const inarrays, const int intype,
		    const int num_outvals,
		    void* const outvals, const int outtype,
		    const reduction* const red)
  {
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (num_outvals>=0);
    assert (outvals || (proc!=-1 && proc!=CCTK_MyProc(cgh)));
    
    assert (num_inarrays>=0);
    assert (inarrays);
    for (int n=0; n<num_inarrays; ++n) {
      assert (inarrays[n]);
    }
    
    assert (num_dims>=0 && num_dims<=dim);
    for (int d=0; d<num_dims; ++d) {
      assert (dims[d]>=0);
    }
    
    int lsize = 1;
    for (int d=0; d<num_dims; ++d) {
      lsize *= dims[d];
    }
    
    const bool do_local_reduction = num_outvals == 1;
    
    if (! do_local_reduction) {
      assert (num_outvals == lsize);
    }
    
    vect<int,dim> mylsh, mynghostzones;
    vect<vect<int,2>,dim> mybbox;
    for (int d=0; d<num_dims; ++d) {
      mylsh[d] = dims[d];
      mybbox[d][0] = 0;
      mybbox[d][1] = 0;
      mynghostzones[d] = 0;
    }
    for (int d=num_dims; d<dim; ++d) {
      mylsh[d] = 1;
      mybbox[d][0] = 0;
      mybbox[d][1] = 0;
      mynghostzones[d] = 0;
    }
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    vector<char> myoutvals (vartypesize * num_inarrays * num_outvals);
    vector<char> mycounts  (vartypesize * num_inarrays * num_outvals);
    
    Initialise (cgh, proc, num_inarrays * num_outvals, &myoutvals[0], outtype,
                &mycounts[0], red);
    if (do_local_reduction) {
      Reduce   (cgh, proc, &mylsh[0], &mybbox[0][0], &mynghostzones[0],
                num_inarrays, inarrays, intype,
                num_inarrays * num_outvals, &myoutvals[0], outtype,
                &mycounts[0], red);
    } else {
      Copy     (cgh, proc, lsize, num_inarrays, inarrays, intype,
                num_inarrays * num_outvals, &myoutvals[0], outtype,
                &mycounts[0]);
    }
    Finalise   (cgh, proc, num_inarrays * num_outvals, outvals, outtype,
                &myoutvals[0], &mycounts[0], red);
    
    return 0;
  }
  
  
  
  int ReduceGVs (const cGH* const cgh, const int proc,
		 const int num_outvals, const int outtype, void* const outvals,
		 const int num_invars, const int* const invars,
		 const reduction* const red)
  {
    int ierr;
    
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (num_outvals>=0);
    assert (num_outvals==1);
    assert (outvals || (proc!=-1 && proc!=CCTK_MyProc(cgh)));
    
    assert (num_invars>=0);
    assert (invars);
    for (int n=0; n<num_invars; ++n) {
      assert (invars[n]>=0 && invars[n]<CCTK_NumVars());
    }
    
    if (num_invars==0) return 0;
    
    assert (num_invars>0);
    const int vi = invars[0];
    assert (vi>=0 && vi<CCTK_NumVars());
    
    const int grpdim = CCTK_GroupDimFromVarI(vi);
    assert (grpdim>=0 && grpdim<=dim);
    for (int n=0; n<num_invars; ++n) {
      assert (CCTK_GroupDimFromVarI(invars[n]) == grpdim);
    }
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    
    
    bool const reduce_arrays = CCTK_GroupTypeFromVarI(vi) != CCTK_GF;
    
    for (int n=0; n<num_invars; ++n) {
      if ((CCTK_GroupTypeFromVarI(invars[n]) != CCTK_GF) != reduce_arrays) {
        CCTK_WARN (0, "Cannot (yet) reduce grid functions and grid arrays/scalars at the same time");
      }
    }
    
    // global mode
    if (! reduce_arrays && is_meta_mode()) {
      CCTK_WARN (0, "Grid function reductions are not possible in meta mode");
    }
    if (! reduce_arrays && is_global_mode()) {
      CCTK_WARN (0, "Grid function reduction operators in global mode are not yet implemented");
    }
    
    
    
    vector<char> myoutvals (vartypesize * num_invars * num_outvals);
    vector<char> mycounts  (vartypesize * num_invars * num_outvals);
    
    Initialise (cgh, proc, num_invars * num_outvals, &myoutvals[0], outtype,
                &mycounts[0], red);
    
    // Multiple maps are not supported
    // (because we don't know how to select a map)
    assert (maps == 1);
    const int m = 0;
    
    int const minrl = reduce_arrays ? 0 : reflevel==-1 ? 0                      : reflevel;
    int const maxrl = reduce_arrays ? 1 : reflevel==-1 ? vhh.at(m)->reflevels() : reflevel+1;
    
    BEGIN_GLOBAL_MODE(cgh) {
      for (int rl=minrl; rl<maxrl; ++rl) {
        enter_level_mode (const_cast<cGH*>(cgh), rl);
        BEGIN_MAP_LOOP(cgh, reduce_arrays ? CCTK_ARRAY : CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cgh, reduce_arrays ? CCTK_ARRAY : CCTK_GF) {
              
            assert (grpdim<=dim);
            int lsh[dim], bbox[2*dim], nghostzones[dim];
            ierr = CCTK_GrouplshVI(cgh, grpdim, lsh, vi);
            assert (!ierr);
            ierr = CCTK_GroupbboxVI(cgh, 2*grpdim, bbox, vi);
            assert (!ierr);
            ierr = CCTK_GroupnghostzonesVI(cgh, grpdim, nghostzones, vi);
            assert (!ierr);
            for (int d=0; d<grpdim; ++d) {
              assert (lsh[d]>=0);
              assert (nghostzones[d]>=0 && 2*nghostzones[d]<=lsh[d]);
            }
              
            vector<const void*> inarrays (num_invars);
            for (int n=0; n<num_invars; ++n) {
              inarrays[n] = CCTK_VarDataPtrI(cgh, 0, invars[n]);
              assert (inarrays[n]);
            }
              
            const int intype = CCTK_VarTypeI(vi);
            for (int n=0; n<num_invars; ++n) {
              assert (CCTK_VarTypeI(invars[n]) == intype);
            }
              
            vect<int,dim> mylsh, mynghostzones;
            vect<vect<int,2>,dim> mybbox;
            for (int d=0; d<grpdim; ++d) {
              mylsh[d] = lsh[d];
              mybbox[d][0] = bbox[2*d  ];
              mybbox[d][1] = bbox[2*d+1];
              mynghostzones[d] = nghostzones[d];
            }
            for (int d=grpdim; d<dim; ++d) {
              mylsh[d] = 1;
              mybbox[d][0] = 0;
              mybbox[d][1] = 0;
              mynghostzones[d] = 0;
            }
              
            Reduce (cgh, proc, &mylsh[0], &mybbox[0][0], &mynghostzones[0],
                    num_invars, &inarrays[0], intype,
                    num_invars * num_outvals, &myoutvals[0], outtype,
                    &mycounts[0], red);
              
          } END_LOCAL_COMPONENT_LOOP;
        } END_MAP_LOOP;
        leave_level_mode (const_cast<cGH*>(cgh));
      } // for rl
    } END_GLOBAL_MODE;
    
    Finalise (cgh, proc, num_invars * num_outvals, outvals, outtype,
              &myoutvals[0], &mycounts[0], red);
    
    return 0;
  }
  
  
  
#define REDUCTION(OP)                                                     \
  int OP##_arrays (const cGH * const cgh, const int proc,                 \
		   const int num_dims, const int * const dims,            \
		   const int num_inarrays,                                \
		   const void * const * const inarrays, const int intype, \
		   const int num_outvals,                                 \
		   void * const outvals, const int outtype)               \
  {                                                                       \
    const OP red;                                                         \
    return ReduceArrays                                                   \
      (cgh, proc, num_dims, dims,                                         \
       num_inarrays, inarrays, intype, num_outvals, outvals, outtype,     \
       &red);                                                             \
  }                                                                       \
                                                                          \
  int OP##_GVs (const cGH * const cgh, const int proc,                    \
	        const int num_outvals,                                    \
	        const int outtype, void * const outvals,                  \
	        const int num_invars, const int * const invars)           \
  {                                                                       \
    const OP red;                                                         \
    return ReduceGVs (cgh, proc,                                          \
		      num_outvals, outtype, outvals, num_invars, invars,  \
		      &red);                                              \
  }
  
  REDUCTION(count);
  REDUCTION(origin);
  REDUCTION(minimum);
  REDUCTION(maximum);
  REDUCTION(product);
  REDUCTION(sum);
  REDUCTION(sum_abs);
  REDUCTION(sum_squared);
  REDUCTION(average);
  REDUCTION(norm1);
  REDUCTION(norm2);
  REDUCTION(norm_inf);
  
#undef REDUCTION
  
  
  
  void CarpetReduceStartup ()
  {
    CCTK_RegisterReductionOperator (count_GVs,       "count");
    CCTK_RegisterReductionOperator (origin_GVs,      "origin");
    CCTK_RegisterReductionOperator (minimum_GVs,     "minimum");
    CCTK_RegisterReductionOperator (maximum_GVs,     "maximum");
    CCTK_RegisterReductionOperator (product_GVs,     "product");
    CCTK_RegisterReductionOperator (sum_GVs,         "sum");
    CCTK_RegisterReductionOperator (sum_abs_GVs,     "sum_abs");
    CCTK_RegisterReductionOperator (sum_squared_GVs, "sum_squared");
    CCTK_RegisterReductionOperator (average_GVs,     "average");
    CCTK_RegisterReductionOperator (norm1_GVs,       "norm1");
    CCTK_RegisterReductionOperator (norm2_GVs,       "norm2");
    CCTK_RegisterReductionOperator (norm_inf_GVs,    "norm_inf");
    
    CCTK_RegisterReductionArrayOperator (count_arrays,       "count");
    CCTK_RegisterReductionArrayOperator (origin_arrays,      "origin");
    CCTK_RegisterReductionArrayOperator (minimum_arrays,     "minimum");
    CCTK_RegisterReductionArrayOperator (maximum_arrays,     "maximum");
    CCTK_RegisterReductionArrayOperator (product_arrays,     "product");
    CCTK_RegisterReductionArrayOperator (sum_arrays,         "sum");
    CCTK_RegisterReductionArrayOperator (sum_abs_arrays,     "sum_abs");
    CCTK_RegisterReductionArrayOperator (sum_squared_arrays, "sum_squared");
    CCTK_RegisterReductionArrayOperator (average_arrays,     "average");
    CCTK_RegisterReductionArrayOperator (norm1_arrays,       "norm1");
    CCTK_RegisterReductionArrayOperator (norm2_arrays,       "norm2");
    CCTK_RegisterReductionArrayOperator (norm_inf_arrays,    "norm_inf");
  }
  
} // namespace CarpetReduce
