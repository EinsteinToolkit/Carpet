// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/reduce.cc,v 1.30 2003/11/12 17:29:30 schnetter Exp $

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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/reduce.cc,v 1.30 2003/11/12 17:29:30 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetReduce_reduce_cc);
}



namespace CarpetReduce {
  
  using namespace Carpet;
  
  
  
  // Helper functions
  template<class T> inline T mymin (const T x, const T y) { return min(x, y); }
  template<> inline complex<float> mymin (const complex<float> x, const complex<float> y) { return complex<float> (min(x.real(), y.real()), min(x.imag(), y.imag())); }
  template<> inline complex<double> mymin (const complex<double> x, const complex<double> y) { return complex<double> (min(x.real(), y.real()), min(x.imag(), y.imag())); }
#ifdef LDBL_MAX
  template<> inline complex<long double> mymin (const complex<long double> x, const complex<long double> y) { return complex<long double> (min(x.real(), y.real()), min(x.imag(), y.imag())); }
#endif
  
  template<class T> inline T mymax (const T x, const T y) { return max(x, y); }
  template<> inline complex<float> mymax (const complex<float> x, const complex<float> y) { return complex<float> (max(x.real(), y.real()), max(x.imag(), y.imag())); }
  template<> inline complex<double> mymax (const complex<double> x, const complex<double> y) { return complex<double> (max(x.real(), y.real()), max(x.imag(), y.imag())); }
#ifdef LDBL_MAX
  template<> inline complex<long double> mymax (const complex<long double> x, const complex<long double> y) { return complex<long double> (max(x.real(), y.real()), max(x.imag(), y.imag())); }
#endif
  
  template<class T> inline T mymin () { return mymin(numeric_limits<T>::min(),-numeric_limits<T>::max()); }
  template<> inline complex<float> mymin () { return complex<float> (mymin<float>(),mymin<float>()); }
  template<> inline complex<double> mymin () { return complex<double> (mymin<double>(),mymin<double>()); }
#ifdef LDBL_MAX
  template<> inline complex<long double> mymin () { return complex<long double> (mymin<long double>(),mymin<long double>()); }
#endif
  
  template<class T> inline T mymax () { return mymax(numeric_limits<T>::max(),-numeric_limits<T>::min()); }
  template<> inline complex<float> mymax () { return complex<float> (mymax<float>(),mymax<float>()); }
  template<> inline complex<double> mymax () { return complex<double> (mymax<double>(),mymax<double>()); }
#ifdef LDBL_MAX
  template<> inline complex<long double> mymax () { return complex<long double> (mymax<long double>(),mymax<long double>()); }
#endif
  
  
  
  // Poor man's RTTI
  enum ared { do_count, do_minimum, do_maximum, do_product, do_sum, do_sum_abs,
	      do_sum_squared, do_average, do_norm1, do_norm2, do_norm_inf };
  
  
  
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
      static inline void initialise (T& accum) { accum = 0; }
      static inline void reduce (T& accum, const T& val) { accum += 1; }
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
      static inline void initialise (T& accum) { accum = mymax<T>(); }
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
      static inline void initialise (T& accum) { accum = mymin<T>(); }
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
      static inline void initialise (T& accum) { accum = 1; }
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
      static inline void initialise (T& accum) { accum = 0; }
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
      static inline void initialise (T& accum) { accum = 0; }
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
      static inline void initialise (T& accum) { accum = 0; }
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
      static inline void initialise (T& accum) { accum = 0; }
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
      static inline void initialise (T& accum) { accum = 0; }
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
      static inline void initialise (T& accum) { accum = 0; }
      static inline void reduce (T& accum, const T& val) { accum += abs(val)*abs(val); }
      static inline void finalise (T& accum, const T& cnt) { accum = sqrt(accum / cnt); }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  template<> inline void norm2::op<char>         ::finalise (char         & accum, const char         & cnt) { accum = (char         )sqrt((double)accum / cnt); }
  template<> inline void norm2::op<signed char>  ::finalise (signed char  & accum, const signed char  & cnt) { accum = (signed char  )sqrt((double)accum / cnt); }
  template<> inline void norm2::op<unsigned char>::finalise (unsigned char& accum, const unsigned char& cnt) { accum = (unsigned char)sqrt((double)accum / cnt); }
  template<> inline void norm2::op<short>        ::finalise (short        & accum, const short        & cnt) { accum = (short        )sqrt((double)accum / cnt); }
  template<> inline void norm2::op<int>          ::finalise (int          & accum, const int          & cnt) { accum = (int          )sqrt((double)accum / cnt); }
  template<> inline void norm2::op<long>         ::finalise (long         & accum, const long         & cnt) { accum = (long         )sqrt((double)accum / cnt); }
  template<> inline void norm2::op<long long>    ::finalise (long long    & accum, const long long    & cnt) { accum = (long long    )sqrt((double)accum / cnt); }
  
  struct norm_inf : reduction {
    norm_inf () { }
    ared thered () const { return do_norm_inf; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = 0; }
      static inline void reduce (T& accum, const T& val) { accum = mymax(accum, (T)abs(val)); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  
  
  template<class T,class OP>
  void initialise (void* const outval, void* const cnt)
  {
    OP::initialise (*(T*)outval);
    *(T*)cnt = 0;
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
#define INITIALISE(OP,T)                                                \
      case do_##OP:                                                     \
	initialise<T,OP::op<T> > (&((char*)myoutvals)[vartypesize*n],   \
				  &((char*)mycounts )[vartypesize*n]);  \
	break;
#define TYPECASE(N,T)				\
      case N: {					\
	switch (red->thered()) {		\
	  INITIALISE(count,T);			\
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
#define COPY(T)                                                         \
          ((T*)myoutvals)[n+lsize*m] = ((const T*)inarrays[m])[n];      \
          ((T*)mycounts )[n+lsize*m] = (T)1;
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
#define REDUCE(OP,T)                                                    \
	case do_##OP:                                                   \
	  reduce<T,OP::op<T> > (mylsh, mybbox, mynghostzones, inarrays[n], \
				&((char*)myoutvals)[vartypesize*n],     \
				&((char*)mycounts )[vartypesize*n]);    \
	  break;
#define TYPECASE(N,T)				\
      case N: {					\
	switch (red->thered()) {		\
	  REDUCE(count,T);			\
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
#define FINALISE(OP,T)                                                  \
	case do_##OP:                                                   \
	  finalise<T,OP::op<T> > (&((char*)outvals)[vartypesize*n],     \
                                  &        counts  [vartypesize*n]);    \
	  break;
#define TYPECASE(N,T)				\
	case N: {				\
	  switch (red->thered()) {		\
	    FINALISE(count,T);			\
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
    if (! reduce_arrays && reflevel == -1) {
      CCTK_WARN (0, "Grid function reduction operators in global mode are not yet implemented");
    }
    
    
    
    vector<char> myoutvals (vartypesize * num_invars * num_outvals);
    vector<char> mycounts  (vartypesize * num_invars * num_outvals);
    
    Initialise (cgh, proc, num_invars * num_outvals, &myoutvals[0], outtype,
                &mycounts[0], red);
    
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
