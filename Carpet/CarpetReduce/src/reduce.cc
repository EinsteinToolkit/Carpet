// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/reduce.cc,v 1.13 2002/10/24 10:51:13 schnetter Exp $

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex>

#include <mpi.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/dist.hh"

#include "Carpet/Carpet/src/carpet.hh"

#include "reduce.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetReduce/src/reduce.cc,v 1.13 2002/10/24 10:51:13 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetReduce_reduce_cc);
}



#ifndef LLONG_MAX
// #  warning "no long long int support"
#endif
#ifndef LDBL_MAX
// #  warning "no long double support"
#endif



namespace CarpetReduce {
  
  using namespace Carpet;
  
  
  
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
    ared thered () const { return do_count; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = 0; }
      static inline void reduce (T& accum, const T& val) { ++accum; }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_SUM; }
  };
  
  struct minimum : reduction {
    ared thered () const { return do_minimum; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { assert (0); }
      static inline void reduce (T& accum, const T& val) { accum = min(accum, val); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_MIN; }
  };
  template<> inline void minimum::op<char>         ::initialise (char         & accum) { accum = CHAR_MAX; }
  template<> inline void minimum::op<signed char>  ::initialise (signed char  & accum) { accum = SCHAR_MAX; }
  template<> inline void minimum::op<unsigned char>::initialise (unsigned char& accum) { accum = UCHAR_MAX; }
  template<> inline void minimum::op<short>        ::initialise (short        & accum) { accum = SHRT_MAX; }
  template<> inline void minimum::op<int>          ::initialise (int          & accum) { accum = INT_MAX; }
  template<> inline void minimum::op<long>         ::initialise (long         & accum) { accum = LONG_MAX; }
#ifdef LLONG_MAX
  template<> inline void minimum::op<long long>    ::initialise (long long    & accum) { accum = LLONG_MAX; }
#endif
  template<> inline void minimum::op<float>        ::initialise (float        & accum) { accum = FLT_MAX; }
  template<> inline void minimum::op<double>       ::initialise (double       & accum) { accum = DBL_MAX; }
#ifdef LDBL_MAX
  template<> inline void minimum::op<long double>  ::initialise (long double  & accum) { accum = LDBL_MAX; }
#endif
  template<> inline void minimum::op<complex<float> >        ::initialise (complex<float>        & accum) { accum = complex<float>(FLT_MAX,FLT_MAX); }
  template<> inline void minimum::op<complex<double> >       ::initialise (complex<double>       & accum) { accum = complex<double>(DBL_MAX,DBL_MAX); }
#ifdef LDBL_MAX
  template<> inline void minimum::op<complex<long double> >  ::initialise (complex<long double>  & accum) { accum = complex<long double>(LDBL_MAX,LDBL_MAX); }
#endif
  template<> inline void minimum::op<complex<float> >        ::reduce (complex<float>& accum, const complex<float>& val) { accum = complex<float>(min(accum.real(), val.real()), min(accum.imag(), val.imag())); }
  template<> inline void minimum::op<complex<double> >       ::reduce (complex<double>& accum, const complex<double>& val) { accum = complex<double>(min(accum.real(), val.real()), min(accum.imag(), val.imag())); }
#ifdef LDBL_MAX
  template<> inline void minimum::op<complex<long double> >  ::reduce (complex<long double>& accum, const complex<long double>& val) { accum = complex<long double>(min(accum.real(), val.real()), min(accum.imag(), val.imag())); }
#endif
  
  struct maximum : reduction {
    ared thered () const { return do_maximum; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { assert (0); }
      static inline void reduce (T& accum, const T& val) { accum = max(accum, val); }
      static inline void finalise (T& accum, const T& cnt) { }
    };
    MPI_Op mpi_op () const { return MPI_MAX; }
  };
  template<> inline void maximum::op<char>         ::initialise (char         & accum) { accum = CHAR_MIN; }
  template<> inline void maximum::op<signed char>  ::initialise (signed char  & accum) { accum = SCHAR_MIN; }
  template<> inline void maximum::op<unsigned char>::initialise (unsigned char& accum) { accum = 0; }
  template<> inline void maximum::op<short>        ::initialise (short        & accum) { accum = SHRT_MIN; }
  template<> inline void maximum::op<int>          ::initialise (int          & accum) { accum = INT_MIN; }
  template<> inline void maximum::op<long>         ::initialise (long         & accum) { accum = LONG_MIN; }
#ifdef LLONG_MIN
  template<> inline void maximum::op<long long>    ::initialise (long long    & accum) { accum = LLONG_MIN; }
#endif
  template<> inline void maximum::op<float>        ::initialise (float        & accum) { accum = -FLT_MAX; }
  template<> inline void maximum::op<double>       ::initialise (double       & accum) { accum = -DBL_MAX; }
#ifdef LDBL_MAX
  template<> inline void maximum::op<long double>  ::initialise (long double  & accum) { accum = -LDBL_MAX; }
#endif
  template<> inline void maximum::op<complex<float> >        ::initialise (complex<float>        & accum) { accum = complex<float>(-FLT_MAX,-FLT_MAX); }
  template<> inline void maximum::op<complex<double> >       ::initialise (complex<double>       & accum) { accum = complex<double>(-DBL_MAX,-DBL_MAX); }
#ifdef LDBL_MAX
  template<> inline void maximum::op<complex<long double> >  ::initialise (complex<long double>  & accum) { accum = complex<long double>(-LDBL_MAX,-LDBL_MAX); }
#endif
  template<> inline void maximum::op<complex<float> >        ::reduce (complex<float>& accum, const complex<float>& val) { accum = complex<float>(max(accum.real(), val.real()), max(accum.imag(), val.imag())); }
  template<> inline void maximum::op<complex<double> >       ::reduce (complex<double>& accum, const complex<double>& val) { accum = complex<double>(max(accum.real(), val.real()), max(accum.imag(), val.imag())); }
#ifdef LDBL_MAX
  template<> inline void maximum::op<complex<long double> >  ::reduce (complex<long double>& accum, const complex<long double>& val) { accum = complex<long double>(max(accum.real(), val.real()), max(accum.imag(), val.imag())); }
#endif
  
  struct product : reduction {
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
  template<> inline void norm2::op<char>         ::finalise (char         & accum, const char         & cnt) { accum = sqrt((double)accum / cnt); }
  template<> inline void norm2::op<signed char>  ::finalise (signed char  & accum, const signed char  & cnt) { accum = sqrt((double)accum / cnt); }
  template<> inline void norm2::op<unsigned char>::finalise (unsigned char& accum, const unsigned char& cnt) { accum = sqrt((double)accum / cnt); }
  template<> inline void norm2::op<short>        ::finalise (short        & accum, const short        & cnt) { accum = sqrt((double)accum / cnt); }
  template<> inline void norm2::op<int>          ::finalise (int          & accum, const int          & cnt) { accum = sqrt((double)accum / cnt); }
  template<> inline void norm2::op<long>         ::finalise (long         & accum, const long         & cnt) { accum = sqrt((double)accum / cnt); }
  template<> inline void norm2::op<long long>    ::finalise (long long    & accum, const long long    & cnt) { accum = sqrt((double)accum / cnt); }
  
  struct norm_inf : reduction {
    ared thered () const { return do_norm_inf; }
    bool uses_cnt () const { return false; }
    template<class T>
    struct op {
      static inline void initialise (T& accum) { accum = 0; }
      static inline void reduce (T& accum, const T& val) { accum = max(accum, abs(val)); }
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
  void reduce (const int* const dims, const int* const nghostzones,
	       const void* const inarray, void* const outval, void* const cnt)
  {
    const T *myinarray = (const T*)inarray;
    T myoutval = *(T*)outval;
    T mycnt = *(T*)cnt;
    for (int k=nghostzones[2]; k<dims[2]-nghostzones[2]; ++k) {
      for (int j=nghostzones[1]; j<dims[1]-nghostzones[1]; ++j) {
	for (int i=nghostzones[0]; i<dims[0]-nghostzones[0]; ++i) {
	  const int index = i + dims[0] * (j + dims[1] * k);
	  OP::reduce (myoutval, myinarray[index]);
	  ++mycnt;
	}
      }
    }
    *(T*)outval = myoutval;
    *(T*)cnt    = mycnt;
  }
  
  template<class T,class OP>
  void finalise (void* const outval, const void* const cnt)
  {
    OP::finalise (*(T*)outval, *(T*)cnt);
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
#define INITIALISE(OP,T)						\
      case do_##OP:							\
	initialise<T,OP::op<T> > (&((char*)myoutvals)[vartypesize*n],	\
				  &((char*)mycounts)[vartypesize*n]);	\
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
  
  
  
  void Reduce (const cGH* const cgh, const int proc,
	       const int* const mydims, const int* const mynghostzones,
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
      assert (mydims[d]>=0);
      assert (mynghostzones[d]>=0 && 2*mynghostzones[d]<=mydims[d]);
    }
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    assert (myoutvals);
    assert (mycounts);
    
    for (int n=0; n<num_outvals; ++n) {
      
      switch (outtype) {
#define REDUCE(OP,T)							\
	case do_##OP:							\
	  reduce<T,OP::op<T> > (mydims, mynghostzones, inarrays[n],	\
				&((char*)myoutvals)[vartypesize*n],	\
				&((char*)mycounts)[vartypesize*n]);	\
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
    
    char *counts = 0;
    if (proc==-1 || proc==CCTK_MyProc(cgh)) {
      counts = new char [vartypesize * num_outvals];
    }
    
    if (proc == -1) {
      MPI_Allreduce ((void*)myoutvals, outvals, num_outvals,
		     CarpetMPIDatatype(outtype), red->mpi_op(),
		     CarpetMPIComm());
      if (red->uses_cnt()) {
	MPI_Allreduce ((void*)mycounts, counts, num_outvals,
		       CarpetMPIDatatype(outtype), MPI_SUM,
		       CarpetMPIComm());
      }
    } else {
      MPI_Reduce ((void*)myoutvals, outvals, num_outvals,
		  CarpetMPIDatatype(outtype), red->mpi_op(),
		  proc, CarpetMPIComm());
      if (red->uses_cnt()) {
	MPI_Reduce ((void*)mycounts, counts, num_outvals,
		    CarpetMPIDatatype(outtype), MPI_SUM,
		    proc, CarpetMPIComm());
      }
    }
    
    if (proc==-1 || proc==CCTK_MyProc(cgh)) {
      
      for (int n=0; n<num_outvals; ++n) {
	
	assert (outvals);
	assert (counts);
	
	switch (outtype) {
#define FINALISE(OP,T)							\
	case do_##OP:							\
	  finalise<T,OP::op<T> >					\
	    (&((char*)outvals)[vartypesize*n], &counts[vartypesize*n]);	\
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
      
      delete [] counts;
      
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
    assert (num_inarrays == num_outvals);
    assert (inarrays);
    for (int n=0; n<num_inarrays; ++n) {
      assert (inarrays[n]);
    }
    
    assert (num_dims>=0 && num_dims<=dim);
    for (int d=0; d<num_dims; ++d) {
      assert (dims[d]>=0);
    }
    
    int mydims[dim], mynghostzones[dim];
    for (int d=0; d<num_dims; ++d) {
      mydims[d] = dims[d];
      mynghostzones[d] = 0;
    }
    for (int d=num_dims; d<dim; ++d) {
      mydims[d] = 1;
      mynghostzones[d] = 0;
    }
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    char *myoutvals = new char [vartypesize * num_outvals];
    char *mycounts  = new char [vartypesize * num_outvals];
    
    Initialise (cgh, proc, num_outvals, myoutvals, outtype, mycounts, red);
    Reduce     (cgh, proc, mydims, mynghostzones,
		num_inarrays, inarrays, intype,
		num_outvals, myoutvals, outtype, mycounts, red);
    Finalise   (cgh, proc, num_outvals, outvals, outtype, myoutvals, mycounts,
		red);
    
    delete [] myoutvals;
    delete [] mycounts;
    
    return 0;
  }
  
  
  
  int ReduceGVs (const cGH* const cgh, const int proc,
		 const int num_outvals, const int outtype, void* const outvals,
		 const int num_invars, const int* const invars,
		 const reduction* const red)
  {
    int ierr;
    
    // global mode
    if (component != -1) {
      CCTK_WARN (0, "It is not possible to use a grid variable reduction operator in local mode");
    }
    assert (component == -1);
    
    assert (cgh);
    
    assert (proc == -1 || (proc>=0 && proc<CCTK_nProcs(cgh)));
    
    assert (num_outvals>=0);
    assert (outvals || (proc!=-1 && proc!=CCTK_MyProc(cgh)));
    
    assert (num_invars>=0);
    assert (num_invars == num_outvals);
    assert (invars);
    for (int n=0; n<num_invars; ++n) {
      assert (invars[n]>=0 && invars[n]<CCTK_NumVars());
    }
    
    if (num_outvals==0) return 0;
    
    const int vi = invars[0];
    assert (vi>=0 && vi<CCTK_NumVars());
    
    const int num_dims = CCTK_GroupDimFromVarI(vi);
    assert (num_dims>=0 && num_dims<=dim);
    for (int n=0; n<num_invars; ++n) {
      assert (CCTK_GroupDimFromVarI(invars[n]) == num_dims);
    }
    
    const int vartypesize = CCTK_VarTypeSize(outtype);
    assert (vartypesize>=0);
    
    char *myoutvals = new char [vartypesize * num_outvals];
    char *mycounts  = new char [vartypesize * num_outvals];
    
    Initialise (cgh, proc, num_outvals, myoutvals, outtype, mycounts, red);
    
    BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
      
      int dims[dim], nghostzones[dim];
      ierr = CCTK_GrouplshVI(cgh, num_dims, dims, vi);
      assert (!ierr);
      ierr = CCTK_GroupnghostzonesVI(cgh, num_dims, nghostzones, vi);
      assert (!ierr);
      for (int d=0; d<num_dims; ++d) {
	assert (dims[d]>=0);
	assert (nghostzones[d]>=0 && 2*nghostzones[d]<=dims[d]);
      }
      
      const void** inarrays = new const void * [num_invars];
      for (int n=0; n<num_invars; ++n) {
	inarrays[n] = CCTK_VarDataPtrI(cgh, 0, invars[n]);
	assert (inarrays[n]);
      }
      
      const int intype = CCTK_VarTypeI(vi);
      for (int n=0; n<num_invars; ++n) {
	assert (CCTK_VarTypeI(invars[n]) == intype);
      }
      
      int mydims[dim], mynghostzones[dim];
      for (int d=0; d<num_dims; ++d) {
	mydims[d] = dims[d];
	mynghostzones[d] = nghostzones[d];
      }
      for (int d=num_dims; d<dim; ++d) {
	mydims[d] = 1;
	mynghostzones[d] = 0;
      }
      
      Reduce (cgh, proc, mydims, mynghostzones, num_invars, inarrays, intype,
	      num_outvals, myoutvals, outtype, mycounts, red);
      
      delete [] inarrays;
      
    } END_LOCAL_COMPONENT_LOOP(cgh);
    
    Finalise (cgh, proc, num_outvals, outvals, outtype, myoutvals, mycounts,
	      red);
    
    delete [] myoutvals;
    delete [] mycounts;
    
    return 0;
  }
  
  
  
#define REDUCTION(OP)							 \
  int OP##_arrays (cGH * const cgh, const int proc,			 \
		   const int num_dims, int * const dims,		 \
		   const int num_inarrays,				 \
		   void ** const inarrays, const int intype,		 \
		   const int num_outvals,				 \
		   void * const outvals, const int outtype)		 \
  {									 \
    const OP red;							 \
    return ReduceArrays							 \
      (cgh, proc, num_dims, dims,					 \
       num_inarrays, inarrays, intype, num_outvals, outvals, outtype,	 \
       &red);								 \
  }									 \
  									 \
  int OP##_GVs (cGH * const cgh, const int proc,			 \
	        const int num_outvals,					 \
	        const int outtype, void * const outvals,		 \
	        const int num_invars, int * const invars)		 \
  {									 \
    const OP red;							 \
    return ReduceGVs (cgh, proc,					 \
		      num_outvals, outtype, outvals, num_invars, invars, \
		      &red);						 \
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
  
  
  
  int CarpetReduceStartup ()
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
    
    return 0;
  }
  
} // namespace CarpetReduce
