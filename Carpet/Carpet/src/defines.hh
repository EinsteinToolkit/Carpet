#ifndef DEFINES_HH
#define DEFINES_HH

#include "cctk.h"

#include <defs.hh>



namespace Carpet {
  
  typedef vect<CCTK_INT,dim> jvect;
  typedef vect<CCTK_REAL,dim> rvect;
  typedef bbox<CCTK_INT,dim> jbbox;
  typedef bbox<CCTK_REAL,dim> rbbox;
  
  typedef vect<vect<CCTK_INT,2>,dim> jjvect;
  
} // namespace Carpet

#endif // !defined(DEFINES_HH)
