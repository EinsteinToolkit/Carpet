// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/defines.hh,v 1.2 2004/03/02 10:29:54 schnetter Exp $

#ifndef DEFINES_HH
#define DEFINES_HH

#include "cctk.h"

#include <bbox.hh>
#include <bboxset.hh>
#include <vect.hh>



namespace Carpet {
  
  const int dim = 3;
  
  typedef vect<bool,dim> bvect;
  typedef vect<int,dim> ivect;
  typedef vect<CCTK_INT,dim> jvect;
  typedef vect<CCTK_REAL,dim> rvect;
  typedef bbox<int,dim> ibbox;
  typedef bbox<CCTK_INT,dim> jbbox;
  typedef bbox<CCTK_REAL,dim> rbbox;
  typedef bboxset<int,dim> ibset;
  
  typedef vect<vect<bool,2>,dim> bbvect;
  typedef vect<vect<int,2>,dim> iivect;
  typedef vect<vect<CCTK_INT,2>,dim> jjvect;
  
} // namespace Carpet

#endif // !defined(DEFINES_HH)
