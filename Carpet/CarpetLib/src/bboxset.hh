#ifndef BBOXSET_HH
#define BBOXSET_HH

#include "bboxset1.hh"
#ifndef CARPET_NO_BBOXSET2
#  include "bboxset2.hh"
#endif

#ifdef CARPET_BBOXSET2
using namespace bboxset2;
#else
using namespace bboxset1;
#endif

#endif  // #ifndef BBOXSET_HH
