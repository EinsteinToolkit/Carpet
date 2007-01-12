#include <iostream>

#include "defs.hh"
#include "region.hh"

using namespace std;



bool
operator== (region_t const & a, region_t const & b)
{
  return
    a.extent == b.extent and
    all (all (a.outer_boundaries == b.outer_boundaries)) and
    all (all (a.refinement_boundaries == b.refinement_boundaries)) and
    a.map == b.map and
    a.processor == b.processor;
}



istream &
operator>> (istream & is, region_t & reg)
{
  skipws (is);
  consume (is, "region_t");
  skipws (is);
  consume (is, '(');
  
  skipws (is);
  consume (is, "extent");
  skipws (is);
  consume (is, '=');
  is >> reg.extent;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "outer_boundaries");
  skipws (is);
  consume (is, '=');
  is >> reg.outer_boundaries;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "refinement_boundaries");
  skipws (is);
  consume (is, '=');
  is >> reg.refinement_boundaries;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "map");
  skipws (is);
  consume (is, '=');
  is >> reg.map;
  skipws (is);
  consume (is, ',');
  
  skipws (is);
  consume (is, "processor");
  skipws (is);
  consume (is, '=');
  is >> reg.processor;
  skipws (is);
  consume (is, ')');
  
  return is;
}



ostream &
operator<< (ostream & os, region_t const & reg)
{
  os << "region_t("
     << "extent=" << reg.extent << ","
     << "outer_boundaries=" << reg.outer_boundaries << ","
     << "refinement_boundaries=" << reg.refinement_boundaries << ","
     << "map=" << reg.map << ","
     << "processor=" << reg.processor << ")";
  return os;
}
