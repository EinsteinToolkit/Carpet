// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gf.hh,v 1.12 2004/03/23 12:40:27 schnetter Exp $

#ifndef GF_HH
#define GF_HH

#include <assert.h>
#include <math.h>

#include <iostream>
#include <string>

#include "bbox.hh"
#include "bboxset.hh"
#include "data.hh"
#include "defs.hh"
#include "dh.hh"
#include "ggf.hh"
#include "th.hh"
#include "vect.hh"

using namespace std;



// A real grid function
template<class T,int D>
class gf: public ggf<D> {
  
  // Types
  typedef vect<int,D>    ivect;
  typedef bbox<int,D>    ibbox;
  typedef bboxset<int,D> ibset;
  typedef list<ibbox>    iblist;
  typedef vector<iblist> iblistvect;
  
  typedef data<T,D>*    tdata;	        // data ...
  typedef vector<tdata> mdata;	        // ... for each multigrid level
  typedef vector<mdata> cdata;	        // ... for each component
  typedef vector<cdata> rdata;	        // ... for each refinement level
  typedef vector<rdata> fdata;          // ... for each time level

public:
  
  // Constructors
  gf (const int varindex, const operator_type transport_operator,
      th<D>& t, dh<D>& d,
      const int tmin, const int tmax, const int prolongation_order_time,
      const int vectorlength, const int vectorindex,
      gf* const vectorleader);
  
  // Destructors
  virtual ~gf ();
  
  
  
  // Helpers
  
protected:
  
  virtual gdata<D>* typed_data (int tl, int rl, int c, int ml)
  {
    return new data<T,D>(this->varindex, this->transport_operator,
                         this->vectorlength, this->vectorindex,
                         this->vectorleader
                         ? (data<T,D>*)(*this->vectorleader)(tl,rl,c,ml)
                         : NULL);
  }
  
  
  
  // Access to the data
  
public:
  
  virtual const data<T,D>* operator() (int tl, int rl, int c, int ml) const;
  
  virtual data<T,D>* operator() (int tl, int rl, int c, int ml);
  
  
  
  // Output
  virtual ostream& output (ostream& os) const;
};



#endif // GF_HH
