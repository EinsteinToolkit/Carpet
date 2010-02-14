#ifndef GF_HH
#define GF_HH

#include <cassert>
#include <cmath>
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
template<typename T>
class gf: public ggf {
  
  // Types
  typedef data<T>*      tdata;	        // data ...
  typedef vector<tdata> mdata;	        // ... for each multigrid level
  typedef vector<mdata> cdata;	        // ... for each component
  typedef vector<cdata> rdata;	        // ... for each refinement level
  typedef vector<rdata> fdata;          // ... for each time level

public:
  
  // Constructors
  gf (const int varindex, const operator_type transport_operator,
      th& t, dh& d,
      const int prolongation_order_time,
      const int vectorlength, const int vectorindex,
      gf* const vectorleader);
  
  // Destructors
  virtual ~gf ();
  
  
  
  // Helpers
  
  virtual gdata* typed_data (int tl, int rl, int lc, int ml) const
  {
    data<T>* const vl =
      this->vectorleader
      ? (data<T>*)(*this->vectorleader)(tl,rl,lc,ml)
      : NULL;
    return new data<T>(this->varindex,
                       h.refcent, this->transport_operator,
                       this->vectorlength, this->vectorindex,
                       vl);
  }
  
  virtual gdata* new_typed_data () const
  {
    return new data<T>(this->varindex,
                       h.refcent, this->transport_operator,
                       1, 0, NULL);
  }
  
  
  
  // Access to the data
  
  virtual const data<T>* operator() (int tl, int rl, int lc, int ml) const;
  
  virtual data<T>* operator() (int tl, int rl, int lc, int ml);
  
  
  
  // Output
  virtual size_t memory () const;
  virtual ostream& output (ostream& os) const;
private:
  gf ();                        // canonical default construtor
  gf (const gf &);              // canonical copy construtor
  gf & operator= (const gf &);  // canonical copy

};



#endif // GF_HH
