/***************************************************************************
                          gdata.cc  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.cc,v 1.2 2001/03/05 14:31:03 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <cassert>
#include <fstream>

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GDATA_HH)
#  include "gdata.hh"
#endif



// Constructors
template<int D>
generic_data<D>::generic_data ()
  : _has_storage(false)
{ }

// Destructors
template<int D>
generic_data<D>::~generic_data () { }



// Output
template<int D>
template<int DD>
void generic_data<D>::write_ascii (const string name, const int time,
				   const vect<int,D>& org,
				   const vect<int,DD>& dirs,
				   const int tl, const int rl,
				   const int c, const int ml)
  const
{
  assert (_has_storage);
  
  if (_proc==0) {
    // output on processor 0
    
    int rank;
    MPI_Comm_rank (dist::comm, &rank);
    if (rank == 0) {
      
      ofstream file(name.c_str(), ios::app);
      assert (file.good());
      
      file << "# iteration " << time << endl
	   << "# time level " << tl << "   refinement level " << rl
	   << "   component " << c << "   multigrid level " << ml << endl
	   << "# column format: it tl rl c ml";
      assert (D<=3);
      for (int d=0; d<D; ++d) file << " " << "xyz"[d];
      file << " data" << endl;
      
      const vect<int,DD> lo = extent().lower()[dirs];
      const vect<int,DD> up = extent().upper()[dirs];
      const vect<int,DD> str = extent().stride()[dirs];
      const bbox<int,DD> ext(lo,up,str);
      
      for (bbox<int,DD>::iterator it=ext.begin(); it!=ext.end(); ++it) {
        ivect index(org);
        for (int d=0; d<DD; ++d) index[dirs[d]] = (*it)[d];
        file << time << " " << tl << " " << rl << " " << c << " " << ml << " ";
        for (int d=0; d<D; ++d) file << index[d] << " ";
	write_ascii_output_element (file, index);
        for (int d=0; d<D; ++d) {
          if (index[d]!=extent().upper()[d]) break;
          file << endl;
        }
      }
      
      file.close();
      assert (file.good());
      
    }
    
  } else {
    // copy to processor 0 and output there
    
    generic_data* tmp = make_typed(_extent, 0);
    tmp->copy_from (this, _extent);
    tmp->write_ascii (name, time, org, dirs, tl, rl, c, ml);
    delete tmp;
    
  }
}



template<int D>
ostream& operator<< (ostream& os, const generic_data<D>& f) {
  return f.out(os);
}



#if defined(TMPL_EXPLICIT)
template class generic_data<1>;
template ostream& operator<< (ostream& os, const generic_data<1>& d);

template class generic_data<2>;
template ostream& operator<< (ostream& os, const generic_data<2>& d);

template class generic_data<3>;
template ostream& operator<< (ostream& os, const generic_data<3>& d);

template void generic_data<1>::write_ascii
(const string name, const int time,
 const vect<int,1>& org, const vect<int,1>& dirs,
 const int tl, const int rl, const int c, const int ml) const;

template void generic_data<2>::write_ascii
(const string name, const int time,
 const vect<int,2>& org, const vect<int,1>& dirs,
 const int tl, const int rl, const int c, const int ml) const;
template void generic_data<2>::write_ascii
(const string name, const int time,
 const vect<int,2>& org, const vect<int,2>& dirs,
 const int tl, const int rl, const int c, const int ml) const;

template void generic_data<3>::write_ascii
 (const string name, const int time,
  const vect<int,3>& org, const vect<int,1>& dirs,
  const int tl, const int rl, const int c, const int ml) const;
template void generic_data<3>::write_ascii
(const string name, const int time,
 const vect<int,3>& org, const vect<int,2>& dirs,
 const int tl, const int rl, const int c, const int ml) const;
template void generic_data<3>::write_ascii
(const string name, const int time,
 const vect<int,3>& org, const vect<int,3>& dirs,
 const int tl, const int rl, const int c, const int ml) const;
#endif
