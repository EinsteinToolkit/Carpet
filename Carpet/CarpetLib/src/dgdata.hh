/***************************************************************************
                          dgdata.hh  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dgdata.hh,v 1.1 2001/06/12 14:56:57 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DGDATA_HH
#define DGDATA_HH

#include <assert.h>

#include <iostream>

#include "defs.hh"

using namespace std;



// Forward declaration
class dimgeneric_data;

// Output
ostream& operator<< (ostream& os, const dimgeneric_data& d);



// A generic data storage without type information
class dimgeneric_data {

protected:                      // should be readonly

  // Fields
  bool _has_storage;		// has storage associated (on some processor)
  bool _owns_storage;		// owns the storage
  int _size;			// size

  int _proc;			// stored on processor

public:

  // Constructors
  dimgeneric_data ();

  // Destructors
  virtual ~dimgeneric_data ();
  
  // Processor management
  virtual void change_processor (const int newproc, void* const mem=0) = 0;

  // Accessors
  bool has_storage () const {
    return _has_storage;
  }
  bool owns_storage () const {
    assert (_has_storage);
    return _owns_storage;
  }
  
  virtual const void* storage () const = 0;
  
  virtual void* storage () = 0;
  
  int size () const {
    assert (_has_storage);
    return _size;
  }

  int proc () const {
    assert (_has_storage);
    return _proc;
  }
  
  // Output
  virtual ostream& output (ostream& os) const = 0;
};



inline ostream& operator<< (ostream& os, const dimgeneric_data& d) {
  return d.output(os);
}



#if defined(TMPL_IMPLICIT)
#  include "dgdata.cc"
#endif

#endif // DGDATA_HH
