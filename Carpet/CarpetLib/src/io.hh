/***************************************************************************
                          io.hh  -  I/O routines
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/io.hh,v 1.1 2001/03/01 13:40:10 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef IO_HH
#define IO_HH

#include <cassert>
#include <string>

#include <AMRwriter.hh>
#include <H5IO.hh>
#include <HDFIO.hh>
#include <IEEEIO.hh>
#include <IO.hh>

#include "data.hh"
#include "defs.hh"

template<class T, int D>
class io {
  
  enum iotype { ieee, hdf, h5 };
  
//   template<int DD>
//   void write_ascii (const data<T,D>& d, const string name,
//                     const vect<int,DD>& dirs);
  
//   static void write_ASCII_1D (const data<T,D>& d, const string name,
// 			      const int dir);
//   static void write_ASCII_2D (const data<T,D>& d, const string name,
// 			      const int dir1, const int dir2);
//   static void write_ASCII_3D (const data<T,D>& d, const string name);
  
//   static void write_FlexIO (const data<T,D>& d, const string name,
// 			    const iotype iot);
};



#if defined(TMPL_IMPLICIT)
#  include "io.cc"
#endif

#endif // IO_HH
