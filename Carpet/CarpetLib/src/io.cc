/***************************************************************************
                          io.cc  -  I/O routines
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/io.cc,v 1.1 2001/03/01 13:40:10 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "io.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GF_HH)
#  include "io.hh"
#endif

template<class TT>
IObase::DataType datatype(TT dummy)
{ assert(false); return IObase::DataType(); }

template <>
IObase::DataType datatype(char* dummy) { return IObase::uInt8; }
template <>
IObase::DataType datatype(short* dummy) { return IObase::Int16; }
template <>
IObase::DataType datatype(int* dummy) { return IObase::Int32; }
template <>
IObase::DataType datatype(long long* dummy) { return IObase::Int64; }
template <>
IObase::DataType datatype(float* dummy) { return IObase::Float32; }
template <>
IObase::DataType datatype(double* dummy) { return IObase::Float64; }
template <>
IObase::DataType datatype(long double* dummy) { return IObase::Float64; }
