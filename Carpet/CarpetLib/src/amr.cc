/***************************************************************************
                          amr.cc  -  Wavetoy example driver
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/amr.cc,v 1.1 2001/03/01 13:40:10 eschnett Exp $

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
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "dist.hh"
#include "vect.hh"

#include "wave.hh"

int main (int argc, char *argv[]) {

  dist::init (argc, argv);

  typedef vect<int,wave::D> ivect;

  if (argc != 2) {
    cerr << "Error: wrong arguments" << endl
	 << "Synopsis: " << argv[0] << " <number of refinement levels>"
	 << endl;
    exit(1);
  }

  const int ref_fact = 2;
  const int ref_levs = atoi(argv[1]);
  assert (ref_levs>0);
  const int mg_fact = 2;
  const int mg_levs = 1;
  const ivect lb(-24), ub(24), str((int)rint(pow(ref_fact, ref_levs-1)));
  const int tstr = (int)rint(pow(ref_fact, ref_levs-1));
  const int nsteps = 16 / tstr;

  wave w(ref_fact, ref_levs, mg_fact, mg_levs, lb, ub, str, tstr);
  w.init();
  for (int i=0; i<nsteps; ++i) {
    w.update();
  }

  dist::finalize ();

  return 0;
}
