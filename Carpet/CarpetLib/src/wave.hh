/***************************************************************************
                          wave.hh  -  Scalar wave equation application
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/wave.hh,v 1.1 2001/03/01 13:40:10 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef WAVE_HH
#define WAVE_HH

#include <iostream>

#include "bbox.hh"
#include "data.hh"
#include "defs.hh"
#include "dh.hh"
#include "gf.hh"
#include "gh.hh"
#include "th.hh"
#include "vect.hh"

class wave {
  
public:
  // Dimensions
  static const int D=3;
  
private:
  // Types
  typedef vect<int,D> ivect;
  typedef bbox<int,D> ibbox;
  
  // Refinement
  int ref_fact;
  centering ref_cent;
  int ref_levs;
  
  // Multigrid
  int mg_fact;
  centering mg_cent;
  int mg_levs;
  
  // Extent
  ivect global_lb, global_ub, global_str;
  ibbox global_ext;
  
  // Ghosts
  ivect global_lg, global_ug;
  
  // Time
  int time_lb, time_ub, time_str;
  
  // Data
  gh<D> hh;
  th<D> tt;
  dh<D> dd;
  gf<double,D> phi, anaphi;
  
  // World
  typedef vect<double,D> dvect;
  dvect  world_xmin, world_xmax, world_dx;
  double cfl_fact;
  double world_tmin, world_tmax, world_dt;
  
  const char* output_dir;
  
public:
  wave (const int ref_fact, const int ref_levs,
	const int mg_fact, const int mg_levs,
	const ivect& lb, const ivect& ub, const ivect& str, const int tstr);
  void init   ();
  void update ();
  
  void output (const int tl, const int rl, const int ml, const char* const dir)
    const;
  void output_var (const gf<double,D>& var,
		   const int tl, const int rl, const int ml,
		   const char* const dir)
    const;
  template<int DD>
  void output_var_all_dirs (const gf<double,D>& var,
			    const int tl, const int rl, const int ml,
			    const char* const dir)
    const;
  template<int DD>
  void output_var_one_dir (const gf<double,D>& var,
			   const vect<int,DD>& dirs,
			   const int tl, const int rl, const int ml,
			   const char* const dir)
    const;
  
protected:
  void amr_init   (const int rl);
  void amr_update (const int rl);
  
  void sync (const int tl, const int rl, const int ml);
  
  void physics_analytic  (const int tl, const int rl, const int c,
			  const int ml);
  void physics_update    (const int tl, const int rl, const int c,
			  const int ml);
  void physics_boundary  (const int tl, const int rl, const int c,
			  const int ml);
  void physics_setphi    (const int tl, const int rl, const int c,
			  const int ml);
  
};

#endif // WAVE_HH
