/***************************************************************************
                          wave.cc  -  Scalar wave equation application
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/wave.cc,v 1.1 2001/03/01 13:40:10 eschnett Exp $

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
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <strstream>

#include <mpi.h>

#include "vect.hh"

#include "wave.hh"



typedef vect<int,wave::D> ivect;



wave::wave (const int ref_fact, const int ref_levs,
	    const int mg_fact, const int mg_levs,
	    const ivect& lb, const ivect& ub, const ivect& str, const int tstr)
  : ref_fact(ref_fact), ref_cent(vertex_centered), ref_levs(ref_levs),
    mg_fact(mg_fact), mg_cent(vertex_centered), mg_levs(mg_levs),
    global_lb(lb), global_ub(ub), global_str(str),
    global_ext(global_lb,global_ub,global_str),
    global_lg(1,1,1), global_ug(1,1,1),
    time_lb(-2), time_ub(0), time_str(tstr),
    hh(ref_fact,ref_cent,mg_fact,mg_cent,global_ext),
    tt(hh,time_str),
    dd(hh,global_lg,global_ug),
    phi("phi",tt,dd,time_lb,time_ub),
    anaphi("anaphi",tt,dd,time_lb,time_ub),
    output_dir("output")
{
  clog << "WAVE" << endl;
  clog << "hh: " << hh << endl;
  clog << "tt: " << tt << endl;
  clog << "dd: " << dd << endl;
  clog << "phi: " << phi << endl;
  clog << "anaphi: " << anaphi << endl;
}



void wave::init () {
  int size;
  MPI_Comm_size (dist::comm, &size);

  vector<vector<ibbox> > bbss(ref_levs);
  ivect str = global_str;
  ivect lb  = global_lb;
  ivect ub  = global_ub;
  for (int rl=0; rl<ref_levs; ++rl) {
    if (rl>0) {
      // refined boxes have smaller stride
      assert (all(str%ref_fact==0));
      str /= ref_fact;
      // refine (arbitrarily) the center only
      lb /= ref_fact;
      ub /= ref_fact;
    }
    vector<ibbox> bbs(size);
    for (int c=0; c<size; ++c) {
      ivect my_lb = lb;
      ivect my_ub = ub;
      const int zstep = ((ub[D-1] - lb[D-1]) / str[D-1] + size)
                        / size * str[D-1];
      my_lb[D-1] = lb[D-1] + zstep *  c;
      my_ub[D-1] = lb[D-1] + zstep * (c+1) - str[D-1];
      if (c==size-1) my_ub[D-1] = ub[D-1];
      bbs[c] = ibbox(my_lb,my_ub,str);
    }
    bbss[rl] = bbs;
  }
  vector<vector<vector<ibbox> > > bbsss
    = hh.make_multigrid_boxes(bbss, mg_levs);
  vector<vector<int> > pss(bbss.size());
  for (int rl=0; rl<(int)bbss.size(); ++rl) {
    pss[rl] = vector<int>(bbss[rl].size());
    for (int c=0; c<(int)bbss[rl].size(); ++c) {
      pss[rl][c] = c % size;		// distribute among processors
    }
  }
  hh.recompose(bbsss, pss);
  
  clog << "WAVE-INIT" << endl;
  clog << "hh: " << hh << endl;
  clog << "tt: " << tt << endl;
  clog << "dd: " << dd << endl;
  clog << "phi: " << phi << endl;
  clog << "anaphi: " << anaphi << endl;
  
  world_xmin = -5;
  world_xmax = 5;
  world_dx = (world_xmax - world_xmin) / dvect(global_ub - global_lb);
  world_tmin = 0;
  world_tmax = 2;
  cfl_fact = 0.4;
  world_dt = cfl_fact * maxval(world_dx);

  amr_init (0);
}

void wave::amr_init (const int rl) {
  const int tl=0;
  const int ml=0;

  for (int c=0; c<hh.components(rl); ++c) {
    physics_analytic (tl-1,rl,c,ml);
    physics_setphi   (tl-1,rl,c,ml);
  } // for c
  sync (tl-1,rl,ml);

  for (int c=0; c<hh.components(rl); ++c) {
    physics_analytic (tl,rl,c,ml);
    physics_setphi   (tl,rl,c,ml);
  } // for c
  sync (tl,rl,ml);

  // Recursive initialisation of finer grids
  if (rl<hh.reflevels()-1) {
    amr_init (rl+1);
  }
  
  // Check time stepping
  assert (rl==hh.reflevels()-1 || tt.get_time(rl+1,ml)==tt.get_time(rl,ml));

  // Output
  output (tl,rl,ml,output_dir);
}



void wave::update () {
  clog << "WAVE-UPDATE at t=" << tt.get_time(0,0) << endl;
  amr_update (0);
}

void wave::amr_update (const int rl) {
  const int tl=0;
  
  const int ml=0;
  
  // Copy data to old time level and increment time
  for (int c=0; c<hh.components(rl); ++c) {
    anaphi.copy(tl-2,rl,c,ml);
    anaphi.copy(tl-1,rl,c,ml);
    phi   .copy(tl-2,rl,c,ml);
    phi   .copy(tl-1,rl,c,ml);
  }
  tt.advance_time(rl,ml);
  clog << "rl=" << rl << " time=" << tt.get_time(rl,ml) << endl;
  
  // Calculate current time level
  for (int c=0; c<hh.components(rl); ++c) {
    physics_analytic (tl,rl,c,ml);
    physics_update   (tl,rl,c,ml);
    physics_boundary (tl,rl,c,ml);
  }
  sync (tl,rl,ml);
  
  // Recursive integration of finer refinement levels
  if (rl<hh.reflevels()-1) {
    for (int i=0; i<ref_fact; ++i) {
      amr_update (rl+1);
    }
  }
  
  // Restrict from finer levels
  if (rl<hh.reflevels()-1) {
    for (int c=0; c<hh.components(rl); ++c) {
      phi.ref_restrict (tl,rl,c,ml);
    }
    sync (tl,rl,ml);
  }
  
  // Check time stepping
  assert (rl==hh.reflevels()-1 || tt.get_time(rl+1,ml)==tt.get_time(rl,ml));
  
  // Output
  output(tl,rl,ml,output_dir);
}



void wave::sync (const int tl, const int rl, const int ml) {
  if (rl>0) {
    for (int c=0; c<hh.components(rl); ++c) {
      phi.ref_bnd_prolongate (tl,rl,c,ml);
    }
  }
  for (int c=0; c<hh.components(rl); ++c) {
    phi.sync (tl,rl,c,ml);
  }
}



// Output all variables
void wave::output (const int tl, const int rl, const int ml,
		               const char* const dir)
  const
{
  assert (dir!=0);
  
  output_var (phi   , tl, rl, ml, dir);
  output_var (anaphi, tl, rl, ml, dir);
}

// Output one variable
void wave::output_var (const gf<double,D>& var,
		       const int tl, const int rl, const int ml,
		       const char* const dir)
  const
{
  for (int d=0; d<D; ++d) {
    switch (d) {
    case 0: output_var_all_dirs<1> (var, tl, rl, ml, dir); break;
    case 1: output_var_all_dirs<2> (var, tl, rl, ml, dir); break;
    case 2: output_var_all_dirs<3> (var, tl, rl, ml, dir); break;
    default: abort();
    }
  }
}

// Output data of rank DD
template<int DD>
void wave::output_var_all_dirs (const gf<double,D>& var,
				const int tl, const int rl, const int ml,
				const char* const dir)
  const
{
  bbox<int,DD> dirbox(0, D-1, 1);
  for (bbox<int,DD>::iterator bi=dirbox.begin(); bi!=dirbox.end(); ++bi) {
    const vect<int,DD> dirs = *bi;
    bool valid=true;
    for (int d=0; d<DD; ++d)
      for (int dd=d+1; dd<DD; ++dd)
	valid &= dirs[dd] > dirs[d];
    if (valid)
      output_var_one_dir (var, dirs, tl, rl, ml, dir);
  }
}

// Output one rank-DD-slice of data
template<int DD>
void wave::output_var_one_dir (const gf<double,D>& var,
			       const vect<int,DD>& dirs,
			       const int tl, const int rl, const int ml,
			       const char* const dir)
  const
{
  assert (all(dirs>=0 && dirs<D));
  
  const int t = tt.time(tl,rl,ml);
  
  // do not output too often
  if (DD>1 && (t/tt.get_delta(rl,ml)) % 10 != 0) return;
  
  ostrstream name;
  name << dir << "/" << var.name << ".";
  for (int d=0; d<DD; ++d) name << "xyz"[dirs[d]];
  name << ".data" << ends;

  for (int c=0; c<hh.components(rl); ++c) {
    var(tl,rl,c,ml)->write_ascii(name.str(), dirs, t, tl, rl, c, ml);
  }

  name.freeze(0);
}



void wave::physics_analytic (const int tl, const int rl, const int c,
			     const int ml)
{
  if (! hh.is_local(rl,c)) return;

  clog << "analytic tl=" << tl << ", rl=" << rl << ", c=" << c << ", "
       << "ml=" << ml << endl;
  
  const double A      = 1;	// amplitude
  const double radius = 0;	// radius
  const double sigma  = 1;	// width
  
  const double t  = tt.time(tl,rl,ml)   * world_dt + world_tmin;
  const double dt = tt.get_delta(rl,ml) * world_dt;
  clog << "t=" << t << " dt=" << dt << endl;
  
  data<double,D>& ap = *anaphi(tl,rl,c,ml);
  
  const ibbox ext = ap.extent();
  clog << "   bbox=" << ext << endl;
  for (ibbox::iterator it=ext.begin(); it!=ext.end(); ++it) {
    const ivect index = *it;
    
    // equation in cartesian coordinates (t and vector x):
    //    (d2/dt2 - d2/dx2) phi = 0
    // one solution:
    //    phi_k(x,t) = f(t + k x) + g(t - k x)
    // constraint:
    //    abs(k) = 1
    // test:
    //    (d2/dt2 + d2/dx2) phi_k = (1 - k^2) f'' + (1 - k^2) g'' = 0
    
    // equation in D-dim spherical coordinates (only t and r, no Omega):
    //    (d2/dt2 - (D-1)/r d/dr - d2/dr2) phi = 0
    // solution:
    //    phi(r,t) = A/r f(t+r) + B/r f(t-r)
    
    const dvect x = dvect(index - global_lb) * world_dx + world_xmin;
    const double r = hypot(x);
    ap[index] = A * exp(- square((r - t - radius) / sigma));
  }
}



void wave::physics_update (const int tl, const int rl, const int c,
			   const int ml)
{
  if (! hh.is_local(rl,c)) return;

  clog << "updating tl=" << tl << ", rl=" << rl << ", c=" << c << ", "
       << "ml=" << ml << endl;
  
  const data<double,D>& pp = *phi(tl-2,rl,c,ml);
  const data<double,D>& p  = *phi(tl-1,rl,c,ml);
  data<double,D>&       np = *phi(tl  ,rl,c,ml);
  
  const ivect str = p.extent().stride();
  
  const double dt = tt.get_delta(rl,ml) * world_dt;
  const dvect  dx = dvect(str) * world_dx;
  const dvect  dtdx2 = square(dvect(dt) / dx);
  clog << "dt=" << dt << " dx=" << dx << endl;
  
  const ivect lb = p.extent().lower() + str;
  const ivect ub = p.extent().upper() - str;
  const ibbox ext(lb,ub,str);
  clog << "   bbox=" << ext << endl;
  for (ibbox::iterator it=ext.begin(); it!=ext.end(); ++it) {
    const ivect index = *it;
    
    double z = 2 * p[index] - pp[index];
    for (int d=0; d<D; ++d) {
      const ivect offset = ivect::dir(d) * str;
      z += dtdx2[d] * (p[index-offset] - 2*p[index] + p[index+offset]);
    }
    np [index] = z;
  }
}



void wave::physics_boundary (const int tl, const int rl, const int c,
			     const int ml)
{
  if (! hh.is_local(rl,c)) return;

  clog << "boundary conditions tl=" << tl << ", rl=" << rl << ", "
       << "c=" << c << ", ml=" << ml << endl;
  
  const double var0 = 0;	// neutral boundary value
  const double v0   = 1;	// wave speed
  
  const data<double,D>& p   = *phi(tl-1,rl,c,ml);
  data<double,D>&       np  = *phi(tl  ,rl,c,ml);
  
  const ivect str = p.extent().stride();
  
  const double dt = tt.get_delta(rl,ml) * world_dt;
  const dvect  dx = dvect(str) * world_dx;
  
  // Loop over all boundaries
  for (int d=0; d<D; ++d) {
    for (int sign=-1; sign<=+1; sign+=2) {
      
      // Is this an outer boundary?
      if ((sign==-1 && p.extent().lower()[d] < hh.baseextent.lower()[d])
 	  || (sign==+1 && p.extent().upper()[d] > hh.baseextent.upper()[d])) {
	
	const ivect dir = ivect(sign) * ivect::dir(d); // points outwards
	
	// Find bounding box for this boundary
	ivect lb = p.extent().lower();
	ivect ub = p.extent().upper();
	switch (sign) {
	case -1: ub[d] = lb[d]; break;
	case  1: lb[d] = ub[d]; break;
	}
	for (int dd=d+1; dd<D; ++dd) {
	  lb[dd] += str[dd];
	  ub[dd] -= str[dd];
	}
	const ibbox ext(lb,ub,str);
	
	// Walk bounding box
 	clog << "   bbox=" << ext << endl;
	for (ibbox::iterator it=ext.begin(); it!=ext.end(); ++it) {
	  const ivect index0 = *it;
	  const ivect index1 = index0 - dir * str;
	  
	  const dvect x0 = dvect(index0 - global_lb) * world_dx + world_xmin;
	  const dvect x1 = dvect(index1 - global_lb) * world_dx + world_xmin;
	  const double r0 = hypot(x0);
	  const double r1 = hypot(x1);
	  const double varp0 = p [index0];
	  const double varp1 = p [index1];
	  const double varn1 = np[index1];
	  //
	  // f = f0 + u(r - v0*t) / r
	  //
	  // (xi/r) df/dt + v0*df/dxi + v0*xi*(f-f0)/r^2 = 0
	  //
	  np[index0]
	    = ((v0*dt*var0 * (x0[d]/square(r0) + x1[d]/square(r1))
		+ varn1 * (sign*v0*dt/dx[d] - x1[d]/r1 * (1+0.5*v0*dt/r1))
		+ varp1 * (sign*v0*dt/dx[d] + x1[d]/r1 * (1-0.5*v0*dt/r1))
		- varp0 * (sign*v0*dt/dx[d] - x0[d]/r0 * (1-0.5*v0*dt/r0)))
	       / (sign*v0*dt/dx[d] + x0[d]/r0 * (1+0.5*v0*dt/r0)));
	}
	
      } // if is-outer-boundary
      
    } // for sign
  } // for d
}



void wave::physics_setphi (const int tl, const int rl, const int c,
			   const int ml)
{
  if (! hh.is_local(rl,c)) return;

  const data<double,D>& ap = *anaphi(tl,rl,c,ml);
  data<double,D>&       p  = *phi   (tl,rl,c,ml);
  
  const ibbox ext = p.extent();
  for (ibbox::iterator it=ext.begin(); it!=ext.end(); ++it) {
    const ivect index = *it;
    p[index] = ap[index];
  }
}
