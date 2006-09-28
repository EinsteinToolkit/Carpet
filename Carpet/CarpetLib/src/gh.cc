#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dh.hh"
#include "th.hh"
#include "vect.hh"

#include "gh.hh"

using namespace std;



  // Constructors
gh::gh (const vector<ivect> & reffacts_, const centering refcent_,
        const int mgfact_, const centering mgcent_,
        const ibbox baseextent_)
  : reffacts(reffacts_), refcent(refcent_),
    mgfact(mgfact_), mgcent(mgcent_),
    baseextent(baseextent_)
{
  assert (reffacts.size() >= 1);
  assert (all (reffacts.front() == 1));
  for (size_t n = 1; n < reffacts.size(); ++ n) {
    assert (all (reffacts.at(n) >= reffacts.at(n-1)));
    assert (all (reffacts.at(n) % reffacts.at(n-1) == 0));
  }
}

// Destructors
gh::~gh () { }

// Modifiers
void gh::regrid (const mexts& exts,
                 const rbnds& outer_bounds,
                 const rprocs& procs)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Save the old grid hierarchy
  _oldextents = _extents;
  _oldouter_boundaries = _outer_boundaries;
  _oldprocessors = _processors;
  
  _extents = exts;
  _outer_boundaries = outer_bounds;
  _processors = procs;
  
  // Consistency checks
  
  // nota bene: there might be 0 refinement levels.
  
  check_processor_number_consistency ();
  check_multigrid_consistency ();
  check_component_consistency ();
  check_base_grid_extent ();
  check_refinement_levels ();
  
  calculate_base_extents_of_all_levels ();

  if (output_bboxes) {
    do_output_bboxes (cout);
    do_output_bases (cout);
  }
  
  // Recompose the other hierarchies
  for (list<th*>::iterator t=ths.begin(); t!=ths.end(); ++t) {
    (*t)->regrid();
  }
  
  for (list<dh*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
    (*d)->regrid();
  }
}

bool gh::recompose (const int rl,
                    const bool do_prolongate)
{
  // Handle changes in number of mglevels
  if (_oldextents.size() != _extents.size()) {
    _oldextents.resize (_extents.size());
  }
  
  bool const do_recompose = level_did_change(rl);
  
  if (do_recompose) {
    
    // Recompose the other hierarchies
    for (list<dh*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
      (*d)->recompose (rl, do_prolongate);
    }
    
    // Overwrite old with new grid hierarchy
    for (int ml=0; ml<mglevels(); ++ml) {
      _oldextents.at(ml).resize (_extents.at(ml).size());
      _oldextents.at(ml).at(rl) = _extents.at(ml).at(rl);
    }
    _oldouter_boundaries.resize (_outer_boundaries.size());
    _oldouter_boundaries.at(rl) = _outer_boundaries.at(rl);
    _oldprocessors.resize (_processors.size());
    _oldprocessors.at(rl) = _processors.at(rl);
  }
  
  return do_recompose;
}

bool gh::level_did_change (const int rl) const
{
  // Find out whether this level changed
  if (_extents.size() != _oldextents.size()) return true;
  for (int ml=0; ml<mglevels(); ++ml) {
    assert (rl>=0 and rl<reflevels());
    if (rl >= (int)_oldextents.at(ml).size()) return true;
    if (_extents.at(ml).at(rl).size() != _oldextents.at(ml).at(rl).size()) {
      return true;
    }
    for (int c=0; c<components(rl); ++c) {
      if (_extents.at(ml).at(rl).at(c) != _oldextents.at(ml).at(rl).at(c)) {
        return true;
      }
      if (_processors.at(rl).at(c) != _oldprocessors.at(rl).at(c)) return true;
    } // for c
  } // for ml
  return false;
}

void gh::check_processor_number_consistency ()
{
  for (int rl=0; rl<reflevels(); ++rl) {
    assert (processors().size() == extents().at(0).size());
    assert (outer_boundaries().size() == extents().at(0).size());
    for (int c=0; c<components(rl); ++c) {
      assert (processors().at(rl).size() == extents().at(0).at(rl).size());
      assert (outer_boundaries().at(rl).size() == extents().at(0).at(rl).size());
    }
  }
}
  
void gh::check_multigrid_consistency ()
{
  assert (mglevels()>0);
  for (int ml=1; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
	assert (all(extents().at(ml).at(rl).at(c).stride()
		    == ivect(mgfact) * extents().at(ml-1).at(rl).at(c).stride()));
        // TODO: put the check back in, taking outer boundaries into
        // account
#if 0
	assert (extents().at(ml).at(rl).at(c)
		.contracted_for(extents().at(ml-1).at(rl).at(c))
		.is_contained_in(extents().at(ml-1).at(rl).at(c)));
#endif
      }
    }
  }
}
  
void gh::check_component_consistency ()
{
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      assert (components(rl)>0);
      for (int c=0; c<components(rl); ++c) {
        const ibbox &b  = extents().at(ml).at(rl).at(c);
        const ibbox &b0 = extents().at(ml).at(rl).at(0);
	assert (all(b.stride() == b0.stride()));
	assert (b.is_aligned_with(b0));
        for (int cc=c+1; cc<components(rl); ++cc) {
          assert ((b & extents().at(ml).at(rl).at(cc)).empty());
        }
      }
    }
  }
}

void gh::check_base_grid_extent ()
{
  if (reflevels()>0) {
    for (int c=0; c<components(0); ++c) {
      // TODO: put the check back in, taking outer boundaries into
      // account
#if 0
      assert (extents().at(0).at(c).at(0).is_contained_in(baseextent));
#endif
    }
  }
}
  
// Check proper nesting, i.e., whether the fine grids are contained in
// the coarser grids
void gh::check_refinement_levels ()
{
  bool have_error = false;
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=1; rl<reflevels(); ++rl) {
      assert (all(extents().at(ml).at(rl-1).at(0).stride()
                  == ((reffacts.at(rl) / reffacts.at(rl-1))
                      * extents().at(ml).at(rl).at(0).stride())));
      // Check contained-ness:
      // first take all coarse grids ...
      ibset all;
      for (int c=0; c<components(rl-1); ++c) {
        all |= extents().at(ml).at(rl-1).at(c);
      }
#if 0
      // ... remember their size ...
      const int sz = all.size();
      // ... then add the coarsified fine grids ...
      for (int c=0; c<components(rl); ++c) {
        all |= extents().at(ml).at(rl).at(c).contracted_for(extents().at(ml).at(rl-1).at(0));
      }
      // ... and then check the sizes:
      assert (all.size() == sz);
#endif
      for (int c=0; c<components(rl); ++c) {
        ibset const finebox
          = (extents().at(ml).at(rl).at(c).contracted_for
             (extents().at(ml).at(rl-1).at(0)));
        if (! (finebox <= all)) {
          if (! have_error) {
            cout << "The following components are not properly nested, i.e., they are not" << endl
                 << "contained within their next coarser level's components:" << endl;
          }
          cout << "   ml " << ml << " rl " << rl << " c " << c << endl;
          have_error = true;
        }
      }      
    }
  }
  if (have_error) {
    cout << "The grid hierarchy looks like:" << endl;
    for (int ml=0; ml<mglevels(); ++ml) {
      for (int rl=0; rl<reflevels(); ++rl) {
        for (int c=0; c<components(rl); ++c) {
          cout << "   ml " << ml << " rl " << rl << " c " << c << ":   "
               << extents().at(ml).at(rl).at(c) << endl;
        }
      }
    }
    CCTK_WARN (0, "The refinement hierarchy is not properly nested.");
  }
}

void gh::calculate_base_extents_of_all_levels ()
{
  _bases.resize(mglevels());
  for (int ml=0; ml<mglevels(); ++ml) {
    _bases.at(ml).resize(reflevels());
    for (int rl=0; rl<reflevels(); ++rl) {
      _bases.at(ml).at(rl) = ibbox();
      ibbox &bb = _bases.at(ml).at(rl);
      for (int c=0; c<components(rl); ++c) {
        bb = bb.expanded_containing(extents().at(ml).at(rl).at(c));
      }
    }
  }
}

// Accessors
int gh::local_components (const int rl) const
{
  int lc = 0;
  for (int c=0; c<components(rl); ++c) {
    if (is_local(rl,c)) ++lc;
  }
  return lc;
}



// Time hierarchy management
void gh::add (th* t)
{
  ths.push_back(t);
}

void gh::remove (th* t)
{
  ths.remove(t);
}



// Data hierarchy management
void gh::add (dh* d)
{
  dhs.push_back(d);
}

void gh::remove (dh* d)
{
  dhs.remove(d);
}


void gh::do_output_bboxes (ostream& os) const
{
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
        os << endl;
        os << "gh bboxes:" << endl;
        os << "ml=" << ml << " rl=" << rl << " c=" << c << endl;
        os << "extent=" << extents().at(ml).at(rl).at(c) << endl;
        os << "outer_boundary=" << outer_boundaries().at(rl).at(c) << endl;
        os << "processor=" << processors().at(rl).at(c) << endl;
      }
    }
  }
}

void gh::do_output_bases (ostream& os) const
{
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      if (components(rl)>0) {
        os << endl;
        os << "gh bases:" << endl;
        os << "ml=" << ml << " rl=" << rl << endl;
        os << "base=" << bases().at(ml).at(rl) << endl;
      }
    }
  }
}

ostream& gh::output (ostream& os) const
{
  os << "gh:"
     << "reffacts=" << reffacts << ",refcentering=" << refcent << ","
     << "mgfactor=" << mgfact << ",mgcentering=" << mgcent << ","
     << "extents=" << extents() << ","
     << "outer_boundaries=" << outer_boundaries() << ","
     << "processors=" << processors() << ","
     << "dhs={";
  const char * sep = "";
  for (list<dh*>::const_iterator d = dhs.begin();
       d != dhs.end(); ++d)
  {
    os << sep;
    (*d)->output(os);
    sep = ",";
  }
  os << "}";
  return os;
}
