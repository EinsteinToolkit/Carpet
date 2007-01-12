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
  assert (all (baseextent.stride() % reffacts.at(reffacts.size()-1) == 0));
}

// Destructors
gh::~gh () { }

// Modifiers
void gh::regrid (mregs const & regs)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Save the old grid hierarchy
  _oldregions = _regions;
  _regions = regs;
  
  // Consistency checks
  
  // nota bene: there might be 0 refinement levels.
  
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
  if (_oldregions.size() != _regions.size()) {
    _oldregions.resize (_regions.size());
  }
  
  bool const do_recompose = level_did_change(rl);
  
  if (do_recompose) {
    
    // Recompose the other hierarchies
    for (list<dh*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
      (*d)->recompose (rl, do_prolongate);
    }
    
    // Overwrite old with new grid hierarchy
    for (int ml=0; ml<mglevels(); ++ml) {
      _oldregions.at(ml).resize (_regions.at(ml).size());
      _oldregions.at(ml).at(rl) = _regions.at(ml).at(rl);
    }
    
  }
  
  return do_recompose;
}

bool gh::level_did_change (const int rl) const
{
  // Find out whether this level changed
  if (_regions.size() != _oldregions.size()) return true;
  for (int ml=0; ml<mglevels(); ++ml) {
    assert (rl>=0 and rl<reflevels());
    if (rl >= (int)_oldregions.at(ml).size()) return true;
    if (_regions.at(ml).at(rl).size() != _oldregions.at(ml).at(rl).size()) {
      return true;
    }
    for (int c=0; c<components(rl); ++c) {
      if (_regions.at(ml).at(rl).at(c) != _oldregions.at(ml).at(rl).at(c)) {
        return true;
      }
    } // for c
  } // for ml
  return false;
}
  
void gh::check_multigrid_consistency ()
{
  assert (mglevels()>0);
  for (int ml=1; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
	assert (all(extent(ml,rl,c).stride()
		    == ivect(mgfact) * extent(ml-1,rl,c).stride()));
        // TODO: put the check back in, taking outer boundaries into
        // account
#if 0
	assert (extent(ml,rl,c)
		.contracted_for(extent(ml-1,rl,c))
		.is_contained_in(extent(ml-1,rl,c)));
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
        const ibbox &b  = extent(ml,rl,c);
        const ibbox &b0 = extent(ml,rl,0);
	assert (all(b.stride() == b0.stride()));
	assert (b.is_aligned_with(b0));
        for (int cc=c+1; cc<components(rl); ++cc) {
          assert ((b & extent(ml,rl,cc)).empty());
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
      assert (extent(0,c,0).is_contained_in(baseextent));
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
      assert (all(extent(ml,rl-1,0).stride()
                  == ((reffacts.at(rl) / reffacts.at(rl-1))
                      * extent(ml,rl,0).stride())));
      // Check contained-ness:
      // first take all coarse grids ...
      ibset all;
      for (int c=0; c<components(rl-1); ++c) {
        all |= extent(ml,rl-1,c);
      }
#if 0
      // ... remember their size ...
      const int sz = all.size();
      // ... then add the coarsified fine grids ...
      for (int c=0; c<components(rl); ++c) {
        all |= extent(ml,rl,c).contracted_for(extent(ml,rl-1,0));
      }
      // ... and then check the sizes:
      assert (all.size() == sz);
#endif
      for (int c=0; c<components(rl); ++c) {
        ibset const finebox
          = (extent(ml,rl,c).contracted_for (extent(ml,rl-1,0)));
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
               << extent(ml,rl,c) << endl;
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
        bb = bb.expanded_containing(extent(ml,rl,c));
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
        os << "extent=" << extent(ml,rl,c) << endl;
        os << "outer_boundaries=" << outer_boundaries(rl,c) << endl;
        os << "refinement_boundaries=" << refinement_boundaries(rl,c) << endl;
        os << "processor=" << processor(rl,c) << endl;
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
     << "regions=" << regions() << ","
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
