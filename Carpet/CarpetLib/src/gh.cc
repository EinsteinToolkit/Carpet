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

using namespace CarpetLib;



  // Constructors
gh::
gh (vector<ivect> const & reffacts_, centering const refcent_,
    int const mgfact_,  centering const mgcent_,
    vector<vector<ibbox> > const & baseextents_,
    i2vect const & boundary_width_)
  : reffacts(reffacts_), refcent(refcent_),
    mgfact(mgfact_), mgcent(mgcent_),
    baseextents(baseextents_),
    boundary_width(boundary_width_)
{
  assert (reffacts.size() >= 1);
  assert (all (reffacts.front() == 1));
  for (int rl = 1; rl < (int)reffacts.size(); ++ rl) {
    assert (all (reffacts.AT(rl) >= reffacts.AT(rl-1)));
    assert (all (reffacts.AT(rl) % reffacts.AT(rl-1) == 0));
  }
  assert (refcent == vertex_centered or refcent == cell_centered);
  
  assert (mgfact >= 1);
  assert (mgcent == vertex_centered or mgcent == cell_centered);
  
  assert (baseextents.size() >= 1);
  assert (baseextents.AT(0).size() >= 1);
  assert (baseextents.AT(0).size() == reffacts.size());
  for (int ml = 1; ml < (int)baseextents.size(); ++ ml) {
    assert (baseextents.AT(ml).size() == baseextents.AT(ml-1).size());
  }
  for (int ml = 0; ml < (int)baseextents.size(); ++ ml) {
    for (int rl = 1; rl < (int)baseextents.AT(ml).size(); ++ rl) {
      ibbox const & cbox = baseextents.AT(ml).AT(rl-1);
      ibbox const & fbox = baseextents.AT(ml).AT(rl);
      assert (all (cbox.stride() * reffacts.AT(rl-1) == 
                   fbox.stride() * reffacts.AT(rl)));
    }
  }
  
  assert (all (all (boundary_width >= 0)));
  for (int ml = 0; ml < (int)baseextents.size(); ++ ml) {
    for (int rl = 0; rl < (int)baseextents.AT(ml).size(); ++ rl) {
      ibbox const & box = baseextents.AT(ml).AT(rl);
      // This condition must hold even for zero-sized grid arrays
      assert (all (box.shape() / box.stride() >=
                   boundary_width[0] + boundary_width[1]));
    }
  }
}

// Destructors
gh::
~gh ()
{
}



// Modifiers
void
gh::
regrid (rregs const & superregs, mregs const & regs)
{
  DECLARE_CCTK_PARAMETERS;
  
  superregions = superregs;
  
  // Save the grid hierarchy
  oldregions.clear ();
  swap (oldregions, regions);
  regions = regs;
  
  
  
  // Consistency checks
  
  // Note: there may be zero refinement levels
  
  // Check multigrid consistency
  assert (mglevels()>0);
  for (int ml=1; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
	assert (all(extent(ml,rl,c).stride() ==
                    mgfact * extent(ml-1,rl,c).stride()));
	assert (extent(ml,rl,c)
		.contracted_for(extent(ml-1,rl,c))
		.is_contained_in(extent(ml-1,rl,c)));
      }
    }
  }
  
  // Check component consistency
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      assert (components(rl)>=0);
      for (int c=0; c<components(rl); ++c) {
        ibbox const & b  = extent(ml,rl,c);
        ibbox const & b0 = extent(ml,rl,0);
	assert (all(b.stride() == b0.stride()));
	assert (b.is_aligned_with(b0));
        for (int cc=c+1; cc<components(rl); ++cc) {
          assert ((b & extent(ml,rl,cc)).empty());
        }
      }
    }
  }
  
  // Check base grid extent
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
        assert (extent(ml,rl,c).is_contained_in(baseextents.at(ml).at(rl)));
      }
    }
  }

  // Check proper nesting, i.e., whether the fine grids are contained
  // in the coarser grids
  {
    bool have_error = false;
    for (int ml=0; ml<mglevels(); ++ml) {
      for (int rl=1; rl<reflevels(); ++rl) {
        if (components(rl)>0) {
          assert (components(rl-1)>0);
          assert (all (extent(ml,rl,0).stride() * reffacts.AT(rl) ==
                       extent(ml,rl-1,0).stride() * reffacts.AT(rl-1)));
          // Check contained-ness:
          // first take all coarse grids
          ibset coarse;
          for (int c=0; c<components(rl-1); ++c) {
            coarse += extent(ml,rl-1,c);
          }
          coarse.normalize();
          // then check all fine grids
          for (int c=0; c<components(rl); ++c) {
            ibbox const & fine =
              extent(ml,rl,c).contracted_for(extent(ml,rl-1,0));
            if (not (fine <= coarse)) {
              if (not have_error) {
                cout << "The following components are not properly nested, i.e.," << endl
                     << "they are not contained within the next coarser level's components:" << endl;
                have_error = true;
              }
              cout << "   ml " << ml << " rl " << rl << " c " << c << ":   "
                   << fine << endl;
            }
          } // for c
        }   // if c
      }     // for rl
    }       // for ml
    if (have_error) {
      cout << "The grid hierarchy is:" << endl;
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
  
  // Output
  if (output_bboxes) {
    do_output_bboxes (cout);
    do_output_bases (cout);
  }
  
  // Regrid the other hierarchies
  for (list<th*>::iterator t=ths.begin(); t!=ths.end(); ++t) {
    (*t)->regrid();
  }
  
  for (list<dh*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
    (*d)->regrid();
  }
}



bool
gh::
recompose (int const rl,
           bool const do_prolongate)
{
  // Handle changes in number of mglevels
  if (oldregions.size() != regions.size()) {
    oldregions.resize (regions.size());
  }
  
  bool const do_recompose = level_did_change(rl);
  
  if (do_recompose) {
    
    // Recompose the other hierarchies
    for (list<dh*>::iterator d=dhs.begin(); d!=dhs.end(); ++d) {
      (*d)->recompose (rl, do_prolongate);
    }
    
    // Overwrite old with new grid hierarchy
    for (int ml=0; ml<mglevels(); ++ml) {
      oldregions.AT(ml).resize (regions.AT(ml).size());
      oldregions.AT(ml).AT(rl) = regions.AT(ml).AT(rl);
    }
    
  }
  
  return do_recompose;
}



bool
gh::
level_did_change (int const rl)
  const
{
  // Find out whether this level changed
  if (regions.size() != oldregions.size()) return true;
  for (int ml=0; ml<mglevels(); ++ml) {
    assert (rl>=0 and rl<reflevels());
    if (rl >= (int)oldregions.AT(ml).size()) return true;
    if (regions.AT(ml).AT(rl).size() != oldregions.AT(ml).AT(rl).size()) {
      return true;
    }
    for (int c=0; c<components(rl); ++c) {
      if (regions.AT(ml).AT(rl).AT(c) != oldregions.AT(ml).AT(rl).AT(c)) {
        return true;
      }
    } // for c
  } // for ml
  return false;
}



// Accessors

int
gh::
local_components (int const rl)
  const
{
  int lc = 0;
  for (int c=0; c<components(rl); ++c) {
    lc +=  is_local(rl,c);
  }
  return lc;
}



// Find the refinement level and component to which a grid point
// belongs.  This uses a tree search over the superregions in the grid
// struction, which should scale reasonably (O(n log n)) better with
// the number of componets components.
void
gh::
locate_position (rvect const & rpos,
                 int const ml,
                 int const minrl, int const maxrl,
                 int & rl, int & c, ivect & aligned_ipos) const
{
  assert (ml>=0 and ml<mglevels());
  assert (minrl>=0 and minrl<=maxrl and maxrl<=reflevels());
  
  // Try finer levels first
  for (rl = maxrl-1; rl >= minrl; --rl) {
    
    // Align (round) the position to the nearest existing grid point
    // on this refinement level
    ivect const str = baseextent(ml,rl).stride();
    aligned_ipos = ivect(floor(rpos / rvect(str) + rvect(0.5))) * str;
    
    gh::cregs const & regs = superregions.AT(rl);
    
    // Search all superregions linearly.  Each superregion corresponds
    // to a "refined region", and the number of superregions is thus
    // presumably independent of the number of processors.
    for (size_t r = 0; r < regs.size(); ++r) {
      region_t const & reg = regs.AT(r);
      if (reg.extent.contains(aligned_ipos)) {
        // We found the superregion to which this grid point belongs
        
        // Search the superregion hierarchically
        pseudoregion_t const * const preg =
          reg.processors->search(aligned_ipos);
        assert (preg);
        
        // We now know the refinement level, component, and index to
        // which this grid point belongs
        c = preg->component;
        return;
      }
    }
  } // for rl
  
  // The point does not belong to any component on any refinement
  // level
  rl = -1;
  c = -1;
}

void
gh::
locate_position (ivect const & ipos,
                 int const ml,
                 int const minrl, int const maxrl,
                 int & rl, int & c, ivect & aligned_ipos) const
{
  assert (ml>=0 and ml<mglevels());
  assert (minrl>=0 and minrl<=maxrl and maxrl<=reflevels());
  
  // Try finer levels first
  for (rl = maxrl-1; rl >= minrl; --rl) {
    
    // Align (round) the position to the nearest existing grid point
    // on this refinement level
    ivect const str = baseextent(ml, rl).stride();
    aligned_ipos = ivect(floor(rvect(ipos) / rvect(str) + rvect(0.5))) * str;
    
    gh::cregs const & regs = superregions.AT(rl);
    
    // Search all superregions linearly.  Each superregion corresponds
    // to a "refined region", and the number of superregions is thus
    // presumably independent of the number of processors.
    for (size_t r = 0; r < regs.size(); ++r) {
      region_t const & reg = regs.AT(r);
      if (reg.extent.contains(aligned_ipos)) {
        // We found the superregion to which this grid point belongs
        
        // Search the superregion hierarchically
        pseudoregion_t const * const preg =
          reg.processors->search(aligned_ipos);
        assert (preg);
        
        // We now know the refinement level, component, and index to
        // which this grid point belongs
        c = preg->component;
        return;
      }
    }
  } // for rl
  
  // The point does not belong to any component on any refinement
  // level
  rl = -1;
  c = -1;
}



// Time hierarchy management

void
gh::
add (th * const t)
{
  ths.push_back (t);
}

void
gh::
remove (th * const t)
{
  ths.remove (t);
}



// Data hierarchy management

void
gh::
add (dh * const d)
{
  dhs.push_back (d);
}

void
gh::
remove (dh * const d)
{
  dhs.remove (d);
}



// Memory usage

size_t
gh::
memory ()
  const
{
  return
    memoryof (reffacts) +
    memoryof (refcent) +
    memoryof (mgfact) +
    memoryof (mgcent) +
    memoryof (baseextents) +
    memoryof (boundary_width) +
    memoryof (regions) +
    memoryof (oldregions) +
    memoryof (ths) +
    memoryof (dhs);
}



// Output

void
gh::
do_output_bboxes (ostream & os)
  const
{
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      for (int c=0; c<components(rl); ++c) {
        os << eol;
        os << "gh bboxes:" << eol;
        os << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
        os << "extent=" << extent(ml,rl,c) << eol;
        os << "outer_boundaries=" << outer_boundaries(rl,c) << eol;
        os << "processor=" << processor(rl,c) << eol;
      }
    }
  }
}

void
gh::
do_output_bases (ostream & os)
  const
{
  for (int ml=0; ml<mglevels(); ++ml) {
    for (int rl=0; rl<reflevels(); ++rl) {
      if (components(rl)>0) {
        os << eol;
        os << "gh bases:" << eol;
        os << "ml=" << ml << " rl=" << rl << eol;
        os << "base=" << baseextents.AT(ml).AT(rl) << eol;
      }
    }
  }
}

ostream &
gh::
output (ostream & os)
  const
{
  os << "gh:"
     << "reffacts=" << reffacts << ",refcentering=" << refcent << ","
     << "mgfactor=" << mgfact << ",mgcentering=" << mgcent << ","
     << "baseextents=" << baseextents << ","
     << "boundary_width=" << boundary_width << ","
     << "regions=" << regions << ","
     << "dhs={";
  {
    bool isfirst = true;
    for (list<dh*>::const_iterator
           d = dhs.begin(); d != dhs.end(); ++ d, isfirst = false)
    {
      if (not isfirst) os << ",";
      os << *d;
    }
  }
  os << "}";
  return os;
}
