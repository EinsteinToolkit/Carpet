// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dh.cc,v 1.39 2003/07/17 15:40:28 schnetter Exp $

#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "vect.hh"

#include "dh.hh"

using namespace std;



// Constructors
template<int D>
dh<D>::dh (gh<D>& h, const ivect& lghosts, const ivect& ughosts,
	   int prolongation_order_space)
  : h(h),
    lghosts(lghosts), ughosts(ughosts),
    prolongation_order_space(prolongation_order_space)
{
  assert (all(lghosts>=0 && ughosts>=0));
  h.add(this);
  CHECKPOINT;
  recompose(false);
}

// Destructors
template<int D>
dh<D>::~dh ()
{
  CHECKPOINT;
  h.remove(this);
}

// Helpers
template<int D>
int dh<D>::prolongation_stencil_size () const {
  assert (prolongation_order_space>=0);
  return prolongation_order_space/2;
}

// Modifiers
template<int D>
void dh<D>::recompose (const int initialise_upto) {
  DECLARE_CCTK_PARAMETERS;
  
  CHECKPOINT;
  
  boxes.clear();
  
  boxes.resize(h.reflevels());
  for (int rl=0; rl<h.reflevels(); ++rl) {
    boxes[rl].resize(h.components(rl));
    for (int c=0; c<h.components(rl); ++c) {
      boxes[rl][c].resize(h.mglevels(rl,c));
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
       	const ibbox intr = h.extents[rl][c][ml];
	
       	// Interior
       	// (the interior of the grid has the extent as specified by
       	// the user)
       	boxes[rl][c][ml].interior = intr;
	
       	// Exterior (add ghost zones)
       	// (the content of the exterior is completely determined by
       	// the interior of this or other components; the content of
       	// the exterior is redundant)
	ivect ldist(lghosts), udist(ughosts);
	for (int d=0; d<D; ++d) {
	  if (h.outer_boundaries[rl][c][d][0]) ldist[d] = 0;
	  if (h.outer_boundaries[rl][c][d][1]) udist[d] = 0;
	}
        boxes[rl][c][ml].exterior = intr.expand(ldist, udist);
	
       	// Boundaries (ghost zones only)
       	// (interior + boundaries = exterior)
       	boxes[rl][c][ml].boundaries = boxes[rl][c][ml].exterior - intr;
	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	
      	// Sync boxes
      	const int cs = h.components(rl);
      	boxes[rl][c][ml].send_sync.resize(cs);
      	boxes[rl][c][ml].recv_sync.resize(cs);
      	
      	// Refinement boxes
      	if (rl>0) {
      	  const int csm1 = h.components(rl-1);
      	  boxes[rl][c][ml].send_ref_coarse.resize(csm1);
      	  boxes[rl][c][ml].recv_ref_coarse.resize(csm1);
      	  boxes[rl][c][ml].recv_ref_bnd_coarse.resize(csm1);
      	}
      	if (rl<h.reflevels()-1) {
      	  const int csp1 = h.components(rl+1);
      	  boxes[rl][c][ml].recv_ref_fine.resize(csp1);
      	  boxes[rl][c][ml].send_ref_fine.resize(csp1);
      	  boxes[rl][c][ml].send_ref_bnd_fine.resize(csp1);
      	}
	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	const ibset& bnds = boxes[rl][c][ml].boundaries;
	
      	// Sync boxes
      	for (int cc=0; cc<h.components(rl); ++cc) {
      	  assert (ml<h.mglevels(rl,cc));
      	  // intersect boundaries with interior of that component
      	  const ibset ovlp = bnds & boxes[rl][cc][ml].interior;
      	  for (typename ibset::const_iterator b=ovlp.begin();
	       b!=ovlp.end(); ++b) {
	    boxes[rl][c ][ml].recv_sync[cc].push_back(*b);
	    boxes[rl][cc][ml].send_sync[c ].push_back(*b);
	  }
      	} // for cc
      	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	const ibbox& intr = boxes[rl][c][ml].interior;
      	const ibbox& extr = boxes[rl][c][ml].exterior;
	
      	// Multigrid boxes
      	if (ml>0) {
      	  const ibbox intrf = boxes[rl][c][ml-1].interior;
      	  const ibbox extrf = boxes[rl][c][ml-1].exterior;
      	  // Restriction (interior)
      	  {
      	    // (the restriction must fill all of the interior of the
      	    // coarse grid, and may use the exterior of the fine grid)
      	    const ibbox recv = intr;
            assert (! recv.empty());
      	    const ibbox send = recv.expanded_for(extrf);
            assert (! send.empty());
      	    assert (send.is_contained_in(extrf));
      	    boxes[rl][c][ml-1].send_mg_coarse.push_back(send);
      	    boxes[rl][c][ml  ].recv_mg_fine  .push_back(recv);
      	  }
      	  // Prolongation (interior)
      	  {
      	    // (the prolongation may use the exterior of the coarse
      	    // grid, and may fill only the interior of the fine grid,
      	    // and the bbox must be as large as possible)
      	    const ibbox recv = extr.contracted_for(intrf) & intrf;
            assert (! recv.empty());
      	    const ibbox send = recv.expanded_for(extr);
            assert (! send.empty());
      	    boxes[rl][c][ml-1].recv_mg_coarse.push_back(recv);
      	    boxes[rl][c][ml  ].send_mg_fine  .push_back(send);
      	  }
      	} // if not finest multigrid level
        
      } // for ml
    } // for c
  } // for rl
  
  // TODO: prefer boxes from the same processor
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
      	const ibbox& intr = boxes[rl][c][ml].interior;
      	const ibbox& extr = boxes[rl][c][ml].exterior;
	
      	// Refinement boxes
      	if (rl<h.reflevels()-1) {
      	  for (int cc=0; cc<h.components(rl+1); ++cc) {
      	    const ibbox intrf = boxes[rl+1][cc][ml].interior;
      	    // Prolongation (interior)
      	    {
      	      // (the prolongation may use the exterior of the coarse
      	      // grid, and must fill all of the interior of the fine
      	      // grid)
              const int pss = prolongation_stencil_size();
              ibset recvs = (extr.expand(-pss,-pss).contracted_for(intrf)
                             & intrf);
              const iblistvect& rrc = boxes[rl+1][cc][ml].recv_ref_coarse;
              for (typename iblistvect::const_iterator lvi=rrc.begin();
                   lvi!=rrc.end(); ++lvi) {
                for (typename iblist::const_iterator li=lvi->begin();
                     li!=lvi->end(); ++li) {
                  recvs -= *li;
                }
              }
              assert (recvs.setsize() <= 1);
              if (recvs.setsize() == 1) {
                const ibbox recv = *recvs.begin();
                const ibbox send = recv.expanded_for(extr);
                assert (! send.empty());
                assert (send.is_contained_in(extr));
                boxes[rl+1][cc][ml].recv_ref_coarse[c ].push_back(recv);
                boxes[rl  ][c ][ml].send_ref_fine  [cc].push_back(send);
              }
      	    }
      	    // Prolongation (boundaries)
      	    {
              const int pss = prolongation_stencil_size();
      	      const ibset& bndsf = boxes[rl+1][cc][ml].boundaries;
      	      // coarsify boundaries of fine component
      	      for (typename ibset::const_iterator bi=bndsf.begin();
		   bi!=bndsf.end(); ++bi) {
		const ibbox& bndf = *bi;
		// (the prolongation may use the exterior of the
		// coarse grid, and must fill all of the boundary of
		// the fine grid)
                ibset recvs = (extr.expand(-pss,-pss).contracted_for(bndf)
                               & bndf);
                const iblistvect& rrbc
                  = boxes[rl+1][cc][ml].recv_ref_bnd_coarse;
                for (typename iblistvect::const_iterator lvi=rrbc.begin();
                     lvi!=rrbc.end(); ++lvi) {
                  for (typename iblist::const_iterator li=lvi->begin();
                       li!=lvi->end(); ++li) {
                    recvs -= *li;
                  }
                }
                const iblistvect& rs = boxes[rl+1][cc][ml].recv_sync;
                for (typename iblistvect::const_iterator lvi=rs.begin();
                     lvi!=rs.end(); ++lvi) {
                  for (typename iblist::const_iterator li=lvi->begin();
                       li!=lvi->end(); ++li) {
                    recvs -= *li;
                  }
                }
                recvs.normalize();
                for (typename ibset::const_iterator si = recvs.begin();
                     si != recvs.end(); ++si) {
                  const ibbox & recv = *si;
                  const ibbox send = recv.expanded_for(extr);
                  assert (! send.empty());
                  assert (send.is_contained_in(extr));
                  boxes[rl+1][cc][ml].recv_ref_bnd_coarse[c ].push_back(recv);
                  boxes[rl  ][c ][ml].send_ref_bnd_fine  [cc].push_back(send);
                }
      	      }
      	    }
      	    // Restriction (interior)
      	    {
      	      // (the restriction may fill the interior of the of the
      	      // coarse grid, and may use the interior of the fine
      	      // grid, and the bbox must be as large as possible)
#if 0
	      // (the restriction must not fill points that are used
	      // to prolongate the boundaries)
              const int pss = prolongation_stencil_size();
              ibset recvs = intrf.contracted_for(intr) & intr;
              for (int ccc=0; ccc<h.components(rl); ++ccc) {
                const iblist& sendlist
                  = boxes[rl][ccc][ml].send_ref_bnd_fine[cc];
                for (typename iblist::const_iterator sli = sendlist.begin();
                     sli != sendlist.end(); ++sli) {
                  recvs -= (*sli).expand(pss+1,pss+1);
                }
              }
              recvs.normalize();
              for (typename ibset::const_iterator si = recvs.begin();
                   si != recvs.end(); ++si) {
                const ibbox & recv = *si;
                assert (! recv.empty());
                const ibbox send = recv.expanded_for(intrf);
                assert (! send.empty());
                boxes[rl+1][cc][ml].send_ref_coarse[c ].push_back(send);
                boxes[rl  ][c ][ml].recv_ref_fine  [cc].push_back(recv);
              }
#else
              ivect buf[2];
              for (int f=0; f<2; ++f) {
                for (int d=0; d<D; ++d) {
                  buf[f][d] = (h.outer_boundaries[rl+1][cc][d][f]
                               ? 0 : buffer_width);
                }
              }
              const ibbox recv = (intrf.contracted_for(intr)
                                  .expand(-buf[0], -buf[1])
                                  & intr);
              const ibbox send = recv.expanded_for(intrf);
              assert (send.empty() == recv.empty());
              if (! send.empty()) {
                boxes[rl+1][cc][ml].send_ref_coarse[c ].push_back(send);
                boxes[rl  ][c ][ml].recv_ref_fine  [cc].push_back(recv);
              }
#endif
      	    }
            
      	  } // for cc
      	} // if not finest refinement level
	
      } // for ml
    } // for c
  } // for rl
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	
	// Boundaries that are not synced, or are neither synced nor
	// prolonged to from coarser grids (outer boundaries)
        ibset& sync_not = boxes[rl][c][ml].sync_not;
        ibset& recv_not = boxes[rl][c][ml].recv_not;
        
        // The whole boundary
        sync_not = boxes[rl][c][ml].boundaries;
        recv_not = boxes[rl][c][ml].boundaries;
        
        // Subtract boxes received during synchronisation
        const iblistvect& recv_sync = boxes[rl][c][ml].recv_sync;
        for (typename iblistvect::const_iterator lvi=recv_sync.begin();
             lvi!=recv_sync.end(); ++lvi) {
          for (typename iblist::const_iterator li=lvi->begin();
               li!=lvi->end(); ++li) {
            sync_not -= *li;
            recv_not -= *li;
          }
        }
        
        // Subtract boxes received during prolongation
        const iblistvect& recv_ref_bnd_coarse
          = boxes[rl][c][ml].recv_ref_bnd_coarse;
        for (typename iblistvect::const_iterator
               lvi=recv_ref_bnd_coarse.begin();
             lvi!=recv_ref_bnd_coarse.end(); ++lvi) {
          for (typename iblist::const_iterator li=lvi->begin();
               li!=lvi->end(); ++li) {
            recv_not -= *li;
          }
        }
	
      } // for ml
    } // for c
  } // for rl
  
  // Calculate bases
  bases.resize(h.reflevels());
  for (int rl=0; rl<h.reflevels(); ++rl) {
    if (h.components(rl)==0) {
      bases[rl].resize(0);
    } else {
      bases[rl].resize(h.mglevels(rl,0));
      for (int ml=0; ml<h.mglevels(rl,0); ++ml) {
	bases[rl][ml].exterior = ibbox();
	bases[rl][ml].interior = ibbox();
	for (int c=0; c<h.components(rl); ++c) {
	  bases[rl][ml].exterior
	    = (bases[rl][ml].exterior
	       .expanded_containing(boxes[rl][c][ml].exterior));
	  bases[rl][ml].interior
	    = (bases[rl][ml].interior
	       .expanded_containing(boxes[rl][c][ml].interior));
	}
	bases[rl][ml].boundaries
	  = bases[rl][ml].exterior - bases[rl][ml].interior;
      }
    }
  }
  
  if (output_bboxes) {
    cout << endl << h << endl;
    for (int rl=0; rl<h.reflevels(); ++rl) {
      for (int c=0; c<h.components(rl); ++c) {
	for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	  cout << endl;
	  cout << "dh bboxes:" << endl;
	  cout << "rl=" << rl << " c=" << c << " ml=" << ml << endl;
	  cout << "exterior=" << boxes[rl][c][ml].exterior << endl;
	  cout << "interior=" << boxes[rl][c][ml].interior << endl;
	  cout << "send_mg_fine=" << boxes[rl][c][ml].send_mg_fine << endl;
	  cout << "send_mg_coarse=" << boxes[rl][c][ml].send_mg_coarse << endl;
	  cout << "recv_mg_fine=" << boxes[rl][c][ml].recv_mg_fine << endl;
	  cout << "recv_mg_coarse=" << boxes[rl][c][ml].recv_mg_coarse << endl;
	  cout << "send_ref_fine=" << boxes[rl][c][ml].send_ref_fine << endl;
	  cout << "send_ref_coarse=" << boxes[rl][c][ml].send_ref_coarse << endl;
	  cout << "recv_ref_fine=" << boxes[rl][c][ml].recv_ref_fine << endl;
	  cout << "recv_ref_coarse=" << boxes[rl][c][ml].recv_ref_coarse << endl;
	  cout << "send_sync=" << boxes[rl][c][ml].send_sync << endl;
	  cout << "send_ref_bnd_fine=" << boxes[rl][c][ml].send_ref_bnd_fine << endl;
	  cout << "boundaries=" << boxes[rl][c][ml].boundaries << endl;
	  cout << "recv_sync=" << boxes[rl][c][ml].recv_sync << endl;
	  cout << "recv_ref_bnd_coarse=" << boxes[rl][c][ml].recv_ref_bnd_coarse << endl;
	  cout << "sync_not=" << boxes[rl][c][ml].sync_not << endl;
	  cout << "recv_not=" << boxes[rl][c][ml].recv_not << endl;
	}
      }
    }
    for (int rl=0; rl<h.reflevels(); ++rl) {
      if (h.components(rl)>0) {
	for (int ml=0; ml<h.mglevels(rl,0); ++ml) {
	  cout << endl;
	  cout << "dh bases:" << endl;
	  cout << "rl=" << rl << " ml=" << ml << endl;
	  cout << "exterior=" << bases[rl][ml].exterior << endl;
	  cout << "interior=" << bases[rl][ml].interior << endl;
	  cout << "boundaries=" << bases[rl][ml].boundaries << endl;
	}
      }
    }
  } // if output_bboxes
  
  for (int rl=0; rl<h.reflevels(); ++rl) {
    for (int c=0; c<h.components(rl); ++c) {
      for (int ml=0; ml<h.mglevels(rl,c); ++ml) {
	
	// Assert that all boundaries are synced or received
        {
          const ibset& sync_not = boxes[rl][c][ml].sync_not;
#if 0
          const ibset& recv_not = boxes[rl][c][ml].recv_not;
#endif
          
          // Check that no boundaries are left over
          if (rl==0) assert (sync_not.empty());
#if 0
          assert (recv_not.empty());
#endif
        }
	
        // Assert that the interior is received exactly once during
        // prolongation, and that nothing else is received
        {
          if (rl==0) {
            const iblistvect& recv_ref_coarse
              = boxes[rl][c][ml].recv_ref_coarse;
            assert (recv_ref_coarse.empty());
          } else {              // rl!=0
            const iblistvect& recv_ref_coarse
              = boxes[rl][c][ml].recv_ref_coarse;
            ibset intr = boxes[rl][c][ml].interior;
            for (typename iblistvect::const_iterator
                   lvi=recv_ref_coarse.begin();
                 lvi!=recv_ref_coarse.end(); ++lvi) {
              for (typename iblist::const_iterator li=lvi->begin();
                   li!=lvi->end(); ++li) {
                const int old_sz = intr.size();
                const int this_sz = li->size();
                intr -= *li;
                const int new_sz = intr.size();
                assert (new_sz + this_sz == old_sz);
              }
            }
#if 0
            assert (intr.empty());
#endif
          }
        }
        
        // Assert that the boundaries are received at most once during
        // prolongation and synchronisation, and that nothing else is
        // received
        {
          const iblistvect& recv_sync = boxes[rl][c][ml].recv_sync;
          const iblistvect& recv_ref_bnd_coarse
            = boxes[rl][c][ml].recv_ref_bnd_coarse;
          ibset bnds = boxes[rl][c][ml].boundaries;
          for (typename iblistvect::const_iterator lvi=recv_sync.begin();
               lvi!=recv_sync.end(); ++lvi) {
            for (typename iblist::const_iterator li=lvi->begin();
                 li!=lvi->end(); ++li) {
              const int old_sz = bnds.size();
              const int this_sz = li->size();
              bnds -= *li;
              const int new_sz = bnds.size();
              assert (new_sz + this_sz == old_sz);
            }
          }
          for (typename iblistvect::const_iterator
                 lvi=recv_ref_bnd_coarse.begin();
               lvi!=recv_ref_bnd_coarse.end(); ++lvi) {
            for (typename iblist::const_iterator li=lvi->begin();
                 li!=lvi->end(); ++li) {
              const int old_sz = bnds.size();
              const int this_sz = li->size();
              bnds -= *li;
              const int new_sz = bnds.size();
              assert (new_sz + this_sz == old_sz);
            }
          }
        }
        
      } // for ml
    } // for c
  } // for rl
  
  for (typename list<ggf<D>*>::iterator f=gfs.begin();
       f!=gfs.end(); ++f) {
    (*f)->recompose(initialise_upto);
  }
}



// Grid function management
template<int D>
void dh<D>::add (ggf<D>* f) {
  CHECKPOINT;
  gfs.push_back(f);
}

template<int D>
void dh<D>::remove (ggf<D>* f) {
  CHECKPOINT;
  gfs.remove(f);
}



// Output
template<int D>
void dh<D>::output (ostream& os) const {
  os << "dh<" << D << ">:"
     << "ghosts=[" << lghosts << "," << ughosts << "],"
     << "gfs={";
  int cnt=0;
  for (typename list<ggf<D>*>::const_iterator f = gfs.begin();
       f != gfs.end(); ++f) {
    if (cnt++) os << ",";
    (*f)->output(os);
  }
  os << "}";
}



template class dh<3>;
