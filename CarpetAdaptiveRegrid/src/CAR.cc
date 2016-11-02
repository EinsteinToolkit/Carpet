#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"
#include "vect.hh"

#include "CAR.hh"
#include "carpet.hh"

//
// For the moment this is going to live as one large file with one
// large routine. However, it should be broken up into pieces later
// when possible.
//

namespace CarpetAdaptiveRegrid {

//
// The following "local_" variables contain copies of the list of
// boxes and the outer boundaries. However, these should be
// PROCESSOR INDEPENDENT; that is, these are never used in calls to
// SplitRegions. All box computation is done with respect to these
// local variables, so should always be independent of the number of
// processors. The results are split across processors and put into
// the "standard" variables before being passed back to Carpet.
//

static gh::mregs local_regsss;

//
// Keep track of the last iteration on which we were called. This
// means that we can return if we've been called on this timestep
// (because we want to ensure proper nesting we regrid all levels
// finer than the one on which we are called. Therefore if we're
// called on the same iteration, as long as we're called in
// coarse->fine order, we only need to regrid once.)
//
// FIXME: Note that this may well cause problems with checkpoint /
// restart; should this be moved to a grid scalar?
//

static int last_iteration = -1;

using namespace std;
using namespace Carpet;

//
// The following Fortran helper routines are defined in
// CAR_utils.F90.
//

extern "C" {
void CCTK_FCALL CCTK_FNAME(copy_mask)(const int &snx, const int &sny,
                                      const int &snz, const int *smask,
                                      const int sbbox[3][3], const int &dnx,
                                      const int &dny, const int &dnz,
                                      int *dmask, const int dbbox[3][3]);
void CCTK_FCALL CCTK_FNAME(check_box)(
    const int &nx, const int &ny, const int &nz, const int *mask, int *sum_x,
    int *sum_y, int *sum_z, int *sig_x, int *sig_y, int *sig_z,
    const int bbox[3][3], int newbbox1[3][3], int newbbox2[3][3],
    const int &min_width, const CCTK_REAL &min_fraction, int &didit);
}

//
// The following helper routines translate between real coordinates
// and integer (Carpet bbox) coordinates on a given refinement
// level. Basically copied from CarpetRegrid, "pos" <-> "position"
// and "int" <-> "integer". Actually defined at the bottom of this
// file.
//

ivect pos2int(const cGH *const cctkGH, const gh &hh, const rvect &rpos,
              const int rl);
rvect int2pos(const cGH *const cctkGH, const gh &hh, const ivect &ipos,
              const int rl);

//
// The real routine that does all the work
//

CCTK_INT CarpetAdaptiveRegrid_Regrid(CCTK_POINTER_TO_CONST const cctkGH_,
                                     CCTK_POINTER const regsss_,
                                     CCTK_INT force) {
  DECLARE_CCTK_PARAMETERS;

  const cGH *const cctkGH = (const cGH *)cctkGH_;

  //
  // The following variables are the important ones that Carpet uses
  // to set up the grids. They contain the bounding boxes (bb), the
  // outer boundary specifications (ob) and the processor numbers
  // (p). I believe the "s" stands for "set"; so the "obss" is a set
  // of sets of outer boundary specifications. However, each "s" is
  // implemented as a vector, not a set. For the ob and p the
  // indices on the vectors are the refinement levels and the maps;
  // for the bounding boxes the order is multigrid levels,
  // refinement levels, maps.
  //

  gh::mregs &regsss = *(gh::mregs *)regsss_;

  gh const &hh = *vhh.at(map);

  assert(is_singlemap_mode());

  if (local_bbsss.empty()) { // It's the first call
    //
    // We set up the "local_" variables using just the base grid
    // which we get from looking at the base extent; of course, the
    // outer boundary set up is trivial.
    //
    // Note that this will need improving (in fact the whole
    // regridding mechanism need minor changes) in order to do mesh
    // refinement with multiple maps.
    //
    const ibbox &baseext = vdd.at(map)->bases.at(mglevel).at(reflevel).exterior;
    vector<ibbox> tmp_bbs;
    tmp_bbs.push_back(baseext);
    vector<bbvect> tmp_obs;
    tmp_obs.push_back(bbvect(true));
    vector<bbvect> tmp_rbs;
    tmp_rbs.push_back(bbvect(false));
    vector<vector<ibbox> > tmp_bbss(1);
    vector<vector<bbvect> > tmp_obss(1);
    vector<vector<bbvect> > tmp_rbss(1);
    tmp_bbss.at(0) = tmp_bbs;
    tmp_obss.at(0) = tmp_obs;
    tmp_rbss.at(0) = tmp_rbs;
    MakeMultigridBoxes(cctkGH, tmp_bbss, tmp_obss, local_bbsss);
    local_obss = tmp_obss;
    local_rbss = tmp_rbss;
    last_iteration = cctkGH->cctk_iteration;
    //
    // Having set up the base grid we then set any finer grids
    // according to the coordinate parameters, as if we were the
    // standard CarpetRegrid.
    //
    int do_recompose =
        ManualCoordinateList(cctkGH, hh, bbsss, obss, rbss, pss, local_bbsss,
                             local_obss, local_rbss);

    if (verbose) {
      ostringstream buf;
      buf << "Done with manual coordinate list. Total list is:" << endl
          << local_bbsss;
      CCTK_INFO(buf.str().c_str());
    }

    return do_recompose;
  }

  // FIXME: We should check that the local reflevel "agrees"
  // with what is passed in.

  // In force mode (force == true) we do not check the
  // CarpetAdaptiveregrid parameters

  if (!force) {

    assert(regrid_every == -1 || regrid_every == 0 ||
           regrid_every % maxmglevelfact == 0);

    // Return if no regridding is desired
    if (regrid_every == -1)
      return 0;

    // Return if we want to regrid during initial data only, and this
    // is not the time for initial data
    if (regrid_every == 0 && cctkGH->cctk_iteration != 0)
      return 0;

    // Return if we want to regrid regularly, but not at this time
    if (regrid_every > 0 && cctkGH->cctk_iteration != 0 &&
        (cctkGH->cctk_iteration - 1) % regrid_every != 0) {
      return 0;
    }
  }

  //
  // Even if force is set, there are times that I'm not going to do
  // any regridding, so be told.
  //

  if (reflevel == maxreflevels - 1)
    return 0;

  // Return if we want to regrid regularly, but not at this time
  if (regrid_every > 0 && cctkGH->cctk_iteration != 0 &&
      (cctkGH->cctk_iteration - 1) % regrid_every != 0) {
    return 0;
  }

  // Return if we have already been called on this iteration
  if (cctkGH->cctk_iteration == last_iteration) {
    return 0;
  } else {
    last_iteration = cctkGH->cctk_iteration;
  }

  //
  // Set up a few local variables for later use. In particular we
  // need to keep track of the level and map on which we were
  // called.
  //

  int do_recompose;
  do_recompose = 1;

  int sum_handle = CCTK_ReductionArrayHandle("sum");

  int called_on_ml = mglevel;
  int called_on_rl = reflevel;
  int called_on_grouptype = mc_grouptype;
  int called_on_map = carpetGH.map;

  int finest_current_rl = local_bbsss.at(0).size();
  finest_current_rl = min(finest_current_rl, (maxreflevels - 1));

  //
  // Loop over all levels finer than this one.
  //

  leave_singlemap_mode(const_cast<cGH *>(cctkGH));
  leave_level_mode(const_cast<cGH *>(cctkGH));
  for (int rl = called_on_rl; rl < finest_current_rl; ++rl) {
    enter_level_mode(const_cast<cGH *>(cctkGH), rl);
    enter_singlemap_mode(const_cast<cGH *>(cctkGH), called_on_map,
                         called_on_grouptype);

    if (verbose) {
      ostringstream buf;
      buf << "Entering level " << rl << " (of " << finest_current_rl
          << "), map " << called_on_map;
      CCTK_INFO(buf.str().c_str());
    }

    vector<ibbox> bbs = local_bbsss.at(mglevel).at(reflevel);

    stack<ibbox> final;

    vector<vector<ibbox> > bbss = bbsss.at(0);
    vector<vector<ibbox> > local_bbss = local_bbsss.at(0);

    ivect const reffact = spacereffacts.at(rl + 1) / spacereffacts.at(rl);

    bool did_regrid = false;

    //
    // For later use get the physical domain specification in real
    // coordinates.
    //

    rvect physical_min, physical_max;
    rvect interior_min, interior_max;
    rvect exterior_min, exterior_max;
    rvect base_spacing;
    int ierr = GetDomainSpecification(
        dim, &physical_min[0], &physical_max[0], &interior_min[0],
        &interior_max[0], &exterior_min[0], &exterior_max[0], &base_spacing[0]);
    assert(!ierr);

    //
    // Loop over all boxes on this level.
    //

    for (vector<ibbox>::const_iterator bbi = bbs.begin(); bbi != bbs.end();
         ++bbi) {

      ivect low = bbi->lower();
      ivect upp = bbi->upper();

      // low and upp now define the starting bbox.

      ibbox bb(low, upp, bbs.at(0).stride());

      if (verbose) {
        ostringstream buf;
        buf << "Found the local size of the box: " << endl << bb;
        CCTK_INFO(buf.str().c_str());
      }

      vector<int> local_mask(prod(bb.shape() / bb.stride()), 0),
          mask(prod(bb.shape() / bb.stride()), 0);

      if (veryverbose) {
        ostringstream buf;
        buf << "Allocated mask size: " << bb.shape() / bb.stride()
            << " (points: " << prod(bb.shape() / bb.stride()) << ")";
        CCTK_INFO(buf.str().c_str());
      }

      //
      // Setup the mask.
      // That is, for any errors that we can locally see that exceed
      // the threshold, set the mask to 1.
      //

      const ibbox &baseext =
          vdd.at(map)->bases.at(mglevel).at(reflevel).exterior;
      ivect imin = (bb.lower() - baseext.lower()) / bb.stride(),
            imax = (bb.upper() - baseext.lower()) / bb.stride();

      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        const CCTK_REAL *error_var_ptr = static_cast<const CCTK_REAL *>(
            CCTK_VarDataPtr(cctkGH, 0, error_var));
        const CCTK_REAL *x_var_ptr = static_cast<const CCTK_REAL *>(
            CCTK_VarDataPtr(cctkGH, 0, "Grid::x"));
        const CCTK_REAL *y_var_ptr = static_cast<const CCTK_REAL *>(
            CCTK_VarDataPtr(cctkGH, 0, "Grid::y"));
        const CCTK_REAL *z_var_ptr = static_cast<const CCTK_REAL *>(
            CCTK_VarDataPtr(cctkGH, 0, "Grid::z"));

        // These can actually be negative if the parent shrinks.
        // Of course, the final grid should still be properly
        // nested...
        //          assert(all(imin >= 0));
        //          assert(all(imax >= 0));
        // FIXME: Why should the following assert be true?
        //      assert(all(imax < ivect::ref(cctkGH->cctk_lsh)));
        assert(all(imin <= imax));

        for (int k = 0; k < cctkGH->cctk_lsh[2]; ++k) {
          for (int j = 0; j < cctkGH->cctk_lsh[1]; ++j) {
            for (int i = 0; i < cctkGH->cctk_lsh[0]; ++i) {
              int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL local_error = abs(error_var_ptr[index]);
              if (local_error > max_error) {
                int ii = i + cctkGH->cctk_lbnd[0] - imin[0];
                int jj = j + cctkGH->cctk_lbnd[1] - imin[1];
                int kk = k + cctkGH->cctk_lbnd[2] - imin[2];
                // Check that this point actually intersects with
                // this box (if this component was actually a
                // different grid on the same processor, it need not)
                if ((ii >= 0) and (jj >= 0) and (kk >= 0) and
                    (ii <= imax[0] - imin[0]) and (jj <= imax[1] - imin[1]) and
                    (kk <= imax[2] - imin[2])) {
                  assert(ii >= 0);
                  assert(jj >= 0);
                  assert(kk >= 0);
                  assert(ii <= imax[0] - imin[0]);
                  assert(jj <= imax[1] - imin[1]);
                  assert(kk <= imax[2] - imin[2]);
                  int mindex = ii +
                               (imax[0] - imin[0] + 1) *
                                   (jj + (imax[1] - imin[1] + 1) * kk);
                  local_mask[mindex] = 1;
                  if (veryverbose) {
                    CCTK_VInfo(CCTK_THORNSTRING,
                               "In error at point"
                               "\n(%g,%g,%g) [%d,%d,%d] [[%d,%d,%d]]",
                               (double)x_var_ptr[index],
                               (double)y_var_ptr[index],
                               (double)z_var_ptr[index], ii, jj, kk, i, j, k);
                  }
                }
              }
            }
          }
        }
      }
      END_LOCAL_COMPONENT_LOOP;

      //
      // Check the error on child level, if such a level exists
      // Also only worry if there's a grandchild level.
      // This should fix the "orphaned grandchild" problem
      //

      if ((int)local_bbss.size() > reflevel + 2) {

        int currentml = mglevel;
        int currentrl = reflevel;
        int currentmap = carpetGH.map;

        leave_singlemap_mode(const_cast<cGH *>(cctkGH));
        leave_level_mode(const_cast<cGH *>(cctkGH));

        enter_level_mode(const_cast<cGH *>(cctkGH), currentrl + 1);
        enter_singlemap_mode(const_cast<cGH *>(cctkGH), currentmap, CCTK_GF);

        const ibbox &child_baseext =
            vdd.at(map)->bases.at(mglevel).at(reflevel).exterior;
        ivect child_levoff = child_baseext.lower() / (bb.stride() / reffact);

        if (verbose) {
          ostringstream buf;
          buf << "Checking for errors on child level " << reflevel << " map "
              << currentmap;
          CCTK_INFO(buf.str().c_str());
        }

        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          const CCTK_REAL *error_var_ptr = static_cast<const CCTK_REAL *>(
              CCTK_VarDataPtr(cctkGH, 0, error_var));
          const CCTK_REAL *x_var_ptr = static_cast<const CCTK_REAL *>(
              CCTK_VarDataPtr(cctkGH, 0, "Grid::x"));
          const CCTK_REAL *y_var_ptr = static_cast<const CCTK_REAL *>(
              CCTK_VarDataPtr(cctkGH, 0, "Grid::y"));
          const CCTK_REAL *z_var_ptr = static_cast<const CCTK_REAL *>(
              CCTK_VarDataPtr(cctkGH, 0, "Grid::z"));

          //            assert(all(imin >= 0));
          //            assert(all(imax >= 0));
          // FIXME: Why should the following assert be true?
          //      assert(all(imax < ivect::ref(cctkGH->cctk_lsh)));
          assert(all(imin <= imax));

          if (verbose) {
            ostringstream buf;
            buf << "Checking for errors on child level " << reflevel << " map "
                << currentmap << " component " << component << endl
                << "Mask: " << imin << imax << endl
                << "Component: " << ivect::ref(cctkGH->cctk_lbnd)
                << ivect::ref(cctkGH->cctk_lsh) << endl
                << "level offset: " << child_levoff;
            CCTK_INFO(buf.str().c_str());
          }

          for (int k = 0; k < cctkGH->cctk_lsh[2]; ++k) {
            for (int j = 0; j < cctkGH->cctk_lsh[1]; ++j) {
              for (int i = 0; i < cctkGH->cctk_lsh[0]; ++i) {
                int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
                CCTK_REAL local_error = abs(error_var_ptr[index]);
                if (local_error > max_error) {
                  // Correct for the change in level !
                  int ii = (i + cctkGH->cctk_lbnd[0] -
                            cctkGH->cctk_nghostzones[0] + child_levoff[0]) /
                               reffact[0] -
                           imin[0];
                  int jj = (j + cctkGH->cctk_lbnd[1] -
                            cctkGH->cctk_nghostzones[1] + child_levoff[1]) /
                               reffact[1] -
                           imin[1];
                  int kk = (k + cctkGH->cctk_lbnd[2] -
                            cctkGH->cctk_nghostzones[2] + child_levoff[1]) /
                               reffact[2] -
                           imin[2];
                  // Check that this point actually intersects with
                  // this box (if this component was actually a
                  // different grid on the same processor, it need not)
                  if ((ii >= 0) and (jj >= 0) and (kk >= 0) and
                      (ii <= imax[0] - imin[0]) and
                      (jj <= imax[1] - imin[1]) and (kk <= imax[2] - imin[2])) {
                    assert(ii >= 0);
                    assert(jj >= 0);
                    assert(kk >= 0);
                    assert(ii <= imax[0] - imin[0]);
                    assert(jj <= imax[1] - imin[1]);
                    assert(kk <= imax[2] - imin[2]);
                    int mindex = ii +
                                 (imax[0] - imin[0] + 1) *
                                     (jj + (imax[1] - imin[1] + 1) * kk);
                    local_mask[mindex] = 1;
                    if (veryverbose) {
                      CCTK_VInfo(CCTK_THORNSTRING,
                                 "In error at point"
                                 "\n(%g,%g,%g) [%d,%d,%d] [[%d,%d,%d]]",
                                 (double)x_var_ptr[index],
                                 (double)y_var_ptr[index],
                                 (double)z_var_ptr[index], ii, jj, kk, i, j, k);
                    }
                  }
                }
              }
            }
          }
        }
        END_LOCAL_COMPONENT_LOOP;

        leave_singlemap_mode(const_cast<cGH *>(cctkGH));
        leave_level_mode(const_cast<cGH *>(cctkGH));

        enter_level_mode(const_cast<cGH *>(cctkGH), currentrl);
        enter_singlemap_mode(const_cast<cGH *>(cctkGH), currentmap, CCTK_GF);
      }

      //
      // Reduce errors (MPI sum).
      // This sets up the mask globally, although the mask is now
      // non-zero where in error instead of just being 1. It's still
      // 0 at points not in error.
      //

      int ierr = CCTK_ReduceLocArrayToArray1D(
          cctkGH, -1, sum_handle, &local_mask.front(), &mask.front(),
          local_mask.size(), CCTK_VARIABLE_INT);

      if (ierr) {
        ostringstream buf;
        buf << "The reduction on level " << reflevel << " failed "
            << "with error " << ierr;
        CCTK_WARN(0, buf.str().c_str());
      }

      //
      // Set errors to 0/1; so, where the mask indicates that a
      // point is in error, set it precisely to 1.
      //

      for (int k = 0; k < imax[2] - imin[2] + 1; k++) {
        for (int j = 0; j < imax[1] - imin[1] + 1; j++) {
          for (int i = 0; i < imax[0] - imin[0] + 1; i++) {
            int index =
                i + (imax[0] - imin[0] + 1) * (j + (imax[1] - imin[1] + 1) * k);
            if (mask[index]) {
              mask[index] = 1;
            }
          }
        }
      }

      //
      // Pad the errors: stage 1 - buffer points marked as 2.
      //

      for (int k = 0; k < imax[2] - imin[2] + 1; k++) {
        for (int j = 0; j < imax[1] - imin[1] + 1; j++) {
          for (int i = 0; i < imax[0] - imin[0] + 1; i++) {
            int index =
                i + (imax[0] - imin[0] + 1) * (j + (imax[1] - imin[1] + 1) * k);
            if (mask[index] == 1) {
              for (int kk = max((int)(k - pad), 0);
                   kk < min((int)(k + pad + 1), imax[2] - imin[2] + 1); ++kk) {
                for (int jj = max((int)(j - pad), 0);
                     jj < min((int)(j + pad + 1), imax[1] - imin[1] + 1);
                     ++jj) {
                  for (int ii = max((int)(i - pad), 0);
                       ii < min((int)(i + pad + 1), imax[0] - imin[0] + 1);
                       ++ii) {
                    int mindex = ii +
                                 (imax[0] - imin[0] + 1) *
                                     (jj + (imax[1] - imin[1] + 1) * kk);
                    if (!mask[mindex])
                      mask[mindex] = 2;
                  }
                }
              }
            }
          }
        }
      }
      //
      // stage 2: all buffer points marked truly in error.
      // Also mark if there are any errors.
      //
      bool should_regrid = false;
      for (int k = 0; k < imax[2] - imin[2] + 1; k++) {
        if (maskpicture) {
          cout << endl << "k = " << k << endl;
        }
        for (int j = 0; j < imax[1] - imin[1] + 1; j++) {
          if (maskpicture) {
            cout << endl;
          }
          for (int i = 0; i < imax[0] - imin[0] + 1; i++) {
            int index =
                i + (imax[0] - imin[0] + 1) * (j + (imax[1] - imin[1] + 1) * k);
            if (mask[index] > 1)
              mask[index] = 1;
            if ((veryverbose) and (mask[index])) {
              CCTK_VInfo(CCTK_THORNSTRING, "Mask set at point"
                                           "\n[%d,%d,%d]",
                         i, j, k);
            }
            if (maskpicture) {
              if (mask[index]) {
                cout << "x";
              } else {
                cout << " ";
              }
            }
            should_regrid |= (mask[index]);
            did_regrid |= should_regrid;
          }
        }
      }

      if (verbose) {
        ostringstream buf;
        buf << "Finished looking for errors on level " << reflevel << endl
            << "should_regrid " << should_regrid << " did_regrid "
            << did_regrid;
        CCTK_INFO(buf.str().c_str());
      }

      //
      // For this box on this level we now have the marked
      // points. If there are any errors then we should actually
      // create the new boxes. Starting from the actual bbox, placed
      // in the "todo" stack, we iterate over the "todo" stack
      // creating new boxes which are either accepted (and hence
      // placed on the "final" stack) or still need work done to
      // them (hence placed on the "todo" stack).
      //

      if (should_regrid) {

        stack<ibbox> todo;

        todo.push(bb);

        // Define vector of masks (contains error mask)

        stack<vector<int> > masklist;

        masklist.push(mask);

        //
        // Loop over all entries in todo:
        //   Setup appropriate 1d array memory
        //   Call fortran routine
        //   If return is:
        //     zero: add bbox to final
        //     one:  add new bbox to todo and assoc mask to masklist
        //     two:  add both new bboxs to todo and assoc masks to
        //           masklist
        //

        while (!todo.empty()) {

          ibbox bb = todo.top();
          todo.pop();
          vector<int> mask = masklist.top();
          masklist.pop();

          int nx = bb.shape()[0] / bb.stride()[0];
          int ny = bb.shape()[1] / bb.stride()[1];
          int nz = bb.shape()[2] / bb.stride()[2];

          if (verbose) {
            ostringstream buf;
            buf << "todo loop. Box: " << endl << bb;
            CCTK_INFO(buf.str().c_str());
          }

          vector<int> sum_x(nx, 0);
          vector<int> sig_x(nx, 0);
          vector<int> sum_y(ny, 0);
          vector<int> sig_y(ny, 0);
          vector<int> sum_z(nz, 0);
          vector<int> sig_z(nz, 0);

          int fbbox[3][3], fbbox1[3][3], fbbox2[3][3];

          for (int d = 0; d < 3; ++d) {
            fbbox[0][d] = bb.lower()[d];
            fbbox[1][d] = bb.upper()[d];
            fbbox[2][d] = bb.stride()[d];
          }

          int didit;

          //
          // This is the actual Fortran routine that does the work
          // of the Berger-Rigoutsos algorithm.
          //

          CCTK_FNAME(check_box)
          (nx, ny, nz, &mask.front(), &sum_x.front(), &sum_y.front(),
           &sum_z.front(), &sig_x.front(), &sig_y.front(), &sig_z.front(),
           fbbox, fbbox1, fbbox2, min_width, min_fraction, didit);

          if (didit == 0) { // Box was accepted

            final.push(bb);

            if (verbose) {
              ostringstream buf;
              buf << "todo loop. Box pushed to final: " << endl << bb;
              CCTK_INFO(buf.str().c_str());
            }
          } else if (didit == 1) { // Box was replaced by a new single
                                   // box
            ibbox newbbox1(ivect::ref(&fbbox1[0][0]), ivect::ref(&fbbox1[1][0]),
                           ivect::ref(&fbbox1[2][0]));
            todo.push(newbbox1);

            int dnx = newbbox1.shape()[0] / newbbox1.stride()[0];
            int dny = newbbox1.shape()[1] / newbbox1.stride()[1];
            int dnz = newbbox1.shape()[2] / newbbox1.stride()[2];

            vector<int> newmask1(prod(newbbox1.shape() / newbbox1.stride()), 0);

            CCTK_FNAME(copy_mask)
            (nx, ny, nz, &mask.front(), fbbox, dnx, dny, dnz, &newmask1.front(),
             fbbox1);
            masklist.push(newmask1);

            if (verbose) {
              ostringstream buf;
              buf << "todo loop. New (single) box created: " << endl
                  << newbbox1;
              CCTK_INFO(buf.str().c_str());
            }
          } else if (didit == 2) { // Box was replaced with two boxes

            ibbox newbbox1(ivect::ref(&fbbox1[0][0]), ivect::ref(&fbbox1[1][0]),
                           ivect::ref(&fbbox1[2][0]));
            todo.push(newbbox1);
            ibbox newbbox2(ivect::ref(&fbbox2[0][0]), ivect::ref(&fbbox2[1][0]),
                           ivect::ref(&fbbox2[2][0]));
            todo.push(newbbox2);

            int dnx = newbbox1.shape()[0] / newbbox1.stride()[0];
            int dny = newbbox1.shape()[1] / newbbox1.stride()[1];
            int dnz = newbbox1.shape()[2] / newbbox1.stride()[2];

            vector<int> newmask1(prod(newbbox1.shape() / newbbox1.stride()), 0);

            CCTK_FNAME(copy_mask)
            (nx, ny, nz, &mask.front(), fbbox, dnx, dny, dnz, &newmask1.front(),
             fbbox1);
            masklist.push(newmask1);

            dnx = newbbox2.shape()[0] / newbbox2.stride()[0];
            dny = newbbox2.shape()[1] / newbbox2.stride()[1];
            dnz = newbbox2.shape()[2] / newbbox2.stride()[2];

            vector<int> newmask2(prod(newbbox2.shape() / newbbox2.stride()), 0);

            CCTK_FNAME(copy_mask)
            (nx, ny, nz, &mask.front(), fbbox, dnx, dny, dnz, &newmask2.front(),
             fbbox2);
            masklist.push(newmask2);

            if (verbose) {
              ostringstream buf;
              buf << "todo loop. New (double) box created. Box 1: " << endl
                  << newbbox1
                  << "                                     Box 2: " << endl
                  << newbbox2;
              CCTK_INFO(buf.str().c_str());
            }
          } else {
            CCTK_WARN(0, "The fortran routine must be confused.");
          }

        } // loop over todo vector (boxes needing to be done).
      }   // should regrid.
    }     // Loop over boxes on the parent grid.

    if (did_regrid) { // If we actually did something, reconvert the
                      // boxes to correct Carpet style, plus correct
                      // the boundaries.
      // Fixup the stride
      vector<ibbox> newbbs;
      vector<bbvect> obs;
      vector<bbvect> rbs;
      while (!final.empty()) {
        ibbox bb = final.top();
        final.pop();

        if (veryverbose) {
          ostringstream buf;
          buf << "Looping over the final list. Box is:" << endl << bb;
          CCTK_INFO(buf.str().c_str());
        }

        ivect ilo = bb.lower();
        ivect ihi = bb.upper();
        rvect lo = int2pos(cctkGH, hh, ilo, reflevel);
        rvect hi = int2pos(cctkGH, hh, ihi, reflevel);
        rvect str = base_spacing * ipow((CCTK_REAL)mgfact, basemglevel) /
                    rvect(spacereffacts.at(reflevel));
        rbbox newbbcoord(lo, hi, str);

        if (veryverbose) {
          ostringstream buf;
          buf << "Dealing with boundaries. Coord box is:" << endl << newbbcoord;
          CCTK_INFO(buf.str().c_str());
        }

        // Set the correct ob here.

        bbvect ob(false);
        bbvect rb(true);
        for (int d = 0; d < dim; ++d) {
          assert(mglevel == 0);

          // Find the size of the physical domain

          rvect const spacing = base_spacing *
                                ipow((CCTK_REAL)mgfact, basemglevel) /
                                rvect(ipow(reffact, reflevel + 1));
          ierr = ConvertFromPhysicalBoundary(dim, &physical_min[0],
                                             &physical_max[0], &interior_min[0],
                                             &interior_max[0], &exterior_min[0],
                                             &exterior_max[0], &spacing[0]);
          assert(!ierr);

          // If need be clip the domain

          rvect lo = newbbcoord.lower();
          if (newbbcoord.lower()[d] < physical_min[d]) {
            lo[d] = exterior_min[d];
          }
          rvect up = newbbcoord.upper();
          if (newbbcoord.upper()[d] > physical_max[d]) {
            up[d] = exterior_max[d];
          }
          rvect str = newbbcoord.stride();

          // Set the ob if outside the physical domain

          ob[d][0] = abs(lo[d] - exterior_min[d]) < 1.0e-6 * spacing[d];
          ob[d][1] = abs(up[d] - exterior_max[d]) < 1.0e-6 * spacing[d];
          rb[d][0] = !ob[d][0];
          rb[d][1] = !ob[d][1];

          if (veryverbose) {
            ostringstream buf;
            buf << "Done clipping domain:" << endl
                << lo << endl
                << up << endl
                << str;
            CCTK_INFO(buf.str().c_str());
          }

          // Check that the striding is correct.

          CCTK_REAL remainder = fmod((up[d] - lo[d]), str[d]) / str[d];

          if (abs(remainder) > 1.e-6) {
            if (ob[d][0]) {
              up[d] += str[d] * (1 - remainder);
            } else if (ob[d][1]) {
              lo[d] -= str[d] * remainder;
            }
          }

          if (veryverbose) {
            ostringstream buf;
            buf << "Corrected coords for striding:" << endl
                << lo << endl
                << up << endl
                << str;
            CCTK_INFO(buf.str().c_str());
          }

          newbbcoord = rbbox(lo, up, str);
        }
        if (verbose) {
          ostringstream buf;
          buf << "Done dealing with boundaries. Coord box is:" << endl
              << newbbcoord << endl
              << "obox is:" << endl
              << ob;
          CCTK_INFO(buf.str().c_str());
        }

        //
        // Convert back to integer coordinates
        // We have to do this on the fine grid to ensure that
        // it is correct for an outer boundary with odd numbers
        // of ghost zones where the bbox does not align with the
        // parent.
        //

        ilo = pos2int(cctkGH, hh, newbbcoord.lower(), reflevel + 1);
        ihi = pos2int(cctkGH, hh, newbbcoord.upper(), reflevel + 1);
        ivect istr = bb.stride() / reffact;

        // Check that the width is sufficient
        // This can only be too small if the domain was clipped
        for (int d = 0; d < dim; ++d) {
          if (ihi[d] - ilo[d] < min_width * istr[d]) {
            if (ob[d][0]) {
              if (ob[d][1]) {
                CCTK_WARN(0, "The domain is too small?!");
              }
              ihi[d] = ilo[d] + min_width * istr[d];
            } else if (ob[d][1]) {
              if (ob[d][0]) {
                CCTK_WARN(0, "The domain is too small?!");
              }
              ilo[d] = ihi[d] - min_width * istr[d];
            } else {
              ostringstream buf;
              buf << "The grid is unclipped and too small?" << endl
                  << ilo << endl
                  << ihi << endl
                  << istr << endl
                  << d;
              CCTK_WARN(0, buf.str().c_str());
            }
          }
        }

        if (veryverbose) {
          ostringstream buf;
          buf << "Corrected integer coords for min_width:" << endl
              << ilo << endl
              << ihi << endl
              << istr;
          CCTK_INFO(buf.str().c_str());
        }

        ibbox newbb(ilo, ihi, istr);

        if (verbose) {
          ostringstream buf;
          buf << "After dealing with boundaries. Final box is:" << endl
              << newbb;
          CCTK_INFO(buf.str().c_str());
        }

        newbbs.push_back(newbb);
        obs.push_back(ob);
        rbs.push_back(rb);
      }

      // FIXME: check if the newbbs is really different
      // from the current bbs
      //        if not, set do_recompose = 0
      bbs = newbbs;

      // Set local bbss

      //
      // This mess ensures that the local bbsss is correctly set for
      // multigrid levels, then splits over regions (remember that
      // the local_bbsss is processor independent), then goes
      // through the whole multigrid procedure with the processor
      // dependent bbsss.
      //

      if ((int)bbss.size() < reflevel + 2) {
        if (verbose) {
          CCTK_INFO("Adding new refinement level");
        }
        local_bbss.resize(reflevel + 2);
        bbss.resize(reflevel + 2);
        local_obss.resize(reflevel + 2);
        obss.resize(reflevel + 2);
        local_rbss.resize(reflevel + 2);
        rbss.resize(reflevel + 2);
        pss.resize(reflevel + 2);
      }
      local_bbss.at(reflevel + 1) = bbs;
      local_obss.at(reflevel + 1) = obs;
      local_rbss.at(reflevel + 1) = rbs;
      MakeMultigridBoxes(cctkGH, local_bbss, local_obss, local_bbsss);

      // make multiprocessor aware
      gh::cprocs ps;
      SplitRegions(cctkGH, bbs, obs, rbs, ps);

      bbss.at(reflevel + 1) = bbs;
      obss.at(reflevel + 1) = obs;
      rbss.at(reflevel + 1) = rbs;
      pss.at(reflevel + 1) = ps;

    } // did_regrid?
    else {
      if ((int)local_bbss.size() > reflevel + 1) {
        if (verbose) {
          CCTK_INFO("Removing refinement level");
        }
      }
      local_bbss.resize(reflevel + 1);
      bbss.resize(reflevel + 1);
      local_obss.resize(reflevel + 1);
      obss.resize(reflevel + 1);
      local_rbss.resize(reflevel + 1);
      rbss.resize(reflevel + 1);
      // Set local bbsss
      MakeMultigridBoxes(cctkGH, local_bbss, local_obss, local_bbsss);

      pss.resize(reflevel + 1);

      do_recompose = 1;
    }

    // make multigrid aware
    MakeMultigridBoxes(cctkGH, bbss, obss, bbsss);

    leave_singlemap_mode(const_cast<cGH *>(cctkGH));
    leave_level_mode(const_cast<cGH *>(cctkGH));
  }

  enter_level_mode(const_cast<cGH *>(cctkGH), called_on_rl);
  enter_singlemap_mode(const_cast<cGH *>(cctkGH), called_on_map,
                       called_on_grouptype);

  if (verbose) {
    ostringstream buf;
    buf << "Done with it all. Iteration " << cctkGH->cctk_iteration << " level "
        << reflevel << endl
        << "Total list is:" << endl
        << local_bbsss;
    CCTK_INFO(buf.str().c_str());
  }

  return do_recompose;
}

ivect pos2int(const cGH *const cctkGH, const gh &hh, const rvect &rpos,
              const int rl) {
  rvect global_lower, global_upper;
  for (int d = 0; d < dim; ++d) {
    const int ierr = CCTK_CoordRange(cctkGH, &global_lower[d], &global_upper[d],
                                     d + 1, 0, "cart3d");
    if (ierr < 0) {
      global_lower[d] = 0;
      global_upper[d] = 1;
    }
  }
  const ivect global_extent(hh.baseextent.upper() - hh.baseextent.lower());

  const rvect scale = rvect(global_extent) / (global_upper - global_lower);
  const ivect levfac = hh.reffacts.at(rl);
  assert(all(hh.baseextent.stride() % levfac == 0));
  const ivect istride = hh.baseextent.stride() / levfac;

  const ivect ipos =
      (ivect(floor((rpos - global_lower) * scale / rvect(istride) +
                   (CCTK_REAL)0.5)) *
       istride);

  return ipos;
}

rvect int2pos(const cGH *const cctkGH, const gh &hh, const ivect &ipos,
              const int rl) {
  rvect global_lower, global_upper;
  for (int d = 0; d < dim; ++d) {
    const int ierr = CCTK_CoordRange(cctkGH, &global_lower[d], &global_upper[d],
                                     d + 1, 0, "cart3d");

    if (ierr < 0) {
      global_lower[d] = 0;
      global_upper[d] = 1;
    }
  }
  const ivect global_extent(hh.baseextent.upper() - hh.baseextent.lower());

  const rvect scale = rvect(global_extent) / (global_upper - global_lower);
  const ivect levfac = hh.reffacts.at(rl);
  assert(all(hh.baseextent.stride() % levfac == 0));
  const ivect istride = hh.baseextent.stride() / levfac;

  const rvect rpos = rvect(ipos) / scale + global_lower;

  return rpos;
}

} // namespace CarpetAdaptiveRegrid
