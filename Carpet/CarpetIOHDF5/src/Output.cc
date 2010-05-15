#include <cassert>
#include <cstdlib>
#include <cstring>
#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "operators.hh"
#include "CarpetIOHDF5.hh"
#include "CactusBase/IOUtil/src/ioGH.h"


namespace CarpetIOHDF5
{

using namespace std;
using namespace Carpet;


// add attributes to an HDF5 dataset
static int AddAttributes (const cGH *const cctkGH, const char *fullname,
                          int vdim, int refinementlevel,
                          const ioRequest* const request,
                          const ibbox& bbox, hid_t dataset, bool is_index = false);


int WriteVarUnchunked (const cGH* const cctkGH,
                       hid_t outfile,
                       CCTK_REAL & io_bytes,
                       const ioRequest* const request,
                       bool called_from_checkpoint)
{
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;
  char *fullname = CCTK_FullName(request->vindex);
  const int gindex = CCTK_GroupIndexFromVarI (request->vindex);
  assert (gindex >= 0 and gindex < (int) Carpet::arrdata.size ());
  const int var = request->vindex - CCTK_FirstVarIndexI (gindex);
  assert (var >= 0 and var < CCTK_NumVars ());
  cGroup group;
  CCTK_GroupData (gindex, &group);


  // Scalars and arrays have only one refinement level 0,
  // regardless of what the current refinement level is.
  // Output for them must be called in global mode.
  int refinementlevel = reflevel;
  if (group.grouptype == CCTK_SCALAR or group.grouptype == CCTK_ARRAY) {
    assert (do_global_mode);
    refinementlevel = 0;
  }

  // HDF5 doesn't like 0-dimensional arrays
  if (group.grouptype == CCTK_SCALAR) group.dim = 1;

  // If the user requested so, select single precision data output
  hid_t memdatatype, filedatatype;
  HDF5_ERROR (memdatatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype, 0));
  HDF5_ERROR (filedatatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype,
                          out_single_precision and not called_from_checkpoint));

  // create a file access property list to use the CORE virtual file driver
  hid_t plist;
  HDF5_ERROR (plist = H5Pcreate (H5P_FILE_ACCESS));
  HDF5_ERROR (H5Pset_fapl_core (plist, 0, 0));

  // Traverse all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype) {
    // Collect the set of all components' bboxes
    ibset bboxes;
    BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype) {
      // Using "interior" removes ghost zones and refinement boundaries.
      bboxes += arrdata.at(gindex).at(Carpet::map).dd->
                boxes.at(mglevel).at(refinementlevel).at(component).interior;
    } END_COMPONENT_LOOP;

    // Normalise the set, i.e., try to represent the set with fewer bboxes.
    //
    // According to Cactus conventions, DISTRIB=CONSTANT arrays
    // (including grid scalars) are assumed to be the same on all
    // processors and are therefore stored only by processor 0.
    if (group.disttype != CCTK_DISTRIB_CONSTANT) bboxes.normalize();

    // Loop over all components in the bbox set
    int bbox_id = 0;
    for (ibset::const_iterator bbox  = bboxes.begin();
                               bbox != bboxes.end();
                               bbox++, bbox_id++) {
      // Get the shape of the HDF5 dataset (in Fortran index order)
      hsize_t shape[dim];
      hsize_t num_elems = 1;
      for (int d = 0; d < group.dim; ++d) {
        assert (group.dim-1-d>=0 and group.dim-1-d<dim);
        shape[group.dim-1-d] = (bbox->shape() / bbox->stride())[d];
        num_elems *= shape[group.dim-1-d];
      }

      // Don't create zero-sized components
      if (num_elems == 0) continue;

      // create the dataset on the I/O processor
      // skip DISTRIB=CONSTANT components from processors other than 0
      hid_t memfile = -1, memdataset = -1;
      hid_t dataspace = -1, dataset = -1;
      if (dist::rank() == 0 and
          (bbox_id == 0 or group.disttype != CCTK_DISTRIB_CONSTANT)) {
        // Construct a file-wide unique HDF5 dataset name
        // (only add parts with varying metadata)
        ostringstream datasetname;
        datasetname << fullname
                    << " it=" << cctkGH->cctk_iteration
                    << " tl=" << request->timelevel;
        if (mglevels > 1) datasetname << " ml=" << mglevel;
        if (group.grouptype == CCTK_GF) {
          if (maps > 1) datasetname << " m="  << Carpet::map;
          datasetname << " rl=" << refinementlevel;
        }
        if (bboxes.setsize () > 1 and group.disttype != CCTK_DISTRIB_CONSTANT) {
          datasetname << " c=" << bbox_id;
        }

        // We create a temporary HDF5 file in memory (using the core VFD)
        // in order to do the recombination of individual processor components
        // for a single dataset.
        // Although this requires some more memory, it should be much faster
        // than recombining an HDF5 dataset on a disk file.
        HDF5_ERROR (memfile = H5Fcreate ("tempfile", H5F_ACC_EXCL, H5P_DEFAULT,
                                         plist));
        HDF5_ERROR (dataspace = H5Screate_simple (group.dim, shape, NULL));
        HDF5_ERROR (memdataset = H5Dcreate (memfile, datasetname.str().c_str(),
                                         filedatatype, dataspace, H5P_DEFAULT));

        // remove an already existing dataset of the same name
        if (request->check_exist) {
          H5E_BEGIN_TRY {
            H5Gunlink (outfile, datasetname.str().c_str());
          } H5E_END_TRY;
        }
        // enable compression if requested
        hid_t plist;
        HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_CREATE));
        const int compression_lvl = request->compression_level >= 0 ?
                                    request->compression_level :
                                    compression_level;
        if (compression_lvl) {
          HDF5_ERROR (H5Pset_chunk (plist, group.dim, shape));
          HDF5_ERROR (H5Pset_deflate (plist, compression_lvl));
        }
        // enable checksums if requested
        if (use_checksums) {
          HDF5_ERROR (H5Pset_chunk (plist, group.dim, shape));
          HDF5_ERROR (H5Pset_filter (plist, H5Z_FILTER_FLETCHER32, 0, 0, NULL));
        }
        HDF5_ERROR (dataset = H5Dcreate (outfile, datasetname.str().c_str(),
                                         filedatatype, dataspace, plist));
        HDF5_ERROR (H5Pclose (plist));
      }

      // Loop over all components
      bool first_time = true;
      BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype) {
        // Get the intersection of the current component with this combination
        // (use either the interior or exterior here, as we did above)
        ibbox const overlap = *bbox &
          arrdata.at(gindex).at(Carpet::map).dd->
          boxes.at(mglevel).at(refinementlevel).at(component).interior;

        // Continue if this component is not part of this combination
        if (overlap.empty()) continue;

        // Copy the overlap to the local processor
        const ggf* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
        const gdata* const data = (*ff) (request->timelevel,
                                         refinementlevel,
                                         component, mglevel);
        gdata* const processor_component =
          data->make_typed (request->vindex, error_centered, op_sync);

        processor_component->allocate (overlap, 0);
        for (comm_state state; not state.done(); state.step()) {
          processor_component->copy_from (state, data, overlap);
        }

        // Write data
        if (dist::rank() == 0) {
          const void *data = (const void *) processor_component->storage();

          // As per Cactus convention, DISTRIB=CONSTANT arrays
          // (including grid scalars) are assumed to be the same on
          // all processors and are therefore stored only by processor
          // 0.
          //
          // Warn the user if this convention is violated.
          if (bbox_id > 0 and group.disttype == CCTK_DISTRIB_CONSTANT) {
            const void *proc0 = CCTK_VarDataPtrI (cctkGH, request->timelevel,
                                                  request->vindex);

            if (memcmp (proc0, data,
                        num_elems * CCTK_VarTypeSize(group.vartype))) {
              CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "values for DISTRIB=CONSTANT grid variable '%s' "
                          "(timelevel %d) differ between processors 0 and %d; "
                          "only the array from processor 0 will be stored",
                          fullname, request->timelevel, bbox_id);
            }
          } else {
            hsize_t overlapshape[dim];

            // before HDF5-1.6.4 the H5Sselect_hyperslab() function expected
            // the 'start' argument to be of type 'hssize_t'
            slice_start_size_t overlaporigin[dim];
            for (int d = 0; d < group.dim; ++d) {
              assert (group.dim-1-d>=0 and group.dim-1-d<dim);
              overlaporigin[group.dim-1-d] =
                ((overlap.lower() - bbox->lower()) / overlap.stride())[d];
              overlapshape[group.dim-1-d]  =
                (overlap.shape() / overlap.stride())[d];
            }

            // Write the processor component as a hyperslab into the recombined
            // dataset.
            hid_t overlap_dataspace;
            HDF5_ERROR (overlap_dataspace =
                        H5Screate_simple (group.dim, overlapshape, NULL));
            HDF5_ERROR (H5Sselect_hyperslab (dataspace, H5S_SELECT_SET,
                                             overlaporigin, NULL,
                                             overlapshape, NULL));
            HDF5_ERROR (H5Dwrite (memdataset, memdatatype, overlap_dataspace,
                                  dataspace, H5P_DEFAULT, data));
            io_bytes +=
              H5Sget_simple_extent_npoints (overlap_dataspace) *
              H5Tget_size (filedatatype);
            HDF5_ERROR (H5Sclose (overlap_dataspace));

            // Add metadata information on the first time through
            // (have to do it inside of the COMPONENT_LOOP so that we have
            //  access to the cGH elements)
            if (first_time) {
              error_count += AddAttributes (cctkGH, fullname, group.dim,
                                            refinementlevel, request, *bbox,
                                            dataset);
              first_time = false;
            }
          }
        }

        // Delete temporary copy of this component
        delete processor_component;

      } END_COMPONENT_LOOP;

      // Finally create the recombined dataset in the real HDF5 file on disk
      // (skip DISTRIB=CONSTANT components from processors other than 0)
      if (dist::rank() == 0 and
          (bbox_id == 0 or group.disttype != CCTK_DISTRIB_CONSTANT)) {
        void *data = malloc (H5Dget_storage_size (memdataset));
        HDF5_ERROR (H5Dread (memdataset, filedatatype, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, data));
        HDF5_ERROR (H5Dclose (memdataset));
        HDF5_ERROR (H5Fclose (memfile));
        HDF5_ERROR (H5Dwrite (dataset, filedatatype, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, data));
        free (data);
        HDF5_ERROR (H5Sclose (dataspace));
        HDF5_ERROR (H5Dclose (dataset));
      }

    } // for bboxes
  } END_MAP_LOOP;

  HDF5_ERROR (H5Pclose (plist));

  free (fullname);

  // return the number of errors that occured during this output
  return error_count;
}


int WriteVarChunkedSequential (const cGH* const cctkGH,
                               hid_t outfile,
                               CCTK_REAL & io_bytes,
                               const ioRequest* const request,
                               bool called_from_checkpoint)
{
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;
  char *fullname = CCTK_FullName(request->vindex);
  const int gindex = CCTK_GroupIndexFromVarI (request->vindex);
  assert (gindex >= 0 and gindex < (int) Carpet::arrdata.size ());
  const int var = request->vindex - CCTK_FirstVarIndexI (gindex);
  assert (var >= 0 and var < CCTK_NumVars ());
  cGroup group;
  CCTK_GroupData (gindex, &group);


  // Scalars and arrays have only one refinement level 0,
  // regardless of what the current refinement level is.
  // Output for them must be called in global mode.
  int refinementlevel = reflevel;
  if (group.grouptype == CCTK_SCALAR or group.grouptype == CCTK_ARRAY) {
    assert (do_global_mode);
    refinementlevel = 0;
  }

  // HDF5 doesn't like 0-dimensional arrays
  if (group.grouptype == CCTK_SCALAR) group.dim = 1;

  // If the user requested so, select single precision data output
  hid_t memdatatype, filedatatype;
  HDF5_ERROR (memdatatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype, 0));
  HDF5_ERROR (filedatatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype,
                         out_single_precision and not  called_from_checkpoint));

  // Traverse all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype) {
    BEGIN_COMPONENT_LOOP (cctkGH, group.grouptype) {
      // Using "exterior" includes ghost zones and refinement boundaries.
      ibbox& bbox = arrdata.at(gindex).at(Carpet::map).dd->
                    boxes.at(mglevel).at(refinementlevel).at(component).exterior;

      // Get the shape of the HDF5 dataset (in Fortran index order)
      hsize_t shape[dim];
      hsize_t num_elems = 1;
      for (int d = 0; d < group.dim; ++d) {
        assert (group.dim-1-d>=0 and group.dim-1-d<dim);
        shape[group.dim-1-d] = (bbox.shape() / bbox.stride())[d];
        num_elems *= shape[group.dim-1-d];
      }

      // Don't create zero-sized components
      if (num_elems == 0) continue;

      // Copy the overlap to the local processor
      const ggf* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
      const gdata* const data = (*ff) (request->timelevel, refinementlevel,
                                       component, mglevel);
      gdata* const processor_component =
        data->make_typed (request->vindex,error_centered, op_sync);

      processor_component->allocate (bbox, 0);
      for (comm_state state; not state.done(); state.step()) {
        processor_component->copy_from (state, data, bbox);
      }

      // Write data on I/O processor 0
      if (dist::rank() == 0) {
        const void *data = (const void *) processor_component->storage();

        // As per Cactus convention, DISTRIB=CONSTANT arrays
        // (including grid scalars) are assumed to be the same on
        // all processors and are therefore stored only by processor 0.
        //
        // Warn the user if this convention is violated.
        if (component > 0 and group.disttype == CCTK_DISTRIB_CONSTANT) {
          const void *proc0 = CCTK_VarDataPtrI (cctkGH, request->timelevel,
                                                request->vindex);

          if (memcmp (proc0, data, num_elems*CCTK_VarTypeSize(group.vartype))) {
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "values for DISTRIB=CONSTANT grid variable '%s' "
                        "(timelevel %d) differ between processors 0 and %d; "
                        "only the array from processor 0 will be stored",
                        fullname, request->timelevel, component);
          }
        } else {
          // Construct a file-wide unique HDF5 dataset name
          // (only add parts with varying metadata)
          ostringstream datasetname;
          datasetname << fullname
                      << " it=" << cctkGH->cctk_iteration
                      << " tl=" << request->timelevel;
          if (mglevels > 1) datasetname << " ml=" << mglevel;
          if (group.grouptype == CCTK_GF) {
            if (maps > 1) datasetname << " m="  << Carpet::map;
            datasetname << " rl=" << refinementlevel;
          }
          if (arrdata.at(gindex).at(Carpet::map).dd->
              boxes.at(mglevel).at(refinementlevel).size () > 1 and
              group.disttype != CCTK_DISTRIB_CONSTANT) {
            datasetname << " c=" << component;
          }

          // remove an already existing dataset of the same name
          if (request->check_exist) {
            H5E_BEGIN_TRY {
              H5Gunlink (outfile, datasetname.str().c_str());
            } H5E_END_TRY;
          }

          hsize_t shape[dim];
          hssize_t origin[dim];
          for (int d = 0; d < group.dim; ++d) {
            assert (group.dim-1-d>=0 and group.dim-1-d<dim);
            origin[group.dim-1-d] = (bbox.lower() / bbox.stride())[d];
            shape[group.dim-1-d]  = (bbox.shape() / bbox.stride())[d];
          }

          // Write the component as an individual dataset
          hid_t plist, dataspace, dataset;
          HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_CREATE));
          // enable compression if requested
          const int compression_lvl = request->compression_level >= 0 ?
                                      request->compression_level :
                                      compression_level;
          if (compression_lvl) {
            HDF5_ERROR (H5Pset_chunk (plist, group.dim, shape));
            HDF5_ERROR (H5Pset_deflate (plist, compression_lvl));
          }
          // enable checksums if requested
          if (use_checksums) {
            HDF5_ERROR (H5Pset_chunk (plist, group.dim, shape));
            HDF5_ERROR (H5Pset_filter (plist, H5Z_FILTER_FLETCHER32, 0, 0, NULL));
          }
          HDF5_ERROR (dataspace = H5Screate_simple (group.dim, shape, NULL));
          HDF5_ERROR (dataset = H5Dcreate (outfile, datasetname.str().c_str(),
                                           filedatatype, dataspace, plist));
          io_bytes +=
            H5Sget_simple_extent_npoints (dataspace) *
            H5Tget_size (filedatatype);
          HDF5_ERROR (H5Pclose (plist));
          HDF5_ERROR (H5Sclose (dataspace));
          HDF5_ERROR (H5Dwrite (dataset, memdatatype, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, data));
          error_count += AddAttributes (cctkGH, fullname, group.dim,
                                        refinementlevel, request, bbox,dataset);
          HDF5_ERROR (H5Dclose (dataset));
        }

      } // if dist::rank() == 0

      // Delete temporary copy of this component
      delete processor_component;

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  free (fullname);

  // return the number of errors that occured during this output
  return error_count;
}


int WriteVarChunkedParallel (const cGH* const cctkGH,
                             hid_t outfile,
                             CCTK_REAL & io_bytes,
                             const ioRequest* const request,
                             bool called_from_checkpoint,
                             hid_t indexfile)
{
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;
  char *fullname = CCTK_FullName(request->vindex);
  const int gindex = CCTK_GroupIndexFromVarI (request->vindex);
  assert (gindex >= 0 and gindex < (int) Carpet::arrdata.size ());
  const int var = request->vindex - CCTK_FirstVarIndexI (gindex);
  assert (var >= 0 and var < CCTK_NumVars ());
  cGroup group;
  CCTK_GroupData (gindex, &group);


  // Scalars and arrays have only one refinement level 0,
  // regardless of what the current refinement level is.
  // Output for them must be called in global mode.
  int refinementlevel = reflevel;
  if (group.grouptype == CCTK_SCALAR or group.grouptype == CCTK_ARRAY) {
    assert (do_global_mode);
    refinementlevel = 0;
  }

  // HDF5 doesn't like 0-dimensional arrays
  if (group.grouptype == CCTK_SCALAR) group.dim = 1;

  // If the user requested so, select single precision data output
  hid_t memdatatype, filedatatype;
  HDF5_ERROR (memdatatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype, 0));
  HDF5_ERROR (filedatatype = CCTKtoHDF5_Datatype (cctkGH, group.vartype,
                          out_single_precision and not called_from_checkpoint));

  // Traverse all maps
  BEGIN_MAP_LOOP (cctkGH, group.grouptype) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, group.grouptype) {

      const ggf* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
      const ibbox& bbox = (*ff) (request->timelevel, refinementlevel,
                                 group.disttype == CCTK_DISTRIB_CONSTANT ?
                                 0 : component, mglevel)->extent();

      // Don't create zero-sized components
      if (bbox.empty()) continue;

      // As per Cactus convention, DISTRIB=CONSTANT arrays
      // (including grid scalars) are assumed to be the same on
      // all processors and are therefore stored only by processor 0.
      void* data = cctkGH->data[request->vindex][request->timelevel];
      const void* mydata = data;
      if (group.disttype == CCTK_DISTRIB_CONSTANT) {

        MPI_Datatype datatype;
        switch (group.vartype) {
#define TYPECASE(N,T)                                                     \
          case  N: { T dummy; datatype = dist::datatype(dummy); } break;
#include "carpet_typecase.hh"
#undef TYPECASE
          default: assert (0 and "invalid datatype");
        }

        const size_t size = bbox.size() * CCTK_VarTypeSize (group.vartype);
        if (dist::rank() > 0) {
          data = malloc (size);
        }
        MPI_Bcast (data, bbox.size(), datatype, 0, MPI_COMM_WORLD);

        if (memcmp (mydata, data, size)) {
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "values for DISTRIB=CONSTANT grid variable '%s' "
                      "(timelevel %d) differ between processors 0 and %d; "
                      "only the array from processor 0 will be stored",
                      fullname, request->timelevel, component);
        }
      }

      // Construct a file-wide unique HDF5 dataset name
      // (only add parts with varying metadata)
      ostringstream datasetname;
      datasetname << fullname
                  << " it=" << cctkGH->cctk_iteration
                  << " tl=" << request->timelevel;
      if (mglevels > 1) datasetname << " ml=" << mglevel;
      if (group.grouptype == CCTK_GF) {
        if (maps > 1) datasetname << " m="  << Carpet::map;
        datasetname << " rl=" << refinementlevel;
      }
      if (arrdata.at(gindex).at(Carpet::map).dd->
          boxes.at(mglevel).at(refinementlevel).size () > 1 and
          group.disttype != CCTK_DISTRIB_CONSTANT) {
        datasetname << " c=" << component;
      }

      // remove an already existing dataset of the same name
      if (request->check_exist) {
        H5E_BEGIN_TRY {
          H5Gunlink (outfile, datasetname.str().c_str());
          if (indexfile != -1)
            H5Gunlink (indexfile, datasetname.str().c_str());
        } H5E_END_TRY;
      }

      // Get the shape of the HDF5 dataset (in Fortran index order)
      hsize_t shape[dim];
      hsize_t index_shape[dim];
      hssize_t origin[dim];
      for (int d = 0; d < group.dim; ++d) {
        assert (group.dim-1-d>=0 and group.dim-1-d<dim);
        origin[group.dim-1-d] = (bbox.lower() / bbox.stride())[d];
        shape[group.dim-1-d]  = (bbox.shape() / bbox.stride())[d];
        index_shape[group.dim-1-d]  = 1;
      }

      // Write the component as an individual dataset
      hid_t plist, dataspace, dataset, index_dataspace, index_dataset;
      HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_CREATE));
      // enable compression if requested
      const int compression_lvl = request->compression_level >= 0 ?
                                  request->compression_level :
                                  compression_level;
      if (compression_lvl) {
        HDF5_ERROR (H5Pset_chunk (plist, group.dim, shape));
        HDF5_ERROR (H5Pset_deflate (plist, compression_lvl));
      }
      // enable checksums if requested
      if (use_checksums) {
        HDF5_ERROR (H5Pset_chunk (plist, group.dim, shape));
        HDF5_ERROR (H5Pset_filter (plist, H5Z_FILTER_FLETCHER32, 0, 0, NULL));
      }
      HDF5_ERROR (dataspace = H5Screate_simple (group.dim, shape, NULL));
      HDF5_ERROR (dataset = H5Dcreate (outfile, datasetname.str().c_str(),
                                       filedatatype, dataspace, plist));

      if (indexfile != -1) {
        HDF5_ERROR (index_dataspace = H5Screate_simple (group.dim,
                                                        index_shape, NULL));
        HDF5_ERROR (index_dataset = H5Dcreate (indexfile, datasetname.str().c_str(),
                                               filedatatype, index_dataspace, H5P_DEFAULT));
      }

      io_bytes +=
        H5Sget_simple_extent_npoints (dataspace) * H5Tget_size (filedatatype);
      HDF5_ERROR (H5Pclose (plist));
      HDF5_ERROR (H5Sclose (dataspace));
      HDF5_ERROR (H5Dwrite (dataset, memdatatype, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, data));
      error_count += AddAttributes (cctkGH, fullname, group.dim,refinementlevel,
                                    request, bbox, dataset);
      HDF5_ERROR (H5Dclose (dataset));

      if (indexfile != -1) {
        HDF5_ERROR (H5Sclose (index_dataspace));
        error_count += AddAttributes (cctkGH, fullname, group.dim,refinementlevel,
                                      request, bbox, index_dataset, true);
        HDF5_ERROR (H5Dclose (index_dataset));
      }

      if (data != mydata) free (data);

    } END_LOCAL_COMPONENT_LOOP;
  } END_MAP_LOOP;

  free (fullname);

  // return the number of errors that occured during this output
  return error_count;
}


// add attributes to an HDF5 dataset
static int AddAttributes (const cGH *const cctkGH, const char *fullname,
                          int vdim, int refinementlevel,
                          const ioRequest* request,
                          const ibbox& bbox, hid_t dataset, bool is_index)
{
  assert (vdim>=0 and vdim<=dim);
  int error_count = 0;

  // Legacy arguments
  hid_t attr, dataspace, datatype;
  HDF5_ERROR (dataspace = H5Screate (H5S_SCALAR));
  HDF5_ERROR (attr = H5Acreate (dataset, "level", H5T_NATIVE_INT,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &refinementlevel));
  HDF5_ERROR (H5Aclose (attr));

  HDF5_ERROR (attr = H5Acreate (dataset, "carpet_mglevel", H5T_NATIVE_INT,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &mglevel));
  HDF5_ERROR (H5Aclose (attr));

  HDF5_ERROR (attr = H5Acreate (dataset, "timestep", H5T_NATIVE_INT,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &cctkGH->cctk_iteration));
  HDF5_ERROR (H5Aclose (attr));

  HDF5_ERROR (attr = H5Acreate (dataset, "group_timelevel", H5T_NATIVE_INT,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &request->timelevel));
  HDF5_ERROR (H5Aclose (attr));

  HDF5_ERROR (attr = H5Acreate (dataset, "time", HDF5_REAL,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, HDF5_REAL, &cctkGH->cctk_time));
  HDF5_ERROR (H5Aclose (attr));

  HDF5_ERROR (datatype = H5Tcopy (H5T_C_S1));
  HDF5_ERROR (H5Tset_size (datatype, strlen (fullname)));
  HDF5_ERROR (attr = H5Acreate (dataset, "name", datatype,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, datatype, fullname));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (H5Tclose (datatype));
  
  // Specify whether the coordinate system is Cartesian or not
  if (CCTK_IsFunctionAliased ("MultiPatch_MapIsCartesian")) {
    int const map_is_cartesian = MultiPatch_MapIsCartesian (Carpet::map);
    HDF5_ERROR (attr = H5Acreate (dataset, "MapIsCartesian", H5T_NATIVE_INT,
                                  dataspace, H5P_DEFAULT));
    HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, & map_is_cartesian));
    HDF5_ERROR (H5Aclose (attr));
  }
  
  HDF5_ERROR (H5Sclose (dataspace));

  // store cctk_bbox and cctk_nghostzones (for grid arrays only)
  if (CCTK_GroupTypeFromVarI (request->vindex) != CCTK_SCALAR) {
    vector<int> cctk_bbox(2*vdim);
    hsize_t size = cctk_bbox.size();
    HDF5_ERROR (dataspace = H5Screate_simple (1, &size, NULL));
    CCTK_GroupbboxVI (cctkGH, size, &cctk_bbox[0], request->vindex);
    HDF5_ERROR (attr = H5Acreate (dataset, "cctk_bbox", H5T_NATIVE_INT,
                                  dataspace, H5P_DEFAULT));
    HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &cctk_bbox[0]));
    HDF5_ERROR (H5Aclose (attr));
    HDF5_ERROR (H5Sclose (dataspace));

    ivect cctk_nghostzones;
    size = vdim;
    HDF5_ERROR (dataspace = H5Screate_simple (1, &size, NULL));
    CCTK_GroupnghostzonesVI (cctkGH, size, &cctk_nghostzones[0],
                             request->vindex);
    HDF5_ERROR (attr = H5Acreate (dataset, "cctk_nghostzones", H5T_NATIVE_INT,
                                  dataspace, H5P_DEFAULT));
    HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &cctk_nghostzones[0]));
    HDF5_ERROR (H5Aclose (attr));
    HDF5_ERROR (H5Sclose (dataspace));
  }

  // write bbox attributes if we have coordinate system info
  CCTK_REAL origin[dim], delta[dim];
  int coord_system_handle = -1;
  if (CCTK_IsFunctionAliased ("Coord_GroupSystem"))
  {
    char *groupname = CCTK_GroupNameFromVarI (request->vindex);
    coord_system_handle = Coord_GroupSystem (cctkGH, groupname);
    free (groupname);
  }

  

  hsize_t size = vdim;
  HDF5_ERROR (dataspace = H5Screate_simple (1, &size, NULL));

  CCTK_INT levoffdenom[dim];
#if 0                           // dh::dbases
    const ibbox& baseext =
      vdd.at(Carpet::map)->bases.at(mglevel).at(reflevel).exterior;
#endif
  const ibbox& baseext =
    vhh.at(Carpet::map)->baseextents.at(mglevel).at(reflevel);

  const ivect pos = (bbox.lower() - baseext.lower()) / bbox.stride();

  rvect global_lower;
  rvect coord_delta;
  const int m = Carpet::map;
  if (CCTK_GroupTypeFromVarI (request->vindex) == CCTK_GF) {
    rvect const cctk_origin_space =
      origin_space.at(m).at(mglevel);
    rvect const cctk_delta_space =
      delta_space.at(m) * rvect (mglevelfact);
    for (int d=0; d<dim; ++d) {
       // lower boundary of Carpet's integer indexing
       global_lower[d] = cctk_origin_space[d];
       // grid spacing of Carpet's integer indexing
       coord_delta[d] =
          cctk_delta_space[d] / cctkGH->cctk_levfac[d];
      levoffdenom[d] = cctkGH->cctk_levoffdenom[d];
    }
  } else {
    for (int d=0; d<dim; ++d) {
      global_lower[d] = 0.0;
      coord_delta[d] = 1.0;
      levoffdenom[d] = 1;
    }
  }
  
  for (int d=0; d<dim; ++d) {
      origin[d] = global_lower[d] + coord_delta[d] * (cctkGH->cctk_levoff[d] / levoffdenom[d] + pos[d]);
      delta[d] = coord_delta[d];
  }

  HDF5_ERROR (attr = H5Acreate (dataset, "origin", HDF5_REAL,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, HDF5_REAL, origin));
  HDF5_ERROR (H5Aclose (attr));
  HDF5_ERROR (attr = H5Acreate (dataset, "delta", HDF5_REAL,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, HDF5_REAL, delta));
  HDF5_ERROR (H5Aclose (attr));

  ivect iorigin = bbox.lower() / bbox.stride();
  HDF5_ERROR (attr = H5Acreate (dataset, "iorigin", H5T_NATIVE_INT,
                                dataspace, H5P_DEFAULT));
  HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &iorigin[0]));
  HDF5_ERROR (H5Aclose (attr));

  hsize_t shape[vdim];
  for (int d = 0; d < vdim; ++d) {
    assert (vdim-1-d>=0 and vdim-1-d<vdim);
    shape[vdim-1-d]  = (bbox.shape() / bbox.stride())[d];
  }

  if (is_index) {
    HDF5_ERROR (attr = H5Acreate (dataset, "h5shape", H5T_NATIVE_HSIZE,
                                  dataspace, H5P_DEFAULT));
    HDF5_ERROR (H5Awrite (attr, H5T_NATIVE_HSIZE, &shape[0]));
    HDF5_ERROR (H5Aclose (attr));
  }

  HDF5_ERROR (H5Sclose (dataspace));

  // return the number of errors that occured during this output
  return error_count;
}

} // namespace CarpetIOHDF5
