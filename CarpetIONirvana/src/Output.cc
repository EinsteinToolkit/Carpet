#include <cassert>
#include <cstdlib>
#include <cstring>
#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CarpetIONirvana.hh"
#include "operators.hh"
#include "typeprops.hh"

namespace CarpetIONirvana {

using namespace std;
using namespace Carpet;
using namespace Nirvana;

static int GetMetaData(const cGH *const cctkGH, const char *name,
                       const char *groupname, int filenum, int vdim,
                       int refinementlevel, const ioRequest *request,
                       const ibbox &bbox, metadata &md);

int WriteVar(const cGH *const cctkGH, const string &filename, const int filenum,
             CCTK_REAL &io_bytes, const ioRequest *const request) {
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;
  char *fullname = CCTK_FullName(request->vindex);
  const int gindex = CCTK_GroupIndexFromVarI(request->vindex);
  assert(gindex >= 0 and gindex < (int)Carpet::arrdata.size());
  char *groupname = CCTK_GroupName(gindex);
  const int var = request->vindex - CCTK_FirstVarIndexI(gindex);
  assert(var >= 0 and var < CCTK_NumVars());
  cGroup group;
  CCTK_GroupData(gindex, &group);

  // Scalars and arrays have only one refinement level 0,
  // regardless of what the current refinement level is.
  // Output for them must be called in global mode.
  int refinementlevel = reflevel;
  if (group.grouptype == CCTK_SCALAR or group.grouptype == CCTK_ARRAY) {
    assert(do_global_mode);
    refinementlevel = 0;
  }

  // HDF5 doesn't like 0-dimensional arrays
  if (group.grouptype == CCTK_SCALAR)
    group.dim = 1;

  // Traverse all maps
  BEGIN_LOCAL_MAP_LOOP(cctkGH, group.grouptype) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, group.grouptype) {

      // const ggf* ff = arrdata.at(gindex).at(Carpet::map).data.at(var);
      // const ibbox& bbox = (*ff) (request->timelevel, refinementlevel,
      //                            group.disttype == CCTK_DISTRIB_CONSTANT ?
      //                            0 : component, mglevel)->extent();
      const dh *dd = arrdata.at(gindex).at(Carpet::map).dd;
      const ibbox &bbox =
          dd->light_boxes.AT(mglevel)
              .AT(refinementlevel)
              .AT(group.disttype == CCTK_DISTRIB_CONSTANT ? 0 : component)
              .exterior;

      // Don't create zero-sized components
      if (bbox.empty())
        continue;

      // As per Cactus convention, DISTRIB=CONSTANT arrays
      // (including grid scalars) are assumed to be the same on
      // all processors and are therefore stored only by processor 0.
      void *data = cctkGH->data[request->vindex][request->timelevel];
      const void *mydata = data;
      if (group.disttype == CCTK_DISTRIB_CONSTANT) {

        MPI_Datatype datatype;
        switch (specific_cactus_type(group.vartype)) {
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T dummy;                                                                   \
    datatype = dist::mpi_datatype(dummy);                                      \
  } break;
#include "typecase.hh"
#undef TYPECASE
        default:
          assert(0 and "invalid datatype");
        }

        const size_t size = bbox.size() * CCTK_VarTypeSize(group.vartype);
        if (dist::rank() > 0) {
          data = malloc(size);
        }
        MPI_Bcast(data, bbox.size(), datatype, 0, MPI_COMM_WORLD);

        if (memcmp(mydata, data, size)) {
          CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "values for DISTRIB=CONSTANT grid variable '%s' "
                     "(timelevel %d) differ between processors 0 and %d; "
                     "only the array from processor 0 will be stored",
                     fullname, request->timelevel, component);
        }
      }

      vector<double> buffer(bbox.size());
      for (int i = 0; i < bbox.size(); ++i)
        buffer[i] = ((double *)data)[i];

      // the metadata describing the mesh
      metadata md;

      error_count += GetMetaData(cctkGH, fullname, groupname, filenum,
                                 group.dim, refinementlevel, request, bbox, md);

      // cout << fullname << " " << cctkGH->cctk_iteration << endl;
      string past_checkpoint = "";
      const bool write_only =
          true; // this prevents the file-format to scan existing files
      CarpetN5 File(filename, write_only, dist::size(), false, past_checkpoint);

      File.WriteMesh(md);
      File.WriteMesh(md, filenum);
      File.WriteField<double>(md, buffer);
      File.WriteField<double>(md, buffer, filenum);

      if (data != mydata)
        free(data);
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_LOCAL_MAP_LOOP;

  free(fullname);
  free(groupname);

  // return the number of errors that occured during this output
  return error_count;
}

// add attributes to an HDF5 dataset
static int GetMetaData(const cGH *const cctkGH, const char *name,
                       const char *groupname, int filenum, int vdim,
                       int refinementlevel, const ioRequest *request,
                       const ibbox &bbox, metadata &md) {
  assert(vdim >= 0 and vdim <= dim);
  int error_count = 0;

  // Get the shape of the dataset
  vector<int> shape(dim, 0);
  vector<int> iorigin(dim, 0);
  for (int d = 0; d < vdim; ++d) {
    iorigin[d] = (bbox.lower() / bbox.stride())[d];
    shape[d] = (bbox.shape() / bbox.stride())[d];
  }

  int size = 3;
  vector<int> cctk_nghostzones(size);
  CCTK_GroupnghostzonesVI(cctkGH, size, &cctk_nghostzones[0], request->vindex);

  attribute<int> acycle("cycle", cctkGH->cctk_iteration);
  attribute<int> areflevel("reflevel", refinementlevel);
  attribute<int> atimelevel("timelevel", request->timelevel);
  attribute<int> amap("map", Carpet::map);
  attribute<int> afilenum("filenum", filenum);
  attribute<double> atime("time", cctkGH->cctk_time);
  attribute<vector<int> > alsh("lsh", shape);
  attribute<vector<int> > aghosts("ghosts", cctk_nghostzones);
  attribute<vector<int> > aiorigin("iorigin", iorigin);
  attribute<int> acomp("comp", Carpet::component);
  attribute<string> avarname("varname", name);
  attribute<string> avargroupname("vargroupname", groupname);
  attribute<string> avargrouptype;
  attribute<string> ameshname;
  attribute<string> acoordinates;
  attribute<vector<double> > aorigin;
  attribute<vector<double> > adelta;

  // Specify whether the coordinate system is Cartesian or not
  if (CCTK_IsFunctionAliased("MultiPatch_MapIsCartesian")) {
    int const map_is_cartesian = MultiPatch_MapIsCartesian(Carpet::map);
    if (!map_is_cartesian) {
      ameshname = attribute<string>("meshname", "Curvilinear");
      acoordinates = attribute<string>(
          "coordinates",
          "grid::coordinates"); // the group-name of the coordinates
    } else {
      ameshname = attribute<string>("meshname", "Carpet-AMR");
      acoordinates = attribute<string>(
          "coordinates",
          "none"); // we don't need external coordinates to represent mesh
    }
  } else {
    ameshname = attribute<string>("meshname", "Carpet-AMR");
    acoordinates = attribute<string>(
        "coordinates",
        "none"); // we don't need external coordinates to represent mesh
  }

  // write bbox attributes if we have coordinate system info
  CCTK_REAL origin[dim], delta[dim];
  int coord_system_handle = -1;
  if (CCTK_IsFunctionAliased("Coord_GroupSystem")) {
    char *groupname = CCTK_GroupNameFromVarI(request->vindex);
    coord_system_handle = Coord_GroupSystem(cctkGH, groupname);
    free(groupname);
  }

  CCTK_INT coord_handles[dim];
  if (coord_system_handle >= 0 and
      Util_TableGetIntArray(coord_system_handle, vdim, coord_handles,
                            "COORDINATES") >= 0) {
#if 0 // dh::dbases
    const ibbox& baseext =
      vdd.at(Carpet::map)->bases.at(mglevel).at(reflevel).exterior;
#endif
    const ibbox &baseext =
        vhh.at(Carpet::map)->baseextents.at(mglevel).at(reflevel);

    const ivect pos = (bbox.lower() - baseext.lower()) / bbox.stride();

    rvect global_lower;
    rvect coord_delta;
    const int m = Carpet::map;
    if (CCTK_GroupTypeFromVarI(request->vindex) == CCTK_GF) {
      rvect const cctk_origin_space = origin_space.at(m).at(mglevel);
      rvect const cctk_delta_space = delta_space.at(m) * rvect(mglevelfact);
      for (int d = 0; d < dim; ++d) {
        // lower boundary of Carpet's integer indexing
        global_lower[d] = cctk_origin_space[d];
        // grid spacing of Carpet's integer indexing
        coord_delta[d] = cctk_delta_space[d] / cctkGH->cctk_levfac[d];
      }
    } else {
      for (int d = 0; d < dim; ++d) {
        global_lower[d] = 0.0;
        coord_delta[d] = 1.0;
      }
    }

    vector<fp> origin(dim, 0);
    vector<fp> delta(dim, 0);
    for (int d = 0; d < dim; ++d) {
      origin[d] =
          global_lower[d] +
          coord_delta[d] *
              (cctkGH->cctk_levoff[d] / cctkGH->cctk_levoffdenom[d] + pos[d]);
      delta[d] = coord_delta[d];
    }

    aorigin = attribute<vector<fp> >("origin", origin);
    adelta = attribute<vector<fp> >("delta", delta);
  }

  // find out tensor type
  char tensortypealias[100];
  int table;
  const int gindex = CCTK_GroupIndexFromVarI(request->vindex);
  const int numvars = CCTK_NumVarsInGroupI(gindex);
  assert(numvars > 0);
  table = CCTK_GroupTagsTableI(gindex);
  assert(table >= 0);

  int ierr = Util_TableGetString(table, sizeof tensortypealias, tensortypealias,
                                 "tensortypealias");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    /* assume a scalar */
    strcpy(tensortypealias, "scalar");
  } else if (ierr < 0) {
    char *groupname = CCTK_GroupName(gindex);
    assert(groupname);
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error in tensor type alias declaration for group \"%s\"",
               groupname);
    free(groupname);
  }

  if (CCTK_EQUALS(tensortypealias, "scalar")) {
    /* scalar */
    avargrouptype = attribute<string>("vargrouptype", "scalar");
  } else if (CCTK_EQUALS(tensortypealias, "4scalar")) {
    /* 4-scalar */
    avargrouptype = attribute<string>("vargrouptype", "4scalar");
  } else if (CCTK_EQUALS(tensortypealias, "u") ||
             CCTK_EQUALS(tensortypealias, "d")) {
    /* vector */
    avargrouptype = attribute<string>("vargrouptype", "vector");
    assert(numvars == dim);
  } else if (CCTK_EQUALS(tensortypealias, "4u") ||
             CCTK_EQUALS(tensortypealias, "4d")) {
    /* 4-vector */
    avargrouptype = attribute<string>("vargrouptype", "4scalar");
    assert(numvars == dim + 1);
  } else if (CCTK_EQUALS(tensortypealias, "du")) {
    /* tensor */
    avargrouptype = attribute<string>("vargrouptype", "tensor");
    assert(numvars == dim * dim);
  } else if (CCTK_EQUALS(tensortypealias, "uu_sym") ||
             CCTK_EQUALS(tensortypealias, "dd_sym")) {
    /* symmetric tensor */
    avargrouptype = attribute<string>("vargrouptype", "symmetric-tensor");
    assert(numvars == dim * (dim + 1) / 2);
  } else if (CCTK_EQUALS(tensortypealias, "ddd_sym")) {
    /* 3rd rank tensor, symmetric in last 2 indices */
    avargrouptype = attribute<string>("vargrouptype", "rank3-tensor-symmetric");
    assert(numvars == dim * dim * (dim + 1) / 2);
  } else if (CCTK_EQUALS(tensortypealias, "4uu_sym") ||
             CCTK_EQUALS(tensortypealias, "4dd_sym")) {
    /* symmetric 4-tensor */
    avargrouptype = attribute<string>("vargrouptype", "4tensor-symmetric");
    assert(numvars == (dim + 1) * (dim + 2) / 2);
  } else if (CCTK_EQUALS(tensortypealias, "weylscalars_real")) {
    /* Weyl scalars, stored as 10 real numbers */
    avargrouptype = attribute<string>("vargrouptype", "scalar");
    assert(numvars == 10);
  } else {
    char *groupname = CCTK_GroupName(gindex);
    assert(groupname);
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Illegal tensor type alias for group \"%s\"", groupname);
    free(groupname);
  }

  // Append attributes to metadata
  md << afilenum;     // the filenumber that contains the data
  md << ameshname;    // the name of the mesh
  md << acycle;       // the current cycle/iteration
  md << atime;        // the time corresponding to cycle
  md << amap;         // the map
  md << areflevel;    // the refinement level
  md << atimelevel;   // the time level
  md << acomp;        // the grid component
  md << alsh;         // the number of points of the local piece of grid data
  md << aiorigin;     // the integer origin of the grid
  md << aorigin;      // the origin in coordinate space
  md << adelta;       // the delta spacing
  md << acoordinates; // the name of the coordinate group
  md << aghosts;      // number of ghostzones

  md << avarname;      // the name of the variable
  md << avargroupname; // the name of the variable's group
  md << avargrouptype; // tensor type

  // return the number of errors that occured during this output
  return error_count;
}

} // namespace CarpetIONirvana
