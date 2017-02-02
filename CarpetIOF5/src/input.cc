#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>

#include <F5/F5F.h>
#include <F5/F5R.h>
#include <F5/F5iterate.h>
#include <F5/F5uniform.h>

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include <bbox.hh>
#include <vect.hh>

#include <carpet.hh>

#include "iof5.hh"

// F5 helper
namespace {
char const *ArrayTypeName(ArrayType const array_type) {
  switch (array_type) {
  case UnknownArrayType:
    return "Unknown";
  case Contiguous:
    return "Contiguous";
  case SeparatedCompound:
    return "SeparatedCompound";
  case Constant:
    return "Constant";
  case FragmentedContiguous:
    return "FragmentedContiguous";
  case FragmentedSeparatedCompound:
    return "FragmentedSeparatedCompound";
  case DirectProduct:
    return "DirectProduct";
  case IndexPermutation:
    return "IndexPermutation";
  case UniformSampling:
    return "UniformSampling";
  case FragmentedUniformSampling:
    return "FragmentedUniformSampling";
  default:
    return "(illegal)";
  }
}
}

namespace CarpetIOF5 {

using namespace std;
using namespace Carpet;

// Use a class for reading data, so that we have an easy way to pass
// variables between the various iterators
class input_iterator_t {
  cGH const *cctkGH;

  vector<bool> const input_var; // whether to input this variable
  bool const input_past_timelevels;
  bool const input_metadata;

  double time;
  char const *gridname;
  char const *topologyname;
  int index_depth; // 0=vertex, 1=cell
  int topological_dimension;
  int reflevel;
  char const *fieldname;
  ArrayType fieldtype; // determines how to interpret the fragment
  int varindex;
  char const *fragmentname;
  int map, component;

  // string chartname;

  scatter_t &scatter;

public:
  input_iterator_t(cGH const *const cctkGH_, vector<bool> const &input_var_,
                   bool const input_past_timelevels_,
                   bool const input_metadata_, scatter_t &scatter_)
      : cctkGH(cctkGH_), input_var(input_var_),
        input_past_timelevels(input_past_timelevels_),
        input_metadata(input_metadata_), time(nan), gridname(NULL),
        topologyname(NULL), index_depth(-1), topological_dimension(-1),
        reflevel(-1), fieldname(NULL), fieldtype(UnknownArrayType),
        varindex(-1), fragmentname(NULL), map(-1), component(-1),
        scatter(scatter_) {}

private:
  void read_timeslice(F5Path *const path) {
    indent_t indent;
    cout << indent << "time=" << time << "\n";

    F5iterate_grids(path, NULL, grid_iterator, this, NULL, NULL);

    // TODO: synchronise all read grid functions
  }

  void read_grid(F5Path *const path) {
    indent_t indent;
    cout << indent << "grid=" << gridname << "\n";

    if (input_metadata) {
      if (reflevel == 0) {
        indent_t indent2;
        cout << indent2 << "reading metadata\n";
        // hid_t const metadata_group = path->Grid_hid;
        ostringstream pathname;
        pathname << FIBER_CONTENT_GRIDS << "/" << gridname;
        hid_t group;
        group = H5Gopen(path->ContentsGroup_hid, pathname.str().c_str(),
                        H5P_DEFAULT);
        assert(group >= 0);
        read_metadata(cctkGH, group);
        herr_t const herr = H5Gclose(group);
        assert(not herr);
      }
    }

    // F5iterate_vertex_fields(path, NULL, field_iterator, this, NULL, NULL);
    F5iterate_topologies(path, NULL, topology_iterator, this);
  }

  void read_topology(F5Path *const path) {
    indent_t indent;

    herr_t herr;

    cout << indent << "topology=" << topologyname << ""
         << " (" << (index_depth == 0 ? "vertex" : "cell") << ")\n"
         << indent << "topological dimension=" << topological_dimension << "\n";

    // Ignore topologies that are only an alias for another topology
    H5G_stat_t stat;
    herr = H5Gget_objinfo(path->Grid_hid, topologyname, false, &stat);
    assert(not herr);
    if (stat.type == H5G_LINK) {
      char linkval[100000];
      herr =
          H5Gget_linkval(path->Grid_hid, topologyname, sizeof linkval, linkval);
      assert(not herr);
      indent_t indent2;
      cout << indent2 << "alias for topology \"" << linkval << "\"\n"
           << indent2 << "ignoring this topology\n";
      return;
    }

    // Determine refinement level from the topology
    // (Could also use topology name instead)
    hsize_t hreffact[FIBER_MAX_RANK];
    int const iret = F5LAget_dimensions(path->Topology_hid,
                                        FIBER_HDF5_REFINEMENT_INFO, hreffact);
    assert(iret == dim);
    hsize_t hreffact2[FIBER_MAX_RANK];
    void *const pret = F5Tpermute_dimensions(path, dim, hreffact2, hreffact);
    assert(pret);
    ivect const reffact = h2v(hreffact2);
    int rl;
    for (rl = 0; rl < reflevels; ++rl) {
      if (all(reffact == Carpet::spacereffacts.AT(rl)))
        break;
    }
    assert(rl < reflevels);
    reflevel = rl;
    cout << indent << "refinement level " << reflevel << "\n";

    // F5iterate_topology_fields
    //   (path, NULL, field_iterator, this, chartname.c_str(), NULL);
    F5iterate_topology_fields(path, NULL, field_iterator, this, NULL, NULL);
  }

  void read_field(F5Path *const path) {
    indent_t indent;
    cout << indent << "field=" << fieldname << "\n";

    interpret_fieldname(cctkGH, fieldname, varindex);
    // TODO: check all variables in the group
    // TODO: loop over all variables in the group
    if (varindex >= 0 and input_var.at(varindex)) {

      int major_version, minor_version, release_version;
      fieldtype = F5Fget_field_enum(path, fieldname, &major_version,
                                    &minor_version, &release_version);
      cout << indent << "field_type=" << ArrayTypeName(fieldtype) << "\n";

      int const iret = F5Fopen(path, fieldname);
      assert(iret);

      // Do we need to iterate over fragments?
      // TODO: Should instead check whether attribute
      // FIBER_HDF5_TYPEID_ATTRIB exists (existence indicates
      // fragmentation)
      int const is_fragmented = F5Fis_fragmented(path, fieldname);
      cout << indent << (is_fragmented ? "fragmented" : "not fragmented")
           << "\n";
      if (is_fragmented) {
        F5iterate_field_fragments(path, NULL, fragment_iterator, this);
      } else {
        read_fragment(path);
      }

      F5Fclose(path);

    } else {
      indent_t indent2;
      cout << indent2 << "ignoring this field\n";
    }
    // TODO: keep track of which fields have been read, and complain
    // about unread ones
    // TODO: keep track of which part of a field has been read, and
    // complain about unread parts
  }

  void read_fragment(F5Path *const path) {
    indent_t indent;

    if (fragmentname) {
      cout << indent << "fragment=" << fragmentname << "\n";
    } else {
      cout << indent << "no fragment\n";
    }

    int major_version, minor_version, release_version;
    fieldtype = F5Fget_field_enum(path, fieldname, &major_version,
                                  &minor_version, &release_version);
    cout << indent << "field_type=" << ArrayTypeName(fieldtype) << "\n";

    // Determine map and component from fragment name
    if (fragmentname) {
      interpret_fragmentname(cctkGH, fragmentname, map, component);
    } else {
      map = 0;
      component = CCTK_MyProc(cctkGH);
    }
    cout << indent << "map " << map << " component " << component << "\n";

    // This routine has a bug, and does not recognise
    // "SeparatedCompound" as separated
    // int const is_separated = F5Fis_separatedcompound(path, fieldname);
    int const is_separated = fieldtype == FragmentedSeparatedCompound or
                             fieldtype == SeparatedCompound;
    cout << indent << (is_separated ? "separated" : "contiguous") << "\n";
    hid_t const type_id = F5Fget_type(path);
    if (F5Tis_convertible(type_id, H5T_NATIVE_DOUBLE)) {
      cout << indent << "compound type: scalar\n";
      read_variable(path, NULL, varindex);
    } else if (F5Tis_convertible(type_id, F5T_VEC3_DOUBLE)) {
      cout << indent << "compound type: vector\n";
      // This assumes separated storage; we don't support contiguous
      // storage yet
      assert(is_separated);
      read_variable(path, "Dx", varindex + 0);
      read_variable(path, "Dy", varindex + 1);
      read_variable(path, "Dz", varindex + 2);
    } else if (F5Tis_convertible(type_id, F5T_METRIC33_DOUBLE)) {
      cout << indent << "compound type: symmetric tensor\n";
      // This assumes separated storage; we don't support contiguous
      // storage yet
      assert(is_separated);
      read_variable(path, "gxx", varindex + 0);
      read_variable(path, "gxy", varindex + 1);
      read_variable(path, "gxz", varindex + 2);
      read_variable(path, "gyy", varindex + 3);
      read_variable(path, "gyz", varindex + 4);
      read_variable(path, "gzz", varindex + 5);
    } else {
      // Unknown tensor type
      assert(0);
    }
  }

  void read_variable(F5Path *const path, char const *const name,
                     int const var) {
    indent_t indent;

    herr_t herr;
    int iret;
    void *pret;

    assert(path);
    assert(var >= 0);

    cout << indent << "dataset=" << (name ? name : "(null)") << "\n";

    assert(var >= 0);
    {
      char *const fullname = CCTK_FullName(var);
      cout << indent << "variable=" << fullname << "\n";
      free(fullname);
    }

    // Determine fragment properties

    // A dataset can either be the same as a field, or can be in a
    // group containing all fragments, and/or can be in a group
    // containing all compound (vector) elements.
    // The possible group hierarchies are (the rightmost element is
    // the dataset):
    //    .../FIELD
    //    .../FIELD/ELEMENT
    //    .../FIELD/FRAGMENT
    //    .../FIELD/FRAGMENT/ELEMENT
    // The field is the Cactus group/variable name, a fragment
    // describes map and component, the element describes the tensor
    // element.
    // We describe this structure with two booleans:
    //    bool fragment_is_group: the fragment is a group
    //    bool field_is_dataset:  the field is a dataset
    // Note that the F5 library will already have openend the field,
    // either as group or as dataset.

    // There is a group for all fragments if (a) we expect a
    // fragment, and (b) there exists a group with this name. (We
    // probably could examine F5's internal state as well, but
    // looking for an HDF5 group is the most robust approach.)

    // Note: F5 has a function F5Fis_group(path) that returns
    // whether the field is a group or a dataset.

    bool fragment_is_group = false;
    if (fragmentname) {
      H5O_info_t info;
      herr = H5Oget_info_by_name(path->Field_hid, fragmentname, &info,
                                 H5P_DEFAULT);
      assert(not herr);
      fragment_is_group = info.type == H5O_TYPE_GROUP;
    }
    cout << indent << "fragment_is_group=" << fragment_is_group << "\n";
    hid_t fragment;
    if (fragment_is_group) {
      fragment = H5Gopen(path->Field_hid, fragmentname, H5P_DEFAULT);
      assert(fragment >= 0);
    } else {
      fragment = path->Field_hid;
    }

    // The field consists of a single dataset if we cannot open the
    // element
    hid_t element = -1;
    char const *datasetname = NULL;
    if (name) {
      datasetname = name;
    } else if (fragmentname and not fragment_is_group) {
      datasetname = fragmentname;
    }
    if (datasetname) {
      H5E_BEGIN_TRY { element = H5Dopen(fragment, datasetname, H5P_DEFAULT); }
      H5E_END_TRY;
      assert(element >= 0);
    }
    bool const field_is_dataset = element < 0;
    cout << indent << "field_is_dataset=" << field_is_dataset << "\n";
    if (field_is_dataset) {
      // Can't open the fragment -- we assume the field consists of
      // a single element and has already been opened
      element = fragment;
    }
    assert(element >= 0);

    // Check index depth
    int index_depth_;
    iret = F5Tget_index_depth(path, &index_depth_);
    assert(iret);
    assert(index_depth_ == index_depth);

    // Read the fragment offset. This is stored with the dataset
    // group.
    ivect foff = ivect(0);
    if (fragmentname) {
      hsize_t hoff[FIBER_MAX_RANK];
      iret = F5LAget_dimensions(fragment_is_group ? fragment : element,
                                FIBER_FRAGMENT_OFFSET_ATTRIBUTE, hoff);
      assert(iret == dim);
      hsize_t hoff2[FIBER_MAX_RANK];
      pret = F5Tpermute_dimensions(path, dim, hoff2, hoff);
      assert(pret);
      foff = h2v(hoff2);
    }
    assert(all(foff >= 0));

#if 0
      // Read the fragment size. This is stored with the field -- why
      // is this different from the offset?
      hsize_t hlen[FIBER_MAX_RANK];
      iret =
        F5LAget_dimensions(path->Field_hid,
                           FIBER_FIELD_DATASPACE_DIMENSIONS_ATTRIBUTE, hlen);
      assert(iret == dim);
      hsize_t hlen2[FIBER_MAX_RANK];
      pret = F5Tpermute_dimensions(path, dim, hlen2, hlen);
      assert(pret);
      ivect const flen = h2v(hlen2);
      assert(all(flen>=0));
#endif
    hid_t const space = H5Dget_space(element);
    assert(space >= 0);
    iret = H5Sget_simple_extent_ndims(space);
    assert(iret == dim);
    hsize_t hlen[dim];
    iret = H5Sget_simple_extent_dims(space, hlen, NULL);
    hsize_t hlen2[dim];
    pret = F5Tpermute_dimensions(path, dim, hlen2, hlen);
    assert(pret);
    ivect const flen = h2v(hlen2);
    assert(all(flen >= 0));
    herr = H5Sclose(space);
    assert(not herr);

    ibbox const fbox(foff, foff + flen - 1, ivect(1));
    {
      indent_t indent2;
      cout << indent2 << "dataset bbox is " << foff << ":" << foff + flen - 1
           << "\n";
    }

    fragdesc_t fragdesc;
    fragdesc.varindex = var;
    fragdesc.reflevel = reflevel;
    fragdesc.map = map;
    fragdesc.component = component;
    // TODO: set timelevel correctly
    fragdesc.timelevel = 0;
    fragdesc.imin = fbox.lower();
    fragdesc.imax = fbox.upper();

    vector<char> data(fragdesc.npoints() * fragdesc.vartypesize());

    int const vartype = CCTK_VarTypeI(var);
    hid_t type;
    switch (vartype) {
    case CCTK_VARIABLE_INT:
      switch (sizeof(CCTK_INT)) {
      case 1:
        type = H5T_NATIVE_INT8;
        break;
      case 2:
        type = H5T_NATIVE_INT16;
        break;
      case 4:
        type = H5T_NATIVE_INT32;
        break;
      case 8:
        type = H5T_NATIVE_INT64;
        break;
      // case 16: type = H5T_NATIVE_INT128; break;
      default:
        CCTK_ERROR("Unsupported CCTK_INT type");
      }
      break;
    case CCTK_VARIABLE_REAL:
      switch (sizeof(CCTK_REAL)) {
      case 4:
        type = H5T_NATIVE_FLOAT;
        break;
      case 8:
        type = H5T_NATIVE_DOUBLE;
        break;
      case 16:
        type = H5T_NATIVE_LDOUBLE;
        break; // ???
      default:
        CCTK_ERROR("Unsupported CCTK_REAL type");
      }
      break;
    default:
      CCTK_ERROR("Unsupported type");
    }
    assert(type >= 0);
    assert(CCTK_VarTypeSize(vartype) == (int)H5Tget_size(type));

    herr = H5Dread(element, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    assert(not herr);

    scatter.send(fragdesc, &data[0]);

    if (not field_is_dataset) {
      herr = H5Dclose(element);
      assert(not herr);
    }

    if (fragment_is_group) {
      herr = H5Gclose(fragment);
      assert(not herr);
    }
  }

public:
  void iterate(hid_t const object) {
    F5iterate_timeslices(object, NULL, timeslice_iterator, this);
  }

  static herr_t timeslice_iterator(F5Path *const path, double const time,
                                   void *const userdata) {
    input_iterator_t *const iterator = (input_iterator_t *)userdata;
    iterator->time = time;
    iterator->read_timeslice(path);

    if (iterator->input_past_timelevels) {
      return 0; // continue
    } else {
      return 1; // done
    }
    return 0;
  }

  static herr_t grid_iterator(F5Path *const path, char const *const gridname,
                              void *const userdata) {
    input_iterator_t *const iterator = (input_iterator_t *)userdata;
    iterator->gridname = gridname;
    iterator->read_grid(path);
    return 0;
  }

  static herr_t topology_iterator(F5Path *const path,
                                  char const *const topologyname,
                                  int const index_depth,
                                  int const topological_dimension,
                                  void *const userdata) {
    input_iterator_t *const iterator = (input_iterator_t *)userdata;
    iterator->topologyname = topologyname;
    iterator->index_depth = index_depth;
    iterator->topological_dimension = topological_dimension;
    iterator->read_topology(path);
    iterator->topologyname = NULL;
    iterator->index_depth = -1;
    iterator->topological_dimension = -1;
    return 0;
  }

  static herr_t field_iterator(F5Path *const path, char const *const fieldname,
                               void *const userdata) {
    input_iterator_t *const iterator = (input_iterator_t *)userdata;
    iterator->fieldname = fieldname;
    iterator->read_field(path);
    iterator->fieldname = NULL;
    return 0;
  }

  static herr_t fragment_iterator(F5Path *const path,
                                  char const *const fragmentname,
                                  void *const userdata) {
    input_iterator_t *const iterator = (input_iterator_t *)userdata;
    iterator->fragmentname = fragmentname;
    iterator->read_fragment(path);
    iterator->fragmentname = NULL;
    return 0;
  }

}; // class input_iterator_t

void read_metadata(cGH const *const cctkGH, hid_t const file) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH);
  assert(file >= 0);

  herr_t herr;

  CCTK_INFO("Reading simulation metadata...");

  // Open a group holding all metadata
  hid_t const group = H5Gopen(file, metadata_group, H5P_DEFAULT);
  assert(group >= 0);

  // General information
  string fullversion;
  ReadAttribute(group, "Cactus version", fullversion);
  cout << "Cactus version: " << fullversion << "\n";

  // Unique identifiers
  string config_id;
  ReadAttribute(group, "config id", config_id);
  cout << "UniqueConfigID:     " << config_id << "\n";
  string build_id;
  ReadAttribute(group, "build id", build_id);
  cout << "UniqueBuildID:      " << build_id << "\n";
  string simulation_id;
  ReadAttribute(group, "simulation id", simulation_id);
  cout << "UniqueSimulationID: " << simulation_id << "\n";
  string run_id;
  ReadAttribute(group, "run id", run_id);
  cout << "UniqueRunID:        " << run_id << "\n";

  // Parameters
  string parameters;
  ReadLargeAttribute(group, all_parameters, parameters);
  IOUtil_SetAllParameters(parameters.c_str());

  // Grid structure
  string gs;
  ReadLargeAttribute(group, grid_structure, gs);
  deserialise_grid_structure(cctkGH, gs);

  herr = H5Gclose(group);
  assert(not herr);
}

void input(cGH const *const cctkGH, hid_t const file,
           vector<bool> const &input_var, bool const input_past_timelevels,
           bool const input_metadata, scatter_t &scatter) {
  // TODO: not yet implemented
  assert(not input_metadata);
  input_iterator_t iterator(cctkGH, input_var, input_past_timelevels,
                            input_metadata, scatter);
  iterator.iterate(file);
}

} // end namespace CarpetIOF5
