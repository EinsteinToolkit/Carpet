#include "input.hpp"
#include "util.hpp"

#include <carpet.hh>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <iomanip>
#include <memory>
#include <string>

namespace CarpetSimulationIO {
using namespace SimulationIO;
using namespace std;

input_file_t::input_file_t(const cGH *cctkGH, io_dir_t io_dir,
                           const string &projectname, int ioproc, int nioprocs)
    : cctkGH(cctkGH) {
  auto filename = generate_filename(cctkGH, io_dir, projectname, "",
                                    cctkGH->cctk_iteration, -1, -1);
  auto file = H5::H5File(filename, H5F_ACC_RDONLY);
  project = readProject(file);
}

void input_file_t::read_vars(const vector<int> &varindices, int reflevel,
                             int timelevel) const {
  DECLARE_CCTK_PARAMETERS;

  auto parameter_iteration = project->parameters().at("iteration");
  auto parameter_timelevel = project->parameters().at("timelevel");

  for (int varindex : varindices) {
    if (verbose)
      CCTK_VINFO("Reading variable \"%s\"",
                 charptr2string(CCTK_FullName(varindex)).c_str());
    const int groupindex = CCTK_GroupIndexFromVarI(varindex);
    assert(groupindex >= 0);
    const int varindex0 = CCTK_FirstVarIndexI(groupindex);
    assert(varindex0 >= 0);
    const int groupvarindex = varindex - varindex0;
    cGroup groupdata;
    int ierr = CCTK_GroupData(groupindex, &groupdata);
    assert(not ierr);

    const string varname(charptr2string(CCTK_VarName(varindex)));
    const string groupname(tolower(charptr2string(CCTK_GroupName(groupindex))));
    const string fullname(tolower(charptr2string(CCTK_FullName(varindex))));

    // Manifold, TangentSpace, and Basis
    string manifoldname;
    string tangentspacename;
    string basisname;
    if (groupdata.grouptype == CCTK_GF) {
      manifoldname = "manifold";
      tangentspacename = "tangentspace";
      basisname = "basis";
    } else {
      manifoldname = groupname;
      tangentspacename = groupname;
      basisname = groupname;
    }
    const int dimension = groupdata.dim;
    assert(dimension <= Carpet::dim);

    // Manifold
    const auto manifold = project->manifolds().at(manifoldname);
    // TangentSpace
    const auto tangentspace = project->tangentspaces().at(tangentspacename);

    // TensorType
    vector<int> tensorindices;
    string tensortypename;
    bool uniqueindices =
        true; // whether these tensor indices are unique within this group
    const int coordinates_group = CCTK_GroupIndex("grid::coordinates");
    if (groupindex == coordinates_group) {
      // Coordinates are a special case
      tensorindices = {};
      tensortypename = tensortypes_scalars.at(dimension);
      uniqueindices = false;
    } else {
      // TODO: use the tensor type instead
      if (groupdata.numvars == 1) {
        tensorindices = {};
        tensortypename = tensortypes_scalars.at(dimension);
      } else if (groupdata.numvars == dimension) {
        tensorindices = {varindex - varindex0};
        tensortypename = tensortypes_vectors.at(dimension);
      } else if (groupdata.numvars == dimension * (dimension + 1) / 2) {
        tensorindices = syminds.at(dimension).at(varindex0);
        tensortypename = tensortypes_symmetric_tensors.at(dimension);
      } else if (groupdata.numvars == dimension * dimension) {
        tensorindices = {groupvarindex / dimension, groupvarindex % dimension};
        tensortypename = tensortypes_tensors.at(dimension);
      } else {
        // Fallback: treat as scalars
        tensorindices = {};
        tensortypename = tensortypes_scalars.at(dimension);
        uniqueindices = false;
      }
    }
    const int tensorrank = tensorindices.size();

    // Get tensor type
    // TODO: look at "tensortypealias" tag
    auto tensortype = project->tensortypes().at(tensortypename);
    assert(tensortype->rank() == tensorrank);
    shared_ptr<TensorComponent> tensorcomponent;
    for (const auto &tc : tensortype->tensorcomponents()) {
      if (tc.second->indexvalues() == tensorindices) {
        tensorcomponent = tc.second;
        break;
      }
    }
    assert(tensorcomponent);

    // Get field
    const string fieldname = uniqueindices ? groupname : fullname;
    if (not project->fields().count(fieldname)) {
      CCTK_VINFO("Variable \"%s\" does not exist in input file",
                 varname.c_str());
      continue;
    }
    auto field = project->fields().at(fieldname);

    const int iteration = cctkGH->cctk_iteration;

    const int min_mapindex = 0;
    const int max_mapindex = Carpet::arrdata.at(groupindex).size();
    for (int mapindex = min_mapindex; mapindex < max_mapindex; ++mapindex) {
      const auto hh = Carpet::arrdata.at(groupindex).at(mapindex).hh;
      const auto dd = Carpet::arrdata.at(groupindex).at(mapindex).dd;

      const int mglevel = Carpet::mglevel;
      assert(mglevel == 0);

      const int min_reflevel = reflevel >= 0 ? reflevel : 0;
      const int max_reflevel = reflevel >= 0 ? reflevel + 1 : hh->reflevels();
      for (int reflevel = min_reflevel; reflevel < max_reflevel; ++reflevel) {

        const int min_timelevel = timelevel >= 0 ? timelevel : 0;
        const int max_timelevel = timelevel >= 0
                                      ? timelevel + 1
                                      : Carpet::groupdata.at(groupindex)
                                            .activetimelevels.at(mglevel)
                                            .at(reflevel);
        for (int timelevel = min_timelevel; timelevel < max_timelevel;
             ++timelevel) {

          // Get configuration
          string value_iteration_name =
              stringify(parameter_iteration->name(), ".", setfill('0'),
                        setw(width_it), iteration);
          string value_timelevel_name =
              stringify(parameter_timelevel->name(), ".", setfill('0'),
                        setw(width_tl), timelevel);
          auto configurationname =
              value_iteration_name + "-" + value_timelevel_name;
          if (not project->configurations().count(configurationname)) {
            CCTK_VINFO(
                "Iteration %d, timelevel %d does not exist in input file",
                iteration, timelevel);
            continue;
          }
          auto configuration = project->configurations().at(configurationname);

          // Get discretization
          const auto create_discretizationname = [&](int rl) {
            ostringstream buf;
            if (groupdata.grouptype != CCTK_GF)
              buf << groupname << "-";
            buf << configuration->name();
            if (max_mapindex > 1)
              buf << "-map." << setfill('0') << setw(width_m) << mapindex;
            if (max_reflevel > 1)
              buf << "-level." << setfill('0') << setw(width_rl) << rl;
            return buf.str();
          };
          auto discretizationname = create_discretizationname(reflevel);
          if (not manifold->discretizations().count(discretizationname)) {
            CCTK_VERROR("Discretization for group \"%s\", map %d, reflevel %d "
                        "does not exist in input file",
                        groupname.c_str(), mapindex, reflevel);
          }
          auto discretization =
              manifold->discretizations().at(discretizationname);

          // Get discrete field
          string discretefieldname =
              stringify(fieldname, "-", discretization->name());
          if (not field->discretefields().count(discretefieldname)) {
            CCTK_VINFO("Discrete Field for variable \"%s\", map %d, reflevel "
                       "%d does not exist in input file",
                       varname.c_str(), mapindex, reflevel);
            continue;
          }
          auto discretefield = field->discretefields().at(discretefieldname);

          // TODO: Keep track of which part of which component has been read

          // Loop over all discrete field blocks in the file
          // TODO: Split this over processes?
          for (const auto &discretefieldblock_pair :
               discretefield->discretefieldblocks()) {
            const auto &discretefieldblock = discretefieldblock_pair.second;

            auto discretizationblock =
                discretefieldblock->discretizationblock();
            auto active = discretizationblock->active();
            // auto will_read = active & need_read;
            // if (will_read.empty())
            //   continue;

            // Get discrete field block component
            // TODO: Access via storage index instead of name
            if (not discretefieldblock->discretefieldblockcomponents().count(
                    tensortype->name()))
              continue;
            auto discretefieldblockcomponent =
                discretefieldblock->discretefieldblockcomponents().at(
                    tensorcomponent->name());

            // Get DataSet
            // TODO: Read this lazily (needs changes in
            // DiscreteFieldBlockComponent)
            auto dataset = discretefieldblockcomponent->copyobj();

            // Read data
            const int min_local_component = 0;
            const int max_local_component = hh->local_components(reflevel);
            for (int local_component = min_local_component;
                 local_component < max_local_component; ++local_component) {
              const auto &local_box =
                  dd->local_boxes.at(mglevel).at(reflevel).at(local_component);
              auto read_box = (active & bboxset2dregion<long long>(
                                            local_box.active, dimension))
                                  .bounding_box();
              if (not read_box.empty()) {
                auto ggf = Carpet::arrdata.at(groupindex)
                               .at(mapindex)
                               .data.at(groupvarindex);
                assert(ggf);
                auto gdata = ggf->data_pointer(timelevel, reflevel,
                                               local_component, mglevel);
                assert(gdata);
                assert(gdata->has_storage());
                auto memdata = gdata->storage();
                H5::DataType memtype = cactustype2hdf5type(groupdata.vartype);
                auto membox = bbox2dbox<long long>(gdata->extent(), dimension);
                assert(read_box <= membox);
                auto memshape = dbox<long long>(
                    membox.lower(),
                    vect2dpoint<long long>(gdata->padded_shape(), dimension));

                // TODO: check dataset's datatype in file
                dataset->readData(memdata, memtype, memshape, read_box);
              } // not read_box.empty()
            }   // local_component

          } // discretefieldblock

        } // timelevel
      }   // reflevel
    }     // mapindex
  }       // varindex
}
}
