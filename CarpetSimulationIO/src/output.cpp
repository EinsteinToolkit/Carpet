#include "output.hpp"
#include "util.hpp"

#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>

#include <Timer.hh>
#include <carpet.hh>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Version.h>

#include <H5Helpers.hpp>
#include <RegionCalculus.hpp>
#include <SimulationIO.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <iostream>

namespace CarpetSimulationIO {
using namespace SimulationIO;
using namespace std;

void add_parametervalue(shared_ptr<Configuration> &configuration,
                        const string &name, const string &value) {
  if (!configuration->project()->parameters().count(name))
    configuration->project()->createParameter(name);
  auto parameter = configuration->project()->parameters().at(name);
  // ostringstream buf;
  // buf << name << "-";
  // for (char ch : value) {
  //   if (ch == '\n')
  //     break;
  //   buf << ch;
  // }
  // auto parametervaluename = buf.str();
  auto checksum = adler32(value);
  auto parametervaluename = stringify(name, "-", hex, setfill('0'),
                                      setw(2 * sizeof(checksum)), checksum);
  if (!parameter->parametervalues().count(parametervaluename)) {
    auto parametervalue = parameter->createParameterValue(parametervaluename);
    parametervalue->setValue(value);
  }
  auto parametervalue = parameter->parametervalues().at(parametervaluename);
  assert(parametervalue->value_type == ParameterValue::type_string);
  assert(parametervalue->value_string == value);
  configuration->insertParameterValue(parametervalue);
}

////////////////////////////////////////////////////////////////////////////////

output_file_t::output_file_t(const cGH *cctkGH, io_dir_t io_dir,
                             const string &projectname, int ioproc,
                             int nioprocs)
    : cctkGH(cctkGH), io_dir(io_dir), projectname(projectname), ioproc(ioproc),
      nioprocs(nioprocs) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Creating project \"%s\"", projectname.c_str());
  project = SimulationIO::createProject(projectname);
}

output_file_t::~output_file_t() {
  assert(not project);
  assert(tasks.empty());
}

void output_file_t::insert_vars(const vector<int> &varindices, int reflevel,
                                int timelevel, bool global) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_INFO("Outputting variables");

  // Parameters
  if (not project->parameters().count("iteration"))
    project->createParameter("iteration");
  auto parameter_iteration = project->parameters().at("iteration");
  if (not project->parameters().count("timelevel"))
    project->createParameter("timelevel");
  auto parameter_timelevel = project->parameters().at("timelevel");
  // Configuration
  if (not project->configurations().count("global")) {
    auto global_configuration = project->createConfiguration("global");
    // Add information to this configuration
    // General information
    add_parametervalue(global_configuration, "Cactus version",
                       CCTK_FullVersion());
    // Unique identifiers
    if (CCTK_IsFunctionAliased("UniqueConfigID"))
      add_parametervalue(global_configuration, "config id",
                         (const char *)UniqueConfigID(cctkGH));
    if (CCTK_IsFunctionAliased("UniqueBuildID"))
      add_parametervalue(global_configuration, "build id",
                         (const char *)UniqueBuildID(cctkGH));
    // Simulation ID and Run ID are not unique for a run; they differ between
    // different processes of the same run. We thus store only the IDs of
    // process 0.
    if (CCTK_MyProc(cctkGH) == 0) {
      if (CCTK_IsFunctionAliased("UniqueSimulationID"))
        add_parametervalue(global_configuration, "simulation id",
                           (const char *)UniqueSimulationID(cctkGH));
      if (CCTK_IsFunctionAliased("UniqueRunID"))
        add_parametervalue(global_configuration, "run id",
                           (const char *)UniqueRunID(cctkGH));
    }
  }
  auto global_configuration = project->configurations().at("global");
  // TensorTypes
  if (not project->tensortypes().count("Scalar0D"))
    project->createStandardTensorTypes();

  for (const int varindex : varindices) {
    if (verbose)
      CCTK_VINFO("Outputting variable \"%s\"",
                 charptr2string(CCTK_FullName(varindex)).c_str());
    const int groupindex = CCTK_GroupIndexFromVarI(varindex);
    assert(groupindex >= 0);
    const int varindex0 = CCTK_FirstVarIndexI(groupindex);
    assert(varindex0 >= 0);
    const int groupvarindex = varindex - varindex0;
    cGroup groupdata;
    int ierr = CCTK_GroupData(groupindex, &groupdata);
    assert(not ierr);
    cGroupDynamicData dynamicdata;
    ierr = CCTK_GroupDynamicData(cctkGH, groupindex, &dynamicdata);
    assert(not ierr);

    // const string varname(charptr2string(CCTK_VarName(varindex)));
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
    if (not project->manifolds().count(manifoldname))
      project->createManifold(manifoldname, global_configuration, dimension);
    const auto manifold = project->manifolds().at(manifoldname);
    // TangentSpace
    if (not project->tangentspaces().count(tangentspacename))
      project->createTangentSpace(tangentspacename, global_configuration,
                                  dimension);
    const auto tangentspace = project->tangentspaces().at(tangentspacename);
    // Basis for TangentSpace
    if (not tangentspace->bases().count(basisname)) {
      auto basis = tangentspace->createBasis(basisname, global_configuration);
      for (int d = 0; d < dimension; ++d)
        basis->createBasisVector(dirnames.at(d), d);
    }
    const auto basis = tangentspace->bases().at(basisname);

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
        tensorindices = syminds.at(dimension).at(groupvarindex);
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
    if (not project->fields().count(fieldname))
      project->createField(fieldname, global_configuration, manifold,
                           tangentspace, tensortype);
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
            auto configuration =
                project->createConfiguration(configurationname);
            if (not parameter_iteration->parametervalues().count(
                    value_iteration_name)) {
              auto value_iteration = parameter_iteration->createParameterValue(
                  value_iteration_name);
              value_iteration->setValue(iteration);
            }
            auto value_iteration =
                parameter_iteration->parametervalues().at(value_iteration_name);
            configuration->insertParameterValue(value_iteration);
            if (not parameter_timelevel->parametervalues().count(
                    value_timelevel_name)) {
              auto value_timelevel = parameter_timelevel->createParameterValue(
                  value_timelevel_name);
              value_timelevel->setValue(timelevel);
            }
            auto value_timelevel =
                parameter_timelevel->parametervalues().at(value_timelevel_name);
            configuration->insertParameterValue(value_timelevel);
            // Cactus parameters
            auto parameters =
                charptr2string(IOUtil_GetAllParameters(cctkGH, 1 /*all*/));
            add_parametervalue(configuration, "All Parameters", parameters);
            // Grid structure
            auto grid_structure = serialize_grid_structure(cctkGH);
            add_parametervalue(configuration, "Grid Structure v6",
                               grid_structure);
          }
          auto configuration = project->configurations().at(configurationname);

          // Get global coordinates
          if (groupindex == coordinates_group and groupvarindex < dimension) {
            string coordinatesystemname =
                stringify(field->name(), "-", configuration->name());
            if (not project->coordinatesystems().count(coordinatesystemname))
              project->createCoordinateSystem(coordinatesystemname,
                                              configuration, manifold);
            auto coordinatesystem =
                project->coordinatesystems().at(coordinatesystemname);
            assert(tensortype->rank() == 0);
            // int direction = tensorcomponent->indexvalues().at(0);
            int direction = groupvarindex;
            string coordinatefieldname =
                stringify(coordinatesystem->name(), "-", direction);
            if (not coordinatesystem->directions().count(direction))
              coordinatesystem->createCoordinateField(coordinatefieldname,
                                                      direction, field);
          }

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
          for (int rl = 0; rl <= reflevel; ++rl) {
            auto discretizationname = create_discretizationname(rl);
            if (not manifold->discretizations().count(discretizationname))
              manifold->createDiscretization(discretizationname, configuration);
          }
          auto discretizationname = create_discretizationname(reflevel);
          auto discretization =
              manifold->discretizations().at(discretizationname);

          // Create local coordinates
          // shared_ptr<CoordinateSystem> coordinatesystem;
          function<void(int, const shared_ptr<DiscretizationBlock> &)>
              create_coordinate_discretefieldblockcomponent;
          // TODO: create coordinates for all group types?
          if (groupdata.grouptype == CCTK_GF) {
            string coordinatesystemname = [&] {
              ostringstream buf;
              buf << "cctkGH.space";
              if (max_mapindex > 1)
                buf << "-map." << setfill('0') << setw(width_m) << mapindex;
              return buf.str();
            }();
            if (not project->coordinatesystems().count(coordinatesystemname)) {
              project->createCoordinateSystem(coordinatesystemname,
                                              global_configuration, manifold);
            }
            auto coordinatesystem =
                project->coordinatesystems().at(coordinatesystemname);
            auto tensortype =
                project->tensortypes().at(tensortypes_scalars.at(dimension));
            assert(tensortype->tensorcomponents().size() == 1);
            auto tensorcomponent =
                tensortype->tensorcomponents().begin()->second;
            vector<
                function<void(int, const shared_ptr<DiscretizationBlock> &)> >
                tasks;
            for (int direction = 0; direction < tangentspace->dimension();
                 ++direction) {
              string coordinatefieldname =
                  stringify(coordinatesystemname, "[", direction, "]");
              string fieldname = coordinatefieldname;
              if (not project->fields().count(fieldname)) {
                project->createField(fieldname, global_configuration, manifold,
                                     tangentspace, tensortype);
              }
              auto field = project->fields().at(fieldname);
              if (not coordinatesystem->coordinatefields().count(
                      coordinatefieldname)) {
                coordinatesystem->createCoordinateField(coordinatefieldname,
                                                        direction, field);
              }
              auto discretefieldname =
                  stringify(fieldname, "-", discretization->name());
              if (not field->discretefields().count(discretefieldname)) {
                field->createDiscreteField(discretefieldname,
                                           global_configuration, discretization,
                                           basis);
              }
              auto discretefield =
                  field->discretefields().at(discretefieldname);

              tasks.push_back([=](
                  int component,
                  const shared_ptr<DiscretizationBlock> &discretizationblock) {
                auto discretefieldblockname = discretizationblock->name();
                if (not discretefield->discretefieldblocks().count(
                        discretefieldblockname)) {
                  discretefield->createDiscreteFieldBlock(
                      discretizationblock->name(), discretizationblock);
                }
                auto discretefieldblock =
                    discretefield->discretefieldblocks().at(
                        discretefieldblockname);
                auto discretefieldblockcomponentname = tensorcomponent->name();
                if (not discretefieldblock->discretefieldblockcomponents()
                            .count(discretefieldblockcomponentname)) {
                  auto discretefieldblockcomponent =
                      discretefieldblock->createDiscreteFieldBlockComponent(
                          discretefieldblockcomponentname, tensorcomponent);
                  dpoint<double> origin(dimension), delta(dimension);
                  if (groupdata.grouptype == CCTK_GF) {
                    const auto &cctk_origin_space =
                        Carpet::origin_space.at(mapindex).at(mglevel);
                    const auto &cctk_delta_space =
                        Carpet::delta_space.at(mapindex);
                    assert(dimension <= size(cctk_origin_space));
                    assert(dimension <= size(cctk_delta_space));
                    for (int d = 0; d < dimension; ++d) {
                      origin[d] = cctk_origin_space[d];
                      delta[d] = cctk_delta_space[d] * Carpet::mglevelfact;
                    }
                  } else {
                    for (int d = 0; d < dimension; ++d) {
                      origin[d] = 0.0;
                      delta[d] = 1.0;
                    }
                  }
                  const auto &light_box =
                      dd->light_boxes.at(mglevel).at(reflevel).at(component);
                  const auto &exterior = light_box.exterior;
                  double data_origin =
                      origin[direction] +
                      exterior.lower()[direction] * delta[direction];
                  vector<double> data_delta(manifold->dimension(), 0.0);
                  data_delta.at(direction) = delta[direction];
                  discretefieldblockcomponent->createDataRange(data_origin,
                                                               data_delta);
                }
              });
            } // direction
            create_coordinate_discretefieldblockcomponent = [=](
                int component,
                const shared_ptr<DiscretizationBlock> &discretizationblock) {
              for (const auto &task : tasks)
                task(component, discretizationblock);
            };
          } // if grouptype == CCTK_GF

          // Get discrete field
          string discretefieldname =
              stringify(fieldname, "-", discretization->name());
          if (not field->discretefields().count(discretefieldname))
            field->createDiscreteField(discretefieldname, configuration,
                                       discretization, basis);
          auto discretefield = field->discretefields().at(discretefieldname);

          if (global) {

            const int min_component = 0;
            const int max_component = hh->components(reflevel);
            for (int component = min_component; component < max_component;
                 ++component) {
              int proc = hh->processor(reflevel, component);
              int other_ioproc = proc;

              // The name choices here must be consistent with the local
              // output below

              // Get discretization block
              string discretizationblockname =
                  stringify("c.", setfill('0'), setw(width_c), component);
              if (not discretization->discretizationblocks().count(
                      discretizationblockname)) {
                auto discretizationblock =
                    discretization->createDiscretizationBlock(
                        discretizationblockname);
                const auto &level_box =
                    dd->level_boxes.at(mglevel).at(reflevel);
                const auto &light_box =
                    dd->light_boxes.at(mglevel).at(reflevel).at(component);
                // const auto &local_box =
                //     dd->local_boxes.at(mglevel).at(reflevel).at(
                //         local_component);
                // TODO: Cut off ghosts?
                discretizationblock->setBox(
                    bbox2dbox<long long>(light_box.exterior, dimension));
                // cout << "box=" << discretizationblock->box() << "\n";
                discretizationblock->setActive(bboxset2dregion<long long>(
                    level_box.active & light_box.owned, dimension));
                // cout << "active=" << discretizationblock->active() << "\n";
              }
              auto discretizationblock =
                  discretization->discretizationblocks().at(
                      discretizationblockname);

              // Get discrete field block
              string discretefieldblockname = discretizationblockname;
              if (not discretefield->discretefieldblocks().count(
                      discretefieldblockname))
                discretefield->createDiscreteFieldBlock(
                    discretizationblock->name(), discretizationblock);
              auto discretefieldblock = discretefield->discretefieldblocks().at(
                  discretefieldblockname);

              // Get discrete field block data
              if (not discretefieldblock->discretefieldblockcomponents().count(
                      tensorcomponent->name()))
                discretefieldblock->createDiscreteFieldBlockComponent(
                    tensorcomponent->name(), tensorcomponent);
              auto discretefieldblockcomponent =
                  discretefieldblock->discretefieldblockcomponents().at(
                      tensorcomponent->name());

              // Create external links
              if (verbose)
                CCTK_VINFO("Creating external link \"%s\"",
                           discretefieldname.c_str());
              auto otherfilename = generate_filename(
                  cctkGH, io_dir_none, projectname, "", cctkGH->cctk_iteration,
                  other_ioproc, nioprocs);
              auto dataset = discretefieldblockcomponent->createExtLink(
                  otherfilename,
                  stringify(discretefieldblockcomponent->getPath(), "/",
                            discretefieldblockcomponent->getName()));

#if 0
                tasks.push_back([=]() {
                  // TODO: Use SimulationIO mechanism for this
                  H5::createExternalLink(
                      h5file,
                      "manifolds/" + manifold->name() + "/discretizations/" +
                          discretization->name() + "/discretizationblocks/" +
                          discretizationblockname,
                      otherfilename,
                      "manifolds/" + manifold->name() + "/discretizations/" +
                          discretization->name() + "/discretizationblocks/" +
                          discretizationblockname);
                  H5::createExternalLink(
                      h5file,
                      "fields/" + field->name() + "/discretefields/" +
                          discretefield->name() + "/discretefieldblocks/" +
                          discretefieldblockname,
                      otherfilename,
                      "fields/" + field->name() + "/discretefields/" +
                          discretefield->name() + "/discretefieldblocks/" +
                          discretefieldblockname);
                });
#endif

              if (bool(create_coordinate_discretefieldblockcomponent))
                create_coordinate_discretefieldblockcomponent(
                    component, discretizationblock);

            } // component

          } else { // not global

            const int min_local_component = 0;
            const int max_local_component = hh->local_components(reflevel);
            for (int local_component = min_local_component;
                 local_component < max_local_component; ++local_component) {
              const int component =
                  hh->get_component(reflevel, local_component);

              // Get discretization block
              string discretizationblockname =
                  stringify("c.", setfill('0'), setw(width_c), component);
              if (not discretization->discretizationblocks().count(
                      discretizationblockname)) {
                auto discretizationblock =
                    discretization->createDiscretizationBlock(
                        discretizationblockname);
                const auto &level_box =
                    dd->level_boxes.at(mglevel).at(reflevel);
                const auto &light_box =
                    dd->light_boxes.at(mglevel).at(reflevel).at(component);
                const auto &local_box =
                    dd->local_boxes.at(mglevel).at(reflevel).at(
                        local_component);
                // TODO: Cut off ghosts?
                discretizationblock->setBox(
                    bbox2dbox<long long>(light_box.exterior, dimension));
                // cout << "box=" << discretizationblock->box() << "\n";
                assert((level_box.active & light_box.owned) ==
                       local_box.active);
                discretizationblock->setActive(
                    bboxset2dregion<long long>(local_box.active, dimension));
                // cout << "active=" << discretizationblock->active() << "\n";
              }
              auto discretizationblock =
                  discretization->discretizationblocks().at(
                      discretizationblockname);

              // Get discrete field block
              string discretefieldblockname = discretizationblockname;
              if (not discretefield->discretefieldblocks().count(
                      discretefieldblockname))
                discretefield->createDiscreteFieldBlock(
                    discretizationblock->name(), discretizationblock);
              auto discretefieldblock = discretefield->discretefieldblocks().at(
                  discretefieldblockname);

              // Get discrete field block data
              if (not discretefieldblock->discretefieldblockcomponents().count(
                      tensorcomponent->name()))
                discretefieldblock->createDiscreteFieldBlockComponent(
                    tensorcomponent->name(), tensorcomponent);
              auto discretefieldblockcomponent =
                  discretefieldblock->discretefieldblockcomponents().at(
                      tensorcomponent->name());

              // Write data
              if (verbose)
                CCTK_VINFO("Creating dataset \"%s\"",
                           discretefieldname.c_str());
              H5::DataType type = cactustype2hdf5type(groupdata.vartype);
              auto dataset = discretefieldblockcomponent->createDataSet(type);
              tasks.push_back([=]() {
                if (verbose)
                  CCTK_VINFO("Writing dataset \"%s\"",
                             discretefieldname.c_str());
                auto ggf = Carpet::arrdata.at(groupindex)
                               .at(mapindex)
                               .data.at(groupvarindex);
                assert(ggf);
                auto gdata = ggf->data_pointer(timelevel, reflevel,
                                               local_component, mglevel);
                assert(gdata);
                assert(gdata->has_storage());
                auto memdata = gdata->storage();
                auto memtype = dataset->datatype(); // by construction
#if 0
                vector<hsize_t> padded_shape(dimension);
                assert(dimension <= size(gdata->padded_shape()));
                for (int d = 0; d < dimension; ++d)
                  padded_shape.at(d) = gdata->padded_shape()[d];
                vector<hsize_t> offset(dimension);
                for (int d = 0; d < dimension; ++d)
                  offset.at(d) = 0;
                vector<hsize_t> shape(dimension);
                assert(dimension <= size(gdata->shape()));
                for (int d = 0; d < dimension; ++d)
                  shape.at(d) = gdata->shape()[d];
                auto memspace =
                    H5::DataSpace(dimension, reversed(padded_shape).data());
                if (dimension > 0)
                  memspace.selectHyperslab(H5S_SELECT_SET,
                                           reversed(shape).data(),
                                           reversed(offset).data());
                dataset->dataset().write(data, memtype, memspace);
#endif
                auto membox = bbox2dbox<long long>(gdata->extent(), dimension);
                auto memshape = dbox<long long>(
                    membox.lower(),
                    membox.lower() + vect2dpoint<long long>(
                                         gdata->padded_shape(), dimension));
                dataset->writeData(memdata, memtype, memshape, membox);
              });

              if (bool(create_coordinate_discretefieldblockcomponent))
                create_coordinate_discretefieldblockcomponent(
                    component, discretizationblock);

            } // local_component
          }   // if not global
        }     // timelevel
      }       // reflevel
    }         // mapindex
  }           // varindex
}

void output_file_t::write() {
  DECLARE_CCTK_PARAMETERS;
  auto tmpname = generate_filename(cctkGH, io_dir, projectname, ".tmp",
                                   cctkGH->cctk_iteration, ioproc, nioprocs);
  const bool create_dirs = ioproc >= 0;
  auto filename =
      generate_filename(cctkGH, io_dir, projectname, "", cctkGH->cctk_iteration,
                        ioproc, nioprocs, create_dirs);
  if (verbose)
    CCTK_VINFO("Creating file \"%s\"", filename.c_str());
  auto fapl = H5::FileAccPropList();
  fapl.setFcloseDegree(H5F_CLOSE_STRONG);
  fapl.setLibverBounds(H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
  static Timers::Timer timer1("SimulationIO::H5File");
  timer1.start();
  // H5F_ACC_EXCL or H5F_ACC_TRUNC,
  auto file =
      H5::H5File(tmpname, H5F_ACC_EXCL, H5::FileCreatPropList::DEFAULT, fapl);
  timer1.stop();
  if (verbose)
    CCTK_VINFO("Writing project \"%s\"", filename.c_str());
  static Timers::Timer timer2("SimulationIO::write");
  timer2.start();
  project->write(file);
  timer2.stop();
  if (verbose)
    CCTK_VINFO("Writing data for project \"%s\"", filename.c_str());
  static Timers::Timer timer3("SimulationIO::tasks");
  timer3.start();
  for (auto &task : tasks)
    std::move(task)();
  timer3.stop();
  if (verbose)
    CCTK_VINFO("Deleting project \"%s\"", filename.c_str());
  static Timers::Timer timer4("SimulationIO::close");
  timer4.start();
  file.close();
  timer4.stop();
  static Timers::Timer timer5("SimulationIO::cleanup");
  timer5.start();
  tasks.clear();
  project.reset();
  H5garbage_collect();
  timer5.stop();
  if (verbose)
    CCTK_VINFO("Done deleting project \"%s\"", filename.c_str());
  if (ioproc >= 0) {
    CCTK_Barrier(cctkGH);
    if (verbose)
      CCTK_INFO("Done waiting for other processes");
  }
  int ierr = rename(tmpname.c_str(), filename.c_str());
  if (ierr)
    CCTK_VERROR("Could not rename output file \"%s\"", filename.c_str());
  if (verbose)
    CCTK_VINFO("Done renaming file \"%s\"", filename.c_str());
}
}
