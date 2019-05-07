#include <iostream>

#include "output.hpp"
#include "util.hpp"

#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>

#include <HighResTimer.hh>
#include <Timer.hh>
#include <carpet.hh>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Version.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <SimulationIO/H5Helpers.hpp>
#include <SimulationIO/RegionCalculus.hpp>
#include <SimulationIO/SimulationIO.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace CarpetSimulationIO {
using namespace SimulationIO;
using namespace std;

void add_parametervalue(shared_ptr<Configuration> &configuration,
                        const string &name, double value) {
  if (!configuration->project()->parameters().count(name))
    configuration->project()->createParameter(name);
  auto parameter = configuration->project()->parameters().at(name);
  auto checksum = adler32(
      static_cast<const unsigned char *>(static_cast<const void *>(&value)),
      sizeof value);
  auto parametervaluename = stringify(name, "-", hex, setfill('0'),
                                      setw(2 * sizeof(checksum)), checksum);
  if (!parameter->parametervalues().count(parametervaluename)) {
    auto parametervalue = parameter->createParameterValue(parametervaluename);
    parametervalue->setValue(value);
  }
  auto parametervalue = parameter->parametervalues().at(parametervaluename);
  assert(parametervalue->value_type == ParameterValue::type_double);
  // assert(parametervalue->value_double == value);
  static_assert(sizeof parametervalue->value_double == sizeof value, "");
  assert(memcmp(&parametervalue->value_double, &value, sizeof value) == 0);
  configuration->insertParameterValue(parametervalue);
}

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
                             const string &projectname,
                             file_format output_format, file_type output_type,
                             bool do_checkpoint, int myioproc, int ioproc_every)
    : cctkGH(cctkGH), io_dir(io_dir), projectname(projectname),
      output_format(output_format), output_type(output_type),
      do_checkpoint(do_checkpoint), iteration(cctkGH->cctk_iteration),
      myioproc(myioproc), ioproc_every(ioproc_every) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Creating project \"%s\"", projectname.c_str());

  static Timers::Timer timer_local_hdf5(
      "SimulationIO::output_file_t[local,hdf5]");
  static Timers::Timer timer_local_asdf(
      "SimulationIO::output_file_t[local,asdf]");
  static Timers::Timer timer_global_hdf5(
      "SimulationIO::output_file_t[global,hdf5]");
  static Timers::Timer timer_global_asdf(
      "SimulationIO::output_file_t[global,asdf]");
  Timers::Timer &timer =
      output_type == file_type::local
          ? (output_format == file_format::hdf5 ? timer_local_hdf5
                                                : timer_local_asdf)
          : (output_format == file_format::hdf5 ? timer_global_hdf5
                                                : timer_local_asdf);
  timer.start();
  assert(do_checkpoint == (io_dir == io_dir_t::checkpoint));
  const int myproc = CCTK_MyProc(cctkGH);
  assert(ioproc_every > 0);
  switch (output_type) {
  case file_type::global:
    assert(myioproc == 0); // process 0 outputs the global file
    break;
  case file_type::local:
    assert(myioproc >= 0 and myioproc <= myproc and
           myioproc % ioproc_every == 0);
    break;
  }
  project = SimulationIO::createProject(projectname);
  timer.stop();
}

output_file_t::~output_file_t() {
  assert(not project);
  assert(tasks.empty());
}

void output_file_t::insert_vars(const vector<int> &varindices,
                                const int reflevel, const int timelevel,
                                const data_handling handle_data) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_INFO("Outputting variables");

  static Timers::Timer timer_local_hdf5(
      "SimulationIO::insert_vars[local,hdf5]");
  static Timers::Timer timer_local_asdf(
      "SimulationIO::insert_vars[local,asdf]");
  static Timers::Timer timer_global_hdf5(
      "SimulationIO::insert_vars[global,hdf5]");
  static Timers::Timer timer_global_asdf(
      "SimulationIO::insert_vars[global,asdf]");
  Timers::Timer &timer =
      output_type == file_type::local
          ? (output_format == file_format::hdf5 ? timer_local_hdf5
                                                : timer_local_asdf)
          : (output_format == file_format::hdf5 ? timer_global_hdf5
                                                : timer_local_asdf);
  timer.start();

  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  // Write options
  WriteOptions write_options;
  write_options.chunk = true;
  write_options.compress = compression_level > 0;
  write_options.compression_method = WriteOptions::compression_method_t::zlib;
  write_options.compression_level = compression_level;
  write_options.shuffle = true;
  write_options.checksum = use_checksums;

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
      // Each checkpointing segment has a different Run ID. Do not
      // store the Run ID, as this prevents combining configurations
      // from different checkpointing segments.
      // if (CCTK_IsFunctionAliased("UniqueRunID"))
      //   add_parametervalue(global_configuration, "run id",
      //                      (const char *)UniqueRunID(cctkGH));
    }
  }

  // Cache some expensive values
  string parameters;
  string grid_structure;

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
    const int grouptable = CCTK_GroupTagsTableI(groupindex);
    assert(grouptable >= 0);

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
    // TODO: Turn this decoding into a separate function
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
      string tensortypealias;
      char buf[100];
      ierr =
          Util_TableGetString(grouptable, sizeof buf, buf, "tensortypealias");
      if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY)
        tensortypealias = "";
      else if (ierr < 0)
        CCTK_VERROR("Error in tensor type alias declaration for group \"%s\"",
                    groupname.c_str());
      else
        tensortypealias = buf;
      if (CCTK_EQUALS(tensortypealias.c_str(), "scalar")) {
        tensorindices = {};
        tensortypename = tensortypes_scalars.at(3);
        uniqueindices = groupdata.numvars == 1;
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "4scalar")) {
        tensorindices = {};
        tensortypename = tensortypes_scalars.at(4);
        uniqueindices = groupdata.numvars == 1;
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "u") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "d")) {
        if (groupdata.numvars != 3)
          CCTK_VERROR(
              "Group \"%s\" has tensor type alias \"%s\" which requires "
              "%d variables, but group has %d variables",
              groupname.c_str(), tensortypealias.c_str(), 3, groupdata.numvars);
        tensorindices = {groupvarindex};
        tensortypename = tensortypes_vectors.at(3);
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "4u") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "4d")) {
        if (groupdata.numvars != 4)
          CCTK_VERROR(
              "Group \"%s\" has tensor type alias \"%s\" which requires "
              "%d variables, but group has %d variables",
              groupname.c_str(), tensortypealias.c_str(), 4, groupdata.numvars);
        tensorindices = {groupvarindex};
        tensortypename = tensortypes_vectors.at(4);
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "uu_sym") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "dd_sym")) {
        if (groupdata.numvars != 6)
          CCTK_VERROR(
              "Group \"%s\" has tensor type alias \"%s\" which requires "
              "%d variables, but group has %d variables",
              groupname.c_str(), tensortypealias.c_str(), 6, groupdata.numvars);
        tensorindices = syminds.at(3).at(groupvarindex);
        tensortypename = tensortypes_symmetric_tensors.at(3);
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "4uu_sym") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "4dd_sym")) {
        if (groupdata.numvars != 10)
          CCTK_VERROR(
              "Group \"%s\" has tensor type alias \"%s\" which requires "
              "%d variables, but group has %d variables",
              groupname.c_str(), tensortypealias.c_str(), 10,
              groupdata.numvars);
        tensorindices = syminds.at(4).at(groupvarindex);
        tensortypename = tensortypes_symmetric_tensors.at(4);
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "uu") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "ud") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "du") ||
                 CCTK_EQUALS(tensortypealias.c_str(), "dd")) {
        if (groupdata.numvars != 9)
          CCTK_VERROR(
              "Group \"%s\" has tensor type alias \"%s\" which requires "
              "%d variables, but group has %d variables",
              groupname.c_str(), tensortypealias.c_str(), 6, groupdata.numvars);
        tensorindices = {groupvarindex / 3, groupvarindex % 3};
        tensortypename = tensortypes_tensors.at(3);
      } else if (CCTK_EQUALS(tensortypealias.c_str(), "dd_sym_d")) {
        if (groupdata.numvars != 18)
          CCTK_VERROR(
              "Group \"%s\" has tensor type alias \"%s\" which requires "
              "%d variables, but group has %d variables",
              groupname.c_str(), tensortypealias.c_str(), 18,
              groupdata.numvars);
        tensorindices = syminds_rank3_symmetric12.at(3).at(groupvarindex);
        tensortypename = tensortypes_tensors_rank3_symmetric12.at(3);
        if (not project->tensortypes().count(tensortypename)) {
          // Create tensor type
          assert(dimension == 3);
          int rank = 3;
          const auto &tensortype =
              project->createTensorType(tensortypename, dimension, rank);
          const auto &tensor_indexvalues =
              syminds_rank3_symmetric12.at(dimension);
          for (int storage_index = 0;
               storage_index < int(tensor_indexvalues.size());
               ++storage_index) {
            const auto &indexvalues = tensor_indexvalues.at(storage_index);
            assert(int(indexvalues.size()) == rank);
            ostringstream buf;
            for (int r = 0; r < rank; ++r)
              buf << indexvalues.at(r);
            auto name = buf.str();
            tensortype->createTensorComponent(name, storage_index, indexvalues);
          }
        }
      } else if (tensortypealias == "") {
        // Use number of variables in group as fallback
        if (groupdata.numvars == 1) {
          tensorindices = {};
          tensortypename = tensortypes_scalars.at(dimension);
        } else if (groupdata.numvars == dimension) {
          tensorindices = {groupvarindex};
          tensortypename = tensortypes_vectors.at(dimension);
        } else if (groupdata.numvars == dimension * (dimension + 1) / 2) {
          tensorindices = syminds.at(dimension).at(groupvarindex);
          tensortypename = tensortypes_symmetric_tensors.at(dimension);
        } else if (groupdata.numvars == dimension * dimension) {
          tensorindices = {groupvarindex / dimension,
                           groupvarindex % dimension};
          tensortypename = tensortypes_tensors.at(dimension);
        } else {
          // Fallback: treat as scalars
          tensorindices = {};
          tensortypename = tensortypes_scalars.at(dimension);
          uniqueindices = false;
        }
      } else {
        CCTK_VERROR("Group \"%s\"' has unknown tensor type alias \"%s\"",
                    groupname.c_str(), tensortypealias.c_str());
      }
    }
    const int tensorrank = tensorindices.size();

    // Get tensor type
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

        // Create efficient index for reverse lookups of component-to-process
        // mapping (should this be put into Carpet?)
        vector<vector<int> > proc2components(nprocs);
        if (groupdata.disttype == CCTK_DISTRIB_CONSTANT) {
          // For distrib=constant groups, output only component 0
          int c = 0;
          int p = hh->processor(reflevel, c);
          assert(p == 0);
          proc2components.at(p).push_back(c);
        } else {
          assert(groupdata.grouptype != CCTK_SCALAR);
          for (int component = 0; component < hh->components(reflevel);
               ++component)
            proc2components.at(hh->processor(reflevel, component))
                .push_back(component);
        }

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
            // Time
            double time = cctkGH->cctk_time;
            add_parametervalue(configuration, "cctk_time", time);
            double delta_time = cctkGH->cctk_delta_time;
            add_parametervalue(configuration, "cctk_delta_time", delta_time);
            // Cactus parameters
            if (parameters.empty())
              parameters =
                  charptr2string(IOUtil_GetAllParameters(cctkGH, 1 /*all*/));
            add_parametervalue(configuration, "All Parameters", parameters);
            // Grid structure
            if (grid_structure.empty())
              grid_structure = serialise_grid_structure(cctkGH);
            add_parametervalue(configuration, "Grid Structure v7",
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
            if (groupdata.grouptype == CCTK_GF &&
                GetMaxRefinementLevels(nullptr) > 1)
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
            if (not project->coordinatesystems().count(coordinatesystemname))
              project->createCoordinateSystem(coordinatesystemname,
                                              global_configuration, manifold);
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
              if (not project->fields().count(fieldname))
                project->createField(fieldname, global_configuration, manifold,
                                     tangentspace, tensortype);
              auto field = project->fields().at(fieldname);
              if (not coordinatesystem->coordinatefields().count(
                      coordinatefieldname))
                coordinatesystem->createCoordinateField(coordinatefieldname,
                                                        direction, field);
              auto discretefieldname =
                  stringify(fieldname, "-", discretization->name());
              if (not field->discretefields().count(discretefieldname))
                field->createDiscreteField(discretefieldname,
                                           global_configuration, discretization,
                                           basis);
              auto discretefield =
                  field->discretefields().at(discretefieldname);

              tasks.push_back([=](int component,
                                  const shared_ptr<DiscretizationBlock>
                                      &discretizationblock) {
                auto discretefieldblockname = discretizationblock->name();
                if (not discretefield->discretefieldblocks().count(
                        discretefieldblockname))
                  discretefield->createDiscreteFieldBlock(
                      discretizationblock->name(), discretizationblock);
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
                  const auto &box = discretizationblock->getBox();
                  double data_origin =
                      origin[direction] +
                      box.lower()[direction] * delta[direction];
                  vector<double> data_delta(manifold->dimension(), 0.0);
                  data_delta.at(direction) = delta[direction];
                  discretefieldblockcomponent->createDataRange(
                      write_options, data_origin, data_delta);
                }
              });
            } // direction
            create_coordinate_discretefieldblockcomponent =
                [=](int component, const shared_ptr<DiscretizationBlock>
                                       &discretizationblock) {
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

          int min_dataproc, max_dataproc;
          switch (output_type) {
          case file_type::global:
            min_dataproc = 0;
            max_dataproc = nprocs;
            break;
          case file_type::local:
            min_dataproc = myproc == myioproc ? myioproc : myproc;
            max_dataproc = myproc == myioproc
                               ? min(myioproc + ioproc_every, nprocs)
                               : myproc + 1;
            break;
          default:
            assert(0);
          }
          for (int dataproc = min_dataproc; dataproc < max_dataproc;
               ++dataproc) {

            for (int component : proc2components[dataproc]) {

              // Get discretization block
              string discretizationblockname =
                  groupdata.disttype == CCTK_DISTRIB_DEFAULT
                      ? stringify("c.", setfill('0'), setw(width_c), component)
                      : "replicated";
              if (not discretization->discretizationblocks().count(
                      discretizationblockname)) {
                auto discretizationblock =
                    discretization->createDiscretizationBlock(
                        discretizationblockname);
                const auto &level_box =
                    dd->level_boxes.at(mglevel).at(reflevel);
                const auto &light_box =
                    dd->light_boxes.at(mglevel).at(reflevel).at(component);
                // TODO: We could cut off inter-process ghost points that can be
                // synchronized. However, we can not cut off prolongation ghost
                // or buffer points on past time levels, since these cannot be
                // recalculated.
                const auto &overall =
                    groupdata.grouptype == CCTK_GF and not output_ghost_zones
                        ? light_box.interior
                        : light_box.exterior;
                ibbox box;
                if (groupdata.grouptype == CCTK_GF and
                    not(do_checkpoint or output_symmetry_zones)) {
                  const int size = 2 * dimension;
                  vector<CCTK_INT> symbnd(size);
                  ierr = GetSymmetryBoundaries(cctkGH, size, symbnd.data());
                  assert(not ierr);
                  const auto &owned = light_box.owned;
                  ivect lo, hi;
                  for (int d = 0; d < Carpet::dim; ++d) {
                    lo[d] =
                        symbnd[2 * d] ? owned.lower()[d] : overall.lower()[d];
                    hi[d] = symbnd[2 * d + 1] ? owned.upper()[d]
                                              : overall.upper()[d];
                  }
                  box = ibbox(lo, hi, overall.stride());
                  assert(box >= owned && box <= overall);
                } else {
                  box = overall;
                }
                discretizationblock->setBox(bbox2box_t(box, dimension));
                auto outer_boundaries = light_box.interior - light_box.owned;
                const ibset active =
                    ((level_box.active & light_box.owned) | outer_boundaries) &
                    box;
                discretizationblock->setActive(
                    bboxset2region_t(active, dimension));
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

              switch (output_type) {
              case file_type::global: {

                int proc = hh->processor(reflevel, component);
                int other_ioproc = proc;

                // Create external links
                switch (output_format) {
                case file_format::hdf5: {
                  // Write HDF5 data
                  if (verbose)
                    CCTK_VINFO("Creating external link \"%s\"",
                               discretefieldname.c_str());
                  auto otherfilename = generate_filename(
                      io_dir_t::none, projectname, "", iteration, output_format,
                      file_type::local, other_ioproc, ioproc_every);
                  discretefieldblockcomponent->createExtLink(
                      write_options, otherfilename,
                      stringify(discretefieldblockcomponent->getPath(), "/",
                                discretefieldblockcomponent->getName()));
                  break;
                }

                case file_format::asdf: {
                  // Write ASDF data
                  if (verbose)
                    CCTK_VINFO("Creating external link \"%s\"",
                               discretefieldname.c_str());
                  auto otherfilename = generate_filename(
                      io_dir_t::none, projectname, "", iteration, output_format,
                      file_type::local, other_ioproc, ioproc_every);
                  const auto concat{
                      [](vector<string> xs, const vector<string> &ys) {
                        auto rs(move(xs));
                        for (const auto &y : ys)
                          rs.push_back(y);
                        return rs;
                      }};
                  discretefieldblockcomponent->createASDFRef(
                      write_options, otherfilename,
                      concat(discretefieldblockcomponent->yaml_path(),
                             {discretefieldblockcomponent->getName()}));
                  break;
                }

                default:
                  assert(0);
                }

                break;
              } // file_type::global

              case file_type::local: {

                switch (output_format) {
                case file_format::hdf5: {
                  // Write HDF5 data
                  if (verbose)
                    CCTK_VINFO("Creating dataset \"%s\"",
                               discretefieldname.c_str());
                  int cactustype = groupdata.vartype;
                  H5::DataType type = cactustype2hdf5type(cactustype);
                  auto dataset = discretefieldblockcomponent->createDataSet(
                      write_options, type);

                  auto task = [=]() {
                    auto memtype = dataset->datatype(); // by construction
                    auto membox0 = bbox2box_t(dd->light_boxes.at(mglevel)
                                                  .at(reflevel)
                                                  .at(component)
                                                  .exterior,
                                              dimension);
                    auto membox = discretizationblock->getBox();
                    assert(membox <= membox0);

                    if (myproc == dataproc) {
                      // We have the data, we need to either write or send them
                      const int local_component =
                          hh->get_local_component(reflevel, component);
                      auto ggf = Carpet::arrdata.at(groupindex)
                                     .at(mapindex)
                                     .data.at(groupvarindex);
                      assert(ggf);
                      auto gdata = ggf->data_pointer(timelevel, reflevel,
                                                     local_component, mglevel);
                      assert(gdata);
                      auto memlayout = box_t(
                          membox0.lower(),
                          membox0.lower() +
                              vect2point_t(gdata->padded_shape(), dimension));
                      assert(membox0 == bbox2box_t(gdata->extent(), dimension));
                      assert(gdata->has_storage());
                      auto memdata = gdata->storage();

                      if (myproc == myioproc) {
                        switch (handle_data) {
                        case data_handling::attach:
                          // Attach data
                          if (verbose)
                            CCTK_VINFO("Attaching dataset \"%s\"",
                                       discretefieldname.c_str());
                          dataset->attachData(memdata, memtype, memlayout,
                                              membox);
                          break;
                        case data_handling::write:
                          // Write data to file
                          if (verbose)
                            CCTK_VINFO("Writing dataset \"%s\"",
                                       discretefieldname.c_str());
                          dataset->writeData(memdata, memtype, memlayout,
                                             membox);
                          break;
                        default:
                          assert(0);
                        }
                      } else {
                        // Send data to I/O process
                        if (verbose)
                          CCTK_VINFO("Sending dataset \"%s\"",
                                     discretefieldname.c_str());
                        send_data(myioproc, memdata, cactustype, memlayout,
                                  membox);
                      }

                    } else if (myproc == myioproc) {
                      // We don't have the data; we need to receive and write
                      // them

                      if (verbose)
                        CCTK_VINFO("Receiving dataset \"%s\"",
                                   discretefieldname.c_str());
                      const vector<char> buf =
                          recv_data(dataproc, cactustype, membox);
                      assert(ptrdiff_t(buf.size()) ==
                             membox.size() * CCTK_VarTypeSize(cactustype));
                      switch (handle_data) {
                      case data_handling::attach:
                        if (verbose)
                          CCTK_VINFO("Attaching dataset \"%s\"",
                                     discretefieldname.c_str());
                        dataset->attachData(buf.data(), memtype, membox,
                                            membox);
                        break;
                      case data_handling::write:
                        if (verbose)
                          CCTK_VINFO("Writing dataset \"%s\"",
                                     discretefieldname.c_str());
                        dataset->writeData(buf.data(), memtype, membox, membox);
                        break;
                      default:
                        assert(0);
                      }

                    } else {
                      assert(0);
                    }
                  };

                  switch (handle_data) {
                  case data_handling::attach:
                    std::move(task)();
                    break;
                  case data_handling::write:
                    tasks.push_back(std::move(task));
                    break;
                  default:
                    assert(0);
                  }

                  break;
                }

                case file_format::asdf: {
                  // Write ASDF data

                  if (verbose)
                    CCTK_VINFO("Creating dataset \"%s\"",
                               discretefieldname.c_str());

                  const auto &box = discretizationblock->box();
                  int cactustype = groupdata.vartype;
                  auto datatype = make_shared<ASDF::datatype_t>(
                      cactustype2asdftype(cactustype));

                  if (myproc == dataproc) {
                    // We have the data, we need to either write or send them

                    auto ggf = Carpet::arrdata.at(groupindex)
                                   .at(mapindex)
                                   .data.at(groupvarindex);
                    assert(ggf);
                    int local_component =
                        hh->get_local_component(reflevel, component);
                    auto gdata = ggf->data_pointer(timelevel, reflevel,
                                                   local_component, mglevel);
                    assert(gdata);
                    const box_t memlayout(
                        box.lower(),
                        box.lower() +
                            vect2point_t(gdata->padded_shape(), dimension));
                    assert(gdata->has_storage());
                    auto memdata = gdata->storage();
                    assert(memdata);
                    size_t mempoints = gdata->size();

                    if (myproc == myioproc) {
                      // Write the data

                      if (verbose)
                        CCTK_VINFO("Writing dataset \"%s\"",
                                   discretefieldname.c_str());
                      discretefieldblockcomponent->createASDFData(
                          write_options, memdata, mempoints, memlayout,
                          datatype);

                    } else {
                      // Send the data

                      if (verbose)
                        CCTK_VINFO("Sending dataset \"%s\"",
                                   discretefieldname.c_str());
                      send_data(myioproc, memdata, cactustype, memlayout, box);
                    }

                  } else if (myproc == myioproc) {
                    // Receive and write the data

                    if (verbose)
                      CCTK_VINFO("Receiving dataset \"%s\"",
                                 discretefieldname.c_str());
                    const vector<char> buf =
                        recv_data(dataproc, cactustype, box);

                    if (verbose)
                      CCTK_VINFO("Writing dataset \"%s\"",
                                 discretefieldname.c_str());
                    discretefieldblockcomponent->createASDFData(
                        write_options,
                        ASDF::make_constant_memoized(shared_ptr<ASDF::block_t>(
                            make_shared<ASDF::typed_block_t<char> >(
                                move(buf)))),
                        datatype);
                  }

                  break;
                }

                default:
                  assert(0);
                }

                break;
              } // file_type::local

              default:
                assert(0);
              } // switch output_type

              if (bool(create_coordinate_discretefieldblockcomponent))
                create_coordinate_discretefieldblockcomponent(
                    component, discretizationblock);

            } // component
          }   // dataproc

        } // timelevel
      }   // reflevel
    }     // mapindex
  }       // varindex
  timer.stop();
}

void output_file_t::write_hdf5() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Writing HDF5 output...");

  static HighResTimer::HighResTimer timer("SimulationIO::output_hdf5");
  auto timer_clock = timer.start();

  const int myproc = CCTK_MyProc(nullptr);
  const string filename =
      generate_filename(io_dir, projectname, "", iteration, output_format,
                        output_type, myioproc, ioproc_every);
  const string tmpname =
      generate_filename(io_dir, projectname, ".tmp", iteration, output_format,
                        output_type, myioproc, ioproc_every);

  if (myproc == myioproc) {
    if (verbose)
      CCTK_VINFO("Writing project \"%s\"...", filename.c_str());
    static HighResTimer::HighResTimer timer1("SimulationIO::write_hdf5");
    auto timer1_clock = timer1.start();
    const int mode = 0777;
    CCTK_CreateDirectory(mode, basename(filename).c_str());
    project->writeHDF5(tmpname);
    timer1_clock.stop(0);
  }

  if (verbose)
    CCTK_VINFO("Writing data for project \"%s\"...", filename.c_str());
  static HighResTimer::HighResTimer timer2("SimulationIO::write_hdf5_data");
  auto timer2_clock = timer2.start();
  for (auto &task : tasks)
    move(task)();
  tasks.clear();
  timer2_clock.stop(0);

  if (output_type == file_type::local) {
    if (verbose)
      CCTK_INFO("Waiting for other processes...");
    CCTK_Barrier(nullptr);
  }

  if (myproc == myioproc) {
    if (verbose)
      CCTK_VINFO("Finalising output \"%s\"...", filename.c_str());
    int ierr = rename(tmpname.c_str(), filename.c_str());
    if (ierr)
      CCTK_VERROR("Could not rename output file \"%s\"", filename.c_str());
  }

  if (verbose)
    CCTK_INFO("Cleaning up...");
  assert(tasks.empty());
  project.reset();

  timer_clock.stop(0);

  if (verbose)
    CCTK_VINFO("Done.");
}

void output_file_t::write_asdf() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Writing ASDF output...");

  static HighResTimer::HighResTimer timer("SimulationIO::output_asdf");
  auto timer_clock = timer.start();

  const int myproc = CCTK_MyProc(nullptr);
  const string filename =
      generate_filename(io_dir, projectname, "", iteration, output_format,
                        output_type, myioproc, ioproc_every);
  const string tmpname =
      generate_filename(io_dir, projectname, ".tmp", iteration, output_format,
                        output_type, myioproc, ioproc_every);

  if (myproc == myioproc) {
    if (verbose)
      CCTK_VINFO("Writing ASDF document \"%s\"...", filename.c_str());
    static HighResTimer::HighResTimer timer1("SimulationIO::write_asdf");
    auto timer1_clock = timer1.start();
    const int mode = 0777;
    CCTK_CreateDirectory(mode, basename(filename).c_str());
    project->writeASDF(tmpname);
    timer1_clock.stop(0);
  }

  if (output_type == file_type::local) {
    if (verbose)
      CCTK_INFO("Waiting for other processes");
    CCTK_Barrier(nullptr);
  }

  if (myproc == myioproc) {
    if (verbose)
      CCTK_VINFO("Finalising output \"%s\"...", filename.c_str());
    int ierr = rename(tmpname.c_str(), filename.c_str());
    if (ierr)
      CCTK_VERROR("Could not rename output file \"%s\"", filename.c_str());
  }

  if (verbose)
    CCTK_INFO("Cleaning up...");
  assert(tasks.empty());
  project.reset();

  timer_clock.stop(0);

  if (verbose)
    CCTK_INFO("Done.");
}

void output_file_t::write() {
  // static Timers::Timer timer_local_hdf5("SimulationIO::write[local,hdf5]");
  // static Timers::Timer timer_local_asdf("SimulationIO::write[local,asdf]");
  // static Timers::Timer timer_global_hdf5("SimulationIO::write[global,hdf5]");
  // static Timers::Timer timer_global_asdf("SimulationIO::write[global,asdf]");
  // Timers::Timer &timer =
  //     output_type == file_type::local
  //         ? (output_format == file_format::hdf5 ? timer_local_hdf5
  //                                               : timer_local_asdf)
  //         : (output_format == file_format::hdf5 ? timer_global_hdf5
  //                                               : timer_local_asdf);
  // timer.start();

  switch (output_format) {
  case file_format::hdf5:
    return write_hdf5();
  case file_format::asdf:
    return write_asdf();
  default:
    assert(0);
  }

  // timer.stop();
}

} // namespace CarpetSimulationIO
