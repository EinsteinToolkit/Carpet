#include "input.hpp"
#include "util.hpp"

#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>

#include <carpet.hh>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace CarpetSimulationIO {
namespace RC = RegionCalculus;
using namespace SimulationIO;
using namespace std;

input_file_t::input_file_t(io_dir_t io_dir, const string &projectname,
                           file_format input_format, int iteration, int ioproc,
                           int nioprocs)
    : input_format(input_format), iteration(iteration) {
  DECLARE_CCTK_PARAMETERS;

  auto filename = generate_filename(io_dir, projectname, "", iteration,
                                    input_format, file_type::global, -1, -1);
  if (verbose)
    CCTK_VINFO("Reading file \"%s\"", filename.c_str());
  switch (input_format) {
  case file_format::hdf5:
    project = readProjectHDF5(filename);
    break;
  case file_format::asdf:
    project = readProjectASDF(filename);
    break;
  default:
    assert(0);
  }
}

void input_file_t::read_params() const {
  // DECLARE_CCTK_PARAMETERS;

  auto parameter_all_params = project->parameters().at("All Parameters");

  // TODO: If there are multiple parameter values, then maybe select
  // by iteration number?
  assert(parameter_all_params->parametervalues().size() == 1);
  auto value_all_params =
      parameter_all_params->parametervalues().begin()->second;
  assert(value_all_params->value_type == ParameterValue::type_string);
  string all_params = value_all_params->value_string;

  IOUtil_SetAllParameters(all_params.c_str());
}

void input_file_t::read_grid_structure(cGH *cctkGH) const {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_INFO("Recovering grid structure");

  auto parameter_grid_structure = project->parameters().at("Grid Structure v7");

  // TODO: If there are multiple parameter values, then maybe select
  // by iteration number?
  assert(parameter_grid_structure->parametervalues().size() == 1);
  auto value_grid_structure =
      parameter_grid_structure->parametervalues().begin()->second;
  assert(value_grid_structure->value_type == ParameterValue::type_string);
  string grid_structure_string = value_grid_structure->value_string;

  grid_structure_v7 grid_structure =
      deserialise_grid_structure(grid_structure_string);

  vector<vector<vector<CarpetLib::region_t> > > superregsss(Carpet::maps);
  for (int m = 0; m < Carpet::maps; ++m)
    superregsss.at(m) = move(grid_structure.maps.at(m).superregions);
  vector<vector<vector<vector<CarpetLib::region_t> > > > regssss(Carpet::maps);

  int type;
  const void *ptr = CCTK_ParameterGet("regrid_in_level_mode", "Carpet", &type);
  assert(ptr != 0);
  assert(type == PARAMETER_BOOLEAN);
  const CCTK_INT regrid_in_level_mode = *static_cast<const CCTK_INT *>(ptr);

  if (not regrid_in_level_mode) {
    // Distribute each map independently

    for (int m = 0; m < Carpet::maps; ++m) {
      vector<vector<CarpetLib::region_t> > &superregss = superregsss.at(m);

      // Make multiprocessor aware
      vector<vector<CarpetLib::region_t> > regss(superregss.size());
      for (size_t rl = 0; rl < superregss.size(); ++rl)
        Carpet::SplitRegions(cctkGH, superregss.at(rl), regss.at(rl));

      // Make multigrid aware
      Carpet::MakeMultigridBoxes(cctkGH, m, regss, regssss.at(m));
    } // for m

  } else { // if regrid_in_level_mode
    // Distribute all maps at the same time

    vector<vector<vector<CarpetLib::region_t> > > regsss(Carpet::maps);

    // Count levels
    vector<int> rls(Carpet::maps);
    for (int m = 0; m < Carpet::maps; ++m)
      rls.at(m) = superregsss.at(m).size();
    int maxrl = 0;
    for (int m = 0; m < Carpet::maps; ++m)
      maxrl = max(maxrl, rls.at(m));
    // All maps must have the same number of levels
    for (int m = 0; m < Carpet::maps; ++m) {
      superregsss.at(m).resize(maxrl);
      regsss.at(m).resize(maxrl);
    }

    // Make multiprocessor aware
    for (int rl = 0; rl < maxrl; ++rl) {
      vector<vector<CarpetLib::region_t> > superregss(Carpet::maps);
      for (int m = 0; m < Carpet::maps; ++m)
        superregss.at(m) = move(superregsss.at(m).at(rl));
      vector<vector<CarpetLib::region_t> > regss(Carpet::maps);
      Carpet::SplitRegionsMaps(cctkGH, superregss, regss);
      for (int m = 0; m < Carpet::maps; ++m) {
        superregsss.at(m).at(rl) = move(superregss.at(m));
        regsss.at(m).at(rl) = move(regss.at(m));
      }
    } // for rl

    // Make multigrid aware
    Carpet::MakeMultigridBoxesMaps(cctkGH, regsss, regssss);

  } // if regrid_in_level_mode

  // Regrid
  for (int m = 0; m < Carpet::maps; ++m)
    Carpet::RegridMap(cctkGH, m, superregsss.at(m), regssss.at(m), false);

  // Set time hierarchy correctly after RegridMap created it
  for (int rl = 0; rl < Carpet::vhh.at(0)->reflevels(); ++rl) {
    for (int tl = 0; tl < Carpet::tt->timelevels; ++tl)
      Carpet::tt->set_time(Carpet::mglevel, rl, tl,
                           grid_structure.times.at(rl).at(tl));
    Carpet::tt->set_delta(Carpet::mglevel, rl, grid_structure.deltas.at(rl));
  }

  // We are in level mode here (we probably shouldn't be, but that's how Cactus
  // I/O was designed), so we need to fix up certain things. See Carpet's
  // function leave_level_mode to see what happens when we leave level mode.
  {
    // TODO: Switch to global mode instead of fixing up things manually
    assert(Carpet::reflevel != -1);
    const int tl = 0;
    Carpet::global_time = grid_structure.global_time;
    cctkGH->cctk_time =
        Carpet::tt->get_time(Carpet::mglevel, Carpet::reflevel, tl);
    Carpet::delta_time = grid_structure.delta_time;
    cctkGH->cctk_delta_time = Carpet::delta_time;
  }

  Carpet::PostRegrid(cctkGH);

  for (int rl = 0; rl < Carpet::reflevels; ++rl)
    Carpet::Recompose(cctkGH, rl, false);

  Carpet::RegridFree(cctkGH, false);

  // Check ghost and buffer widths and prolongation orders
  // TODO: Instead of only checking them, set them during regridding.
  // (This requires setting them during regridding instead of during setup.)
  for (int m = 0; m < Carpet::maps; ++m) {
    const auto &dd = *Carpet::vdd.at(m);
    const auto &map_structure = grid_structure.maps.at(m);
    assert(map_structure.ghost_widths.size() == dd.ghost_widths.size());
    for (size_t rl = 0; rl < map_structure.ghost_widths.size(); ++rl)
      assert(all(
          all(map_structure.ghost_widths.at(rl) == dd.ghost_widths.at(rl))));
    assert(map_structure.buffer_widths.size() == dd.buffer_widths.size());
    for (size_t rl = 0; rl < map_structure.buffer_widths.size(); ++rl)
      assert(all(
          all(map_structure.buffer_widths.at(rl) == dd.buffer_widths.at(rl))));
    assert(map_structure.prolongation_orders_space.size() ==
           dd.prolongation_orders_space.size());
    for (size_t rl = 0; rl < map_structure.prolongation_orders_space.size();
         ++rl)
      assert(map_structure.prolongation_orders_space.at(rl) ==
             dd.prolongation_orders_space.at(rl));
  }
}

void input_file_t::read_vars(const vector<int> &varindices, const int reflevel1,
                             const int timelevel1) const {
  DECLARE_CCTK_PARAMETERS;

  auto parameter_iteration = project->parameters().at("iteration");
  auto parameter_timelevel = project->parameters().at("timelevel");

  // A schedule deciding which discretization blocks should be read for which
  // local components, minimizing the number of discretization blocks that need
  // to be read for each local component.
  typedef vector<shared_ptr<DiscretizationBlock> > read_schedule_t;
  struct component_read_schedules_t {
    int mapindex;
    int reflevel;
    int proc;
    map<int, read_schedule_t> read_schedules; // [local_component]
  };
  map<shared_ptr<Discretization>, component_read_schedules_t>
      discretization_read_schedules;

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
        tensorindices = {groupvarindex};
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
    if (not project->fields().count(fieldname)) {
      CCTK_VINFO("Variable \"%s\" does not exist in input file",
                 varname.c_str());
      continue;
    }
    const auto field = project->fields().at(fieldname);

    // Manifold
    const auto manifold = field->manifold();
    assert(manifold->dimension() == dimension);

    // TangentSpace
    const auto tangentspace = field->tangentspace();
    assert(tangentspace->dimension() == dimension);

    const int min_mapindex = 0;
    const int max_mapindex = Carpet::arrdata.at(groupindex).size();
    for (int mapindex = min_mapindex; mapindex < max_mapindex; ++mapindex) {
      const auto hh = Carpet::arrdata.at(groupindex).at(mapindex).hh;
      const auto dd = Carpet::arrdata.at(groupindex).at(mapindex).dd;

      const int mglevel = Carpet::mglevel;
      assert(mglevel == 0);

      const int min_reflevel = reflevel1 >= 0 ? reflevel1 : 0;
      const int max_reflevel = reflevel1 >= 0 ? reflevel1 + 1 : hh->reflevels();
      for (int reflevel = min_reflevel; reflevel < max_reflevel; ++reflevel) {

        const int min_timelevel = timelevel1 >= 0 ? timelevel1 : 0;
        const int max_timelevel = timelevel1 >= 0
                                      ? timelevel1 + 1
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

          // TODO: traverse all discretizations, then find matching components
          // instead of the other way around. This would ensure that every
          // dataset is read at most once, so that it can be flushed from memory
          // afterwards.

          // Get discretization
          // TODO: traverse all discretizations, choose by semantics instead of
          // by name
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
          auto discretizationname = create_discretizationname(reflevel);
          if (not manifold->discretizations().count(discretizationname))
            CCTK_VERROR("Discretization for group \"%s\", map %d, reflevel %d "
                        "does not exist in input file",
                        groupname.c_str(), mapindex, reflevel);
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
          // auto discretization = discretefield->discretization();

          // Determine which discretization blocks should be read for which
          // local components, minimizing the number of discretization blocks
          // that need to be read.
          const int proc = CCTK_MyProc(nullptr);
          const int min_local_component = 0;
          const int max_local_component = hh->local_components(reflevel);
          if (not discretization_read_schedules.count(discretization)) {
            component_read_schedules_t component_read_schedules{mapindex,
                                                                reflevel, proc};
            // We need to read everything that cannot be recalculated. We can
            // skip inter-process ghost points on all time levels, and we can
            // skip prolongation ghost, buffer points, and restricted points on
            // the current time level. However, prolongated and restricted
            // points on past time levels can not be recalculated and need to be
            // read. To simplify things, we read everything.

            for (int local_component = min_local_component;
                 local_component < max_local_component; ++local_component) {
              const int component =
                  hh->get_component(reflevel, local_component);
              const auto &light_box =
                  dd->light_boxes.at(mglevel).at(reflevel).at(component);
              const auto need_read =
                  bboxset2region_t(ibset(light_box.exterior), dimension);
              // We use a greedy algorithm.
              // 1. Find overlaps for local components and discretization blocks
              struct overlap_t {
                shared_ptr<DiscretizationBlock> discretizationblock;
                RC::region_t overlap;
              };
              vector<overlap_t> overlaps;
              for (const auto &discretefieldblock_pair :
                   discretefield->discretefieldblocks()) {
                const auto &discretefieldblock = discretefieldblock_pair.second;
                auto discretizationblock =
                    discretefieldblock->discretizationblock();
                auto box = discretizationblock->box();
                if (not need_read.isdisjoint(box))
                  overlaps.emplace_back(
                      overlap_t{discretizationblock, need_read & box});
              }
              // 2. Sort overlaps by size (largest first)
              sort(overlaps.begin(), overlaps.end(),
                   [&](const overlap_t &a, const overlap_t &b) {
                     // TODO: overlap.size is not O(1)
                     return a.overlap.size() > b.overlap.size();
                   });
              // 3. Walk overlaps and create schedule
              read_schedule_t read_schedule;
              auto still_need_read = need_read;
              for (const auto &overlap : overlaps) {
                if (not still_need_read.isdisjoint(overlap.overlap)) {
                  read_schedule.push_back(overlap.discretizationblock);
                  still_need_read -= overlap.overlap;
                }
              }
              if (not still_need_read.empty()) {
                CCTK_VWARN(
                    CCTK_WARN_ALERT,
                    "Read schedule for component [m=%d,rl=%d,c=%d,p=%d,lc=%d]:",
                    mapindex, reflevel, component, proc, local_component);
                for (const auto &db : read_schedule)
                  CCTK_VWARN(CCTK_WARN_ALERT, "  %s", db->name().c_str());
                cerr << "need_read=" << need_read << "\n";
                cerr << "still_need_read=" << still_need_read << "\n";
                CCTK_ERROR("");
              }
              assert(still_need_read.empty());
              if (verbose) {
                CCTK_VINFO(
                    "Read schedule for component [m=%d,rl=%d,c=%d,p=%d,lc=%d]:",
                    mapindex, reflevel, component, proc, local_component);
                for (const auto &db : read_schedule)
                  CCTK_VINFO("  %s", db->name().c_str());
              }
              component_read_schedules.read_schedules[local_component] =
                  move(read_schedule);
            }
            discretization_read_schedules[discretization] =
                move(component_read_schedules);
          }
          const auto &component_read_schedules =
              discretization_read_schedules.at(discretization);
          assert(component_read_schedules.mapindex == mapindex);
          assert(component_read_schedules.reflevel == reflevel);
          assert(component_read_schedules.proc == proc);

          // Cache ASDF arrays
          // TODO: Do this in SimulationIO
          map<shared_ptr<DiscreteFieldBlockComponent>,
              shared_ptr<ASDF::ndarray> >
              ndarrays;

          for (int local_component = min_local_component;
               local_component < max_local_component; ++local_component) {
            const int component = hh->get_component(reflevel, local_component);
            // TODO: We could cut off inter-process ghost points that can be
            // synchronized (same as for writing).
            const auto &light_box =
                dd->light_boxes.at(mglevel).at(reflevel).at(component);
            const auto need_read =
                bboxset2region_t(ibset(light_box.exterior), dimension);
            auto still_need_read = need_read;
            // auto did_read = RC::region_t(dimension);

            // Loop over the discretization blocks that need to be read
            const auto &discretizationblocks =
                component_read_schedules.read_schedules.at(local_component);
            for (const auto &discretizationblock : discretizationblocks) {

              // Find respective discrete field block
              // TODO: Use a more efficient lookup method
              shared_ptr<DiscreteFieldBlock> discretefieldblock;
              for (const auto &discretefieldblock_pair :
                   discretefield->discretefieldblocks()) {
                const auto &dfb = discretefieldblock_pair.second;
                if (dfb->discretizationblock() == discretizationblock) {
                  discretefieldblock = dfb;
                  break;
                }
              }
              assert(discretefieldblock);

              // Get discrete field block component
              // TODO: Access via storage index instead of name
              if (not discretefieldblock->discretefieldblockcomponents().count(
                      tensorcomponent->name()))
                continue;
              auto discretefieldblockcomponent =
                  discretefieldblock->discretefieldblockcomponents().at(
                      tensorcomponent->name());

              auto box = discretizationblock->box();
              // auto active = discretizationblock->active();
              assert(not box.empty());
              auto read_box =
                  (RC::region_t(box) & still_need_read).bounding_box();
              assert(not read_box.empty());

              // Get DataSet
              // TODO: Read this lazily (needs changes in
              // DiscreteFieldBlockComponent)

              auto ggf = Carpet::arrdata.at(groupindex)
                             .at(mapindex)
                             .data.at(groupvarindex);
              assert(ggf);
              auto gdata = ggf->data_pointer(timelevel, reflevel,
                                             local_component, mglevel);
              assert(gdata);
              assert(gdata->has_storage());
              auto memdata = gdata->storage();
              assert(memdata);
              H5::DataType memtype = cactustype2hdf5type(groupdata.vartype);
              auto membytes =
                  gdata->size() * CCTK_VarTypeSize(groupdata.vartype);
              auto membox = bbox2box_t(gdata->extent(), dimension);
              assert(read_box <= membox);
              auto memlayout =
                  RC::box_t(membox.lower(),
                            membox.lower() +
                                vect2point_t(gdata->padded_shape(), dimension));

              // TODO: check dataset's datatype in file
              if (verbose)
                CCTK_VINFO("Reading dataset \"%s\"", discretefieldname.c_str());
              switch (input_format) {
              case file_format::hdf5: {
                auto dataset = discretefieldblockcomponent->copyobj();
                dataset->readData(memdata, memtype, memlayout, read_box);
                break;
              }
              case file_format::asdf: {
                // TODO: Move most of this into a new "readASDF" function
                if (not ndarrays.count(discretefieldblockcomponent)) {
                  shared_ptr<ASDF::ndarray> ndarray;
                  auto dataset_asdf = discretefieldblockcomponent->asdfdata();
                  auto dataset_ref = discretefieldblockcomponent->asdfref();
                  if (dataset_asdf) {
                    ndarray = dataset_asdf->ndarray();
                  } else if (dataset_ref) {
                    auto ref = dataset_ref->reference();
                    auto rsn = ref->resolve();
                    // TODO: Allow multiple layers of references, i.e.
                    // references to other references
                    ndarray = make_shared<ASDF::ndarray>(rsn.first, rsn.second);
                  } else {
                    assert(0);
                  }
                  ndarrays[discretefieldblockcomponent] = ndarray;
                }
                const auto &ndarray = ndarrays.at(discretefieldblockcomponent);
                const auto block = ndarray->get_data();
                const auto type_size = ndarray->get_datatype()->type_size();
                assert(int(type_size) == CCTK_VarTypeSize(groupdata.vartype));
                assert(block->nbytes() % type_size == 0);
                const auto mem_off_str =
                    HyperSlab::layout2strides(memlayout, read_box, type_size);
                const auto offset = HyperSlab::layout2offset(
                    ndarray->get_offset(), RC::point_t(ndarray->get_strides()),
                    box, read_box);
                HyperSlab::copy(memdata, membytes, mem_off_str.first,
                                mem_off_str.second, block->ptr(),
                                block->nbytes(), offset,
                                RC::point_t(ndarray->get_strides()),
                                read_box.shape(), type_size);
                break;
              }
              default:
                assert(0);
              }

              still_need_read -= read_box;
              // did_read |= read_box;

            } // for discretizationblock

            assert(still_need_read.empty());

#if 0
            bool have_error = false;
            for (const auto &kv : local_components_did_read)
              assert(kv.first >= min_local_component and
                     kv.first < max_local_component);
            for (int local_component = min_local_component;
                 local_component < max_local_component; ++local_component) {
              int component = hh->get_component(reflevel, local_component);
              const auto &need_read =
                  local_components_need_read.at(local_component);
              const auto &did_read =
                  local_components_did_read.at(local_component);
              // We allow reading more
              if (did_read < need_read) {
                have_error = true;
                ostringstream buf;
                buf << "For variable " << varindex << " \"" << varname
                    << "\", time level " << timelevel << ", map " << mapindex
                    << ", refinement level " << reflevel << ", local component "
                    << local_component << ", component " << component
                    << ": did not read all active+boundary points: read "
                    << did_read << " instead of " << need_read
                    << "; unnecessary points are " << (did_read - need_read)
                    << ", missing points are " << (need_read - did_read);
                CCTK_WARN(CCTK_WARN_ALERT, buf.str().c_str());
              }
            }
            if (have_error)
              CCTK_ERROR("Aborting");
#endif

          } // for local_component

          // Free ASDF data blocks to reduce memory usage
          for (const auto &kv : ndarrays) {
            const auto &ndarray = kv.second;
            ndarray->get_data().forget();
          }

        } // timelevel
      }   // reflevel
    }     // mapindex
  }       // varindex
}
} // namespace CarpetSimulationIO
