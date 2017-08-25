#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <F5/F5F.h>
#include <F5/F5R.h>
#include <F5/F5iterate.h>
#include <F5/F5uniform.h>
#include <hdf5.h>

#include <bbox.hh>
#include <defs.hh>
#include <vect.hh>

#include <carpet.hh>

#include "iof5.hh"

std::ostream &operator<<(std::ostream &os, CarpetIOF5::indent_t const &indent) {
  return indent.output(os);
}

namespace CarpetIOF5 {

using namespace std;
using namespace Carpet;

// Indentation
int indent_t::level = 0;

indent_t::indent_t() {
  if (debug)
    cout << *this << "{{{\n";
  ++level;
}

indent_t::~indent_t() {
  --level;
  if (debug)
    cout << *this << "}}}\n";
}

ostream &indent_t::output(ostream &os) const {
  return os << string(width * level, ' ');
}

// Generate a good file name ("alias") for a variable
string generate_basename(cGH const *const cctkGH, int const variable) {
  DECLARE_CCTK_PARAMETERS;

  ostringstream filename_buf;

  if (CCTK_EQUALS(file_content, "variable")) {
    assert(variable >= 0);
    char *const varname = CCTK_FullName(variable);
    assert(varname);
    for (char *p = varname; *p; ++p) {
      *p = tolower(*p);
    }
    filename_buf << varname;
    free(varname);
  } else if (CCTK_EQUALS(file_content, "group")) {
    assert(variable >= 0);
    char *const groupname = CCTK_GroupNameFromVarI(variable);
    assert(groupname);
    for (char *p = groupname; *p; ++p) {
      *p = tolower(*p);
    }
    filename_buf << groupname;
    free(groupname);
  } else if (CCTK_EQUALS(file_content, "thorn")) {
    assert(variable >= 0);
    char const *const impname = CCTK_ImpFromVarI(variable);
    char *const thornname = strdup(impname);
    assert(thornname);
    char *const colon = strchr(thornname, ':');
    assert(colon);
    *colon = '\0';
    for (char *p = thornname; *p; ++p) {
      *p = tolower(*p);
    }
    filename_buf << thornname;
    free(thornname);
  } else if (CCTK_EQUALS(file_content, "everything")) {
    if (CCTK_EQUALS(out_filename, "")) {
      // Obtain the parameter file name
      char buf[10000];
      int const ilen = CCTK_ParameterFilename(sizeof buf, buf);
      assert(ilen < int(sizeof buf) - 1);
      char *parfilename = buf;
      // Remove directory prefix, if any
      char *const slash = strrchr(parfilename, '/');
      if (slash) {
        parfilename = &slash[1];
      }
      // Remove suffix, if it is there
      char *const suffix = strrchr(parfilename, '.');
      if (suffix and strcmp(suffix, ".par") == 0) {
        suffix[0] = '\0';
      }
      filename_buf << parfilename;
    } else {
      filename_buf << out_filename;
    }
  } else {
    assert(0);
  }

  if (out_timesteps_per_file > 0) {
    int const iteration = (cctkGH->cctk_iteration / out_timesteps_per_file *
                           out_timesteps_per_file);
    filename_buf << ".it" << setw(iteration_digits) << setfill('0')
                 << iteration;
  }

  return filename_buf.str();
}

// Create the final file name on a particular processor
string create_filename(cGH const *const cctkGH, string const basename,
                       int const iteration, int const proc,
                       io_dir_t const io_dir, bool const create_directories) {
  DECLARE_CCTK_PARAMETERS;

  char const *IO_dir = NULL;
  char const *F5_dir = NULL;
  switch (io_dir) {
  case io_dir_input:
    IO_dir = IO_filereader_ID_dir;
    F5_dir = filereader_ID_dir;
    break;
  case io_dir_output:
    IO_dir = IO_out_dir;
    F5_dir = out_dir;
    break;
  case io_dir_recover:
    IO_dir = IO_recover_dir;
    F5_dir = recover_dir;
    break;
  case io_dir_checkpoint:
    IO_dir = IO_checkpoint_dir;
    F5_dir = checkpoint_dir;
    break;
  default:
    assert(0);
  }

  bool const use_IO_dir = strcmp(F5_dir, "") == 0;
  string path = use_IO_dir ? IO_dir : F5_dir;

  if (create_subdirs) {
    {
      ostringstream buf;
      buf << path << "/" << basename << ".p"
          << setw(max(0, int(processor_digits) - 4)) << setfill('0')
          << proc / 10000 << "nnnn/";
      path = buf.str();
      if (create_directories) {
        if (proc % 10000 == 0) {
          int ierr = CCTK_CreateDirectory(mode, path.c_str());
          if (ierr < 0)
            CCTK_VERROR("Could not create output directory \"%s\"",
                        path.c_str());
        }
        CCTK_Barrier(cctkGH);
      }
    }
    {
      ostringstream buf;
      buf << path << "/" << basename << ".p"
          << setw(max(0, int(processor_digits) - 2)) << setfill('0')
          << proc / 100 << "nn/";
      path = buf.str();
      if (create_directories) {
        if (proc % 100 == 0) {
          int ierr = CCTK_CreateDirectory(mode, path.c_str());
          if (ierr < 0)
            CCTK_VERROR("Could not create output directory \"%s\"",
                        path.c_str());
        }
        CCTK_Barrier(cctkGH);
      }
    }
  }
  if (one_dir_per_file) {
    ostringstream buf;
    buf << path << basename << ".p" << setw(processor_digits) << setfill('0')
        << proc << "/";
    path = buf.str();
    if (create_directories) {
      int ierr = CCTK_CreateDirectory(mode, path.c_str());
      if (ierr < 0)
        CCTK_VERROR("Could not create output directory \"%s\"", path.c_str());
    }
  }
  ostringstream buf;
  buf << path << "/" << basename;
  if (io_dir == io_dir_recover or io_dir == io_dir_checkpoint) {
    buf << ".i" << setw(iteration_digits) << setfill('0') << iteration;
  }
  buf << ".p" << setw(processor_digits) << setfill('0') << proc;
  buf << out_extension;
  path = buf.str();
  return path;
}

// Generate a good grid name (simulation name)
string generate_gridname(cGH const *const cctkGH) {
#if 0
    char const *const gridname = (char const*) UniqueSimulationID(cctkGH);
    assert(gridname);
    return gridname;
#endif

  // Use the parameter file name, since the simulation ID is too
  // long and doesn't look nice
  char parfilename[10000];
  CCTK_ParameterFilename(sizeof parfilename, parfilename);
  char const *const suffix = ".par";
  if (strlen(parfilename) >= strlen(suffix)) {
    char *const p = parfilename + strlen(parfilename) - strlen(suffix);
    if (strcmp(p, suffix) == 0) {
      *p = '\0';
    }
  }
  char *const s = strrchr(parfilename, '/');
  char *const p = s ? s + 1 : parfilename;

  return p;
}

// Generate a good topology name (refinement level name)
string generate_topologyname(cGH const *const cctkGH, int const gi,
                             ivect const &reffact, ivect const &slice_ipos) {
  ostringstream buf;
  if (gi == -1) {
    // grid function
    // bit d of centering indicates whether direction d is centered
    int const centering =
        vhh.at(0)->refcent == vertex_centered ? 0 : (1 << dim) - 1;
    char name[1000];
    // TODO: teach TopologyName to not always assume 3 dimensions
    TopologyName(name, sizeof name, &v2h(reffact)[0], centering, dim);
    buf << name;
#if 0
      if (centering) {
        buf << "CellCentering";
        for (int d=0; d<slice_dim; ++d) {
          if (centering & (1<<d)) {
            assert(d<4);
            buf << "XYZW"[d];
          } else {
            buf << "_";
          }
        }
      } else {
        buf << "VertexLevel";
      }
      for (int d=0; d<slice_dim; ++d) {
        buf << (d==0 ? "_" : "x") << slice_reffact[d];
      }
#endif
  } else {
    // grid array
    char *const groupname = CCTK_GroupName(gi);
    for (char *p = groupname; *p; ++p) {
      *p = tolower(*p);
    }
    buf << "Group_" << groupname;
    free(groupname);
  }
  // append slice locations
  if (any(slice_ipos >= 0)) {
    buf << "_Slice" << count(slice_ipos < 0) << "D";
    if (any(slice_ipos < 0)) {
      buf << "_";
      for (int d = 0; d < dim; ++d) {
        if (slice_ipos[d] < 0) {
          assert(d < 4);
          buf << "XYZW"[d];
        }
      }
    }
  }
  return buf.str();
}

// Generate a good chart name (tensor basis name)
string generate_chartname(cGH const *const cctkGH) {
  // Assume a global tensor basis, where all maps use the same chart
  return FIBER_HDF5_DEFAULT_CHART;

#if 0
    // TODO: Add some interesting information here, e.g. the name of
    // the multi-patch system, or a good name of the patch?
    ostringstream buf;
    buf << "Map_" << m;
    return buf.str();
#endif
}

// Generate a good fragment name (map and component name)
// (We assume that fragment names need to be unique only within a
// topology, not across topologies.)
string generate_fragmentname(cGH const *const cctkGH, int const m,
                             int const c) {
  ostringstream buf;
  buf << "Map_" << m << "_Component_" << c;
  return buf.str();
}

void interpret_fragmentname(cGH const *const cctkGH,
                            char const *const fragmentname, int &m, int &c) {
  m = -1;
  c = -1;
  sscanf(fragmentname, "Map_%d_Component_%d", &m, &c);
  assert(m >= 0 and m < Carpet::maps);
  assert(c >= 0);
}

// Generate a good field name (group or variable name)
string generate_fieldname(cGH const *const cctkGH, int const vi,
                          tensortype_t const tt) {
  int const gi = CCTK_GroupIndexFromVarI(vi);
  int const numvars = CCTK_NumVarsInGroupI(gi);
  string name;
  // Use the variable name instead of the group name if we may
  // output several variables per group
  if (tt == tt_scalar and numvars > 1) {
    char *const fullname = CCTK_FullName(vi);
    name = fullname;
    free(fullname);
  } else {
    char *const groupname = CCTK_GroupName(gi);
    name = groupname;
    free(groupname);
  }
  transform(name.begin(), name.end(), name.begin(), ::tolower);
  string const sep = "::";
  size_t const pos = name.find(sep);
  assert(pos != string::npos);
  name.replace(pos, sep.size(), ".");
  return name;
}

void interpret_fieldname(cGH const *const cctkGH, string fieldname, int &vi) {
  string const sep = ".";
  size_t const pos = fieldname.find(sep);
  if (pos == string::npos) {
    // The field name is not a Cactus group or variable
    vi = -1;
    return;
  }
  fieldname.replace(pos, sep.size(), "::");

  vi = CCTK_VarIndex(fieldname.c_str());
  if (vi < 0) {
    int const gi = CCTK_GroupIndex(fieldname.c_str());
    if (gi < 0) {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Unknown group variable name \"%s\"", fieldname.c_str());
      return;
    }
    vi = CCTK_FirstVarIndexI(gi);
    assert(vi >= 0);
  }
}

char const *const grid_structure = "Grid Structure v6";

string serialise_grid_structure(cGH const *const cctkGH) {
  ostringstream buf;
  buf << setprecision(17);
  buf << grid_structure << "{";
  buf << "maps=" << maps << ",";
  buf << "[";
  for (int m = 0; m < maps; ++m) {
    buf << "[" << m << "]={";
    buf << "superregions:" << vhh.AT(m)->superregions << ",";
    // We could omit the regions
    buf << "regions:" << vhh.AT(m)->regions.AT(0) << ",";
    buf << "ghost_widths:" << vdd.AT(m)->ghost_widths << ",";
    buf << "buffer_widths:" << vdd.AT(m)->buffer_widths << ",";
    buf << "overlap_widths:" << vdd.AT(m)->overlap_widths << ",";
    buf << "prolongation_orders_space:" << vdd.AT(m)->prolongation_orders_space
        << ",";
    buf << "},";
  }
  buf << "],";
  buf << "global_time:" << global_time << ",";
  buf << "delta_time:" << delta_time << ",";
  buf << "times:" << *tt << ",";
  buf << "}.";
  return buf.str();
}

void deserialise_grid_structure(cGH const *const cctkGH, string const str) {
  // Define variables for the metadata that will be read in
  int const my_maps = maps;
  vector<gh::rregs> my_superregions(my_maps);
  vector<gh::rregs> my_regions(my_maps);
  vector<vector<i2vect> > my_ghost_widths(my_maps);
  vector<vector<i2vect> > my_buffer_widths(my_maps);
  vector<vector<i2vect> > my_overlap_widths(my_maps);
  vector<int> my_prolongation_orders_space(my_maps);
  CCTK_REAL my_global_time, my_delta_time;
  th my_tt(tt->h, tt->time_interpolation_during_regridding);

  // Read in the metadata
  istringstream is(str);
  skipws(is);
  consume(is, grid_structure);
  skipws(is);
  consume(is, "{");
  skipws(is);
  consume(is, "maps");
  skipws(is);
  consume(is, '=');
  skipws(is);
  consume(is, '[');
  for (int m = 0; m < my_maps; ++m) {
    skipws(is);
    consume(is, "[");
    int my_m;
    is >> my_m;
    assert(my_m == m);
    skipws(is);
    consume(is, ']');
    skipws(is);
    consume(is, '=');
    skipws(is);
    consume(is, '{');
    skipws(is);
    consume(is, "superregions:");
    skipws(is);
    is >> my_superregions.AT(m);
    skipws(is);
    consume(is, ',');
    skipws(is);
    consume(is, "regions:");
    skipws(is);
    // ignore region specification
    // is >> my_regions.AT(m);
    gh::rregs tmp_regions;
    is >> tmp_regions;
    skipws(is);
    consume(is, ',');
    skipws(is);
    consume(is, "ghost_widths:");
    skipws(is);
    is >> my_ghost_widths.AT(m);
    skipws(is);
    consume(is, ',');
    skipws(is);
    consume(is, "buffer_widths:");
    skipws(is);
    is >> my_buffer_widths.AT(m);
    skipws(is);
    consume(is, ',');
    skipws(is);
    consume(is, "overlap_widths:");
    skipws(is);
    is >> my_overlap_widths.AT(m);
    skipws(is);
    consume(is, ',');
    skipws(is);
    consume(is, "prolongation_orders_space:");
    skipws(is);
    is >> my_prolongation_orders_space.AT(m);
    skipws(is);
    consume(is, ',');
    skipws(is);
    consume(is, '}');
    skipws(is);
    consume(is, ',');
  }
  skipws(is);
  consume(is, ']');
  skipws(is);
  consume(is, ',');
  skipws(is);
  consume(is, "global_time:");
  is >> my_global_time;
  skipws(is);
  consume(is, ',');
  skipws(is);
  consume(is, "delta_time:");
  is >> my_delta_time;
  skipws(is);
  consume(is, ',');
  skipws(is);
  consume(is, "times:");
  is >> my_tt;
  skipws(is);
  consume(is, ',');
  skipws(is);
  consume(is, '}');
  skipws(is);
  consume(is, '.');

  // Change Carpet's grid structure according to what has just been
  // read
  {
    int type;
    void const *const ptr =
        CCTK_ParameterGet("regrid_in_level_mode", "Carpet", &type);
    assert(ptr);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const regrid_in_level_mode = *static_cast<CCTK_INT const *>(ptr);
    assert(regrid_in_level_mode);
  }

  // Count refinement levels
  vector<int> rls(maps);
  for (int m = 0; m < maps; ++m) {
    rls.AT(m) = my_superregions.AT(m).size();
  }
  int maxrl = 0;
  for (int m = 0; m < maps; ++m) {
    maxrl = max(maxrl, rls.AT(m));
  }
  // All maps must have the same number of refinement levels
  for (int m = 0; m < maps; ++m) {
    my_superregions.AT(m).resize(maxrl);
    my_regions.AT(m).resize(maxrl);
  }

  // Make multiprocessor aware
  for (int rl = 0; rl < maxrl; ++rl) {
    vector<gh::cregs> superregss(maps);
    for (int m = 0; m < maps; ++m) {
      superregss.AT(m) = my_superregions.AT(m).AT(rl);
    }
    vector<gh::cregs> regss(maps);
    Carpet::SplitRegionsMaps(cctkGH, superregss, regss);
    for (int m = 0; m < maps; ++m) {
      my_superregions.AT(m).AT(rl) = superregss.AT(m);
      my_regions.AT(m).AT(rl) = regss.AT(m);
    }
  } // for rl

  // Make multigrid aware
  vector<gh::mregs> my_multigrid_regions(maps);
  Carpet::MakeMultigridBoxesMaps(cctkGH, my_regions, my_multigrid_regions);

  // Regrid (ignoring the previous content of the grid hierarchy)
  for (int m = 0; m < maps; ++m) {
    RegridMap(cctkGH, m, my_superregions.AT(m), my_multigrid_regions.AT(m),
              false);
  } // for m

  // Set up time hierarchy after RegridMap created it
  for (int ml = 0; ml < vhh.AT(0)->mglevels(); ++ml) {
    for (int rl = 0; rl < vhh.AT(0)->reflevels(); ++rl) {
      for (int tl = 0; tl < tt->timelevels; ++tl) {
        tt->set_time(ml, rl, tl, my_tt.get_time(ml, rl, tl));
      }
      tt->set_delta(ml, rl, my_tt.get_delta(ml, rl));
    }
  }

  PostRegrid(cctkGH);

  // Ensure we are in global mode
  assert(is_global_mode());
  global_time = my_global_time;
  (const_cast<cGH *>(cctkGH))->cctk_time = global_time;
  delta_time = my_delta_time;
  timereflevelfact = timereffacts.AT(reflevels - 1);
  spacereflevelfact = ivect(-get_deadbeef());
  ivect::ref(cctkGH->cctk_levfac) = spacereflevelfact;
  (const_cast<cGH *>(cctkGH))->cctk_timefac = timereflevelfact;

  // Recompose the grid structure (ignoring the previous content of
  // the grid hierarchy)
  for (int rl = 0; rl < reflevels; ++rl) {
    Recompose(cctkGH, rl, false);
  }

  // Free old grid structure
  RegridFree(cctkGH, false);

  // Check ghost and buffer widths and prolongation orders
  // TODO: Instead of only checking them, set them during
  // regridding; this requires setting them during regridding
  // instead of during setup.
  for (int m = 0; m < maps; ++m) {
    assert(vdd.AT(m)->ghost_widths.size() == my_ghost_widths.AT(m).size());
    for (int rl = 0; rl < int(my_ghost_widths.AT(m).size()); ++rl) {
      assert(all(
          all(vdd.AT(m)->ghost_widths.AT(rl) == my_ghost_widths.AT(m).AT(rl))));
    }
    assert(vdd.AT(m)->buffer_widths.size() == my_buffer_widths.AT(m).size());
    for (int rl = 0; rl < int(my_buffer_widths.AT(m).size()); ++rl) {
      assert(all(all(vdd.AT(m)->buffer_widths.AT(rl) ==
                     my_buffer_widths.AT(m).AT(rl))));
    }
    assert(vdd.AT(m)->overlap_widths.size() == my_overlap_widths.AT(m).size());
    for (int rl = 0; rl < int(my_overlap_widths.AT(m).size()); ++rl) {
      assert(all(all(vdd.AT(m)->overlap_widths.AT(rl) ==
                     my_overlap_widths.AT(m).AT(rl))));
    }
    for (int rl = 0; rl < int(my_prolongation_orders_space.size()); ++rl) {
      assert(vdd.AT(m)->prolongation_orders_space.AT(rl) ==
             my_prolongation_orders_space.AT(rl));
    }
  }
}

} // end namespace CarpetIOF5
