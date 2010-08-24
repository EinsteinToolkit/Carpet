#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdio>
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

#include <bbox.hh>
#include <defs.hh>
#include <vect.hh>

#include <carpet.hh>

#include "iof5.hh"



ostream& operator<<(ostream& os, CarpetIOF5::indent_t const& indent)
{ return indent.output(os); }



namespace CarpetIOF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Indentation
  int indent_t::level = 0;
  
  ostream& indent_t::output(ostream& os) const
  {
    return os << string(width*level, ' ');
  }
  
  
  
  // Generate a good file name ("alias") for a variable
  string
  generate_basename (cGH const *const cctkGH,
                     int const variable)
  {
    DECLARE_CCTK_PARAMETERS;
    
    ostringstream filename_buf;
    
    if (CCTK_EQUALS (file_content, "variable")) {
      assert (variable >= 0);
      char *const varname = CCTK_FullName(variable);
      assert (varname);
      for (char *p = varname; *p; ++p) {
        *p = tolower(*p);
      }
      filename_buf << varname;
      free (varname);
    }
    else if (CCTK_EQUALS (file_content, "group")) {
      assert (variable >= 0);
      char *const groupname = CCTK_GroupNameFromVarI(variable);
      assert (groupname);
      for (char *p = groupname; *p; ++p) {
        *p = tolower(*p);
      }
      filename_buf << groupname;
      free (groupname);
    }
    else if (CCTK_EQUALS (file_content, "thorn")) {
      assert (variable >= 0);
      char const *const impname = CCTK_ImpFromVarI(variable);
      char *const thornname = strdup(impname);
      assert (thornname);
      char *const colon = strchr(thornname, ':');
      assert (colon);
      *colon = '\0';
      for (char *p = thornname; *p; ++p) {
        *p = tolower(*p);
      }
      filename_buf << thornname;
      free (thornname);
    }
    else if (CCTK_EQUALS (file_content, "everything")) {
      filename_buf << out_filename;
    }
    else {
      assert (0);
    }
    
    if (out_timesteps_per_file > 0) {
      int const iteration = (cctkGH->cctk_iteration
                             / out_timesteps_per_file * out_timesteps_per_file);
      filename_buf << ".it"
                   << setw (iteration_digits) << setfill ('0') << iteration;
    }
    
    return filename_buf.str();
  }
  
  
  
  // Create the final file name on a particular processor
  string
  create_filename (cGH const *const cctkGH,
                   string const basename,
                   int const proc,
                   bool const create_directories)
  {
    DECLARE_CCTK_PARAMETERS;
    
    bool const use_IO_out_dir = strcmp(out_dir, "") == 0;
    string path = use_IO_out_dir ? IO_out_dir : out_dir;
    
    if (create_subdirs) {
      {
        ostringstream buf;
        buf << path + "/"
            << basename
            << ".p" << setw (max (0, processor_digits - 4)) << setfill ('0')
            << proc / 10000
            << "nnnn/";
        path = buf.str();
        if (create_directories) {
          check (CCTK_CreateDirectory (mode, path.c_str()) >= 0);
        }
      }
      {
        ostringstream buf;
        buf << path + "/"
            << basename
            << ".p" << setw (max (0, processor_digits - 2)) << setfill ('0')
            << proc / 100
            << "nn/";
        path = buf.str();
        if (create_directories) {
          check (CCTK_CreateDirectory (mode, path.c_str()) >= 0);
        }
      }
    }
    if (one_dir_per_file) {
      ostringstream buf;
      buf << path
          << basename
          << ".p" << setw (processor_digits) << setfill ('0')
          << proc
          << "/";
      path = buf.str();
      if (create_directories) {
        check (CCTK_CreateDirectory (mode, path.c_str()) >= 0);
      }
    }
    ostringstream buf;
    buf << path + "/"
        << basename
        << ".p" << setw (processor_digits) << setfill ('0') << proc
        << out_extension;
    return buf.str();
  }
  
  
  
  // Generate a good grid name (simulation name)
  string
  generate_gridname (cGH const *const cctkGH)
  {
#if 0
    char const *const gridname = (char const*) UniqueSimulationID(cctkGH);
    assert (gridname);
    return gridname;
#endif
    
    // Use the parameter file name, since the simulation ID is too
    // long and doesn't look nice
    char parfilename[10000];
    CCTK_ParameterFilename (sizeof parfilename, parfilename);
    char const *const suffix = ".par";
    if (strlen(parfilename) >= strlen(suffix)) {
      char *const p = parfilename + strlen(parfilename) - strlen(suffix);
      if (strcmp(p, suffix) == 0) {
        *p = '\0';
      }
    }
    char *const s = strrchr(parfilename, '/');
    char *const p = s ? s+1 : parfilename;
    
    return p;
  }
  
  // Generate a good tensor basis name (chart name)
  string
  generate_chartname (cGH const *const cctkGH)
  {
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
  
  // Generate a good fragment name (map name)
  string
  generate_fragmentname (cGH const *const cctkGH, int const m)
  {
    ostringstream buf;
    buf << "Map_" << m;
    return buf.str();
  }
  
  int
  interpret_fragmentname (cGH const *const cctkGH,
                          char const *const fragmentname)
  {
    int m = -1;
    sscanf (fragmentname, "Map_%d", &m);
    assert (m>=0 and m<Carpet::maps);
    return m;
  }
  
  // Generate a good topology name (map and refinement level name)
  string
  generate_topologyname (cGH const *const cctkGH,
                         int const m, ivect const& reffact)
  {
    ostringstream buf;
    // buf << "Map_" << m << "_VertexLevel_";
    buf << "VertexLevel_";
    for (int d=0; d<dim; ++d) {
      if (d>0) buf << "x";
      buf << reffact[d];
    }
    return buf.str();
  }
  
  
  
  char const *const grid_structure = "Grid Structure v5";
  
  string
  serialise_grid_structure (cGH const *const cctkGH)
  {
    ostringstream buf;
    buf << setprecision(17);
    buf << grid_structure << "{";
    buf << "maps=" << maps << ",";
    buf << "[";
    for (int m=0; m<maps; ++m) {
      buf << "[" << m << "]={";
      buf << "superregions:" << vhh.AT(m)->superregions << ",";
      // We could omit the regions
      buf << "regions:" << vhh.AT(m)->regions.AT(0) << ",";
      buf << "ghost_widths:" << vdd.AT(m)->ghost_widths << ",";
      buf << "buffer_widths:" << vdd.AT(m)->buffer_widths << ",";
      buf << "prolongation_orders_space:"
          << vdd.AT(m)->prolongation_orders_space << ",";
      buf << "},";
    }
    buf << "],";
    buf << "times:" << *tt << ",";
    buf << "}.";
    return buf.str();
  }
  
} // end namespace CarpetIOF5
