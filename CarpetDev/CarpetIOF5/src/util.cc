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
  
  indent_t::indent_t()
  {
    if (debug) cout << *this << "{{{\n";
    ++level;
  }
  
  indent_t::~indent_t()
  {
    --level;
    if (debug) cout << *this << "}}}\n";
  }
  
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
      if (CCTK_EQUALS(out_filename, "")) {
        // Obtain the parameter file name
        char buf[10000];
        int const ilen = CCTK_ParameterFilename (sizeof buf, buf);
        assert (ilen < int(sizeof buf) - 1);
        char* parfilename = buf;
        // Remove directory prefix, if any
        char* const slash = strrchr(parfilename, '/');
        if (slash) {
          parfilename = &slash[1];
        }
        // Remove suffix, if it is there
        char* const suffix = strrchr(parfilename, '.');
        if (suffix and strcmp(suffix, ".par")==0) {
          suffix[0] = '\0';
        }
        filename_buf << parfilename;
      } else {
        filename_buf << out_filename;
      }
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
                   io_dir_t const io_dir,
                   bool const create_directories)
  {
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
  
  // Generate a good topology name (refinement level name)
  string
  generate_topologyname (cGH const *const cctkGH,
                         int const gi,
                         ivect const& reffact)
  {
    ostringstream buf;
    if (gi == -1) {
      // grid function
      buf << "VertexLevel_";
      for (int d=0; d<dim; ++d) {
        if (d>0) buf << "x";
        buf << reffact[d];
      }
    } else {
      // grid array
      char *const groupname = CCTK_GroupName(gi);
      for (char *p=groupname; *p; ++p) {
        *p = tolower(*p);
      }
      buf << "Group_" << groupname;
      free (groupname);
    }
    return buf.str();
  }
  
  // Generate a good chart name (tensor basis name)
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
  
  // Generate a good fragment name (map and component name)
  // (We assume that fragment names need to be unique only within a
  // topology, not across topologies.)
  string
  generate_fragmentname (cGH const *const cctkGH, int const m, int const c)
  {
    ostringstream buf;
    buf << "Map_" << m << "_Component_" << c;
    return buf.str();
  }
  
  void
  interpret_fragmentname (cGH const *const cctkGH,
                          char const *const fragmentname,
                          int& m, int& c)
  {
    m = -1;
    c = -1;
    sscanf (fragmentname, "Map_%d_Component_%d", &m, &c);
    assert (m>=0 and m<Carpet::maps);
    assert (c>=0);
  }
  
  // Generate a good field name (group or variable name)
  string
  generate_fieldname (cGH const *const cctkGH,
                      int const vi, tensortype_t const tt)
  {
    int const gi = CCTK_GroupIndexFromVarI(vi);
    int const numvars = CCTK_NumVarsInGroupI(gi);
    string name;
    // Use the variable name instead of the group name if we may
    // output several variables per group
    if (tt == tt_scalar and numvars >1) {
      char *const fullname = CCTK_FullName(vi);
      name = fullname;
      free (fullname);
    } else {
      char *const groupname = CCTK_GroupName(gi);
      name = groupname;
      free (groupname);
    }
    transform (name.begin(), name.end(), name.begin(), ::tolower);
    string const sep = "::";
    size_t const pos = name.find(sep);
    assert (pos != string::npos);
    name.replace (pos, sep.size(), ".");
    return name;
  }
  
  void
  interpret_fieldname(cGH const *const cctkGH, string fieldname, int& vi)
  {
    string const sep = ".";
    size_t const pos = fieldname.find(sep);
    if (pos == string::npos) {
      // The field name is not a Cactus group or variable
      vi = -1;
      return;
    }
    fieldname.replace (pos, sep.size(), "::");
    
    vi = CCTK_VarIndex(fieldname.c_str());
    if (vi < 0) {
      int const gi = CCTK_GroupIndex(fieldname.c_str());
      if (gi < 0) {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unknown group variable name \"%s\"", fieldname.c_str());
        return;
      }
      vi = CCTK_FirstVarIndexI(gi);
      assert (vi>=0);
    }
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
