#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// force HDF5 1.8.x installations to use the new API
#define H5Dcreate_vers 2

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "carpet.hh"

#include "defs.hh"

#include "extending.hh"
#include "file.hh"
#include "utils.hh"
#include "writer.hh"

namespace CarpetIOF5 {

int const Error_none = 0; // no error

int const Error_illegal_varname = -1;
int const Error_group_has_no_storage = -2;

int const Error_lonely_option_string = -3;
int const Error_unterminated_option_string = -4;
int const Error_garbage_after_option_string = -5;
int const Error_invalid_variable_name = -6;

extern "C" int CarpetIOF5_Startup();

static void *Setup(tFleshConfig *const fleshconfig, int const convlevel,
                   cGH *const cctkGH);

extern "C" void CarpetIOF5_Init(CCTK_ARGUMENTS);

static int OutputGH(cGH const *cctkGH);

static void mark_variables(int variable, char const *options, void *ptr);

static int TimeToOutput(cGH const *cctkGH, int variable);

static int TriggerOutput(cGH const *cctkGH, int variable);

static int OutputVarAs(cGH const *cctkGH, char const *varname,
                       char const *alias);

static void WriteParameters(F5::file_t &file);

static string generate_filename(cGH const *cctkGH, int variable);

int CarpetIOF5_Startup() {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("Startup");
  return extending_t::create(Setup);
}

void *Setup(tFleshConfig *const fleshconfig, int const convlevel,
            cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(fleshconfig != 0);
  if (verbose)
    CCTK_INFO("Setup");
  return extending_t::setup(cctkGH, OutputGH, TimeToOutput, TriggerOutput,
                            OutputVarAs);
}

void CarpetIOF5_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("Init");

  *next_output_iteration = 0;
  *next_output_time = cctk_time;
  *this_iteration = -1;
}

int OutputGH(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH != 0);

  if (verbose)
    CCTK_INFO("OutputGH");

  int ierr;

  vector<bool> want_variables(CCTK_NumVars());
  ierr = CCTK_TraverseString(out_vars, mark_variables,
                             static_cast<void *>(&want_variables),
                             CCTK_GROUP_OR_VAR);
  switch (ierr) {
  case -2:
    return Error_lonely_option_string;
  case -3:
    return Error_unterminated_option_string;
  case -4:
    return Error_garbage_after_option_string;
  case -5:
    return Error_invalid_variable_name;
  }
  assert(ierr >= 0);

  ierr = Error_none;
  for (int variable = 0; variable < CCTK_NumVars(); ++variable) {
    if (want_variables.at(variable)) {
      if (TimeToOutput(cctkGH, variable)) {
        ierr = TriggerOutput(cctkGH, variable);
      }
    }
  }

  return ierr;
}

void mark_variables(int const variable, char const *const options,
                    void *const ptr) {
  vector<bool> &want_variables = *static_cast<vector<bool> *>(ptr);
  want_variables.at(variable) = true;
}

int TimeToOutput(cGH const *const cctkGH, int const variable) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH != 0);
  assert(variable >= 0 and variable < CCTK_NumVars());

  assert(Carpet::is_level_mode());

  if (verbose) {
    char *const fullname = CCTK_FullName(variable);
    CCTK_VInfo(CCTK_THORNSTRING, "TimeToOutput \"%s\"", fullname);
    free(fullname);
  }

  bool should_output;

  char const *const my_out_criterion =
      (CCTK_EQUALS(out_criterion, "default") ? IO_out_criterion
                                             : out_criterion);

  if (CCTK_EQUALS(my_out_criterion, "always")) {
    should_output = true;
  } else if (CCTK_EQUALS(my_out_criterion, "never")) {
    should_output = false;
  } else if (CCTK_EQUALS(my_out_criterion, "iteration")) {
    int const my_out_every = out_every == -2 ? IO_out_every : out_every;
    switch (my_out_every) {
    case 0:
      should_output = true;
      break;
    case -1:
      should_output = false;
      break;
    default:
      if (*this_iteration == cctk_iteration) {
        // we already decided to output this iteration
        should_output = true;
      } else if (cctk_iteration >= *next_output_iteration) {
        // it is time for the next output
        should_output = true;
        *this_iteration = cctk_iteration;
        *next_output_iteration = cctk_iteration + my_out_every;
      } else {
        should_output = false;
      }
      break;
    }
  } else if (CCTK_EQUALS(my_out_criterion, "time")) {
    CCTK_REAL const my_out_dt = out_dt == -2 ? IO_out_dt : out_dt;
    if (out_dt == 0) {
      should_output = true;
    } else if (out_dt == -1) {
      should_output = false;
    } else {
      if (*this_iteration == cctk_iteration) {
        // we already decided to output this iteration
        should_output = true;
      } else if (cctk_time / cctk_delta_time >=
                 *next_output_time / cctk_delta_time - dt_fudge) {
        // it is time for the next output
        should_output = true;
        *this_iteration = cctk_iteration;
        *next_output_time = cctk_time + my_out_dt;
      } else {
        should_output = false;
      }
    }
  } else {
    CCTK_WARN(1, "internal error");
    should_output = false;
  }

  if (should_output) {
    extending_t extending(cctkGH);
    int const last_output_iteration = (extending.get_last_output_iteration(
        Carpet::mglevel, Carpet::reflevel, variable));
    assert(last_output_iteration <= cctk_iteration);
    if (last_output_iteration == cctk_iteration) {
      // Skipping output for variable, because this variable has
      // already been output during the current iteration --
      // probably via a trigger during the analysis stage
      should_output = false;
    } else {
      extending.set_last_output_iteration(Carpet::mglevel, Carpet::reflevel,
                                          variable, cctk_iteration);
    }
  }

  return should_output;
}

int TriggerOutput(cGH const *const cctkGH, int const variable) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH != 0);
  assert(variable >= 0 and variable < CCTK_NumVars());

  char *const fullname = CCTK_FullName(variable);
  assert(fullname);

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "TriggerOutput \"%s\"", fullname);
  }

  string const alias = generate_filename(cctkGH, variable);

  int const ierr = OutputVarAs(cctkGH, fullname, alias.c_str());

  free(fullname);

  return ierr;
}

int OutputVarAs(cGH const *const cctkGH, char const *const varname,
                char const *const alias) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH != 0);
  assert(varname != 0);
  assert(alias != 0);

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs \"%s\" \"%s\"", varname, alias);
  }

  int const variable = CCTK_VarIndex(varname);
  if (variable < 0) {
    return Error_illegal_varname;
  }
  assert(variable >= 0 and variable < CCTK_NumVars());

  int const group = CCTK_GroupIndexFromVarI(variable);
  assert(group >= 0 and group < CCTK_NumGroups());

  if (CCTK_ActiveTimeLevelsGI(cctkGH, group) == 0) {
    return Error_group_has_no_storage;
  }

  extending_t extending(cctkGH);

  bool const use_IO_out_dir = strcmp(out_dir, "") == 0;
  string const path = use_IO_out_dir ? IO_out_dir : out_dir;
  string const basename = alias;

  bool const did_truncate = extending.get_did_truncate(basename);
  bool const do_truncate = not did_truncate and IO_TruncateOutputFiles(cctkGH);
  extending.set_did_truncate(basename);

  int const proc = CCTK_MyProc(cctkGH);
  bool have_metafile;     // whether there is a metadata file
  int metadata_processor; // the processor which outputs the metadata file
  int output_processor;   // the processor which outputs our data
  if (CCTK_EQUALS(out_mode, "proc")) {
    have_metafile = true;
    metadata_processor = 0;
    output_processor = proc;
  } else if (CCTK_EQUALS(out_mode, "np")) {
    have_metafile = true;
    metadata_processor = 0;
    output_processor = proc / out_proc_every * out_proc_every;
  } else if (CCTK_EQUALS(out_mode, "onefile")) {
    have_metafile = false;
    metadata_processor = 0;
    output_processor = 0;
  } else {
    assert(0);
  }

  F5::file_t *metafile = NULL;
  if (have_metafile and proc == metadata_processor) {
    metafile = new F5::file_t(cctkGH, path, basename, string(out_extension),
                              do_truncate, true, false);
  }

  F5::file_t *file = NULL;
  if (proc == output_processor) {
    file = new F5::file_t(cctkGH, path, basename, string(out_extension),
                          do_truncate, not have_metafile, true);
  }

  if (do_truncate) {
    // Output parameters once after the output file has been created
    if (proc == metadata_processor) {
      if (CCTK_EQUALS(out_save_parameters, "all") or
          CCTK_EQUALS(out_save_parameters, "only set")) {
        WriteParameters(have_metafile ? *metafile : *file);
      } else if (CCTK_EQUALS(out_save_parameters, "no")) {
        // do nothing
      } else {
        assert(0);
      }
    }
  }

  if (metafile) {
    writer_t writer(cctkGH, variable);
    writer.write(*metafile);
    delete metafile;
    metafile = NULL;
  }
  {
    writer_t writer(cctkGH, variable);
    // TODO: handle the case where not all processors are writing to
    // their own file
    assert(proc == output_processor);
    assert(file);
    writer.write(*file);
    delete file;
    file = NULL;
  }

  return Error_none;
}

void WriteParameters(F5::file_t &file) {
  DECLARE_CCTK_PARAMETERS;

  cGH const *const cctkGH = file.get_cctkGH();

  hid_t const hdf5_file = file.get_hdf5_file();

  hid_t const attribute_group =
      F5::open_or_create_group(hdf5_file, "Parameters and Global Attributes");
  assert(attribute_group >= 0);

  // unique configuration identifier
  if (CCTK_IsFunctionAliased("UniqueConfigID")) {
    F5::write_or_check_attribute(
        attribute_group, "config id",
        static_cast<char const *>(UniqueConfigID(cctkGH)));
  }

  // unique build identifier
  if (CCTK_IsFunctionAliased("UniqueBuildID")) {
    F5::write_or_check_attribute(
        attribute_group, "build id",
        static_cast<char const *>(UniqueBuildID(cctkGH)));
  }

  // unique simulation identifier
  if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
    F5::write_or_check_attribute(
        attribute_group, "simulation id",
        static_cast<char const *>(UniqueSimulationID(cctkGH)));
  }

  // unique run identifier
  if (CCTK_IsFunctionAliased("UniqueRunID")) {
    F5::write_or_check_attribute(
        attribute_group, "run id",
        static_cast<char const *>(UniqueRunID(cctkGH)));
  }

  // Output Cactus parameters as single string
  {
    char *const parameters = IOUtil_GetAllParameters(cctkGH, true);
    assert(parameters);
    // Create a dataset, since the data may not fit into an attribute
    hsize_t const size = strlen(parameters) + 1;
    hid_t const dataspace = H5Screate_simple(1, &size, NULL);
    assert(dataspace >= 0);
    hid_t properties = H5Pcreate(H5P_DATASET_CREATE);
    assert(properties >= 0);
    check(not H5Pset_chunk(properties, 1, &size));
    if (compression_level > 0) {
      check(not H5Pset_shuffle(properties));
      check(not H5Pset_deflate(properties, compression_level));
    }
    if (write_checksum) {
      check(not H5Pset_fletcher32(properties));
    }
    hid_t const dataset =
        H5Dcreate(attribute_group, "All Parameters", H5T_NATIVE_CHAR, dataspace,
                  H5P_DEFAULT, properties, H5P_DEFAULT);
    assert(dataset >= 0);
    check(not H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       parameters));
    check(not H5Dclose(dataset));
    check(not H5Pclose(properties));
    check(not H5Sclose(dataspace));
    free(parameters);
  }

  check(not H5Gclose(attribute_group));

// This is far too slow to be useful
#if 0
    hid_t const parameter_group
      = F5::open_or_create_group (hdf5_file, "Cactus parameters");
    assert (parameter_group >= 0);
    
    int first = 1;
    for (;;)
    {
      cParamData const * parameter_data;
      char * parameter_fullname;
      
      int const ierr
        = CCTK_ParameterWalk (first, 0,
                              & parameter_fullname, & parameter_data);
      if (ierr > 0) break;
      assert (ierr >= 0);
      
      int type;
      void const * const parameter_value
        = CCTK_ParameterGet (parameter_data->name, parameter_data->thorn,
                             & type);
      assert (type == parameter_data->type);
      assert (parameter_value != 0);
      
      switch (type)
      {
      case PARAMETER_BOOLEAN:
      case PARAMETER_INT:
        {
          CCTK_INT const value
            = * static_cast<CCTK_INT const *> (parameter_value);
          F5::write_or_check_attribute
            (parameter_group, parameter_fullname, value);
        }
        break;
      case PARAMETER_REAL:
        {
          CCTK_REAL const value
            = * static_cast<CCTK_REAL const *> (parameter_value);
          F5::write_or_check_attribute
            (parameter_group, parameter_fullname, value);
        }
        break;
      case PARAMETER_KEYWORD:
      case PARAMETER_STRING:
        {
          char const * const value
            = * static_cast<char const * const *> (parameter_value);
          F5::write_or_check_attribute
            (parameter_group, parameter_fullname, value);
        }
        break;
      default:
        assert (0);
      }
      
      free (parameter_fullname);
      
      first = 0;
    }
    
    check (not H5Gclose (parameter_group));
#endif
}

string generate_filename(cGH const *const cctkGH, int const variable) {
  DECLARE_CCTK_PARAMETERS;

  assert(variable >= 0);

  ostringstream filename_buf;

  if (CCTK_EQUALS(file_content, "variable")) {
    char *const varname = CCTK_FullName(variable);
    assert(varname);
    for (char *p = varname; *p; ++p) {
      *p = tolower(*p);
    }
    filename_buf << varname;
    free(varname);
  } else if (CCTK_EQUALS(file_content, "group")) {
    char *const groupname = CCTK_GroupNameFromVarI(variable);
    assert(groupname);
    for (char *p = groupname; *p; ++p) {
      *p = tolower(*p);
    }
    filename_buf << groupname;
    free(groupname);
  } else if (CCTK_EQUALS(file_content, "thorn")) {
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
    filename_buf << out_filename;
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

} // namespace CarpetIOF5
