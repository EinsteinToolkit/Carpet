#include <cassert>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <map>
#include <string>

#include "cctk.h"
#include "util_Table.h"

#include "CactusBase/IOUtil/src/ioGH.h"

#include "CarpetTimers.hh"

#include "CarpetIOHDF5.hh"



// That's a hack
namespace Carpet {
  void UnsupportedVarType (const int vindex);
}


#define GetParameter(parameter) \
  outdim == 0 ? out0D_##parameter : \
  outdim == 1 ? out1D_##parameter : out2D_##parameter

namespace CarpetIOHDF5 {

  using namespace std;
  using namespace Carpet;



  // IO processor
  const int ioproc = 0;
  const int nioprocs = 1;



  // Global configuration parameters
  bool stop_on_parse_errors = false;



  // Definition of static members
  template<int outdim> char* IOHDF5<outdim>::my_out_slice_dir;
  template<int outdim> char* IOHDF5<outdim>::my_out_slice_vars;
  template<int outdim> vector<ioRequest*> IOHDF5<outdim>::slice_requests;



  template<int outdim>
  int IOHDF5<outdim>::Startup()
  {
    ostringstream msg;
    msg << "AMR " << outdim << "D HDF5 I/O provided by CarpetIOHDF5";
    CCTK_RegisterBanner (msg.str().c_str());

    ostringstream extension_name;
    extension_name << "CarpetIOHDF5_" << outdim << "D";
    const int GHExtension =
      CCTK_RegisterGHExtension(extension_name.str().c_str());
    CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);

    ostringstream method_name;
    method_name << "IOHDF5_" << outdim << "D";
    const int IOMethod = CCTK_RegisterIOMethod (method_name.str().c_str());
    CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
    CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
    CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
    CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);

    return 0;
  }



  template<int outdim>
  void* IOHDF5<outdim>::SetupGH (tFleshConfig* const fc,
                                 const int convLevel, cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    const void *dummy;

    dummy = &fc;
    dummy = &convLevel;
    dummy = &cctkGH;
    dummy = &dummy;

    if (not CCTK_Equals (verbose, "none")) {
      CCTK_VInfo (CCTK_THORNSTRING, "I/O Method 'IOHDF5_%dD' registered: "
                  "%dD AMR output of grid variables to HDF5 files",
                  outdim, outdim);
    }

    const int numvars = CCTK_NumVars();
    slice_requests.resize (numvars);

    // initial I/O parameter check
    my_out_slice_dir = 0;
    my_out_slice_vars = strdup ("");
    stop_on_parse_errors = strict_io_parameter_check != 0;
    CheckSteerableParameters (cctkGH);
    stop_on_parse_errors = false;

    // We register only once, ergo we get only one handle.  We store
    // that statically, so there is no need to pass anything to
    // Cactus.
    return NULL;
  }



  template<int outdim>
  void IOHDF5<outdim>::CheckSteerableParameters (const cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;

    // re-parse the 'IOHDF5::out%dD_dir' parameter if it has changed
    const char* the_out_dir = GetParameter(dir);
    if (CCTK_EQUALS (the_out_dir, "")) {
      the_out_dir = io_out_dir;
    }

    if (not my_out_slice_dir or strcmp (the_out_dir, my_out_slice_dir)) {
      free (my_out_slice_dir);
      my_out_slice_dir = strdup (the_out_dir);

      // create the output directory
      const int result = IOUtil_CreateDirectory (cctkGH, my_out_slice_dir, 0, 0);
      if (result < 0) {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Problem creating %dD-output directory '%s'",
                    outdim, my_out_slice_dir);
      } else if (result > 0 and CCTK_Equals (verbose, "full")) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "%dD-output directory '%s' already exists",
                    outdim, my_out_slice_dir);
      }
    }

    // re-parse the 'IOHDF5::out%d_vars' parameter if it has changed
    const char* const out_slice_vars = GetParameter(vars);
    if (strcmp (out_slice_vars, my_out_slice_vars)) {
      ostringstream parameter_name;
      parameter_name << "IOHDF5::out" << outdim << "D_vars";
#ifdef IOUTIL_PARSER_HAS_OUT_DT
      IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                                 parameter_name.str().c_str(),
                                 stop_on_parse_errors, out_slice_vars,
                                 -1, -1.0, &slice_requests[0]);
#else
      IOUtil_ParseVarsForOutput (cctkGH, CCTK_THORNSTRING,
                                 parameter_name.str().c_str(),
                                 stop_on_parse_errors, out_slice_vars,
                                 -1, &slice_requests[0]);
#endif

      // notify the user about the new setting
      if (not CCTK_Equals (verbose, "none")) {
        int count = 0;
        ostringstream msg;
        msg << "Periodic " << outdim << "D AMR output requested for:";
        for (int vi=0; vi< CCTK_NumVars(); ++vi) {
          if (slice_requests.at(vi)) {
            ++count;
            char* const fullname = CCTK_FullName(vi);
            msg << "\n" << "   " << fullname;
            free (fullname);
          }
        }
        if (count > 0) {
          CCTK_INFO (msg.str().c_str());
        }
      }

      // save the last setting of 'IOHDF5::out%d_vars' parameter
      free (my_out_slice_vars);
      my_out_slice_vars = strdup (out_slice_vars);
    }
  }



  template<int outdim>
  int IOHDF5<outdim>::OutputGH (const cGH* const cctkGH)
  {
    static Carpet::Timer * timer = NULL;
    if (not timer) {
      ostringstream timer_name;
      timer_name << "CarpetIOHDF5<" << outdim << ">::OutputGH";
      timer = new Carpet::Timer (timer_name.str().c_str());
    }

    timer->start();
    for (int vi=0; vi<CCTK_NumVars(); ++vi) {
      if (TimeToOutput(cctkGH, vi)) {
        TriggerOutput(cctkGH, vi);
      }
    }
    timer->stop();

    return 0;
  }



  template<int outdim>
  int IOHDF5<outdim>::TimeToOutput (const cGH* const cctkGH, const int vindex)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert (vindex >= 0 and vindex < CCTK_NumVars());

    if (CCTK_GroupTypeFromVarI(vindex) != CCTK_GF and not do_global_mode) {
      return 0;
    }

    CheckSteerableParameters (cctkGH);

    // check if output for this variable was requested
    if (not slice_requests.at(vindex)) {
      return 0;
    }

    // check whether this refinement level should be output
    if (not (slice_requests.at(vindex)->refinement_levels & (1 << reflevel))) {
      return 0;
    }

    // check if output for this variable was requested individually by
    // a "<varname>{ out_every = <number> }" option string
    // this will overwrite the output criterion setting
    const char* myoutcriterion = GetParameter(criterion);
    if (CCTK_EQUALS(myoutcriterion, "default")) {
      myoutcriterion = io_out_criterion;
    }
    if (slice_requests.at(vindex)->out_every >= 0) {
      myoutcriterion = "divisor";
    }

    if (CCTK_EQUALS (myoutcriterion, "never")) {
      return 0;
    }

    // check whether to output at this iteration
    bool output_this_iteration = false;

    if (CCTK_EQUALS (myoutcriterion, "iteration")) {
      int myoutevery = GetParameter(every);
      if (myoutevery == -2) {
        myoutevery = io_out_every;
      }
      if (myoutevery > 0) {
        if (cctk_iteration == this_iteration_slice[outdim]) {
          // we already decided to output this iteration
          output_this_iteration = true;
        } else if (cctk_iteration >=
                   last_output_iteration_slice[outdim] + myoutevery) {
          // it is time for the next output
          output_this_iteration = true;
          last_output_iteration_slice[outdim] = cctk_iteration;
          this_iteration_slice[outdim] = cctk_iteration;
        }
      }
    } else if (CCTK_EQUALS (myoutcriterion, "divisor")) {
      int myoutevery = GetParameter(every);
      if (myoutevery == -2) {
        myoutevery = io_out_every;
      }
      if (slice_requests[vindex]->out_every >= 0) {
        myoutevery = slice_requests[vindex]->out_every;
      }
      if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
        // we already decided to output this iteration
        output_this_iteration = true;
      }
    } else if (CCTK_EQUALS (myoutcriterion, "time")) {
      CCTK_REAL myoutdt = GetParameter(dt);
      if (myoutdt == -2) {
        myoutdt = io_out_dt;
      }
      if (myoutdt == 0 or cctk_iteration == this_iteration_slice[outdim]) {
        output_this_iteration = true;
      } else if (myoutdt > 0) {
        int do_output =
          cctk_time / cctk_delta_time >=
          (last_output_time_slice[outdim] + myoutdt) / cctk_delta_time - 1.0e-12;
        MPI_Bcast (&do_output, 1, MPI_INT, 0, dist::comm());
        if (do_output) {
          // it is time for the next output
          output_this_iteration = true;
          last_output_time_slice[outdim] = cctk_time;
          this_iteration_slice[outdim] = cctk_iteration;
        }
      }
    } // select output criterion

    return output_this_iteration ? 1 : 0;
  }



  template<int outdim>
  int IOHDF5<outdim>::TriggerOutput (const cGH* const cctkGH, const int vindex)
  {
    DECLARE_CCTK_PARAMETERS;

    assert (vindex >= 0 and vindex < CCTK_NumVars());

    char* const fullname = CCTK_FullName(vindex);

    int retval;
    if (one_file_per_group) {
      char* const alias = CCTK_GroupNameFromVarI (vindex);
      for (char* p = alias; *p; ++p) *p = (char) tolower (*p);
      retval = OutputVarAs (cctkGH, fullname, alias);
      free (alias);
    } else {
      const char* const alias = CCTK_VarName (vindex);
      retval = OutputVarAs (cctkGH, fullname, alias);
    }

    free (fullname);

    return retval;
  }



  static void GetVarIndex (const int vindex, const char* const optstring,
                           void* const arg)
  {
    if (optstring) {
      char* const fullname = CCTK_FullName(vindex);
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Option string '%s' will be ignored for HDF5 output of "
                  "variable '%s'", optstring, fullname);
      free (fullname);
    }

    *static_cast<int*>(arg) = vindex;
  }

  template<int outdim>
  int IOHDF5<outdim>::OutputVarAs (const cGH* const cctkGH,
                                   const char* const varname,
                                   const char* const alias)
  {
    DECLARE_CCTK_PARAMETERS;

    int vindex = -1;

    if (CCTK_TraverseString (varname, GetVarIndex, &vindex, CCTK_VAR) < 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "error while parsing variable name '%s' (alias name '%s')",
                  varname, alias);
      return -1;
    }

    if (vindex < 0) {
      return -1;
    }

    if (not (is_level_mode() or
             (is_singlemap_mode() and maps == 1) or
             (is_local_mode() and maps == 1 and
              vhh.at(Carpet::map)->local_components(reflevel) == 1)))
    {
      CCTK_WARN (1, "OutputVarAs must be called in level mode");
      return -1;
    }

    BEGIN_LEVEL_MODE (cctkGH) {

      // Get information
      const int group = CCTK_GroupIndexFromVarI (vindex);
      assert (group >= 0);
      cGroup groupdata;
      {
        int const ierr = CCTK_GroupData (group, & groupdata);
        assert (not ierr);
      }

      // Check information
      if (groupdata.grouptype != CCTK_GF) {
        assert (do_global_mode);
      }

      if (outdim >= groupdata.dim) {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot produce %dD slice HDF5 output file '%s' for variable '%s' "
                    "because it has only %d dimensions",
                    outdim, alias, varname, groupdata.dim);
        return -1;
      }

      // Check for storage
      if (not CCTK_QueryGroupStorageI (cctkGH, group)) {
        // This may be okay if storage is e.g. scheduled only in the
        // analysis bin
        CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot output variable '%s' because it has no storage",
                    varname);
        return 0;
      }

      ostringstream basefilenamebuf;
      basefilenamebuf << my_out_slice_dir << "/" << alias;
      const string basefilename = basefilenamebuf.str();

      // Check if the file has been created already
      bool is_new_file, truncate_file;
      const bool did_output =
        DidOutput (cctkGH, vindex, basefilename, is_new_file, truncate_file);
      if (did_output) {
        return 0;
      }

      // Loop over all direction combinations
      vect<int,outdim> dirs (0);
      bool done;
      do {

        // Output each combination only once
        bool ascending = true;
        for (int d1=0; d1<outdim; ++d1) {
          for (int d2=d1+1; d2<outdim; ++d2) {
            ascending = ascending and dirs[d1] < dirs[d2];
          }
        }

        // Skip output if the dimensions are not ascending
        if (ascending) {

          // Skip output if not requested
          if (DirectionIsRequested(dirs)) {
            OutputDirection (cctkGH, vindex, alias, basefilename, dirs,
                             is_new_file, truncate_file);
          }

        } // if ascending

        // Next direction combination
        done = true;
        for (int d=0; d<outdim; ++d) {
          ++dirs[d];
          if (dirs[d]<groupdata.dim + (outdim == 1 ? 1 : 0)) {
            done = false;
            break;
          }
          dirs[d] = 0;
        }

      } while (not done);       // output all directions

    } END_LEVEL_MODE;

    return 0;
  }



  // Traverse all maps and components on this refinement and multigrid
  // level
  template<int outdim>
  void IOHDF5<outdim>::OutputDirection (const cGH* const cctkGH,
                                        const int vindex,
                                        const string alias,
                                        const string basefilename,
                                        const vect<int,outdim>& dirs,
                                        const bool is_new_file,
                                        const bool truncate_file)
  {
    DECLARE_CCTK_PARAMETERS;

    // Get information
    const int group = CCTK_GroupIndexFromVarI (vindex);
    assert (group >= 0);
    const int vindex0 = CCTK_FirstVarIndexI (group);
    assert (vindex0 >= 0 and vindex >= vindex0);
    const int var = vindex - vindex0;
    cGroup groupdata;
    {
      int const ierr = CCTK_GroupData (group, & groupdata);
      assert (not ierr);
    }

    const int ml = groupdata.grouptype == CCTK_GF ? mglevel : 0;
    const int rl = groupdata.grouptype == CCTK_GF ? reflevel : 0;

    const int num_tl = CCTK_NumTimeLevelsFromVarI (vindex);
    assert (num_tl >= 1);

    const int numvars = CCTK_NumVarsInGroupI(group);


    // Loop over all maps
    const int m_min = 0;
    const int m_max = groupdata.grouptype == CCTK_GF ? Carpet::maps : 1;
    for (int m = m_min; m < m_max; ++ m) {

      hid_t file = -1;
      int error_count = 0;
      error_count += OpenFile (cctkGH, m, vindex, numvars, alias, basefilename,
                               dirs, is_new_file, truncate_file, file);

      // Find the output offset
      const ivect offset =
        groupdata.grouptype == CCTK_GF ? GetOutputOffset (cctkGH, m, dirs) : 0;

      // Traverse all components on this multigrid level, refinement
      // level, and map
      const int c_min = 0;
      const int c_max =
        groupdata.grouptype == CCTK_GF ?
        vhh.at(m)->components(reflevel) :
        groupdata.disttype != CCTK_DISTRIB_CONSTANT ?
        CCTK_nProcs(cctkGH) :
        1;
      for (int c = c_min; c < c_max; ++ c) {

        const ggf* const ff = arrdata.at(group).at(m).data.at(var);

        const int tl_min = 0;
        const int tl_max = output_all_timelevels ? num_tl : 1;
        for (int tl = tl_min; tl < tl_max; ++tl) {

          const gdata* const data = (*ff) (tl, rl, c, ml);
          const ibbox ext = GetOutputBBox (cctkGH, group, m, c,
                                           data->extent());

          rvect coord_lower, coord_upper;
          GetCoordinates (cctkGH, m, groupdata, ext,
                          coord_lower, coord_upper);

          // Apply offset
          ivect offset1 = offset * ext.stride();
          if (groupdata.grouptype == CCTK_GF) {
            const ibbox& baseext = vhh.at(m)->baseextents.at(ml).at(rl);
            offset1 += baseext.lower();
          }
          for (int d=0; d<outdim; ++d) {
            if (dirs[d] < 3) {
              offset1[dirs[d]] = ext.lower()[dirs[d]];
            }
          }

          vector<const gdata*> datas;
          if (one_file_per_group) {
            datas.resize (numvars);
            for (int n=0; n<numvars; ++n) {
              const ggf* const ff1 = arrdata.at(group).at(m).data.at(n);
              datas.at(n) = (*ff1) (tl, rl, c, ml);
            }
          } else {
            datas.resize (1);
            datas.at(0) = data;
          }

          error_count += WriteHDF5(cctkGH, file, datas, ext, vindex, groupdata,
                                   offset1, dirs,
                                   rl, ml, m, c, tl,
                                   coord_lower, coord_upper);

        } // for tl

      } // for c

      error_count += CloseFile (cctkGH, file);
      if (error_count > 0 and abort_on_io_errors) {
        CCTK_WARN (0, "Aborting simulation due to previous I/O errors");
      }

    } // for m
  }



  template<int outdim>
  bool IOHDF5<outdim>::DidOutput (const cGH* const cctkGH,
                                  const int vindex,
                                  const string basefilename,
                                  bool& is_new_file, bool& truncate_file)
  {
    DECLARE_CCTK_PARAMETERS;

    typedef std::map<string, vector<vector<vector<int> > > > filelist;
    static filelist created_files;

    filelist::iterator thisfile = created_files.find (basefilename);
    is_new_file = thisfile == created_files.end();
    truncate_file = is_new_file and IO_TruncateOutputFiles (cctkGH);

    if (is_new_file) {
      const int numelems =
        one_file_per_group ? CCTK_NumGroups() : CCTK_NumVars();
      vector<vector<vector<int> > > last_outputs; // [ml][rl][var]
      last_outputs.resize(mglevels);
      for (int ml=0; ml<mglevels; ++ml) {
        last_outputs.at(ml).resize (maxreflevels);
        for (int rl=0; rl<maxreflevels; ++rl) {
          last_outputs.at(ml).at(rl).resize
            (numelems, cctkGH->cctk_iteration - 1);
        }
      }
      // TODO: this makes a copy of last_outputs, which is expensive;
      // change this to use a reference instead
      thisfile = created_files.insert
        (thisfile, filelist::value_type (basefilename, last_outputs));
      assert (thisfile != created_files.end());
    }

    // Check if this variable has been output already during this
    // iteration
    const int elem =
      one_file_per_group ? CCTK_GroupIndexFromVarI(vindex) : vindex;
    int& last_output = thisfile->second.at(mglevel).at(reflevel).at(elem);
    if (last_output == cctkGH->cctk_iteration) {
      // has already been output during this iteration
      char* const fullname = CCTK_FullName (vindex);
      CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Skipping output for variable '%s', because this variable "
                  "has already been output during the current iteration -- "
                  "probably via a trigger during the analysis stage",
                  fullname);
      free (fullname);
      return true;
    }
    assert (last_output < cctkGH->cctk_iteration);
    last_output = cctkGH->cctk_iteration;

    return false;
  }



  CCTK_REAL io_files;
  CCTK_REAL io_bytes;

  template<int outdim>
  int IOHDF5<outdim>::OpenFile (const cGH* const cctkGH,
                                const int m,
                                const int vindex,
                                const int numvars,
                                const string alias,
                                const string basefilename,
                                const vect<int,outdim>& dirs,
                                const bool is_new_file,
                                const bool truncate_file,
                                hid_t& file)
  {
    DECLARE_CCTK_PARAMETERS;

    int error_count = 0;

    BeginTimingIO (cctkGH);
    io_files = 0;
    io_bytes = 0;

    if (dist::rank() == ioproc) {

      const int grouptype = CCTK_GroupTypeFromVarI(vindex);
      assert (grouptype >= 0);

      // Invent a file name
      ostringstream filenamebuf;
      filenamebuf << basefilename;
      if (maps > 1 and grouptype == CCTK_GF) {
        filenamebuf << "." << m;
      }
      filenamebuf << ".";
      for (int d=0; d<outdim; ++d) {
        const char* const coords = "xyzd";
        filenamebuf << coords[dirs[d]];
      }
      filenamebuf << ".h5";
      // we need a persistent temporary here
      const string filenamestr = filenamebuf.str();
      const char* const filename = filenamestr.c_str();

      // Open the file
      bool file_exists = false;
      if (not truncate_file) {
        H5E_BEGIN_TRY {
          file_exists = H5Fis_hdf5(filename) > 0;
        } H5E_END_TRY;
      }

      if (truncate_file or not file_exists) {
        HDF5_ERROR(file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT,
                                    H5P_DEFAULT));
        // write metadata information
        error_count +=
          WriteMetadata(cctkGH, nioprocs, vindex, numvars, false, file);
      } else {
        HDF5_ERROR (file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT));
      }
      io_files += 1;

    } // if on the I/O processor

    return error_count;
  }

  template<int outdim>
  int IOHDF5<outdim>::CloseFile (const cGH* const cctkGH,
                                 hid_t& file)
  {
    DECLARE_CCTK_PARAMETERS;

    int error_count = 0;

    if (dist::rank() == ioproc) {
      if (file >= 0) {
        HDF5_ERROR(H5Fclose(file));
      }
    }
    if (nioprocs > 1) {
      CCTK_REAL local[2], global[2];
      local[0] = io_files;
      local[1] = io_bytes;
      MPI_Allreduce (local, global, 2, dist::datatype (local[0]), MPI_SUM, dist::comm());
      io_files = global[0];
      io_bytes = global[1];
    }

    EndTimingIO (cctkGH, io_files, io_bytes, true);

    return error_count;
  }



  // Check whether this output direction has been requested
  template<int outdim>
  bool IOHDF5<outdim>::DirectionIsRequested (const vect<int,outdim>& dirs)
  {
    DECLARE_CCTK_PARAMETERS;

    switch (outdim) {

    case 0:
      // Output is always requested (if switched on)
      return true;

    case 1:
      switch (dirs[0]) {
      case 0: return out1D_x;
      case 1: return out1D_y;
      case 2: return out1D_z;
      case 3: return out1D_d;
      default: assert (0);
      }

    case 2:
      if (dirs[0]==0 and dirs[1]==1) return out2D_xy;
      if (dirs[0]==0 and dirs[1]==2) return out2D_xz;
      if (dirs[0]==1 and dirs[1]==2) return out2D_yz;
      assert (0);

//    case 3:
//      // Output is always requested (if switched on)
//      return true;

    default:
      assert (0);
      // Prevent compiler warning about missing return statement
      return false;
    }
  }



  // Get the region that should be output, in terms of grid points;
  // this is the offset perpendicular to the output hyperslab
  template<int outdim>
  ivect IOHDF5<outdim>::GetOutputOffset (const cGH* const cctkGH, const int m,
                                         const vect<int,outdim>& dirs)
  {
    DECLARE_CCTK_PARAMETERS;

    // Default is zero
    ivect offset (0);

    switch (outdim) {

    case 0:
      // 0D output
      offset[0] = GetGridOffset (cctkGH, m, 1,
                                 "out0D_point_xi", /*"out_point_xi"*/ NULL,
                                 "out0D_point_x",  /*"out_point_x"*/  NULL,
                                 /*out_point_x*/ 0.0);
      offset[1] = GetGridOffset (cctkGH, m, 2,
                                 "out0D_point_yi", /*"out_point_yi"*/ NULL,
                                 "out0D_point_y",  /*"out_point_y"*/  NULL,
                                 /*out_point_y*/ 0.0);
      offset[2] = GetGridOffset (cctkGH, m, 3,
                                 "out0D_point_zi", /*"out_point_zi"*/ NULL,
                                 "out0D_point_z",  /*"out_point_z"*/  NULL,
                                 /*out_point_z*/ 0.0);
      break;

    case 1:
      // 1D output
      switch (dirs[0]) {
      case 0:
        offset[1] = GetGridOffset (cctkGH, m, 2,
                                   "out1D_xline_yi", "out_xline_yi",
                                   "out1D_xline_y",  "out_xline_y",
                                   out_xline_y);
        offset[2] = GetGridOffset (cctkGH, m, 3,
                                   "out1D_xline_zi", "out_xline_zi",
                                   "out1D_xline_z",  "out_xline_z",
                                   out_xline_z);
        break;
      case 1:
        offset[0] = GetGridOffset (cctkGH, m, 1,
                                   "out1D_yline_xi", "out_yline_xi",
                                   "out1D_yline_x",  "out_yline_x",
                                   out_yline_x);
        offset[2] = GetGridOffset (cctkGH, m, 3,
                                   "out1D_yline_zi", "out_yline_zi",
                                   "out1D_yline_z",  "out_yline_z",
                                   out_yline_z);
        break;
      case 2:
        offset[0] = GetGridOffset (cctkGH, m, 1,
                                   "out1D_zline_xi", "out_zline_xi",
                                   "out1D_zline_x",  "out_zline_x",
                                   out_zline_x);
        offset[1] = GetGridOffset (cctkGH, m, 2,
                                   "out1D_zline_yi", "out_zline_yi",
                                   "out1D_zline_y",  "out_zline_y",
                                   out_zline_y);
        break;
      case 3:
        // the diagonal: we don't care about the offset
        break;
      default:
        assert (0);
      }
      break;

    case 2:
      // 2D output
      if (dirs[0]==0 and dirs[1]==1) {
        offset[2] = GetGridOffset (cctkGH, m, 3,
                                   "out2D_xyplane_zi", "out_xyplane_zi",
                                   "out2D_xyplane_z",  "out_xyplane_z",
                                   out_xyplane_z);
      } else if (dirs[0]==0 and dirs[1]==2) {
        offset[1] = GetGridOffset (cctkGH, m, 2,
                                   "out2D_xzplane_yi", "out_xzplane_yi",
                                   "out2D_xzplane_y",  "out_xzplane_y",
                                   out_xzplane_y);
      } else if (dirs[0]==1 and dirs[1]==2) {
        offset[0] = GetGridOffset (cctkGH, m, 1,
                                   "out2D_yzplane_xi", "out_yzplane_xi",
                                   "out2D_yzplane_x",  "out_yzplane_x",
                                   out_yzplane_x);
      } else {
        assert (0);
      }
      break;

//    case 3:
//      // 3D output: the offset does not matter
//      break;

    default:
      assert (0);
    }

    return offset;
  }



  // Omit symmetry and ghost zones if requested
  ibbox GetOutputBBox (const cGH* const cctkGH,
                       const int group,
                       const int m, const int c,
                       const ibbox& ext)
  {
    DECLARE_CCTK_PARAMETERS;

    const int groupdim = CCTK_GroupDimI(group);
    assert (groupdim >= 0);
    const int grouptype = CCTK_GroupTypeI(group);
    assert (grouptype >= 0);

    // TODO: This is a bit ad hoc
    CCTK_INT symtable;
    if (grouptype == CCTK_GF and groupdim == cctkGH->cctk_dim) {
      symtable = SymmetryTableHandleForGrid (cctkGH);
      if (symtable < 0) CCTK_WARN (0, "internal error");
    } else {
      symtable = SymmetryTableHandleForGI (cctkGH, group);
      if (symtable < 0) CCTK_WARN (0, "internal error");
    }

    CCTK_INT symbnd[2*dim];
    int const ierr = Util_TableGetIntArray
      (symtable, 2*groupdim, symbnd, "symmetry_handle");
    if (ierr != 2*groupdim) CCTK_WARN (0, "internal error");

    bool is_symbnd[2*dim];
    for (int d=0; d<2*groupdim; ++d) {
      is_symbnd[d] = symbnd[d] >= 0;
    }

    ivect lo = ext.lower();
    ivect hi = ext.upper();
    const ivect str = ext.stride();

    const b2vect obnds       = vhh.at(m)->outer_boundaries(reflevel,c);
    const i2vect ghost_width = arrdata.at(group).at(m).dd->ghost_width;

    for (int d=0; d<groupdim; ++d) {
      bool const output_lower_ghosts =
        obnds[0][d]
        ? (is_symbnd[2*d]
           ? output_symmetry_points
           : out3D_outer_ghosts)
        : out3D_ghosts;
      bool const output_upper_ghosts =
        obnds[1][d]
        ? (is_symbnd[2*d+1]
           ? output_symmetry_points
           : out3D_outer_ghosts)
        : out3D_ghosts;

      if (not output_lower_ghosts) {
        lo[d] += ghost_width[0][d] * str[d];
      }
      if (not output_upper_ghosts) {
        hi[d] -= ghost_width[1][d] * str[d];
      }
    }

    return ibbox(lo,hi,str);
  }



  // Determine coordinates
  void GetCoordinates (const cGH* const cctkGH, const int m,
                       const cGroup& groupdata,
                       const ibbox& ext,
                       rvect& coord_lower, rvect& coord_upper)
  {
    rvect global_lower;
    rvect coord_delta;
    if (groupdata.grouptype == CCTK_GF) {
      rvect const cctk_origin_space = origin_space.at(m).at(mglevel);
      rvect const cctk_delta_space  = delta_space.at(m) * rvect(mglevelfact);
      for (int d=0; d<dim; ++d) {
        // lower boundary of Carpet's integer indexing
        global_lower[d] = cctk_origin_space[d];
        // grid spacing of Carpet's integer indexing
        coord_delta[d]  = (cctk_delta_space[d] /
                           vhh.at(m)->baseextents.at(0).at(0).stride()[d]);
      }
    } else {
      for (int d=0; d<dim; ++d) {
        global_lower[d] = 0.0;
        coord_delta[d]  = 1.0;
      }
    }

    coord_lower = global_lower + coord_delta * rvect(ext.lower());
    coord_upper = global_lower + coord_delta * rvect(ext.upper());
  }



  int GetGridOffset (const cGH* const cctkGH, const int m, const int dir,
                     const char* const iparam, const char* const iglobal,
                     const char* const cparam, const char* const cglobal,
                     const CCTK_REAL cfallback)
  {
    // First choice: explicit coordinate
    const int ncparam = CCTK_ParameterQueryTimesSet (cparam, CCTK_THORNSTRING);
    assert (ncparam >= 0);
    if (ncparam > 0) {
      int ptype;
      const CCTK_REAL* const pcoord
        = ((const CCTK_REAL*)CCTK_ParameterGet
           (cparam, CCTK_THORNSTRING, &ptype));
      assert (pcoord);
      const CCTK_REAL coord = *pcoord;
      assert (ptype == PARAMETER_REAL);
      return CoordToOffset (cctkGH, m, dir, coord, 0);
    }

    // Second choice: explicit index
    const int niparam = CCTK_ParameterQueryTimesSet (iparam, CCTK_THORNSTRING);
    assert (niparam >= 0);
    if (niparam > 0) {
      int ptype;
      const int* const pindex
        = (const int*)CCTK_ParameterGet (iparam, CCTK_THORNSTRING, &ptype);
      assert (pindex);
      const int index = *pindex;
      assert (ptype == PARAMETER_INT);
      return index;
    }

    // Third choice: explicit global coordinate
    const char* iothorn = CCTK_ImplementationThorn ("IO");
    assert (iothorn);
    if (cglobal) {
      const int ncglobal = CCTK_ParameterQueryTimesSet (cglobal, iothorn);
      assert (ncglobal >= 0);
      if (ncglobal > 0) {
        int ptype;
        const CCTK_REAL* const pcoord
          = (const CCTK_REAL*)CCTK_ParameterGet (cglobal, iothorn, &ptype);
        assert (pcoord);
        const CCTK_REAL coord = *pcoord;
        assert (ptype == PARAMETER_REAL);
        return CoordToOffset (cctkGH, m, dir, coord, 0);
      }
    }

    // Fourth choice: explicit global index
    if (iglobal) {
      const int niglobal = CCTK_ParameterQueryTimesSet (iglobal, iothorn);
      assert (niglobal >= 0);
      if (niglobal > 0) {
        int ptype;
        const int* const pindex
          = (const int*)CCTK_ParameterGet (iglobal, iothorn, &ptype);
        assert (pindex);
        const int index = *pindex;
        assert (ptype == PARAMETER_INT);
        return index;
      }
    }

    // Fifth choice: default coordinate
    return CoordToOffset (cctkGH, m, dir, cfallback, 0);
  }



  int CoordToOffset (const cGH* cctkGH, const int m, const int dir,
                     const CCTK_REAL coord, const int ifallback)
  {
    assert (m>=0 and m<Carpet::maps and dir>=1 and dir<=dim);

    assert (mglevel!=-1 and reflevel!=-1 and Carpet::map==-1);

    rvect const cctk_origin_space = origin_space.at(m).at(mglevel);
    rvect const cctk_delta_space  = delta_space.at(m) * rvect (mglevelfact);
    ivect const cctk_levfac = spacereffacts.at (reflevel);
    ibbox const & coarseext = vhh.at(m)->baseextents.at(mglevel).at(0       );
    ibbox const & baseext   = vhh.at(m)->baseextents.at(mglevel).at(reflevel);
    ivect const cctk_levoff      = baseext.lower() - coarseext.lower();
    ivect const cctk_levoffdenom = baseext.stride();

    const CCTK_REAL delta = cctk_delta_space[dir-1] / cctk_levfac[dir-1];
    const CCTK_REAL lower = cctk_origin_space[dir-1] + cctk_delta_space[dir-1] / cctk_levfac[dir-1] * cctk_levoff[dir-1] / cctk_levoffdenom[dir-1];

    const CCTK_REAL rindex = (coord - lower) / delta;
    int cindex = (int)floor(rindex + 0.75);

    return cindex;
  }



  // Output
  template<int outdim>
  int IOHDF5<outdim>::WriteHDF5 (const cGH* cctkGH,
                                 hid_t& file,
                                 vector<const gdata*> const gfdatas,
                                 const bbox<int,dim>& gfext,
                                 const int vi,
                                 const cGroup& groupdata,
                                 const vect<int,dim>& org,
                                 const vect<int,outdim>& dirs,
                                 const int rl,
                                 const int ml,
                                 const int m,
                                 const int c,
                                 const int tl,
                                 const vect<CCTK_REAL,dim>& coord_lower,
                                 const vect<CCTK_REAL,dim>& coord_upper)
  {
    DECLARE_CCTK_PARAMETERS;

    assert (outdim<=dim);

    bool all_on_root = true;
    for (size_t n=0; n<gfdatas.size(); ++n) {
      all_on_root &= gfdatas.at(n)->proc() == ioproc;
    }

    // boolean that says if we are doing 1D-diagonal output
    // This is not beautiful, but works for the moment
    bool const diagonal_output = outdim == 1 and dirs[0] == 3;

    // Check whether the output bbox overlaps
    // with the extent of the data to be output
// FIXME: move this check up in the call stack
    bool output_bbox_overlaps_data_extent;
    if (not diagonal_output) {

      const vect<int,outdim> lo = gfext.lower()[dirs];
      const vect<int,outdim> up = gfext.upper()[dirs];
      const vect<int,outdim> str = gfext.stride()[dirs];
      const bbox<int,outdim> ext(lo,up,str);

      // Check whether the output origin is contained in the extent
      // of the data that should be output
      ivect org1(org);
      for (int d=0; d<outdim; ++d) org1[dirs[d]] = ext.lower()[d];
      output_bbox_overlaps_data_extent = gfext.contains(org1);

    } else {

      gh const & hh = *vhh.at(m);
      ibbox const & base = hh.baseextents.at(mglevel).at(reflevel);

      assert (base.stride()[0] ==  base.stride()[1]
              and base.stride()[0] == base.stride()[2]);

      // Check if any point on the diagonal is in our gf's extent
      output_bbox_overlaps_data_extent = false;
      for (int i=maxval(base.lower());
           i<=minval(base.upper()); i+=base.stride()[0]) {

        ivect const pos = ivect(i,i,i);
        output_bbox_overlaps_data_extent |= gfext.contains(pos);
      }
    }
    // Shortcut if there is nothing to output
    if (not output_bbox_overlaps_data_extent) {
      return 0;
    }

    int error_count = 0;
    if (all_on_root) {
      // output on processor 0

      if (dist::rank() == ioproc) {

        ostringstream datasetname_suffix;
        datasetname_suffix << " it=" << cctkGH->cctk_iteration << " tl=" << tl;
        if (mglevels > 1) datasetname_suffix << " ml=" << ml;
        if (groupdata.grouptype == CCTK_GF) {
          if (maps > 1) datasetname_suffix << " m="  << m;
          datasetname_suffix << " rl=" << rl;
        }
        if (groupdata.grouptype == CCTK_GF or
            groupdata.disttype  != CCTK_DISTRIB_CONSTANT) {
          datasetname_suffix << " c=" << c;
        }

        // enable compression and checksums if requested
        hid_t plist;
        HDF5_ERROR(plist = H5Pcreate(H5P_DATASET_CREATE));
        if (compression_level) {
          HDF5_ERROR(H5Pset_deflate(plist, compression_level));
        }
        if (use_checksums) {
          HDF5_ERROR(H5Pset_filter(plist, H5Z_FILTER_FLETCHER32, 0, 0, 0));
        }

        // enable datatype conversion if requested
        const hid_t mem_type   =
          CCTKtoHDF5_Datatype(cctkGH, groupdata.vartype,false);
        const hid_t slice_type =
          CCTKtoHDF5_Datatype(cctkGH, groupdata.vartype,out_single_precision);

        if (not diagonal_output) { // not outputting the diagonal

          const vect<int,outdim> lo = gfext.lower()[dirs];
          const vect<int,outdim> up = gfext.upper()[dirs];
          const vect<int,outdim> str = gfext.stride()[dirs];
          const bbox<int,outdim> ext(lo,up,str);

          // Check whether the output origin is contained in the extent
          // of the data that should be output
          ivect org1(org);
          for (int d=0; d<outdim; ++d) org1[dirs[d]] = ext.lower()[d];
          assert (gfext.contains(org1));

          // HDF5 wants ranks to be >= 1
          const int rank = outdim > 0 ? outdim : 1;
          vector<hsize_t> mem_shape(dim);
          vector<hsize_t> slice_shape(rank, 1);
          for (int d = 0; d < dim; d++) {
            mem_shape[dim-1-d] = gfext.shape()[d] / gfext.stride()[d];
            if (d < outdim) {
              slice_shape[outdim-1-d] = ext.shape()[d] / ext.stride()[d];
            }
          }

          ivect slice_lower(org - gfext.lower());
          for (int d = 0; d < outdim; d++) {
            slice_lower[dirs[d]] = 0;
          }
          ivect slice_upper(slice_lower);
          for (int d = 0; d < outdim; d++) {
            slice_upper[dirs[d]] = ext.upper()[d] - ext.lower()[d];
          }
          slice_lower /= gfext.stride();
          slice_upper /= gfext.stride();

          slice_start_size_t slice_start[dim];
          hsize_t slice_count[dim];
          for (int d = 0; d < dim; d++) {
            slice_start[dim-1-d] = slice_lower[d];
            slice_count[dim-1-d] = slice_upper[d] - slice_lower[d] + 1;
          }
          if (compression_level or use_checksums) {
            HDF5_ERROR(H5Pset_chunk(plist, slice_shape.size(),&slice_shape[0]));
          }
          hid_t slice_space, mem_space;
          HDF5_ERROR(slice_space = H5Screate_simple(slice_shape.size(),
                                                    &slice_shape[0], NULL));
          HDF5_ERROR(mem_space = H5Screate_simple(mem_shape.size(),
                                                  &mem_shape[0], NULL));
          HDF5_ERROR(H5Sselect_hyperslab(mem_space, H5S_SELECT_SET,
                                         slice_start, NULL, slice_count, NULL));

          vector<int> iorigin(rank, 0);
          vector<double> delta(rank, 0), origin(rank, 0);
          for (int d = 0; d < outdim; d++) {
            assert(gfext.upper()[dirs[d]] - gfext.lower()[dirs[d]] >= 0);
            iorigin[d] = ext.lower()[d];
            delta[d] = 0;
            origin[d] = coord_lower[dirs[d]];
            if (gfext.upper()[dirs[d]] - gfext.lower()[dirs[d]] > 0) {
              delta[d] = (coord_upper[dirs[d]] - coord_lower[dirs[d]]) /
                         (gfext.upper()[dirs[d]] - gfext.lower()[dirs[d]]) * gfext.stride()[dirs[d]];
              origin[d] += (org1[dirs[d]] - gfext.lower()[dirs[d]]) * delta[d];
              iorigin[d] /= gfext.stride()[dirs[d]];
            }
          }

          // now loop over all variables
          for (size_t n = 0; n < gfdatas.size(); n++) {

            // create a unique name for this variable's dataset
            char *fullname = CCTK_FullName(vi + n);
            string datasetname(fullname);
            datasetname.append(datasetname_suffix.str());
  
            // remove an already existing dataset of the same name
            if (slice_requests[vi + n]->check_exist) {
              H5E_BEGIN_TRY {
                H5Gunlink(file, datasetname.c_str());
              } H5E_END_TRY;
            }

            // write the dataset
            hid_t dataset;
            HDF5_ERROR(dataset = H5Dcreate(file, datasetname.c_str(),
                                           slice_type, slice_space, plist));
            HDF5_ERROR(H5Dwrite(dataset, mem_type, mem_space, H5S_ALL,
                                H5P_DEFAULT, gfdatas[n]->storage()));
            error_count += AddSliceAttributes(cctkGH, fullname, rl,
                                              origin, delta, iorigin,
                                              dataset);
            HDF5_ERROR(H5Dclose(dataset));
            free(fullname);
  
            io_bytes += H5Sget_simple_extent_npoints(slice_space) *
                        H5Tget_size(slice_type);
          }

          HDF5_ERROR(H5Sclose(mem_space));
          HDF5_ERROR(H5Sclose(slice_space));

        } else { // taking care of the diagonal

          const ibbox& ext = gfext;

          gh const & hh = *vhh.at(m);
          ibbox const & base = hh.baseextents.at(mglevel).at(reflevel);

          assert (base.stride()[0] == base.stride()[1] and
                  base.stride()[0] == base.stride()[2]);

          // count the number of points on the diagonal
          hsize_t npoints = 0;
          for (int i =  maxval(base.lower());
                   i <= minval(base.upper());
                   i += base.stride()[0]) {
            if(gfext.contains(i)) {
              npoints++;
            }
          }
          assert(npoints > 0);

          // allocate a contiguous buffer for the diagonal points
          char* buffer = new char[CCTK_VarTypeSize(groupdata.vartype) *
                                  npoints * gfdatas.size()];

          // copy diagonal points into contiguous buffer
          hsize_t offset = 0;
          for (int i =  maxval(base.lower());
                   i <= minval(base.upper());
                   i += base.stride()[0]) {
            ivect const pos = ivect(i,i,i);
            if(gfext.contains(pos)) {
              for (size_t n = 0; n < gfdatas.size(); n++) {
                switch (groupdata.vartype) {
#define TYPECASE(N,T)                                                     \
                  case N: { T* typed_buffer = (T*) buffer;                \
                            typed_buffer[offset + n*npoints] =            \
                              (*(const data<T>*)gfdatas[n])[pos];         \
                            break;                                        \
                          }
#include "carpet_typecase.hh"
#undef TYPECASE
                }
              }
              offset++;
            }
          }
          assert(offset == npoints);

          if (compression_level or use_checksums) {
            HDF5_ERROR(H5Pset_chunk(plist, 1, &npoints));
          }
          hid_t slice_space;
          HDF5_ERROR(slice_space = H5Screate_simple(1, &npoints, NULL));

          // loop over all variables and write out diagonals
          for (size_t n = 0; n < gfdatas.size(); n++) {

            // create a unique name for this variable's dataset
            char *fullname = CCTK_FullName(vi + n);
            string datasetname(fullname);
            free(fullname);
            datasetname.append(datasetname_suffix.str());

            // remove an already existing dataset of the same name
            if (slice_requests[vi + n]->check_exist) {
              H5E_BEGIN_TRY {
                H5Gunlink(file, datasetname.c_str());
              } H5E_END_TRY;
            }

            // write the dataset
            hid_t dataset;
            HDF5_ERROR(dataset = H5Dcreate(file, datasetname.c_str(),
                                           slice_type, slice_space, plist));
            HDF5_ERROR(H5Dwrite(dataset, mem_type,
                                H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                buffer + n*npoints*gfdatas.size()));
            HDF5_ERROR(H5Dclose(dataset));

            io_bytes +=
              H5Sget_simple_extent_npoints(slice_space) * H5Tget_size(slice_type);
          }

          HDF5_ERROR(H5Sclose(slice_space));

          // release contiguous buffer
          delete[] buffer;

        } // if(not diagonal_output)

        HDF5_ERROR(H5Pclose(plist));

      } // if(dist::rank() == ioproc)
    } else {
      // copy to processor 0 and output there

      mempool pool;
      vector<gdata*> tmps (gfdatas.size());
      for (size_t n=0; n<gfdatas.size(); ++n) {
        tmps.at(n) = gfdatas.at(n)->make_typed (vi, error_centered, op_sync);
        size_t const memsize =
          tmps.at(n)->allocsize (gfdatas.at(n)->extent(), ioproc);
        void * const memptr = pool.alloc (memsize);
        tmps.at(n)->allocate(gfdatas.at(n)->extent(), ioproc, memptr, memsize);
      }
      for (comm_state state; not state.done(); state.step()) {
        for (size_t n=0; n<gfdatas.size(); ++n) {
          tmps.at(n)->copy_from (state, gfdatas.at(n), gfdatas.at(n)->extent());
        }
      }
      vector<const gdata*> ctmps (gfdatas.size());
      for (size_t n=0; n<gfdatas.size(); ++n) {
        ctmps.at(n) = tmps.at(n);
      }
      WriteHDF5 (cctkGH, file, ctmps, gfext, vi, groupdata, org, dirs,
                 rl, ml, m, c, tl, coord_lower, coord_upper);
      for (size_t n=0; n<gfdatas.size(); ++n) {
        delete tmps.at(n);
      }

    }

    return error_count;
  }



  // Explicit instantiation for all slice output dimensions
  template class IOHDF5<0>;
  template class IOHDF5<1>;
  template class IOHDF5<2>;
//  template class IOHDF5<3>;

} // namespace CarpetIOHDF5
