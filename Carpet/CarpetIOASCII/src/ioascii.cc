#include <alloca.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "Carpet/CarpetLib/src/data.hh"
#include "Carpet/CarpetLib/src/gdata.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/ggf.hh"

#include "carpet.hh"

#include "ioascii.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOASCII/src/ioascii.cc,v 1.18 2001/11/02 16:05:04 schnetter Exp $";



using namespace std;
using namespace Carpet;



static bool CheckForVariable (const cGH* const cgh,
			      const char* const varlist, const int vindex);
static void SetFlag (int index, const char* optstring, void* arg);



int CarpetIOASCIIStartup()
{
  CarpetIOASCII<1>::Startup();
  CarpetIOASCII<2>::Startup();
  CarpetIOASCII<3>::Startup();
  
  return 0;
}



// Definition of static members
template<int outdim> int CarpetIOASCII<outdim>::GHExtension;
template<int outdim> int CarpetIOASCII<outdim>::IOMethod;
template<int outdim> vector<bool> CarpetIOASCII<outdim>::do_truncate;
template<int outdim> vector<vector<int> > CarpetIOASCII<outdim>::last_output;



template<int outdim>
int CarpetIOASCII<outdim>::Startup()
{
  char msg[1000];
  sprintf (msg, "AMR %dD ASCII I/O provided by CarpetIOASCII", outdim);
  CCTK_RegisterBanner (msg);
  
  char extension_name[1000];
  sprintf (extension_name, "CarpetIOASCII_%dD", outdim);
  GHExtension = CCTK_RegisterGHExtension(extension_name);
  CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
  
  char method_name[1000];
  sprintf (method_name, "IOASCII_%dD", outdim);
  IOMethod = CCTK_RegisterIOMethod (method_name);
  CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
  
  return 0;
}



template<int outdim>
void* CarpetIOASCII<outdim>
::SetupGH (tFleshConfig* const fc, const int convLevel, const cGH* const cgh)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Truncate all files if this is not a restart
  do_truncate.resize(CCTK_NumVars(), true);
  
  // No iterations have yet been output
  last_output.resize(maxreflevels);
  for (int rl=0; rl<maxreflevels; ++rl) {
    last_output[rl].resize(CCTK_NumVars(), INT_MIN);
  }
  
  // We register only once, ergo we get only one handle.  We store
  // that statically, so there is no need to pass anything to Cactus.
  return 0;
}



template<int outdim>
int CarpetIOASCII<outdim>
::OutputGH (cGH* const cgh)
{
  for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
    if (TimeToOutput(cgh, vindex)) {
      TriggerOutput(cgh, vindex);
    }
  }
  return 0;
}



template<int outdim>
int CarpetIOASCII<outdim>
::OutputVarAs (cGH* const cgh,
	       const char* const varname, const char* const alias)
{
  DECLARE_CCTK_PARAMETERS;
  
  const int n = CCTK_VarIndex(varname);
  assert (n>=0 && n<CCTK_NumVars());
  const int group = CCTK_GroupIndexFromVarI (n);
  assert (group>=0 && group<(int)Carpet::arrdata.size());
  const int n0 = CCTK_FirstVarIndexI(group);
  assert (n0>=0 && n0<CCTK_NumVars());
  const int var = n - n0;
  assert (var>=0 && var<CCTK_NumVars());
  const int tl = 0;
  
  // Check for storage
  if (! CCTK_QueryGroupStorageI((cGH*)cgh, group)) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Cannot output variable \"%s\" because it has no storage",
		varname);
    return 0;
  }
  
  // Create the output directory
  const char* myoutdir = GetStringParameter("outdir%dD", outdir);
  if (CCTK_MyProc(cgh)==0) {
    CCTK_CreateDirectory (0755, myoutdir);
  }
  
  // Loop over all direction combinations
  vect<int,outdim> dirs;
  for (int d=0; d<outdim; ++d) dirs[d] = 0;
  
  bool done;
  do {
    
    // Output each combination only once
    bool ascending = true;
    for (int d1=0; d1<outdim; ++d1) {
      for (int d2=d1+1; d2<outdim; ++d2) {
	ascending = ascending && dirs[d1] < dirs[d2];
      }
    }
    
    // Skip output if the dimensions are not ascending
    if (ascending) {
      
      // If this output is desired at all
      bool desired;
      switch (outdim) {
      case 1:
	switch (dirs[0]) {
	case 0:
	  desired = out1D_x;
	  break;
	case 1:
	  desired = out1D_y;
	  break;
	case 2:
	  desired = out1D_z;
	  break;
	default:
	  abort();
	}
	break;
      case 2:
	if (dirs[0]==0 && dirs[1]==1) {
	  desired = out2D_xy;
	} else if (dirs[0]==0 && dirs[1]==2) {
	  desired = out2D_xz;
	} else if (dirs[0]==1 && dirs[1]==2) {
	  desired = out2D_yz;
	} else {
	  abort();
	}
	break;
      case 3:
	// Output is always desired (if switched on)
	desired = true;
	break;
      default:
	abort();
      }
      
      // Skip output if not desired
      if (desired) {
	
	// Invent a file name
	char* const filename
	  = (char*)alloca(strlen(myoutdir)+strlen(alias)+100);
	sprintf (filename, "%s/%s.", myoutdir, alias);
	for (int d=0; d<outdim; ++d) {
	  assert (dirs[d]>=0 && dirs[d]<3);
	  const char* const coords = "xyz";
	  sprintf (filename, "%s%c", filename, coords[dirs[d]]);
	}
	const char* const suffixes = "lpv";
	sprintf (filename, "%s%c", filename, suffixes[outdim-1]);
	
	// If this is the first time, then write a nice header on
	// the root processor
	if (do_truncate[n]) {
	  if (CCTK_MyProc(cgh)==0) {
	    struct stat fileinfo;
	    if (! IOUtil_RestartFromRecovery(cgh)
		|| stat(filename, &fileinfo)!=0) {
	      ofstream file(filename, ios::out | ios::trunc);
	      assert (file.good());
	      file << "# " << varname;
	      for (int d=0; d<outdim; ++d) {
		file << " " << "xyz"[dirs[d]];
	      }
	      file << " (" << alias << ")" << endl;
	      file << "#" << endl;
	      file.close();
	      assert (file.good());
	    }
	  }
	}
	
	assert (outdim <= CCTK_GroupDimI(group));
	
	// Find the output offset
	vect<int,dim> offset(0);
	switch (outdim) {
	case 1:
	  switch (dirs[0]) {
	  case 0:
	    offset[1] = GetGridOffset (cgh, 2,
				       "out%dD_xline_yi", "out_xline_yi",
				       "out%dD_xline_y",  "out_xline_y",
				       out_xline_y);
	    offset[2] = GetGridOffset (cgh, 3,
				       "out%dD_xline_zi", "out_xline_zi",
				       "out%dD_xline_z",  "out_xline_z",
				       out_xline_z);
	    break;
	  case 1:
	    offset[0] = GetGridOffset (cgh, 1,
				       "out%dD_yline_xi", "out_yline_xi",
				       "out%dD_yline_x",  "out_yline_x",
				       out_yline_x);
	    offset[2] = GetGridOffset (cgh, 3,
				       "out%dD_yline_zi", "out_yline_zi",
				       "out%dD_yline_z",  "out_yline_z",
				       out_yline_z);
	    break;
	  case 2:
	    offset[0] = GetGridOffset (cgh, 1,
				       "out%dD_zline_xi", "out_zline_xi",
				       "out%dD_zline_x",  "out_zline_x",
				       out_zline_x);
	    offset[1] = GetGridOffset (cgh, 2,
				       "out%dD_zline_yi", "out_zline_yi",
				       "out%dD_zline_y",  "out_zline_y",
				       out_zline_y);
	    break;
	  default:
	    abort();
	  }
	  break;
	case 2:
	  if (dirs[0]==0 && dirs[1]==1) {
	    offset[2] = GetGridOffset
	      (cgh, 3,
	       "out%dD_xyplane_zi", "out_xyplane_zi",
	       "out%dD_xyplane_z",  "out_xyplane_z",
	       out_xyplane_z);
	  } else if (dirs[0]==0 && dirs[1]==2) {
	    offset[1] = GetGridOffset
	      (cgh, 2,
	       "out%dD_xzplane_yi", "out_xzplane_yi",
	       "out%dD_xzplane_y",  "out_xzplane_y",
	       out_xzplane_y);
	  } else if (dirs[0]==1 && dirs[1]==2) {
	    offset[0] = GetGridOffset
	      (cgh, 1,
	       "out%dD_yzplane_xi", "out_yzplane_xi",
	       "out%dD_yzplane_x",  "out_yzplane_x",
	       out_yzplane_x);
	  } else {
	    abort();
	  }
	  break;
	case 3:
	  // The offset doesn't matter in this case
	  break;
	default:
	  abort();
	}
	
	// Traverse all components on this refinement and multigrid
	// level
	BEGIN_COMPONENT_LOOP(cgh) {
	  
	  generic_gf<dim>* ff = 0;
	  
	  assert (var < (int)arrdata[group].data.size());
	  ff = (generic_gf<dim>*)arrdata[group].data[var];
	  
	  const generic_data<dim>* const data
	    = (*ff) (tl, reflevel, component, mglevel);
	  const bbox<int,dim> ext = data->extent();
	  const vect<int,dim> offset1 = offset * ext.stride();
	  
	  data->write_ascii (filename, cgh->cctk_iteration, offset1, dirs,
			     tl, reflevel, component, mglevel);
	  
	} END_COMPONENT_LOOP(cgh);
	
	// Append EOL after every complete set of components
	if (CCTK_MyProc(cgh)==0) {
	  ofstream file(filename, ios::app);
	  assert (file.good());
	  file << endl;
	  file.close();
	  assert (file.good());
	}
	
      } // if (desired)
      
    } // if (ascending)
    
      // Next direction combination
    done = true;
    for (int d=0; d<outdim; ++d) {
      ++dirs[d];
      if (dirs[d]<CCTK_GroupDimI(group)) {
	done = false;
	break;
      }
      dirs[d] = 0;
    }
    
  } while (! done);		// all directions
  
  // Don't truncate again
  do_truncate[n] = false;
  
  return 0;
}



template<int outdim>
int CarpetIOASCII<outdim>
::TimeToOutput (cGH* const cgh, const int vindex)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (vindex>=0 && vindex<CCTK_NumVars());
  
  const int myoutevery = GetIntParameter("out%dD_every", out_every);
  
  if (myoutevery < 0) {
    // Nothing should be output at all
    return 0;
  }
  
  if (cgh->cctk_iteration % myoutevery != 0) {
    // Nothing should be output during this iteration
    return 0;
  }
  
  if (! CheckForVariable(cgh, GetStringParameter("out%dD_vars", ""), vindex)) {
    // This variable should not be output
    return 0;
  }
  
  if (last_output[reflevel][vindex] == cgh->cctk_iteration) {
    // Has already been output during this iteration
    char* varname = CCTK_FullName(vindex);
    CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Skipping output for variable \"%s\", because this variable "
		"has already been output during the current iteration -- "
		"probably via a trigger during the analysis stage",
		varname);
    free (varname);
    return 0;
  }
  
  assert (last_output[reflevel][vindex] < cgh->cctk_iteration);
  
  // Should be output during this iteration
  return 1;
}



template<int outdim>
int CarpetIOASCII<outdim>
::TriggerOutput (cGH* const cgh, const int vindex)
{
  assert (vindex>=0 && vindex<CCTK_NumVars());
  
  char* varname = CCTK_FullName(vindex);
  const int retval = OutputVarAs (cgh, varname, CCTK_VarName(vindex));
  free (varname);
  
  last_output[reflevel][vindex] = cgh->cctk_iteration;
  
  return retval;
}



template<int outdim>
int CarpetIOASCII<outdim>
::GetGridOffset (const cGH* cgh, int dir,
		 const char* itempl, const char* iglobal,
		 const char* ctempl, const char* cglobal,
		 const CCTK_REAL cfallback)
{
  // First choice: explicit coordinate
  char cparam[1000];
  sprintf (cparam, ctempl, outdim);
  if (CCTK_ParameterQueryTimesSet (cparam, CCTK_THORNSTRING) > 0) {
    int ptype;
    const CCTK_REAL* const pcoord
      = (const CCTK_REAL*)CCTK_ParameterGet (cparam, CCTK_THORNSTRING, &ptype);
    assert (pcoord);
    const CCTK_REAL coord = *pcoord;
    assert (ptype == PARAMETER_REAL);
    return CoordToOffset (cgh, dir, coord);
  }
  
  // Second choice: explicit index
  char iparam[1000];
  sprintf (iparam, itempl, outdim);
  if (CCTK_ParameterQueryTimesSet (iparam, CCTK_THORNSTRING) > 0) {
    int ptype;
    const int* const pindex
      = (const int*)CCTK_ParameterGet (iparam, CCTK_THORNSTRING, &ptype);
    assert (pindex);
    const int index = *pindex;
    assert (ptype == PARAMETER_INT);
    return index;
  }
  
  // Third choice: explicit global coordinate
  if (CCTK_ParameterQueryTimesSet (cglobal, "IO") > 0) {
    int ptype;
    const CCTK_REAL* const pcoord
      = (const CCTK_REAL*)CCTK_ParameterGet (cglobal, "IO", &ptype);
    assert (pcoord);
    const CCTK_REAL coord = *pcoord;
    assert (ptype == PARAMETER_REAL);
    return CoordToOffset (cgh, dir, coord);
  }
  
  // Fourth choice: explicit global index
  if (CCTK_ParameterQueryTimesSet (iglobal, "IO") > 0) {
    int ptype;
    const int* const pindex
      = (const int*)CCTK_ParameterGet (iglobal, "IO", &ptype);
    assert (pindex);
    const int index = *pindex;
    assert (ptype == PARAMETER_INT);
    return index;
  }
  
  // Fifth choice: default coordinate
  return CoordToOffset (cgh, dir, cfallback);
}



template<int outdim>
int CarpetIOASCII<outdim>
::CoordToOffset (const cGH* cgh, int dir, CCTK_REAL coord)
{
  CCTK_REAL lower, upper;
  CCTK_CoordRange (cgh, &lower, &upper, dir, 0, "cart3d");
  
  const int npoints = cgh->cctk_gsh[dir-1];
  
  const CCTK_REAL rindex = (coord - lower) / (upper - lower) * (npoints-1);
  const int cindex = (int)floor(rindex + 0.5 + 1e-6);
  assert (cindex>=0 && cindex<npoints);
  
  return cindex;
}



template<int outdim>
const char* CarpetIOASCII<outdim>
::GetStringParameter (const char* const parametertemplate,
		      const char* const fallback)
{
  char parametername[1000];
  sprintf (parametername, parametertemplate, outdim);
  if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
    int ptype;
    const char* const* const ppval
      = (const char* const*)CCTK_ParameterGet (parametername, CCTK_THORNSTRING,
					 &ptype);
    assert (ppval);
    const char* const pval = *ppval;
    assert (ptype == PARAMETER_STRING);
    return pval;
  }
  
  return fallback;
}



template<int outdim>
int CarpetIOASCII<outdim>
::GetIntParameter (const char* const parametertemplate, const int fallback)
{
  char parametername[1000];
  sprintf (parametername, parametertemplate, outdim);
  if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
    int ptype;
    const int* const ppval
      = (const int*)CCTK_ParameterGet (parametername, CCTK_THORNSTRING,
				       &ptype);
    assert (ppval);
    const int pval = *ppval;
    assert (ptype == PARAMETER_INT);
    return pval;
  }
  
  return fallback;
}



bool CheckForVariable (const cGH* const cgh,
		       const char* const varlist, const int vindex)
{
  const int numvars = CCTK_NumVars();
  assert (vindex>=0 && vindex<numvars);
  
  bool* const flags = (bool*)alloca(numvars * sizeof(bool));
  
  for (int i=0; i<numvars; ++i) {
    flags[i] = false;
  }
  
  CCTK_TraverseString (varlist, SetFlag, flags, CCTK_GROUP_OR_VAR);
  
  return flags[vindex];
}

void SetFlag (int index, const char* optstring, void* arg)
{
  ((bool*)arg)[index] = true;
}



// Explicit instantiation for all output dimensions
template class CarpetIOASCII<1>;
template class CarpetIOASCII<2>;
template class CarpetIOASCII<3>;
