// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOASCII/src/ioascii.cc,v 1.1 2001/03/01 13:40:10 eschnett Exp $

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "Carpet/CarpetLib/src/data.hh"
#include "Carpet/CarpetLib/src/gdata.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/Carpet/src/carpet.hh"

#include "ioascii.hh"



namespace Carpet {
  
int CarpetIOASCIIStartup()
{
  IOASCII<1>::Startup();
  IOASCII<2>::Startup();
  IOASCII<3>::Startup();
  
  return 0;
}



// Definition of static members
template<int outdim> int IOASCII<outdim>::GHExtension;
template<int outdim> int IOASCII<outdim>::IOMethod;
template<int outdim> vector<bool> IOASCII<outdim>::do_truncate;
template<int outdim> vector<int> IOASCII<outdim>::last_output;



template<int outdim>
int IOASCII<outdim>::Startup()
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
void* IOASCII<outdim>::SetupGH (tFleshConfig *fc, int convLevel, cGH *cgh)
{
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_CreateDirectory (0755, outdir);
  
  // Truncate all files if this is not a restart
  do_truncate.resize(CCTK_NumVars(), ! IOUtil_RestartFromRecovery(cgh));
  
  // No iterations have yet been output
  last_output.resize(CCTK_NumVars(), -1);
  
  // We register only once, ergo we get only one handle.  We store
  // that statically, so there is no need to pass it to Cactus.
  return 0;
}



template<int outdim>
int IOASCII<outdim>::OutputGH (cGH* cgh) {
  for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
    if (TimeToOutput(cgh, vindex)) {
      TriggerOutput(cgh, vindex);
    }
  }
  return 0;
}

template<int outdim>
int IOASCII<outdim>::OutputVarAs (cGH* cgh, const char* varname,
				const char* alias) {
  DECLARE_CCTK_PARAMETERS;
  
  const int n = CCTK_VarIndex(varname);
  assert (n>=0 && n<CCTK_NumVars());
  const int group = CCTK_GroupIndexFromVarI (n);
  assert (group>=0 && group<(int)Carpet::gfdata.size());
  const int n0 = CCTK_FirstVarIndexI(group);
  assert (n0>=0 && n0<CCTK_NumVars());
  const int var = n - n0;
  assert (var>=0 && var<CCTK_NumVars());
  const int tl = max(0, CCTK_NumTimeLevelsFromVarI(n) - 1);
  
  switch (CCTK_GroupTypeI(group)) {
    
  case CCTK_SCALAR: {
    // Don't output scalars
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Cannout output variable \"%s\" because it is a scalar",
		varname);
    return 0;
  }
    
  case CCTK_ARRAY:
  case CCTK_GF: {
    
    // Check for storage
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot output variable \"%s\" because it has no storage",
		  varname);
      return 0;
    }
    
    // Loop over all direction combinations
    vect<int,outdim> dirs;
    for (int d=0; d<outdim; ++d) dirs[d] = 0;
    
    for (;;) {
      
      // Output each combination only once
      bool ascending = true;
      for (int d1=0; d1<outdim; ++d1) {
	for (int d2=d1+1; d2<outdim; ++d2) {
	  ascending = ascending && dirs[d1] < dirs[d2];
	}
      }
      
      if (ascending) {
	
	// Invent a file name
	char filename[strlen(outdir)+strlen(alias)+100];
	sprintf (filename, "%s/%s.", outdir, alias);
	for (int d=0; d<outdim; ++d) {
	  assert (dirs[d]>=0 && dirs[d]<3);
	  sprintf (filename, "%s%c", filename, "xyz"[dirs[d]]);
	}
	sprintf (filename, "%s%c", filename, "lpv"[outdim-1]);
	
	// If this is the first time, then write a nice header
	if (do_truncate[n]) {
	  ofstream file(filename, ios::trunc);
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
	
	// Traverse all components on all refinement levels
  	assert (mglevel>=0);
	assert (reflevel==0);
	for (reflevel=0; reflevel<hh->reflevels(); ++reflevel) {
	  assert (component==-1);
	  for (component=0; component<hh->components(reflevel); ++component) {
	    
	    switch (CCTK_GroupTypeI(group)) {
	      
	    case CCTK_ARRAY:
	      assert (var < (int)Carpet::arrdata[group].data.size());
	      (*arrdata[group].data[var]) (tl, reflevel, component, mglevel)
		->write_ascii (filename, cgh->cctk_time, dirs,
			       tl, reflevel, component, mglevel);
	      break;
	      
	    case CCTK_GF:
	      assert (var < (int)Carpet::gfdata[group].data.size());
	      (*gfdata[group].data[var]) (tl, reflevel, component, mglevel)
		->write_ascii (filename, cgh->cctk_time, dirs,
			       tl, reflevel, component, mglevel);
	      break;
	      
	    default:
	      abort();
	    }
	    
	  }
	  component = -1;
	}
	reflevel = 0;
	
      } // if (ascending)
      
      // Next direction combination
      for (int d=0; d<outdim; ++d) {
	++dirs[d];
	if (dirs[d]<dim) goto notyetdone;
	dirs[d] = 0;
      }
      break;
      
    notyetdone: ;
      
    } // all directions
    
    break;
  }
    
  default:
    abort();
  }
  
  // Don't truncate again
  do_truncate[n] = false;
  
  return 0;
}



template<int outdim>
int IOASCII<outdim>::TimeToOutput (cGH* cgh, int vindex) {
  assert (vindex>=0 && vindex<(int)last_output.size());
  if (last_output[vindex] < cgh->cctk_iteration) {
    return 1;
  } else if (last_output[vindex] == cgh->cctk_iteration) {
    char* varname = CCTK_FullName(vindex);
    CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Skipping output for variable \"%s\", because this variable "
		"has already been output during the current iteration -- "
		"probably via a trigger during the analysis stage",
		varname);
    free (varname);
    return 0;
  } else {
    abort();
  }
}



template<int outdim>
int IOASCII<outdim>::TriggerOutput (cGH* cgh, int vindex) {
  assert (vindex>=0 && vindex<CCTK_NumVars());
  
  char* varname = CCTK_FullName(vindex);
  int retval = OutputVarAs (cgh, varname, CCTK_VarName(vindex));
  free (varname);
  
  last_output[vindex] = cgh->cctk_iteration;
  
  return retval;
}



// Explicit instantiation for all output dimensions
template IOASCII<1>;
template IOASCII<2>;
template IOASCII<3>;

} // namespace Carpet
