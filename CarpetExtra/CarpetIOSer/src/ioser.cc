#include <alloca.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iomanip>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "Carpet/CarpetLib/src/data.hh"
#include "Carpet/CarpetLib/src/gdata.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

#include "ioser.hh"

#include "bbhutil.h"
extern "C" {
int vsxynt(const char *name, double t, double *x, double *y, int n);
}


static const char* rcsid = "$Header:$";

// That's a hack
namespace Carpet {
  void UnsupportedVarType (const int vindex);
}


using namespace std;
using namespace Carpet;



static bool CheckForVariable (const cGH* const cgh,
			      const char* const varlist, const int vindex);
static void SetFlag (int index, const char* optstring, void* arg);

//void assign_ser_output_element (double& range, const ivect& index);
//void assign_ser_output_element (double& range, const vect<int,dim> index);

template<int D,int DD>
static void write_ser ( const char *gfname,
                        const generic_data<D>* const gfdata,
                        const bbox<int,D>& gfext,
                        const int vi,
                        const int time,
                        const vect<int,D>& org,
                        const vect<int,DD>& dirs,
                        const int tl,
                        const int rl,
                        const int c,
                        const int ml,
                        const CCTK_REAL coord_time,
                        const vect<CCTK_REAL,D>& coord_lower,
                        const vect<CCTK_REAL,D>& coord_upper,
                        const int strip_left, const int strip_right);


int CarpetIOSerStartup()
{
  CarpetIOSer<1>::Startup();
  CarpetIOSer<2>::Startup();
  CarpetIOSer<3>::Startup();
  
  return 0;
}



// Definition of static members
template<int outdim> int CarpetIOSer<outdim>::GHExtension;
template<int outdim> int CarpetIOSer<outdim>::IOMethod;
template<int outdim> vector<bool> CarpetIOSer<outdim>::do_truncate;
template<int outdim> vector<vector<int> > CarpetIOSer<outdim>::last_output;



template<int outdim>
int CarpetIOSer<outdim>::Startup()
{
  char msg[1000];
  sprintf (msg, "AMR %dD Ser I/O provided by CarpetIOSer", outdim);
  CCTK_RegisterBanner (msg);
  
  char extension_name[1000];
  sprintf (extension_name, "CarpetIOSer_%dD", outdim);
  GHExtension = CCTK_RegisterGHExtension(extension_name);
  CCTK_RegisterGHExtensionSetupGH (GHExtension, SetupGH);
  
  char method_name[1000];
  sprintf (method_name, "IOSer_%dD", outdim);
  IOMethod = CCTK_RegisterIOMethod (method_name);
  CCTK_RegisterIOMethodOutputGH (IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput (IOMethod, TriggerOutput);
  
  return 0;
}



template<int outdim>
void* CarpetIOSer<outdim>
::SetupGH (tFleshConfig* const fc, const int convLevel, cGH* const cgh)
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
int CarpetIOSer<outdim>
::OutputGH (const cGH* const cgh)
{
  for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
    if (TimeToOutput(cgh, vindex)) {
      TriggerOutput(cgh, vindex);
    }
  }
  return 0;
}



template<int outdim>
int CarpetIOSer<outdim>
::OutputVarAs (const cGH* const cgh,
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
  if (! CCTK_QueryGroupStorageI(cgh, group)) {
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
	
	ofstream file;
	
	if (CCTK_MyProc(cgh)==0) {
	  // If this is the first time, then write a nice header on
	  // the root processor
	  if (do_truncate[n]) {
	    struct stat fileinfo;
	    if (! IOUtil_RestartFromRecovery(cgh)
		|| stat(filename, &fileinfo)!=0) {
	      file.open (filename, ios::out | ios::trunc);
	      assert (file.good());
	      file << "# " << varname;
	      for (int d=0; d<outdim; ++d) {
		file << " " << "xyz"[dirs[d]];
	      }
	      file << " (" << alias << ")" << endl;
	      file << "#" << endl;
	      assert (file.good());
	    }
	  }
	  if (! file.is_open()) {
	    file.open (filename, ios::app);
	    assert (file.good());
	  }
	  file << setprecision(15);
	  assert (file.good());
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
	default:
	  abort();
	}
	
	// Traverse all components on this refinement and multigrid
	// level
	BEGIN_COMPONENT_LOOP(cgh) {
	  
	  const generic_gf<dim>* ff = 0;
	  
	  assert (var < (int)arrdata[group].data.size());
	  ff = (generic_gf<dim>*)arrdata[group].data[var];
	  
	  const generic_data<dim>* const data
	    = (*ff) (tl, reflevel, component, mglevel);
	  const bbox<int,dim> ext = data->extent();
	  const vect<int,dim> offset1 = offset * ext.stride();
	  
	  // coordinates
	  const double coord_time = cgh->cctk_time;
	  vect<double,dim> global_lower, global_upper;
	  for (int d=0; d<dim; ++d) {
	    const int ierr = CCTK_CoordRange
	      (cgh, &global_lower[d], &global_upper[d], d+1, 0, "cart3d");
	    if (ierr<0) {
	      global_lower[d] = 0;
	      global_upper[d] = 0;
	    }
	  }
	  const vect<double,dim> global_extent (hh->baseextent.upper() - hh->baseextent.lower() + hh->baseextent.stride() * (dd->lghosts + dd->ughosts));
	  vect<double,dim> coord_delta;
	  for (int d=0; d<dim; ++d) {
	    if (global_extent[d] != 0) {
	      coord_delta[d] = (global_upper[d] - global_lower[d]) / global_extent[d];
	    } else {
	      coord_delta[d] = 0;
	    }
	  }
	  // Note: don't permute the "coord_delta" and "data->extent().lower()"
	  // (you'll pick up the integer operator* then)
 	  const vect<double,dim> coord_lower = global_lower + coord_delta * vect<double,dim>(data->extent().lower());
 	  const vect<double,dim> coord_upper = global_lower + coord_delta * vect<double,dim>(data->extent().upper());
	  
          char *gfname = new char[strlen(alias)+strlen(out1D_tag)+1];
          sprintf(gfname,"%s%s",alias,out1D_tag);
	  write_ser (gfname, data, ext, n, cgh->cctk_iteration, offset1, dirs,
			     tl, reflevel, component, mglevel,
			     coord_time, coord_lower, coord_upper,
                             strip_left, strip_right);
          delete [] gfname;
	  
	  // Append EOL after every component
	  if (CCTK_MyProc(cgh)==0) {
	    if (separate_components) {
	      assert (file.good());
	      file << endl;
	    }
	  }
	  assert (file.good());
	  
	} END_COMPONENT_LOOP(cgh);
	
	// Append EOL after every complete set of components
	if (CCTK_MyProc(cgh)==0) {
	  if (separate_grids) {
	    assert (file.good());
	    file << endl;
	  }
	  file.close();
	  assert (file.good());
	}
	
	assert (! file.is_open());
	
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
int CarpetIOSer<outdim>
::TimeToOutput (const cGH* const cgh, const int vindex)
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
int CarpetIOSer<outdim>
::TriggerOutput (const cGH* const cgh, const int vindex)
{
  assert (vindex>=0 && vindex<CCTK_NumVars());
  
  char* varname = CCTK_FullName(vindex);
  const int retval = OutputVarAs (cgh, varname, CCTK_VarName(vindex));
  free (varname);
  
  last_output[reflevel][vindex] = cgh->cctk_iteration;
  
  return retval;
}



template<int outdim>
int CarpetIOSer<outdim>
::GetGridOffset (const cGH* const cgh, const int dir,
		 const char* const itempl, const char* const iglobal,
		 const char* const ctempl, const char* const cglobal,
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
    return CoordToOffset (cgh, dir, coord, 0);
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
    return CoordToOffset (cgh, dir, coord, 0);
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
  return CoordToOffset (cgh, dir, cfallback, 0);
}



template<int outdim>
int CarpetIOSer<outdim>
::CoordToOffset (const cGH* cgh, const int dir, const CCTK_REAL coord,
		 const int ifallback)
{
  assert (dir>=1 && dir<=dim);
  
  CCTK_REAL lower, upper;
  CCTK_CoordRange (cgh, &lower, &upper, dir, 0, "cart3d");
  
  const int npoints = cgh->cctk_gsh[dir-1];
  
  const CCTK_REAL rindex = (coord - lower) / (upper - lower) * (npoints-1);
  const int levfac = cgh->cctk_levfac[dir-1];
  int cindex = (int)floor(rindex / levfac + 0.5 + 1e-6) * levfac;
  
  if (cindex<0 || cindex>=npoints) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"The specified coordinate value %g is not within the grid range [%g,%g]",
		coord, lower, upper);
    cindex = ifallback;
  }
  
  assert (cindex>=0 && cindex<npoints);
  
  return cindex;
}



template<int outdim>
const char* CarpetIOSer<outdim>
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
int CarpetIOSer<outdim>
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


// Output
/*
template<class T, int D>
void data<T,D>::assign_ser_output_element (double& range, const vect<int,dim> index)
{
  range = (*this)[index];
}
*/


// Output
template<int D,int DD>
void write_ser ( const char *gfname,
                 const generic_data<D>* const gfdata,
                 const bbox<int,D>& gfext,
                 const int vi,
                 const int time,
                 const vect<int,D>& org,
                 const vect<int,DD>& dirs,
                 const int tl,
                 const int rl,
                 const int c,
                 const int ml,
                 const CCTK_REAL coord_time,
                 const vect<CCTK_REAL,D>& coord_lower,
                 const vect<CCTK_REAL,D>& coord_upper,
                 const int strip_left, const int strip_right)
 
{
  assert (DD<=D);

  if (gfdata->proc()==0) {
    // output on processor 0

    int rank;
    MPI_Comm_rank (dist::comm, &rank); 
    if (rank == 0) {


      assert (D>=1 && D<=3); 
      assert (DD==1);
      const char* const coords = "xyz";

      const vect<int,DD> lo = gfext.lower()[dirs];
      const vect<int,DD> up = gfext.upper()[dirs];
      const vect<int,DD> str = gfext.stride()[dirs];
      const bbox<int,DD> ext(lo,up,str);

      // Check whether the output origin is contained in the extent of
      // the data
      vect<int,dim> org1(org);
      for (int d=0; d<DD; ++d) org1[dirs[d]] = ext.lower()[d];
      if (gfext.contains(org1)) {

        const int nelems = ext.num_points();
        double *domain = new double[nelems];
        double *range = new double[nelems];

        int elem = 0;
        for (bbox<int,DD>::iterator it=ext.begin(); it!=ext.end(); ++it) {
          vect<int,dim> index(org);
          for (int d=0; d<DD; ++d) index[dirs[d]] = (*it)[d];
          const int d = dirs[0];
          if (gfext.upper()[d] - gfext.lower()[d] != 0) {
            domain[elem] =  coord_lower[d] + (index[d] - gfext.lower()[d]) *
(coord_upper[d] - coord_lower[d]) / (gfext.upper()[d] - gfext.lower()[d]);          } else {
            domain[elem] = coord_lower[d];
          }

//          assign_ser_output_element (range[elem], index);
//          range[elem] = (*(data<T,D>*)gfdata)[index];
//            elem++;

//          cout << "CarpetIOSer, before switch" << endl << flush;

          switch (CCTK_VarTypeI(vi)) {
#define TYPECASE(N,T)                           \
          case N:                               \
            range[elem] = (*(data<T,D>*)gfdata)[index]; \
            elem++;\
            break;
#include "Carpet/Carpet/src/typecase"
#undef TYPECASE
          default:  
            UnsupportedVarType(vi);
          }
 //         cout << "CarpetIOSer, after switch" << endl << flush;

        }
        assert(elem == nelems);

        char *window_name = new char[strlen(gfname)+100];
        sprintf(window_name,"%s_l%d.%c",gfname,rl,"xyz"[dirs[0]]);
        vsxynt(window_name, coord_time, domain+strip_left, range+strip_left, 
               nelems-strip_left-strip_right);

        delete [] window_name;
        delete [] domain;
        delete [] range;

      } // if extent().contains(org1)

    }

  } else {
    // copy to processor 0 and output there

/*
    generic_data* const tmp = make_typed();
    tmp->allocate(extent(), 0);
    tmp->copy_from (this, extent());
    tmp->write_ser (gfname, time, org, dirs, tl, rl, c, ml,
                      coord_time, coord_lower, coord_upper);
    delete tmp;  
*/

    generic_data<D>* const tmp = gfdata->make_typed();
    tmp->allocate(gfdata->extent(), 0);
    tmp->copy_from (gfdata, gfdata->extent());
    write_ser (gfname, tmp, gfext, vi, time, org, dirs, tl, rl, c, ml,
                      coord_time, coord_lower, coord_upper, strip_left, strip_right);
    delete tmp;



  }
}





// Explicit instantiation for all output dimensions
template class CarpetIOSer<1>;
template class CarpetIOSer<2>;
template class CarpetIOSer<3>;

template
void write_ser ( const char *gfname,
                 const generic_data<3>* const gfdata,
                 const bbox<int,3>& gfext,
                 const int vi,
                 const int time,
                 const vect<int,3>& org,
                 const vect<int,1>& dirs,
                 const int tl,
                 const int rl,
                 const int c,
                 const int ml,
                 const CCTK_REAL coord_time,
                 const vect<CCTK_REAL,3>& coord_lower,
                 const vect<CCTK_REAL,3>& coord_upper,
                 const int strip_left, const int strip_right);

