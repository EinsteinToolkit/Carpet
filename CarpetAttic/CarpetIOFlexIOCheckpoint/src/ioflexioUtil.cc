#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Version.h"

#include "AMRwriter.hh"
#include "AmrGridReader.hh"
#ifdef HDF4
#  include "HDFIO.hh"
#endif
#ifdef HDF5
#  include "H5IO.hh"
#endif
#include "IEEEIO.hh"
#include "IO.hh"

// Hack to stop FlexIO type clash
#undef BYTE
#undef CHAR

//#include "CactusBase/IOUtil/src/ioGH.h"
//#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
//#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"
#include "ioflexio.hh"


namespace CarpetIOFlexIOUtil {

  using namespace std;
  using namespace Carpet;
  using namespace CarpetIOFlexIO;

IObase::DataType FlexIODataType (int cctk_type){
  //we need this to have the FlexIO data types on hand
  //for WriteGFAs

    int retval;

    switch (cctk_type)
      {
      case CCTK_VARIABLE_CHAR:   retval = FLEXIO_CHAR; break;
      case CCTK_VARIABLE_INT:    retval = FLEXIO_INT; break;
      case CCTK_VARIABLE_REAL:   retval = FLEXIO_REAL; break;
#ifdef CCTK_INT2
      case CCTK_VARIABLE_INT2:   retval = FLEXIO_INT2; break;
#endif
#ifdef CCTK_INT4
      case CCTK_VARIABLE_INT4:   retval = FLEXIO_INT4; break;
#endif
#ifdef CCTK_INT8
      case CCTK_VARIABLE_INT8:   retval = FLEXIO_INT8; break;
#endif
#ifdef CCTK_REAL4
      case CCTK_VARIABLE_REAL4:  retval = FLEXIO_REAL4; break;
#endif
#ifdef CCTK_REAL8
      case CCTK_VARIABLE_REAL8:  retval = FLEXIO_REAL8; break;
#endif
#ifdef CCTK_REAL16
      case CCTK_VARIABLE_REAL16: retval = FLEXIO_REAL16; break;
#endif
	
      default: CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
			   "Unsupported CCTK variable datatype %d", cctk_type);
	retval = -1;
	break;
      }

    return (IObase::DataType)retval;
  }


void DumpCommonAttributes (const cGH *cgh, IObase* writer, ioRequest* request)
{
  int dim, vdim;
  CCTK_INT attr_int;
  CCTK_REAL *attr_real;
  char coord_system_name[20];
  DECLARE_CCTK_PARAMETERS

  /* attributes describing the variable */
  char* groupname = CCTK_GroupNameFromVarI (request->vindex);


//	
  char *name = CCTK_FullName (request->vindex);
  writer->writeAttribute("name",IObase::Char,strlen(name)+1,name);
  free(name); 
  
  //CCTK_VInfo (CCTK_THORNSTRING, "DUMPATTRIB");
  //fprintf(stderr,"\nattrib %s\n",groupname);
  writer->writeAttribute("groupname",IObase::String,strlen(groupname)+1,groupname);
  free (groupname);

  CCTK_INT attr_int = CCTK_GroupTypeFromVarI (request->vindex);
  writer->writeAttribute("grouptype",FlexIODataType(CCTK_VARIABLE_INT),1,&attr_int);

  writer->writeAttribute("reflevel",FlexIODataType(CCTK_VARIABLE_INT),1,&reflevel);
  writer->writeAttribute("component",FlexIODataType(CCTK_VARIABLE_INT),1,&component);
  writer->writeAttribute("mglevel",FlexIODataType(CCTK_VARIABLE_INT),1,&mglevel);

  attr_int = CCTK_MaxTimeLevelsVI (request->vindex);
  writer->writeAttribute("ntimelevels",FlexIODataType(CCTK_VARIABLE_INT),1,&attr_int);

  writer->writeAttribute("timelevel",FlexIODataType(CCTK_VARIABLE_INT),1,&request->timelevel);

  writer->writeAttribute("global_size",FlexIODataType(CCTK_VARIABLE_INT),request->hdim,request->hsize);

  // already dumped by amrwriter if we are dealing with a grid function or a grid array:
  if (CCTK_GroupTypeFromVarI (request->vindex) == CCTK_SCALAR)
    writer->writeAttribute("time",FlexIODataType(CCTK_VARIABLE_REAL),1,&cgh->cctk_time);


  /* attributes describing the underlying grid
     These are only stored for CCTK_GF variables if they are associated
     with coordinates. */
  /* FIXME: This is hardcoded for cartesian coordinate systems.
            A better solution would be to be able to query the coordinate
            system which is associated with the variable. */
  vdim = CCTK_GroupDimFromVarI (request->vindex);
  sprintf (coord_system_name, "cart%dd", vdim);

  /* this is already dumped by amrwriter or */
#if 0
  /*
  if (CCTK_GroupTypeFromVarI (request->vindex) == CCTK_GF &&
      CCTK_CoordSystemHandle (coord_system_name) >= 0)
  {
    attr_real = (CCTK_REAL*) malloc (3 * vdim * sizeof (CCTK_REAL));
    for (dim = 0; dim < vdim; dim++)
    {
      CCTK_CoordRange (cgh, &attr_real[dim], &attr_real[dim + vdim], dim + 1,
                       NULL, coord_system_name);

      attr_real[dim + 0*vdim] +=
        request->origin[dim] * cgh->cctk_delta_space[dim];
      attr_real[dim + 2*vdim] =
        cgh->cctk_delta_space[dim] * request->downsample[dim];
      attr_real[dim + 1*vdim] = attr_real[dim + 0*vdim] +
        ((request->extent[dim] + request->downsample[dim]-1) /
         request->downsample[dim] - 1) * attr_real[dim + 2*vdim];
    }

    writer->writeAttribute ("origin", FlexIODataType(CCTK_VARIABLE_REAL), vdim, attr_real);
    writer->writeAttribute ("min_ext", FlexIODataType(CCTK_VARIABLE_REAL), vdim, attr_real);
    writer->writeAttribute ("max_ext", FlexIODataType(CCTK_VARIABLE_REAL), vdim, attr_real+vdim);
    writer->writeAttribute ("delta", FlexIODataType(CCTK_VARIABLE_REAL), vdim, attr_real+2*vdim);
    free (attr_real);
  } 
  */
#endif

  }

} // namespace CarpetIOFlexIOUtil









