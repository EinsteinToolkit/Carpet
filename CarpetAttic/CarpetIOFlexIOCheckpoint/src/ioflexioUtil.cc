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
    int tl;
    CCTK_INT attr_int,dimscalar;
    DECLARE_CCTK_PARAMETERS;

    /* attributes describing the variable */
    
    char *name = CCTK_FullName (request->vindex);
    WriteAttribute(writer,"name",name);
    free(name); 
    
    char* groupname = CCTK_GroupNameFromVarI (request->vindex);
    WriteAttribute(writer,"groupname",groupname);
    free (groupname);
  
    WriteAttribute(writer,"grouptype",CCTK_GroupTypeFromVarI (request->vindex));
    WriteAttribute(writer,"reflevel",reflevel);
    WriteAttribute(writer,"component",component);
    WriteAttribute(writer,"mglevel",mglevel);
    

    WriteAttribute (writer,"ntimelevels",CCTK_MaxTimeLevelsVI (request->vindex));

    // lets get the correct Carpet time level (which is the (-1) * timelevel):
    if (request->timelevel==0)
      tl = 0;
    else
    tl = - request->timelevel;
    WriteAttribute (writer, "timelevel", tl);
    
    WriteAttribute (writer, "carpet_flexio_version", 1);
    WriteAttribute (writer, "cctk_dim", cgh->cctk_dim);
    WriteAttribute (writer, "cctk_iteration", cgh->cctk_iteration);
    WriteAttribute (writer, "cctk_gsh", cgh->cctk_gsh, dim);
    WriteAttribute (writer, "cctk_lsh", cgh->cctk_lsh, dim);
    WriteAttribute (writer, "cctk_lbnd", cgh->cctk_lbnd, dim);
    WriteAttribute (writer, "cctk_delta_time", cgh->cctk_delta_time);
    WriteAttribute (writer, "cctk_delta_space", cgh->cctk_delta_space, dim);
    WriteAttribute (writer, "cctk_origin_space", cgh->cctk_origin_space, dim);
    WriteAttribute (writer, "cctk_bbox", cgh->cctk_bbox, 2*dim);
    WriteAttribute (writer, "cctk_levfac", cgh->cctk_levfac, dim);
    WriteAttribute (writer, "cctk_levoff", cgh->cctk_levoff, dim);
    WriteAttribute (writer, "cctk_levoffdenom", cgh->cctk_levoffdenom, dim);
    WriteAttribute (writer, "cctk_timefac", cgh->cctk_timefac);
    WriteAttribute (writer, "cctk_convlevel", cgh->cctk_convlevel);
    WriteAttribute (writer, "cctk_nghostzones", cgh->cctk_nghostzones, dim);
    WriteAttribute (writer, "cctk_time", cgh->cctk_time);
  }

  
  void WriteAttribute (IObase* writer, const char* name, int value)
  {
    WriteAttribute (writer, name, &value, 1);
  }
  
  void WriteAttribute (IObase* writer, const char* name,
                       const int* values, int nvalues)
  {
    assert (writer);
    assert (name);
    assert (values);
    vector<CCTK_INT4> values1(nvalues);
    for (int i=0; i<nvalues; ++i) {
      values1[i] = values[i];
    }
    writer->writeAttribute (name, IObase::Int32, nvalues, &values1[0]);
  }
  
  void WriteAttribute (IObase* writer, const char* name, CCTK_REAL value)
  {
    WriteAttribute (writer, name, &value, 1);
  }
  
  void WriteAttribute (IObase* writer, const char* name,
                       const CCTK_REAL* values, int nvalues)
  {
    assert (writer);
    assert (name);
    assert (values);
    vector<CCTK_REAL8> values1(nvalues);
    for (int i=0; i<nvalues; ++i) {
      values1[i] = values[i];
    }
    writer->writeAttribute (name, IObase::Float64, nvalues, &values1[0]);
  }

  void WriteAttribute (IObase* writer, const char* name,
                       const char* valuestring)
  {
    assert (writer);
    assert (name);
    assert (valuestring);
    writer->writeAttribute (name, IObase::String, strlen(valuestring)+1, valuestring);
  }
  
} // namespace CarpetIOFlexIOUtil









