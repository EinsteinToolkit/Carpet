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

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5utils.cc,v 1.1 2004/03/08 09:43:41 cott Exp $";
  CCTK_FILEVERSION(Carpet_CarpetIOHDF5_iohdf5utils_cc);
}

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "iohdf5.hh"



namespace CarpetIOHDF5 {
  
  using namespace std;
  using namespace Carpet;
  
  const char* GetStringParameter (const char* const parametername,
				  const char* const fallback)
  {
    if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
      int ptype;
      const char* const* const ppval = (const char* const*)CCTK_ParameterGet
	(parametername, CCTK_THORNSTRING, &ptype);
      assert (ppval);
      const char* const pval = *ppval;
      assert (ptype == PARAMETER_STRING);
      return pval;
    }
    
    return fallback;
  }
  
  
  
  int GetIntParameter (const char* const parametername, int fallback)
  {
    if (CCTK_ParameterQueryTimesSet (parametername, CCTK_THORNSTRING) > 0) {
      int ptype;
      const int* const ppval = (const int*)CCTK_ParameterGet
	(parametername, CCTK_THORNSTRING, &ptype);
      assert (ppval);
      const int pval = *ppval;
      assert (ptype == PARAMETER_INT);
      return pval;
    }
    
    return fallback;
  }
  
  
  
  bool CheckForVariable (const cGH* const cctkGH,
			 const char* const varlist, const int vindex)
  {
    const int numvars = CCTK_NumVars();
    assert (vindex>=0 && vindex<numvars);
    
    vector<bool> flags(numvars);
    
    CCTK_TraverseString (varlist, SetFlag, &flags, CCTK_GROUP_OR_VAR);
    
    return flags.at(vindex);
  }
  
  void SetFlag (int index, const char* optstring, void* arg)
  {
    vector<bool>& flags = *(vector<bool>*)arg;
    flags.at(index) = true;
  }
  
  
  
  void WriteAttribute (const hid_t dataset, const char* const name, const int value)
  {
    WriteAttribute (dataset, name, &value, 1);
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const int* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    hsize_t shape[1];
    shape[0] = nvalues;
    const hid_t dataspace = nvalues==1 ? H5Screate (H5S_SCALAR) : H5Screate_simple (1, shape, NULL);
    assert (dataspace>=0);
    
    const hid_t datatype = H5T_NATIVE_INT;
    
    const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
    assert (attribute>=0);
    herr = H5Awrite (attribute, datatype, values);
    assert (!herr);
    herr = H5Aclose (attribute);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
  }
  
  
  
  void WriteAttribute (const hid_t dataset, const char* const name, const double value)
  {
    WriteAttribute (dataset, name, &value, 1);
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const double* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    hsize_t shape[1];
    shape[0] = nvalues;
    const hid_t dataspace = nvalues==1 ? H5Screate (H5S_SCALAR) : H5Screate_simple (1, shape, NULL);
    assert (dataspace>=0);
    
    const hid_t datatype = H5T_NATIVE_DOUBLE;
    
    const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
    assert (attribute>=0);
    herr = H5Awrite (attribute, datatype, values);
    assert (!herr);
    herr = H5Aclose (attribute);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
  }
  
  
  
  void WriteAttribute (const hid_t dataset, const char* const name, const char value)
  {
    WriteAttribute (dataset, name, &value, 1);
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const char* const values)
  {
    WriteAttribute (dataset, name, values, strlen(values));
  }
  
  void WriteAttribute (const hid_t dataset, const char* const name, const char* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t dataspace = H5Screate (H5S_SCALAR);
    assert (dataspace>=0);
    
    const hid_t datatype = H5Tcopy (H5T_C_S1);
    assert (datatype>=0);
    herr = H5Tset_size (datatype, nvalues);
    assert (!herr);
    
    const hid_t attribute = H5Acreate (dataset, name, datatype, dataspace, H5P_DEFAULT);
    assert (attribute>=0);
    herr = H5Awrite (attribute, datatype, values);
    assert (!herr);
    herr = H5Aclose (attribute);
    assert (!herr);
    
    herr = H5Tclose (datatype);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
  }
  
  
  
  int ReadAttribute (const hid_t dataset, const char* const name, int& value)
  {
    return ReadAttribute (dataset, name, &value, 1);
  }
  
  int ReadAttribute (const hid_t dataset, const char* name, int* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    hsize_t shape[1];
    if (rank==0) {
      shape[0] = 1;
    } else if (rank==1) {
      herr = H5Sget_simple_extent_dims (dataspace, shape, NULL);
      assert (!herr);
    } else {
      assert (0);
    }
    const int length = shape[0];
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (datatype != H5T_NATIVE_INT) return -100;
    
    vector<int> values1(length);
    
    herr = H5Aread (attribute, datatype, &values1.at(0));
    assert (!herr);
    
    for (int i=0; i<min(length, nvalues); ++i) {
      values[i] = values1[i];
    }
    
    herr = H5Tclose (datatype);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  
  
  int ReadAttribute (const hid_t dataset, const char* const name, double& value)
  {
    return ReadAttribute (dataset, name, &value, 1);
  }
  
  int ReadAttribute (const hid_t dataset, const char* const name, double* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    hsize_t shape[1];
    if (rank==0) {
      shape[0] = 1;
    } else if (rank==1) {
      herr = H5Sget_simple_extent_dims (dataspace, shape, NULL);
      assert (!herr);
    } else {
      assert (0);
    }
    const int length = shape[0];
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (datatype != H5T_NATIVE_DOUBLE) return -100;
    
    vector<double> values1(length);
    
    herr = H5Aread (attribute, datatype, &values1.at(0));
    assert (!herr);
    
    for (int i=0; i<min(length, nvalues); ++i) {
      values[i] = values1[i];
    }
    
    herr = H5Tclose (datatype);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  
  
  int ReadAttribute (const hid_t dataset, const char* const name, char& value)
  {
    return ReadAttribute (dataset, name, &value, 1);
  }
  
  int ReadAttribute (const hid_t dataset, const char* const name, char*& values)
  {
    assert (dataset>=0);
    assert (name);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    assert (rank==0);
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (H5Tget_class (datatype) != H5T_STRING) return -100;
    const int length = H5Tget_size (datatype);
    assert (length>=0);
    
    values = (char*) malloc (length+1);
    assert (values);
    
    herr = H5Aread (attribute, datatype, values);
    assert (!herr);
    values[length] = '\0';
    
    herr = H5Tclose (datatype);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  int ReadAttribute (const hid_t dataset, const char* const name, char* const values, const int nvalues)
  {
    assert (dataset>=0);
    assert (name);
    assert (values);
    assert (nvalues>=0);
    
    herr_t herr;
    
    const hid_t attribute = H5Aopen_name (dataset, name);
    if (attribute<0) return attribute;
    
    const hid_t dataspace = H5Aget_space (attribute);
    assert (dataspace>=0);
    
    hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
    assert (rank==0);
    
    const hid_t datatype = H5Aget_type (attribute);
    assert (datatype>=0);
    if (H5Tget_class (datatype) != H5T_STRING) return -100;
    const int length = H5Tget_size (datatype);
    assert (length>=0);
    
    vector<char> values1(length);
    
    herr = H5Aread (attribute, datatype, &values1.at(0));
    assert (!herr);
    
    for (int i=0; i<min(length, nvalues); ++i) {
      values[i] = values1[i];
    }
    
    herr = H5Tclose (datatype);
    assert (!herr);
    
    herr = H5Sclose (dataspace);
    assert (!herr);
    
    herr = H5Aclose (attribute);
    assert (!herr);
    
    return length;
  }
  
  
  
} // namespace CarpetIOHDF5
