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
static const char* rcsid = "$Header:$";
CCTK_FILEVERSION(Carpet_CarpetIOHDF5_Utils_cc);
}

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "bbox.hh"
#include "data.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "CarpetIOHDF5.hh"


namespace CarpetIOHDF5 {

using namespace std;
using namespace Carpet;



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
  
  //    cout << "reading int attribute " << name << endl;

  hsize_t rank = H5Sget_simple_extent_ndims (dataspace);
  hsize_t shape[1];
  if (rank==0) {
    shape[0] = 1;
  } else if (rank==1) {
    herr = H5Sget_simple_extent_dims (dataspace, shape, NULL);
    assert (herr >= 0);
  } else {
    assert (0);
  }
  const int length = shape[0];
  
  const hid_t datatype = H5Aget_type (attribute);
  assert (datatype>=0);

  assert(H5Tequal(datatype, H5T_NATIVE_INT));
  
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
    rank = H5Sget_simple_extent_dims (dataspace, shape, NULL);
    assert (rank == 1);
  } else {
    assert (0);
  }
  const int length = shape[0];
  
  const hid_t datatype = H5Aget_type (attribute);
  assert (datatype>=0);
  assert(H5Tequal(datatype, H5T_NATIVE_DOUBLE));
  
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
  assert (H5Tget_class (datatype) == H5T_STRING);
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
  assert(H5Tget_class (datatype) == H5T_STRING); 
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

herr_t DatasetCounter(hid_t group_id, const char *member_name, void *operator_data)
  /* Counts datasets. Used by GetnDatasets; straight from John Shalf's FlexIO library */
{
  int *count = (int*)operator_data;
  H5G_stat_t objinfo;
  // request info about the type of objects in root group
  if(H5Gget_objinfo(group_id,member_name,1 /* follow links */,&objinfo)<0) {
    return 0; 
  }
  // only count objects that are datasets (not subgroups)
  if(objinfo.type==H5G_DATASET) {
    (*count)++;
  }
return 0;
}


int GetnDatasets(const hid_t reader)
{
  //this is straight from John Shalf's FlexIO library

  int count=0;
  int idx=0;
  while(H5Giterate(reader, /* hid_t loc_id, */
		     "/", /*const char *name, */
		     &idx, /* int *idx, */
		     DatasetCounter,
		     &count)<0){}
  return count;
}

struct H5IO_getname_t {
  //this is straight from John Shalf's FlexIO library
  int index,count;
  char *name;
};


herr_t GetName(hid_t group_id, const char *member_name, void *operator_data)
{
  //this is straight from John Shalf's FlexIO library
  H5IO_getname_t *getn = (H5IO_getname_t*)operator_data;
  // check type first (only respond if it is a dataset)
  H5G_stat_t objinfo;
  // request info about the type of objects in root group
  if(H5Gget_objinfo(group_id,
		      member_name,
		    1 /* follow links */,
		      &objinfo)<0) return 0; // error (probably bad symlink)
  // only count objects that are datasets (not subgroups)
  if(objinfo.type!=H5G_DATASET)
    return 0; // do not increment count if it isn't a dataset.
  if(getn->index==getn->count){
    strcpy(getn->name,member_name);
    return 1; // success
  }
  getn->count++;
  return 0;
}


void GetDatasetName(const hid_t reader, const int _index,  char *name) {
  //this is straight from John Shalf's FlexIO library
  H5IO_getname_t getn;
  int idx=_index;
  getn.index=_index; getn.name=name; getn.count=_index;
  while(H5Giterate(reader, /* hid_t loc_id, */
		     "/", /*const char *name, */
		     &idx, /* int *idx, */
		     GetName,
		     &getn)<0){}
}

hid_t h5DataType (const cGH* const cctkGH, int cctk_type, int single_precision)
{
  hid_t retval;  

  CarpetIOHDF5GH *myGH;
  myGH = (CarpetIOHDF5GH *) CCTK_GHExtension (cctkGH, "CarpetIOHDF5");

  // this is adapted from Thomas Radke's IOHDF5Util. Thanks, Thomas!

  switch (cctk_type)
  {
    case CCTK_VARIABLE_CHAR:      retval = HDF5_CHAR; break;

    case CCTK_VARIABLE_INT:       retval = HDF5_INT;
#ifdef CCTK_INT2
                                  if (single_precision)
                                  {
                                    retval = HDF5_INT2;
                                  }
#endif
                                  break;

    case CCTK_VARIABLE_REAL:      retval = HDF5_REAL;
#ifdef CCTK_REAL4
                                  if (single_precision)
                                  {
                                    retval = HDF5_REAL4;
                                  }
#endif
                                  break;

    case CCTK_VARIABLE_COMPLEX:   retval = myGH->HDF5_COMPLEX;
#ifdef CCTK_REAL4
                                  if (single_precision)
                                  {
                                    retval = retval = myGH->HDF5_COMPLEX8;
                                  }
#endif
                                  break;

#ifdef CCTK_INT1
    case CCTK_VARIABLE_INT1:      retval = HDF5_INT1; break;
#endif
#ifdef CCTK_INT2
    case CCTK_VARIABLE_INT2:      retval = HDF5_INT2; break;
#endif
#ifdef CCTK_INT4
    case CCTK_VARIABLE_INT4:      retval = HDF5_INT4; break;
#endif
#ifdef CCTK_INT8
    case CCTK_VARIABLE_INT8:      retval = HDF5_INT8; break;
#endif
#ifdef CCTK_REAL4
    case CCTK_VARIABLE_REAL4:     retval = HDF5_REAL4; break;
    case CCTK_VARIABLE_COMPLEX8:  retval = myGH->HDF5_COMPLEX8; break;
#endif
#ifdef CCTK_REAL8
    case CCTK_VARIABLE_REAL8:     retval = HDF5_REAL8; break;
    case CCTK_VARIABLE_COMPLEX16: retval = myGH->HDF5_COMPLEX16; break;
#endif
#ifdef CCTK_REAL16
    case CCTK_VARIABLE_REAL16:    retval = HDF5_REAL16; break;
    case CCTK_VARIABLE_COMPLEX32: retval = myGH->HDF5_COMPLEX32; break;
#endif

    default: CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Unsupported CCTK variable datatype %d", cctk_type);
             retval = -1;
  }

  return (retval);
}


} // namespace CarpetIOHDF5
