#include <cctk.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>

#include <hdf5.h>

#include <defs.hh>

#include "iof5.hh"

namespace CarpetIOF5 {

using namespace std;

// Write an int attribute
bool WriteAttribute(hid_t const group, char const *const name,
                    int const ivalue) {
  bool error_flag = false;

  hid_t const dataspace = FAILWARN(H5Screate(H5S_SCALAR));
  hid_t const attribute = FAILWARN(H5Acreate(
      group, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Awrite(attribute, H5T_NATIVE_INT, &ivalue));
  FAILWARN(H5Aclose(attribute));
  FAILWARN(H5Sclose(dataspace));

  return error_flag;
}

// Write a double attribute
bool WriteAttribute(hid_t const group, char const *const name,
                    double const dvalue) {
  bool error_flag = false;

  hid_t const dataspace = FAILWARN(H5Screate(H5S_SCALAR));
  hid_t const attribute = FAILWARN(H5Acreate(
      group, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dvalue));
  FAILWARN(H5Aclose(attribute));
  FAILWARN(H5Sclose(dataspace));

  return error_flag;
}

// Write a string attribute
bool WriteAttribute(hid_t const group, char const *const name,
                    char const *const svalue) {
  bool error_flag = false;

  hid_t const datatype = FAILWARN(H5Tcopy(H5T_C_S1));
  FAILWARN(H5Tset_size(datatype, strlen(svalue) + 1));
  hid_t const dataspace = FAILWARN(H5Screate(H5S_SCALAR));
  hid_t const attribute = FAILWARN(
      H5Acreate(group, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Awrite(attribute, datatype, svalue));
  FAILWARN(H5Aclose(attribute));
  FAILWARN(H5Sclose(dataspace));
  FAILWARN(H5Tclose(datatype));

  return error_flag;
}

bool WriteAttribute(hid_t const group, char const *const name,
                    string const svalue) {
  return WriteAttribute(group, name, svalue.c_str());
}

// Write an array of int attributes
bool WriteAttribute(hid_t const group, char const *const name,
                    int const *const ivalues, int const nvalues) {
  bool error_flag = false;

  hsize_t const size = nvalues;
  hid_t const dataspace = FAILWARN(H5Screate_simple(1, &size, NULL));
  hid_t const attribute = FAILWARN(H5Acreate(
      group, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Awrite(attribute, H5T_NATIVE_INT, ivalues));
  FAILWARN(H5Aclose(attribute));
  FAILWARN(H5Sclose(dataspace));

  return error_flag;
}

// Write an array of double attributes
bool WriteAttribute(hid_t const group, char const *const name,
                    double const *const dvalues, int const nvalues) {
  bool error_flag = false;

  hsize_t const size = nvalues;
  hid_t const dataspace = FAILWARN(H5Screate_simple(1, &size, NULL));
  hid_t const attribute = FAILWARN(H5Acreate(
      group, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Awrite(attribute, H5T_NATIVE_DOUBLE, dvalues));
  FAILWARN(H5Aclose(attribute));
  FAILWARN(H5Sclose(dataspace));

  return error_flag;
}

// Write an array of string attributes
bool WriteAttribute(hid_t const group, char const *const name,
                    char const *const *const svalues, int const nvalues) {
  bool error_flag = false;

  size_t maxstrlen = 0;
  for (int i = 0; i < nvalues; ++i) {
    maxstrlen = max(maxstrlen, strlen(svalues[i]));
  }
  vector<char> svalue(nvalues * (maxstrlen + 1));
  for (int i = 0; i < nvalues; ++i) {
    strncpy(&svalue.at(i * maxstrlen), svalues[i], maxstrlen + 1);
  }

  hid_t const datatype = FAILWARN(H5Tcopy(H5T_C_S1));
  FAILWARN(H5Tset_size(datatype, maxstrlen + 1));
  hsize_t const size = nvalues;
  hid_t const dataspace = FAILWARN(H5Screate_simple(1, &size, NULL));
  hid_t const attribute = FAILWARN(
      H5Acreate(group, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Awrite(attribute, datatype, &svalue.front()));
  FAILWARN(H5Aclose(attribute));
  FAILWARN(H5Sclose(dataspace));
  FAILWARN(H5Tclose(datatype));

  return error_flag;
}

// Write a large string attribute
bool WriteLargeAttribute(hid_t const group, char const *const name,
                         char const *const svalue) {
  bool error_flag = false;

  // Create a dataset, since the data may not fit into an attribute
  hsize_t const size = strlen(svalue) + 1;
  hid_t const dataspace = FAILWARN(H5Screate_simple(1, &size, NULL));
  hid_t const dataset =
      FAILWARN(H5Dcreate(group, name, H5T_NATIVE_CHAR, dataspace, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT));
  FAILWARN(H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    svalue));
  FAILWARN(H5Dclose(dataset));
  FAILWARN(H5Sclose(dataspace));

  return error_flag;
}

// Read an int attribute
bool ReadAttribute(hid_t const group, char const *const name, int &ivalue) {
  bool error_flag = false;

  hid_t const attribute = FAILWARN(H5Aopen_name(group, name));
  FAILWARN(H5Aread(attribute, H5T_NATIVE_INT, &ivalue));
  FAILWARN(H5Aclose(attribute));

  return error_flag;
}

// Read a double attribute
bool ReadAttribute(hid_t const group, char const *const name, double &dvalue) {
  bool error_flag = false;

  hid_t const attribute = FAILWARN(H5Aopen_name(group, name));
  FAILWARN(H5Aread(attribute, H5T_NATIVE_DOUBLE, &dvalue));
  FAILWARN(H5Aclose(attribute));

  return error_flag;
}

// Read a string attribute
bool ReadAttribute(hid_t const group, char const *const name, string &svalue) {
  bool error_flag = false;

  hid_t const attribute = FAILWARN(H5Aopen_name(group, name));
  hid_t const datatype = FAILWARN(H5Aget_type(attribute));
  hsize_t const size = FAILWARN0(H5Tget_size(datatype));
  char buf[size + 1];
  hid_t const memdatatype = FAILWARN(H5Tcopy(H5T_C_S1));
  FAILWARN(H5Tset_size(memdatatype, size));
  FAILWARN(H5Aread(attribute, memdatatype, buf));
  buf[size] = '\0'; // just to be safe
  svalue = string(buf);
  FAILWARN(H5Tclose(memdatatype));
  FAILWARN(H5Tclose(datatype));
  FAILWARN(H5Aclose(attribute));

  return error_flag;
}

// Read a large string attribute
bool ReadLargeAttribute(hid_t const group, char const *const name,
                        string &svalue) {
  bool error_flag = false;

  // Read a dataset
  hid_t const dataset = FAILWARN(H5Dopen(group, name, H5P_DEFAULT));
  hid_t const dataspace = FAILWARN(H5Dget_space(dataset));
  int const rank = FAILWARN(H5Sget_simple_extent_ndims(dataspace));
  if (rank != 1)
    error_flag = true;
  hsize_t dims[rank];
  FAILWARN(H5Sget_simple_extent_dims(dataspace, dims, NULL));
  hsize_t size = 1;
  for (int d = 0; d < rank; ++d)
    size *= dims[d];
  char buf[size + 1];
  FAILWARN(
      H5Dread(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf));
  buf[size] = '\0'; // just to be safe
  svalue = string(buf);
  FAILWARN(H5Sclose(dataspace));
  FAILWARN(H5Dclose(dataset));

  return error_flag;
}

} // namespace CarpetIOF5
