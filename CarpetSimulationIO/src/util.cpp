#include "util.hpp"

#include <carpet.hh>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace CarpetSimulationIO {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// where data is the location of the data in physical memory and
// len is the length of the data in bytes
uint32_t adler32(const unsigned char *data, size_t len) {
  constexpr int MOD_ADLER = 65521;

  uint32_t a = 1, b = 0;

  // Process each byte of the data in order
  for (size_t index = 0; index < len; ++index) {
    a = (a + data[index]) % MOD_ADLER;
    b = (b + a) % MOD_ADLER;
  }

  return (b << 16) | a;
}

uint32_t adler32(const string &data) {
  return adler32((const unsigned char *)data.c_str(), data.size());
}

////////////////////////////////////////////////////////////////////////////////

// Convert a string to lower case
using std::tolower;
string tolower(string str) {
  // transform(str.begin(), str.end(), str.begin(), std::tolower);
  for (char &ch : str)
    ch = tolower(ch);
  return str;
}

// Convert a const char* or char* to string
string charptr2string(const char *const &str) {
  string res{str};
  return res;
}
string charptr2string(char *&&str) {
  string res{str};
  // Cactus convention is that returned char* needs to be freed by the caller
  free(str);
  return res;
}

////////////////////////////////////////////////////////////////////////////////

H5::DataType cactustype2hdf5type(int cactustype) {
  switch (cactustype) {
  case CCTK_VARIABLE_BYTE:
    return H5::getType(CCTK_BYTE{});
  case CCTK_VARIABLE_INT:
    return H5::getType(CCTK_INT{});
  case CCTK_VARIABLE_INT1:
    return H5::getType(CCTK_INT1{});
  case CCTK_VARIABLE_INT2:
    return H5::getType(CCTK_INT2{});
  case CCTK_VARIABLE_INT4:
    return H5::getType(CCTK_INT4{});
  case CCTK_VARIABLE_INT8:
    return H5::getType(CCTK_INT8{});
  case CCTK_VARIABLE_REAL:
    return H5::getType(CCTK_REAL{});
  case CCTK_VARIABLE_REAL4:
    return H5::getType(CCTK_REAL4{});
  case CCTK_VARIABLE_REAL8:
    return H5::getType(CCTK_REAL8{});
  case CCTK_VARIABLE_REAL16:
    return H5::getType(CCTK_REAL16{});
  case CCTK_VARIABLE_COMPLEX:
    return H5::getType(CCTK_COMPLEX{});
  case CCTK_VARIABLE_COMPLEX8:
    return H5::getType(CCTK_COMPLEX8{});
  case CCTK_VARIABLE_COMPLEX16:
    return H5::getType(CCTK_COMPLEX16{});
  case CCTK_VARIABLE_COMPLEX32:
    return H5::getType(CCTK_COMPLEX32{});
  case CCTK_VARIABLE_CHAR:
    return H5::getType(CCTK_CHAR{});
  default:
    CCTK_VERROR("Unsupported Cactus datatype %d", cactustype);
  }
}

////////////////////////////////////////////////////////////////////////////////

string serialize_grid_structure(const cGH *cctkGH) {
  using namespace CarpetLib;
  ostringstream buf;
  const auto &tt = *Carpet::tt;
  buf << setprecision(17);
  buf << "Grid Structure v6"
      << "{"
      << "maps=" << Carpet::maps << ","
      << "[";
  for (int m = 0; m < Carpet::maps; ++m) {
    const auto &hh = *Carpet::vhh.AT(m);
    const auto &dd = *Carpet::vdd.AT(m);
    buf << "[" << m << "]={"
        << "superregions:" << hh.superregions << ","
        // We could omit the regions
        << "regions:" << hh.regions.AT(0) << ","
        << "ghost_widths:" << dd.ghost_widths << ","
        << "buffer_widths:" << dd.buffer_widths << ","
        << "overlap_widths:" << dd.overlap_widths << ","
        << "prolongation_orders_space:" << dd.prolongation_orders_space << ","
        << "},";
  }
  buf << "],"
      << "global_time:" << Carpet::global_time << ","
      << "delta_time:" << Carpet::delta_time << ","
      << "times:" << tt << ","
      << "}.";
  return buf.str();
}

////////////////////////////////////////////////////////////////////////////////

const int width_it = 10, width_tl = 1, width_m = 3, width_rl = 2, width_c = 8;

const vector<string> dirnames{"x", "y", "z", "w"};

const vector<string> tensortypes_scalars{"Scalar0D", "Scalar1D", "Scalar2D",
                                         "Scalar3D", "Scalar4D"};
const vector<string> tensortypes_vectors{"Vector0D", "Vector1D", "Vector2D",
                                         "Vector3D", "Vector4D"};
const vector<string> tensortypes_symmetric_tensors{
    "n/a", "n/a", "SymmetricTensor2D", "SymmetricTensor3D",
    "SymmetricTensor4D"};
const vector<string> tensortypes_tensors{"n/a", "Tensor1D", "Tensor2D",
                                         "Tensor3D", "Tensor4D"};

const vector<vector<vector<int> > > syminds{
    {{0, 0}},
    {{0, 0}},
    {{0, 0}, {0, 1}, {1, 1}},
    {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}},
    {{0, 0},
     {0, 1},
     {0, 2},
     {0, 3},
     {1, 1},
     {1, 2},
     {1, 3},
     {2, 2},
     {2, 3},
     {3, 3}}};

////////////////////////////////////////////////////////////////////////////////

// Generate a good file name ("alias") for a variable
string generate_projectname(const cGH *cctkGH, int vindex) {
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0);

  ostringstream filename_buf;

  auto varname = charptr2string(CCTK_FullName(vindex));
  // varname = tolower(varname);
  bool skip_colon = false;
  for (char ch : varname) {
    switch (ch) {
    case ':':
      if (not skip_colon)
        filename_buf << '-';
      skip_colon = true;
      break;
    case '[':
      filename_buf << '-';
      break;
    case ']':
      break;
    default:
      filename_buf << char(tolower(ch));
      break;
    }
  }

  return filename_buf.str();
}

// Generate a good file name
string generate_projectname(const cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  ostringstream filename_buf;

  if (CCTK_EQUALS(out_filename, "")) {
    // Obtain the parameter file name
    vector<char> buf(10000);
    int ilen = CCTK_ParameterFilename(buf.size(), buf.data());
    assert(ilen < int(buf.size() - 1));
    string parfilename(buf.data());
    // Remove directory prefix, if any
    auto slash = parfilename.rfind('/');
    if (slash != string::npos)
      parfilename = parfilename.substr(slash + 1);
    // Remove suffix, if it is there
    auto suffix = parfilename.rfind('.');
    if (suffix != string::npos and parfilename.substr(suffix) == ".par")
      parfilename = parfilename.substr(0, suffix);
    filename_buf << parfilename;
  } else {
    filename_buf << out_filename;
  }

  return filename_buf.str();
}

// Generate the final file name on a particular processor
string generate_filename(const cGH *cctkGH, io_dir_t io_dir,
                         const string &basename, const string &extra_suffix,
                         int iteration, int ioproc, int nioprocs,
                         bool create_dirs) {
  DECLARE_CCTK_PARAMETERS;

  const int mode = 0777;

  string IO_dir;
  string S5_dir;
  switch (io_dir) {
  case io_dir_none:
    IO_dir = ".";
    S5_dir = ".";
    break;
  case io_dir_input:
    IO_dir = IO_filereader_ID_dir;
    S5_dir = filereader_ID_dir;
    break;
  case io_dir_output:
    IO_dir = IO_out_dir;
    S5_dir = out_dir;
    break;
  case io_dir_recover:
    IO_dir = IO_recover_dir;
    S5_dir = recover_dir;
    break;
  case io_dir_checkpoint:
    IO_dir = IO_checkpoint_dir;
    S5_dir = checkpoint_dir;
    break;
  default:
    assert(0);
  }
  string path = S5_dir.empty() ? IO_dir : S5_dir;
  if (create_dirs) {
    assert(ioproc >= 0);
    if (ioproc == 0) {
      int ierr = CCTK_CreateDirectory(mode, path.c_str());
      if (ierr < 0)
        CCTK_VERROR("Could not create directory \"%s\"", path.c_str());
    }
    CCTK_Barrier(cctkGH);
  }

  if (ioproc >= 0) {
    const int procs_per_dir = 100;

    long long max_dir_nioprocs = 1;
    while (max_dir_nioprocs < nioprocs)
      max_dir_nioprocs *= procs_per_dir;

    for (long long dir_nioprocs = max_dir_nioprocs; dir_nioprocs >= 1;
         dir_nioprocs /= procs_per_dir) {
      ostringstream buf;
      long long dir_ioproc = ioproc / dir_nioprocs * dir_nioprocs;
      buf << path << "/" << basename << ".p" << setw(processor_digits)
          << setfill('0') << dir_ioproc;
      path = buf.str();
      if (create_dirs) {
        assert(ioproc >= 0);
        if (ioproc == dir_ioproc) {
          int ierr = CCTK_CreateDirectory(mode, path.c_str());
          if (ierr < 0)
            CCTK_VERROR("Could not create directory \"%s\"", path.c_str());
        }
        CCTK_Barrier(cctkGH);
      }
    }
  }

  ostringstream buf;
  buf << path << "/" << basename;
  // if (io_dir == io_dir_recover or io_dir == io_dir_checkpoint) {
  buf << ".it" << setw(iteration_digits) << setfill('0') << iteration;
  // }
  if (ioproc >= 0)
    buf << ".p" << setw(processor_digits) << setfill('0') << ioproc;
  buf << extra_suffix << out_extension;
  path = buf.str();

  return path;
}
}
