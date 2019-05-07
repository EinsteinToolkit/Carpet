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
#include <utility>

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

string basename(const string &path) {
  auto pos = path.rfind('/');
  if (pos == string::npos)
    return path;
  return path.substr(0, pos);
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

ASDF::datatype_t cactustype2asdftype(int cactustype) {
  switch (cactustype) {
  case CCTK_VARIABLE_BYTE:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<uint8_t>::value);
  case CCTK_VARIABLE_INT:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<CCTK_INT>::value);
  case CCTK_VARIABLE_INT1:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<int8_t>::value);
  case CCTK_VARIABLE_INT2:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<int16_t>::value);
  case CCTK_VARIABLE_INT4:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<int32_t>::value);
  case CCTK_VARIABLE_INT8:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<int64_t>::value);
  case CCTK_VARIABLE_REAL:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<CCTK_REAL>::value);
  case CCTK_VARIABLE_REAL4:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<ASDF::float32_t>::value);
  case CCTK_VARIABLE_REAL8:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<ASDF::float64_t>::value);
  // case CCTK_VARIABLE_REAL16:
  //   return ASDF::datatype_t(ASDF::get_scalar_type_id<CCTK_REAL16>::value);
  case CCTK_VARIABLE_COMPLEX:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<CCTK_COMPLEX>::value);
  case CCTK_VARIABLE_COMPLEX8:
    return ASDF::datatype_t(ASDF::get_scalar_type_id<ASDF::complex64_t>::value);
  case CCTK_VARIABLE_COMPLEX16:
    return ASDF::datatype_t(
        ASDF::get_scalar_type_id<ASDF::complex128_t>::value);
  // case CCTK_VARIABLE_COMPLEX32:
  //   return ASDF::datatype_t(ASDF::get_scalar_type_id<CCTK_COMPLEX32>::value);
  default:
    CCTK_VERROR("Unsupported Cactus datatype %d", cactustype);
  }
}

MPI_Datatype cactustype2mpitype(int cactustype) {
  switch (cactustype) {
  case CCTK_VARIABLE_BYTE:
    return dist::mpi_datatype(CCTK_BYTE{});
  case CCTK_VARIABLE_INT:
    return dist::mpi_datatype(CCTK_INT{});
  case CCTK_VARIABLE_INT1:
    return dist::mpi_datatype(CCTK_INT1{});
  case CCTK_VARIABLE_INT2:
    return dist::mpi_datatype(CCTK_INT2{});
  case CCTK_VARIABLE_INT4:
    return dist::mpi_datatype(CCTK_INT4{});
  case CCTK_VARIABLE_INT8:
    return dist::mpi_datatype(CCTK_INT8{});
  case CCTK_VARIABLE_REAL:
    return dist::mpi_datatype(CCTK_REAL{});
  case CCTK_VARIABLE_REAL4:
    return dist::mpi_datatype(CCTK_REAL4{});
  case CCTK_VARIABLE_REAL8:
    return dist::mpi_datatype(CCTK_REAL8{});
  case CCTK_VARIABLE_REAL16:
    return dist::mpi_datatype(CCTK_REAL16{});
  case CCTK_VARIABLE_COMPLEX:
    return dist::mpi_datatype(CCTK_COMPLEX{});
  case CCTK_VARIABLE_COMPLEX8:
    return dist::mpi_datatype(CCTK_COMPLEX8{});
  case CCTK_VARIABLE_COMPLEX16:
    return dist::mpi_datatype(CCTK_COMPLEX16{});
  case CCTK_VARIABLE_COMPLEX32:
    return dist::mpi_datatype(CCTK_COMPLEX32{});
  case CCTK_VARIABLE_CHAR:
    return dist::mpi_datatype(CCTK_CHAR{});
  default:
    CCTK_VERROR("Unsupported Cactus datatype %d", cactustype);
  }
}

////////////////////////////////////////////////////////////////////////////////

string serialise_grid_structure(const cGH *cctkGH) {
  using namespace CarpetLib;
  const auto &tt = *Carpet::tt;
  vector<vector<CCTK_REAL> > times;
  times.resize(Carpet::reflevels);
  for (int rl = 0; rl < Carpet::reflevels; ++rl) {
    times.AT(rl).resize(tt.timelevels);
    for (int tl = 0; tl < tt.timelevels; ++tl)
      times.AT(rl).AT(tl) = tt.get_time(Carpet::mglevel, rl, tl);
  }
  vector<CCTK_REAL> deltas;
  deltas.resize(Carpet::reflevels);
  for (int rl = 0; rl < Carpet::reflevels; ++rl)
    deltas.AT(rl) = tt.get_delta(Carpet::mglevel, rl);
  ostringstream buf;
  buf << setprecision(17);
  buf << "Grid Structure v7"
      << "{"
      << "maps=" << Carpet::maps << ","
      << "[";
  for (int m = 0; m < Carpet::maps; ++m) {
    const auto &hh = *Carpet::vhh.AT(m);
    const auto &dd = *Carpet::vdd.AT(m);
    buf << "[" << m << "]={"
        << "superregions:" << hh.superregions << ","
        << "ghost_widths:" << dd.ghost_widths << ","
        << "buffer_widths:" << dd.buffer_widths << ","
        << "overlap_widths:" << dd.overlap_widths << ","
        << "prolongation_orders_space:";
    output(buf, dd.prolongation_orders_space);
    buf << ","
        << "},";
  }
  buf << "],"
      << "times:";
  output(buf, times);
  buf << ","
      << "deltas:";
  output(buf, deltas);
  buf << ","
      << "global_time:" << Carpet::global_time << ","
      << "delta_time:" << Carpet::delta_time << ","
      << "}.";
  return buf.str();
}

grid_structure_v7
deserialise_grid_structure(const string &grid_structure_string) {
  grid_structure_v7 grid_structure;
  istringstream buf(grid_structure_string);
  skipws(buf);
  consume(buf, "Grid Structure v7");
  consume(buf, "{");
  consume(buf, "maps=");
  int nmaps;
  buf >> nmaps;
  assert(nmaps == Carpet::maps);
  consume(buf, ",");
  consume(buf, "[");
  for (int m = 0; m < nmaps; ++m) {
    consume(buf, "[");
    int mm;
    buf >> mm;
    assert(mm == m);
    consume(buf, "]={");
    grid_structure.maps.push_back({});
    consume(buf, "superregions:");
    buf >> grid_structure.maps[m].superregions;
    consume(buf, ",");
    consume(buf, "ghost_widths:");
    buf >> grid_structure.maps[m].ghost_widths;
    consume(buf, ",");
    consume(buf, "buffer_widths:");
    buf >> grid_structure.maps[m].buffer_widths;
    consume(buf, ",");
    consume(buf, "overlap_widths:");
    buf >> grid_structure.maps[m].overlap_widths;
    consume(buf, ",");
    consume(buf, "prolongation_orders_space:");
    buf >> grid_structure.maps[m].prolongation_orders_space;
    consume(buf, ",");
    consume(buf, "},");
  }
  consume(buf, "],");
  consume(buf, "times:");
  buf >> grid_structure.times;
  consume(buf, ",");
  consume(buf, "deltas:");
  buf >> grid_structure.deltas;
  consume(buf, ",");
  consume(buf, "global_time:");
  buf >> grid_structure.global_time;
  consume(buf, ",");
  consume(buf, "delta_time:");
  buf >> grid_structure.delta_time;
  consume(buf, ",");
  consume(buf, "}.");
  return grid_structure;
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
const vector<string> tensortypes_tensors_rank3_symmetric12{
    "n/a", "n/a", "Tensor2D_Rank3Symmetric12", "Tensor3D_Rank3Symmetric12",
    "Tensor4D_Rank3Symmetric12"};

#if 0
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

const vector<vector<vector<int> > > syminds_rank3_symmetric12{
    {{0, 0, 0}},
    {{0, 0, 0}},
    {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 1, 0}, {1, 1, 1}},
    {{0, 0, 0},
     {0, 0, 1},
     {0, 0, 2},
     {0, 1, 0},
     {0, 1, 1},
     {0, 1, 2},
     {0, 2, 0},
     {0, 2, 1},
     {0, 2, 2},
     {1, 1, 0},
     {1, 1, 1},
     {1, 1, 2},
     {1, 2, 0},
     {1, 2, 1},
     {1, 2, 2},
     {2, 2, 0},
     {2, 2, 1},
     {2, 2, 2}},
    {{0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 0, 3}, {0, 1, 0}, {0, 1, 1},
     {0, 1, 2}, {0, 1, 3}, {0, 2, 0}, {0, 2, 1}, {0, 2, 2}, {0, 2, 3},
     {0, 3, 0}, {0, 3, 1}, {0, 3, 2}, {0, 3, 3}, {1, 1, 0}, {1, 1, 1},
     {1, 1, 2}, {1, 1, 3}, {1, 2, 0}, {1, 2, 1}, {1, 2, 2}, {1, 2, 3},
     {1, 3, 0}, {1, 3, 1}, {1, 3, 2}, {1, 3, 3}, {2, 2, 0}, {2, 2, 1},
     {2, 2, 2}, {2, 2, 3}, {2, 3, 0}, {2, 3, 1}, {2, 3, 2}, {2, 3, 3},
     {3, 3, 0}, {3, 3, 1}, {3, 3, 2}, {3, 3, 3}}};
#endif

const vector<vector<vector<int> > > syminds = []() {
  vector<vector<vector<int> > > inds3;
  for (int dim = 0; dim <= 4; ++dim) {
    vector<vector<int> > inds2;
    for (int i = 0; i < dim; ++i)
      for (int j = i; j < dim; ++j)
        inds2.push_back({i, j});
    assert(int(inds2.size()) == dim * (dim + 1) / 2);
    inds3.push_back(std::move(inds2));
  }
  return inds3;
}();

const vector<vector<vector<int> > > syminds_rank3_symmetric12 = []() {
  vector<vector<vector<int> > > inds3;
  for (int dim = 0; dim <= 4; ++dim) {
    vector<vector<int> > inds2;
    for (int i = 0; i < dim; ++i)
      for (int j = i; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          inds2.push_back({i, j, k});
    assert(int(inds2.size()) == dim * dim * (dim + 1) / 2);
    inds3.push_back(std::move(inds2));
  }
  return inds3;
}();

////////////////////////////////////////////////////////////////////////////////

// Generate a good file name ("alias") for a variable
string generate_projectname(int vindex) {
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
string generate_projectname() {
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
string generate_filename(io_dir_t io_dir, const string &basename,
                         const string &extra_suffix, int iteration,
                         file_format output_format, file_type output_type,
                         int myioproc, int ioproc_every) {
  DECLARE_CCTK_PARAMETERS;

  string IO_dir;
  string S5_dir;
  switch (io_dir) {
  case io_dir_t::none:
    IO_dir = ".";
    S5_dir = ".";
    break;
  case io_dir_t::input:
    IO_dir = IO_filereader_ID_dir;
    S5_dir = filereader_ID_dir;
    break;
  case io_dir_t::output:
    IO_dir = IO_out_dir;
    S5_dir = out_dir;
    break;
  case io_dir_t::recover:
    IO_dir = IO_recover_dir;
    S5_dir = recover_dir;
    break;
  case io_dir_t::checkpoint:
    IO_dir = IO_checkpoint_dir;
    S5_dir = checkpoint_dir;
    break;
  default:
    assert(0);
  }
  string path = S5_dir.empty() ? IO_dir : S5_dir;

  if (output_type == file_type::local) {
    const int nprocs = CCTK_nProcs(nullptr);
    const int nioprocs = (nprocs + ioproc_every - 1) / ioproc_every;

    // Create additional directory levels until we can handle all I/O processes
    int max_dir_nioprocs = 1; // Number of I/O processes we can handle
    while (max_dir_nioprocs < nioprocs) {
      // Need one more level
      assert(INT_MAX / max_files_per_directory >= max_dir_nioprocs);
      max_dir_nioprocs *= max_files_per_directory;
    }

    for (int dir_nioprocs = max_dir_nioprocs; dir_nioprocs >= 1;
         dir_nioprocs /= max_files_per_directory) {
      ostringstream buf;
      int dir_ioproc = (myioproc / ioproc_every) / dir_nioprocs * dir_nioprocs;
      buf << path << "/" << basename << ".p" << setw(processor_digits)
          << setfill('0') << dir_ioproc;
      path = buf.str();
    }
  }

  ostringstream buf;
  buf << path << "/" << basename;
  // if (io_dir == io_dir_recover or io_dir == io_dir_checkpoint) {
  buf << ".it" << setw(iteration_digits) << setfill('0') << iteration;
  // }
  if (output_type == file_type::local)
    buf << ".p" << setw(processor_digits) << setfill('0')
        << (myioproc / ioproc_every);
  buf << extra_suffix;
  switch (output_format) {
  case file_format::hdf5:
    buf << out_extension;
    break;
  case file_format::asdf:
    buf << out_extension_asdf;
    break;
  default:
    assert(0);
  }
  path = buf.str();

  return path;
}

tuple<string, int> split_filename(const string &filename) {
  string itprefix = ".it";
  auto itpos = filename.find(itprefix);
  assert(itpos != string::npos);
  string basename = filename.substr(0, itpos);
  string itstr = filename.substr(itpos + itprefix.length());
  int iteration;
  try {
    iteration = stoi(itstr);
  } catch (...) {
    assert(0);
  }
  return {basename, iteration};
}

////////////////////////////////////////////////////////////////////////////////

const int tag = 2;
MPI_Comm comm = MPI_COMM_NULL;

void init_comm() {
  assert(comm == MPI_COMM_NULL);
  MPI_Comm_dup(dist::comm(), &comm);
  assert(comm != MPI_COMM_NULL);
}
void finalise_comm() {
  assert(comm != MPI_COMM_NULL);
  MPI_Comm_free(&comm);
  assert(comm == MPI_COMM_NULL);
}

void send_data(int ioproc, const void *data, int cactustype,
               const dbox<long long> &memlayout,
               const dbox<long long> &membox) {
  ptrdiff_t count = membox.size();
  MPI_Datatype mpitype = cactustype2mpitype(cactustype);
  size_t typesize = CCTK_VarTypeSize(cactustype);
  vector<unsigned char> buf(count * typesize);
  HyperSlab::copy(buf.data(), buf.size(), membox, membox, data,
                  memlayout.size(), memlayout, membox, typesize);
  MPI_Send(buf.data(), count, mpitype, ioproc, tag, comm);
}

vector<char> recv_data(int dataproc, int cactustype,
                       const dbox<long long> &membox) {
  ptrdiff_t count = membox.size();
  MPI_Datatype mpitype = cactustype2mpitype(cactustype);
  size_t typesize = CCTK_VarTypeSize(cactustype);
  vector<char> buf(count * typesize);
  MPI_Recv(buf.data(), count, mpitype, dataproc, tag, comm, MPI_STATUS_IGNORE);
  return buf;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CarpetSimulationIO
