#ifndef UTIL_HH
#define UTIL_HH

#include <mpi.h>

#include <carpet.hh>

#include <cctk.h>

#include <SimulationIO/RegionCalculus.hpp>
#include <SimulationIO/SimulationIO.hpp>

#include <array>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace CarpetSimulationIO {
namespace RC = RegionCalculus;
using namespace SimulationIO;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

enum class file_format { hdf5, asdf };
enum class file_type { local, global };

inline ostream &operator<<(ostream &os, file_format format) {
  switch (format) {
  case file_format::hdf5:
    return os << "file_format::hdf5";
  case file_format::asdf:
    return os << "file_format::asdf";
  }
  assert(0);
}

inline ostream &operator<<(ostream &os, file_type type) {
  switch (type) {
  case file_type::local:
    return os << "local";
  case file_type::global:
    return os << "global";
  }
  assert(0);
}

////////////////////////////////////////////////////////////////////////////////

uint32_t adler32(const unsigned char *data, size_t len);
uint32_t adler32(const string &data);

////////////////////////////////////////////////////////////////////////////////

// This might execute in the wrong order
// template <typename... Args> std::string stringify(Args &&... args) {
//   std::ostringstream buf;
//   std::array<int, sizeof...(Args)>{(buf << args, 0)...};
//   return buf.str();
// }

namespace detail {
inline void stringify(std::ostringstream &buf) {}
template <typename Arg, typename... Args>
void stringify(std::ostringstream &buf, Arg &&arg, Args &&... args) {
  buf << arg;
  stringify(buf, std::forward<Args>(args)...);
}
} // namespace detail
template <typename... Args> std::string stringify(Args &&... args) {
  std::ostringstream buf;
  detail::stringify(buf, std::forward<Args>(args)...);
  return buf.str();
}

////////////////////////////////////////////////////////////////////////////////

template <typename R, int DR, typename T, int DT>
point<R, DR> vect2point(const CarpetLib::vect<T, DT> &ipos) {
  assert(DR <= DT);
  point<R, DR> pos;
  for (int d = 0; d < DR; ++d)
    pos[d] = ipos[d];
  return pos;
}

template <typename R, int DR, typename T, int DT>
box<R, DR> bbox2box(const CarpetLib::bbox<T, DT> &ibox) {
  auto ioffset = (ibox.lower() % ibox.stride() + ibox.stride()) % ibox.stride();
  assert(CarpetLib::all((ibox.lower() - ioffset) % ibox.stride() == 0));
  auto iorigin = (ibox.lower() - ioffset) / ibox.stride();
  auto ishape = ibox.shape() / ibox.stride();
  auto offset = vect2point<R, DR>(iorigin);
  auto shape = vect2point<R, DR>(ishape);
  return box<R, DR>(offset, offset + shape);
}

template <typename R, int DR, typename T, int DT>
region<R, DR> bboxset2region(const CarpetLib::bboxset<T, DT> &iset) {
  vector<CarpetLib::bbox<T, DT> > iboxes;
  iset.serialise(iboxes);
  vector<box<R, DR> > boxes;
  for (const auto &ibox : iboxes)
    boxes.push_back(bbox2box<R, DR>(ibox));
  return region<R, DR>(boxes);
}

template <typename R, typename T, int DT>
dpoint<R> vect2dpoint(const CarpetLib::vect<T, DT> &ipos, int DR) {
  assert(DR <= DT);
  switch (DR) {
  case 0:
    return vect2point<R, 0>(ipos);
  case 1:
    return vect2point<R, 1>(ipos);
  case 2:
    return vect2point<R, 2>(ipos);
  case 3:
    return vect2point<R, 3>(ipos);
  case 4:
    return vect2point<R, 4>(ipos);
  default:
    CCTK_ERROR("internal error");
  }
}

template <typename R, typename T, int DT>
dbox<R> bbox2dbox(const CarpetLib::bbox<T, DT> &ibox, int DR) {
  assert(DR <= DT);
  switch (DR) {
  case 0:
    return bbox2box<R, 0>(ibox);
  case 1:
    return bbox2box<R, 1>(ibox);
  case 2:
    return bbox2box<R, 2>(ibox);
  case 3:
    return bbox2box<R, 3>(ibox);
  case 4:
    return bbox2box<R, 4>(ibox);
  default:
    CCTK_ERROR("internal error");
  }
}

template <typename R, typename T, int DT>
dregion<R> bboxset2dregion(const CarpetLib::bboxset<T, DT> &iset, int DR) {
  assert(DR <= DT);
  switch (DR) {
  case 0:
    return bboxset2region<R, 0>(iset);
  case 1:
    return bboxset2region<R, 1>(iset);
  case 2:
    return bboxset2region<R, 2>(iset);
  case 3:
    return bboxset2region<R, 3>(iset);
  case 4:
    return bboxset2region<R, 4>(iset);
  default:
    CCTK_ERROR("internal error");
  }
}

template <typename T, int DT>
RC::point_t vect2point_t(const CarpetLib::vect<T, DT> &ipos, int DR) {
  return vect2dpoint<long long>(ipos, DR);
}

template <typename T, int DT>
RC::box_t bbox2box_t(const CarpetLib::bbox<T, DT> &ibox, int DR) {
  return bbox2dbox<long long>(ibox, DR);
}

template <typename T, int DT>
RC::region_t bboxset2region_t(const CarpetLib::bboxset<T, DT> &iset, int DR) {
  return bboxset2dregion<long long>(iset, DR);
}

////////////////////////////////////////////////////////////////////////////////

H5::DataType cactustype2hdf5type(int cactustype);
ASDF::datatype_t cactustype2asdftype(int cactustype);
MPI_Datatype cactustype2mpitype(int cactustype);

////////////////////////////////////////////////////////////////////////////////

// Field widths for various quantities
// TODO: make these run-time parameters
extern const int width_it, width_tl, width_m, width_rl, width_c;

extern const vector<string> dirnames;
extern const vector<string> tensortypes_scalars;
extern const vector<string> tensortypes_vectors;
extern const vector<string> tensortypes_symmetric_tensors;
extern const vector<string> tensortypes_tensors;
extern const vector<string> tensortypes_tensors_rank3_symmetric12;

extern const vector<vector<vector<int> > > syminds;
extern const vector<vector<vector<int> > > syminds_rank3_symmetric12;

////////////////////////////////////////////////////////////////////////////////

// Create the final file name on a particular processor
enum class io_dir_t {
  none,
  input,
  output,
  recover,
  checkpoint,
};

string tolower(string str);
string basename(const string &path);

string charptr2string(const char *const &str);
string charptr2string(char *&&str);
string charptr2string(char *const &str);

string generate_projectname(int variable);
string generate_projectname();
string generate_filename(io_dir_t io_dir, const string &basename,
                         const string &extra_suffix, int iteration,
                         file_format output_format, file_type output_type,
                         int myioproc, int ioproc_every);
tuple<string, int> split_filename(const string &filename);

////////////////////////////////////////////////////////////////////////////////

struct grid_structure_v7 {
  struct map_structure {
    vector<vector<Carpet::region_t> > superregions;
    // vector<vector<Carpet::region_t> > regions;
    vector<i2vect> ghost_widths;
    vector<i2vect> buffer_widths;
    vector<i2vect> overlap_widths;
    vector<int> prolongation_orders_space;
  };
  vector<map_structure> maps;
  vector<vector<CCTK_REAL> > times;
  vector<CCTK_REAL> deltas;
  CCTK_REAL global_time;
  CCTK_REAL delta_time;
};

string serialise_grid_structure(const cGH *cctkGH);
grid_structure_v7 deserialise_grid_structure(const string &grid_structure);

////////////////////////////////////////////////////////////////////////////////

void init_comm();
void finalise_comm();
void send_data(int ioproc, const void *data, int cactustype,
               const dbox<long long> &memlayout, const dbox<long long> &membox);
vector<char> recv_data(int dataproc, int cactustype,
                       const dbox<long long> &membox);

////////////////////////////////////////////////////////////////////////////////

} // namespace CarpetSimulationIO

#endif // #ifndef UTIL_HH
