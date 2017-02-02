#ifndef UTIL_HH
#define UTIL_HH

#include <carpet.hh>

#include <cctk.h>

#include <RegionCalculus.hpp>
#include <SimulationIO.hpp>

#include <array>
#include <sstream>
#include <string>
#include <utility>

namespace CarpetSimulationIO {
using namespace SimulationIO;
using namespace std;

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
}
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

////////////////////////////////////////////////////////////////////////////////

H5::DataType cactustype2hdf5type(int cactustype);

////////////////////////////////////////////////////////////////////////////////

// Field widths for various quantities
// TODO: make these run-time parameters
extern const int width_it, width_tl, width_m, width_rl, width_c;

extern const vector<string> dirnames;
extern const vector<string> tensortypes_scalars;
extern const vector<string> tensortypes_vectors;
extern const vector<string> tensortypes_symmetric_tensors;
extern const vector<string> tensortypes_tensors;

extern const vector<vector<vector<int> > > syminds;

////////////////////////////////////////////////////////////////////////////////

// Create the final file name on a particular processor
enum io_dir_t {
  io_dir_none,
  io_dir_input,
  io_dir_output,
  io_dir_recover,
  io_dir_checkpoint,
};

string tolower(string str);

string charptr2string(const char *const &str);
string charptr2string(char *&&str);
string charptr2string(char *const &str);

string generate_projectname(const cGH *cctkGH, int variable);
string generate_projectname(const cGH *cctkGH);
string generate_filename(const cGH *cctkGH, io_dir_t io_dir,
                         const string &basename, int iteration, int ioproc,
                         int nioprocs);

string serialize_grid_structure(const cGH *cctkGH);
}

#endif // #ifndef UTIL_HH