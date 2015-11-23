#ifndef DISTRIBUTE_HH
#define DISTRIBUTE_HH

#include <cctk.h>

#include <defs.hh>
#include <vect.hh>

#include <list>

#include <mpi.h>

namespace CarpetIOF5 {

using namespace std;

// Fragment descriptor
struct fragdesc_t {
  // only integer entries
  int varindex;
  static int const mglevel = 0;
  int reflevel, map, src_component, timelevel;
  ivect imin, imax; // upper bound is inclusive

  // destination
  int component, process;
  int tag;

  // meta-messages
  enum state_t { state_normal, state_sent_all, state_all_sent_all };
  int state;

  fragdesc_t() : state(state_normal) {}

  int num_ints() const { return sizeof(fragdesc_t) / sizeof(int); }
  int npoints() const;
  int vartypesize() const;
  MPI_Datatype datatype() const;
};

// Scatter (distribute) Cactus variables
class scatter_t {

  cGH const *const cctkGH;

  static int const max_concurrent_sends = 100;
  enum tags_t {
    tag_desc,
    tag_data_min,
    tag_data_max = tag_data_min + max_concurrent_sends
  };

  class busy_tags_t {
    vector<bool> tags;

  public:
    busy_tags_t();
    ~busy_tags_t();
    int allocate_tag();
    void free_tag(int tag);
  };
  vector<busy_tags_t> busy_tags;

  struct transmission_t {
    fragdesc_t fragdesc;
    vector<char> data;
    MPI_Request request;
    transmission_t() : request(MPI_REQUEST_NULL) {}
    ~transmission_t() { assert(request == MPI_REQUEST_NULL); }
  };

  // Desired number of public receives to keep posted at all times
  static int const num_public_recvs = 10;
  list<transmission_t *> public_recvs, recvs, sends;
  int num_received, num_sent;
  ptrdiff_t bytes_allocated;

  bool did_send_all;
  int did_receive_sent_all;
  bool did_receive_all_sent_all;
  fragdesc_t to_parent, to_children;

public:
  scatter_t(cGH const *cctkGH_);
  ~scatter_t();

  void send(fragdesc_t const &fd, void const *data);

private:
  void post_public_recvs();
  void do_some_work(bool do_wait = false);

  void maybe_send_did_send();
  void send_all_did_send();
  void set_did_send_all();
  void set_did_receive_sent_all();
  void set_did_receive_all_sent_all();
  void handle_state_transition(fragdesc_t const &fd);

  list<transmission_t *> split_for_sending(fragdesc_t const &fd,
                                           void const *data);
  void write_data(transmission_t *const tm);
};

} // namespace CarpetIOF5

#endif // #ifndef DISTRIBUTE_HH
