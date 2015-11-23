#include "distribute.hh"

#include <cctk_Parameters.h>

#include <cacheinfo.hh>
#include <carpet.hh>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace CarpetIOF5 {

using namespace std;

/*** fragdesc_t *************************************************************/

int fragdesc_t::npoints() const {
  int const np = prod(imax - imin + 1);
  assert(np > 0);
  return np;
}

int fragdesc_t::vartypesize() const {
  assert(varindex >= 0);
  int const vartype = CCTK_VarTypeI(varindex);
  assert(vartype >= 0);
  int const size = CCTK_VarTypeSize(vartype);
  assert(size > 0);
  return size;
}

MPI_Datatype fragdesc_t::datatype() const {
  assert(varindex >= 0);
  int const vartype = CCTK_VarTypeI(varindex);
  assert(vartype >= 0);
  switch (vartype) {
  case CCTK_VARIABLE_BYTE:
    return dist::mpi_datatype<CCTK_BYTE>();
  case CCTK_VARIABLE_INT:
    return dist::mpi_datatype<CCTK_INT>();
  case CCTK_VARIABLE_REAL:
    return dist::mpi_datatype<CCTK_REAL>();
  case CCTK_VARIABLE_COMPLEX:
    return dist::mpi_datatype<CCTK_COMPLEX>();
  default:
    assert(0);
  }
}

/*** scatter_t::busy_tags_t *************************************************/

scatter_t::busy_tags_t::busy_tags_t() {
  // do nothing
}

scatter_t::busy_tags_t::~busy_tags_t() {
  // Ensure all tags are free
  for (int tag = 0; tag < (int)tags.size(); ++tag) {
    assert(not tags.at(tag));
  }
}

int scatter_t::busy_tags_t::allocate_tag() {
  if (tags.empty()) {
    tags.resize(max_concurrent_sends, false);
  }
  for (int tag = 0; tag < (int)tags.size(); ++tag) {
    if (not tags.at(tag)) {
      tags.at(tag) = true;
      return tag;
    }
  }
  return -1; // no free tag available
}

void scatter_t::busy_tags_t::free_tag(int const tag) {
  assert(tags.at(tag));
  tags.at(tag) = false;
}

/*** scatter_t **************************************************************/

scatter_t::scatter_t(cGH const *const cctkGH_)
    : cctkGH(cctkGH_), busy_tags(CCTK_nProcs(cctkGH)), num_received(0),
      num_sent(0), bytes_allocated(0), did_send_all(false),
      did_receive_sent_all(0), did_receive_all_sent_all(false) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO("Creating global scatter object");
  }
  post_public_recvs();
}

scatter_t::~scatter_t() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO("Shutting down global scatter object");
  }

  // Wait until everything has been sent
  bool did_send_everything;
  for (;;) {
    if (verbose) {
      CCTK_INFO("Waiting until something has been transmitted...");
    }
    did_send_everything = sends.empty();
    if (did_send_everything)
      break;
    do_some_work(true);
  }
  if (verbose) {
    CCTK_INFO("We sent all our data");
  }

  // Notify others: we sent all our data
  set_did_send_all();

  // Wait until everything has been transmitted
  bool did_transmit_everything;
  for (;;) {
    if (verbose) {
      CCTK_INFO("Waiting until something has been transmitted...");
    }
    did_transmit_everything =
        did_receive_all_sent_all and recvs.empty() and sends.empty();
    if (did_transmit_everything)
      break;
    do_some_work(true);
  }
  if (verbose) {
    CCTK_INFO("Everything has been transmitted");
  }

  int const my_difference = num_sent - num_received;
  int total_difference;
  MPI_Allreduce(const_cast<int *>(&my_difference), &total_difference, 1,
                MPI_INT, MPI_SUM, dist::comm());
  if (total_difference < 0) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "More received messages than sent messages -- impossible!");
  }
  if (total_difference > 0) {
    CCTK_WARN(CCTK_WARN_ABORT, "Not all sent messages have been received");
  }
  if (verbose) {
    CCTK_INFO("Global number of sent and received messages is consistent");
  }

  // Cancel the public receives
  if (verbose) {
    CCTK_INFO("Cancelling all public receives...");
  }
  while (not public_recvs.empty()) {
    list<transmission_t *>::iterator const tmi = public_recvs.begin();
    transmission_t *const tm = *tmi;
    MPI_Cancel(&tm->request);
    public_recvs.erase(tmi);
  }

  if (verbose) {
    CCTK_INFO("Destroying global scatter object");
  }

  assert(public_recvs.empty());
  assert(recvs.empty());
  assert(sends.empty());
}

// Communication tree
int nsiblings() { return 2; }
// Get process id of my first child (the other children are
// consecutive)
int child() { return dist::rank() * nsiblings() + 1; }
// Get process id if my parent
int parent() {
  return dist::rank() == 0 ? -1 : (dist::rank() - 1) / nsiblings();
}

// Check whether we should tell our parent that all our (ours and
// our childrens') messages have been sent
void scatter_t::maybe_send_did_send() {
  DECLARE_CCTK_PARAMETERS;
  int need_receive = 0;
  for (int p = child(); p < child() + nsiblings(); ++p) {
    if (p < dist::size())
      ++need_receive;
  }
  if (did_send_all and did_receive_sent_all == need_receive) {
    int const p = parent();
    if (p < 0) {
      // We are root; now broadcast this to all
      set_did_receive_all_sent_all();
    } else {
      if (verbose) {
        CCTK_INFO(
            "[Telling our parent that we and our children sent all our data]");
      }
      to_parent.state = fragdesc_t::state_sent_all;
      MPI_Request request;
      MPI_Isend(&to_parent, to_parent.num_ints(), MPI_INT, p, tag_desc,
                dist::comm(), &request);
    }
  }
}

// Broadcast "all message have been sent" to all our children
void scatter_t::send_all_did_send() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO("[Telling our children that all data have been sent]");
  }
  to_children.state = fragdesc_t::state_all_sent_all;
  for (int p = child(); p < child() + nsiblings(); ++p) {
    if (p < dist::size()) {
      MPI_Request request;
      MPI_Isend(&to_children, to_children.num_ints(), MPI_INT, p, tag_desc,
                dist::comm(), &request);
    }
  }
}

// We (this process) sent all our messages
void scatter_t::set_did_send_all() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO("[We sent all our data]");
  }
  assert(not did_send_all);
  did_send_all = true;
  maybe_send_did_send();
}

// One of our children (and all its children) sent all their
// messages
void scatter_t::set_did_receive_sent_all() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO("[One of our children sent all their data]");
  }
  assert(did_receive_sent_all < 2);
  ++did_receive_sent_all;
  maybe_send_did_send();
}

// Our parent broadcast that all messages have been sent
void scatter_t::set_did_receive_all_sent_all() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    CCTK_INFO("[Everything has been sent]");
  }
  assert(not did_receive_all_sent_all);
  did_receive_all_sent_all = true;
  send_all_did_send();
}

// Dispatch upon a state change message
void scatter_t::handle_state_transition(fragdesc_t const &fd) {
  switch (fd.state) {
  case fragdesc_t::state_sent_all:
    set_did_receive_sent_all();
    break;
  case fragdesc_t::state_all_sent_all:
    set_did_receive_all_sent_all();
    break;
  default:
    assert(0);
  }
}

void scatter_t::send(fragdesc_t const &fd, void const *const data) {
  DECLARE_CCTK_PARAMETERS;

  assert(data);

  if (verbose) {
    char *const fullname = CCTK_FullName(fd.varindex);
    CCTK_VInfo(CCTK_THORNSTRING, "Sending data read from variable %s, reflevel "
                                 "%d, map %d, component %d, timelevel %d",
               fullname, fd.reflevel, fd.map, fd.src_component, fd.timelevel);
    free(fullname);
  }

  list<transmission_t *> tosend = split_for_sending(fd, data);

  while (not tosend.empty()) {
    list<transmission_t *>::iterator const tmi = tosend.begin();
    transmission_t *const tm = *tmi;
    tosend.erase(tmi);

    // Choose tag
    int const proc = tm->fragdesc.process;
    int tag;
    while (true) {
      tag = busy_tags.at(proc).allocate_tag();
      if (tag >= 0)
        break;
      // No tag was available; wait for some progress before trying
      // again (until the receiver has received some of our
      // messages)
      do_some_work(true);
    }
    tm->fragdesc.tag = tag_data_min + tag;

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "   Sending data to process %d with tag %d...",
                 tm->fragdesc.process, tm->fragdesc.tag);
    }

    // Send descriptor and data
    MPI_Request req;
    MPI_Isend(&tm->fragdesc, tm->fragdesc.num_ints(), MPI_INT,
              tm->fragdesc.process, tag_desc, dist::comm(), &req);
    assert(tm->request == MPI_REQUEST_NULL);
    MPI_Isend(&tm->data[0], tm->fragdesc.npoints(), tm->fragdesc.datatype(),
              tm->fragdesc.process, tm->fragdesc.tag, dist::comm(),
              &tm->request);
    sends.push_back(tm);
  }

  // Do some work (if some is available)
  do_some_work();

  if (verbose) {
    CCTK_INFO("Done sending");
  }
}

void scatter_t::post_public_recvs() {
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("Posting public receives");
  }

  // Post public receives
  while ((int)public_recvs.size() < num_public_recvs) {
    transmission_t *const tm = new transmission_t;
    MPI_Irecv(&tm->fragdesc, tm->fragdesc.num_ints(), MPI_INT, MPI_ANY_SOURCE,
              tag_desc, dist::comm(), &tm->request);
    public_recvs.push_back(tm);
  }
}

void scatter_t::do_some_work(bool const do_wait) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("Checking for progress");
  }

  // Set up an array of all open requests
  vector<MPI_Request> requests;
  vector<list<transmission_t *>::iterator> iterators;
  int const npublic_recvs = public_recvs.size();
  int const nrecvs = recvs.size();
  int const nsends = sends.size();
  int const nrequests = npublic_recvs + nrecvs + nsends;
  requests.reserve(nrequests);
  iterators.reserve(nrequests);
  for (list<transmission_t *>::iterator tmi = public_recvs.begin();
       tmi != public_recvs.end(); ++tmi) {
    transmission_t const *const tm = *tmi;
    requests.push_back(tm->request);
    iterators.push_back(tmi);
  }
  for (list<transmission_t *>::iterator tmi = recvs.begin(); tmi != recvs.end();
       ++tmi) {
    transmission_t const *const tm = *tmi;
    requests.push_back(tm->request);
    iterators.push_back(tmi);
  }
  for (list<transmission_t *>::iterator tmi = sends.begin(); tmi != sends.end();
       ++tmi) {
    transmission_t const *const tm = *tmi;
    requests.push_back(tm->request);
    iterators.push_back(tmi);
  }
  assert((int)requests.size() == nrequests);

  // Wait (or test) for some open requests
  int outcount;
  vector<int> indices(nrequests);
  vector<MPI_Status> statuses(nrequests);
  if (do_wait) {
    if (verbose) {
      CCTK_INFO("Waiting for some progress...");
    }
    MPI_Waitsome(requests.size(), &requests[0], &outcount, &indices[0],
                 &statuses[0]);
  } else {
    if (verbose) {
      CCTK_INFO("Testing for some progress...");
    }
    MPI_Testsome(requests.size(), &requests[0], &outcount, &indices[0],
                 &statuses[0]);
  }
  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Completed %d transmissions", outcount);
  }

  // Process all completed requests
  for (int n = 0; n < outcount; ++n) {
    int const idx = indices.at(n);
    MPI_Status const &status = statuses.at(n);

    if (idx < npublic_recvs) {

      // We received a new descriptor
      int const source = status.MPI_SOURCE;
      assert(source >= 0 and source < dist::size());

      list<transmission_t *>::iterator const tmi = iterators.at(idx);
      transmission_t *const tm = *tmi;
      assert(tm->request != MPI_REQUEST_NULL);
      tm->request = MPI_REQUEST_NULL;
      public_recvs.erase(tmi);

      if (tm->fragdesc.state != fragdesc_t::state_normal) {
        handle_state_transition(tm->fragdesc);
      } else {

        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Receiving data from process %d with tag %d", source,
                     tm->fragdesc.tag);
        }

        // Prepare receiving the data
        assert(tm->fragdesc.process == dist::rank());
        tm->data.resize(tm->fragdesc.npoints() * tm->fragdesc.vartypesize());
        bytes_allocated += tm->data.size();
        assert(tm->request == MPI_REQUEST_NULL);
        MPI_Irecv(&tm->data[0], tm->fragdesc.npoints(), tm->fragdesc.datatype(),
                  source, tm->fragdesc.tag, dist::comm(), &tm->request);
        recvs.push_back(tm);
        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "   Current buffer size: %td bytes",
                     bytes_allocated);
        }

        post_public_recvs();
      }

    } else if (idx < npublic_recvs + nrecvs) {

      // We completed receiving a dataset; process it
      list<transmission_t *>::iterator const tmi = iterators.at(idx);
      transmission_t *const tm = *tmi;
      assert(tm->request != MPI_REQUEST_NULL);
      tm->request = MPI_REQUEST_NULL;

      if (verbose) {
        char *const fullname = CCTK_FullName(tm->fragdesc.varindex);
        CCTK_VInfo(CCTK_THORNSTRING, "Completed receiving data for variable %s",
                   fullname);
        free(fullname);
      }

      write_data(tm);
      bytes_allocated -= tm->data.size();
      delete tm;
      if (verbose) {
        CCTK_VInfo(CCTK_THORNSTRING, "   Current buffer size: %td bytes",
                   bytes_allocated);
      }
      recvs.erase(tmi);
      ++num_received;

    } else {

      // We completed sending a dataset; forget it
      list<transmission_t *>::iterator const tmi = iterators.at(idx);
      transmission_t *const tm = *tmi;
      assert(tm->request != MPI_REQUEST_NULL);
      tm->request = MPI_REQUEST_NULL;

      if (verbose) {
        char *const fullname = CCTK_FullName(tm->fragdesc.varindex);
        CCTK_VInfo(CCTK_THORNSTRING, "Completed sending data for variable %s",
                   fullname);
        free(fullname);
      }

      int const proc = tm->fragdesc.process;
      int const tag = tm->fragdesc.tag - tag_data_min;
      busy_tags.at(proc).free_tag(tag);

      bytes_allocated -= tm->data.size();
      delete tm;
      if (verbose) {
        CCTK_VInfo(CCTK_THORNSTRING, "   Current buffer size: %td bytes",
                   bytes_allocated);
      }
      sends.erase(tmi);
      ++num_sent;
    }
  }
}

list<scatter_t::transmission_t *>
scatter_t::split_for_sending(fragdesc_t const &fd, void const *const data) {
  DECLARE_CCTK_PARAMETERS;

  int const groupindex = CCTK_GroupIndexFromVarI(fd.varindex);
  assert(groupindex >= 0);
  int const varoffset = fd.varindex - CCTK_FirstVarIndexI(groupindex);
  assert(varoffset >= 0 and varoffset <= fd.varindex);

  gh const &hh = *Carpet::arrdata.at(groupindex).at(fd.map).hh;
  dh const &dd = *Carpet::arrdata.at(groupindex).at(fd.map).dd;

  ibbox const &baseext = hh.baseextents.AT(fd.mglevel).AT(fd.reflevel);

  ibbox const mybox(baseext.lower() + fd.imin * baseext.stride(),
                    baseext.lower() + fd.imax * baseext.stride(),
                    baseext.stride());

  // Distribute load: don't start with sending to process 0, instead
  // start with sending to (self+1)
  int coffset;
  {
    int const nlc = hh.local_components(fd.reflevel);
    if (nlc > 0) {
      // find component of next processor
      int c0 = hh.get_component(fd.reflevel, 0);
      coffset = c0 + nlc;
    } else {
      // estimate component of next processor
      int const nc = hh.components(fd.reflevel);
      coffset = nc * dist::rank() / dist::size();
    }
  }

  ibset done;
  list<transmission_t *> tosend;
  dh::light_cboxes const &light_cbox =
      dd.light_boxes.at(fd.mglevel).at(fd.reflevel);
  for (int ci = 0; ci < hh.components(fd.reflevel); ++ci) {
    int const c = (ci + coffset) % hh.components(fd.reflevel);
    dh::light_dboxes const &light_box = light_cbox.at(c);
    ibbox const &intr = light_box.interior;
    ibbox const ovlp = mybox & intr;
    assert((ovlp & done).empty());
    done += ovlp;

    if (not ovlp.empty()) {
      ibbox const &box = ovlp;
      transmission_t *const tm = new transmission_t;
      tm->fragdesc = fd;
      tm->fragdesc.imin = (box.lower() - baseext.lower()) / baseext.stride();
      tm->fragdesc.imax = (box.upper() - baseext.lower()) / baseext.stride();
      tm->fragdesc.component = c;
      tm->fragdesc.process = hh.processor(fd.reflevel, c);

      ptrdiff_t const vartypesize = tm->fragdesc.vartypesize();
      ptrdiff_t const ni = tm->fragdesc.imax[0] - tm->fragdesc.imin[0] + 1;
      ptrdiff_t const nj = tm->fragdesc.imax[1] - tm->fragdesc.imin[1] + 1;
      ptrdiff_t const nk = tm->fragdesc.imax[2] - tm->fragdesc.imin[2] + 1;
      ptrdiff_t const di = vartypesize;
      ptrdiff_t const dj = di * ni;
      ptrdiff_t const dk = dj * nj;
      ptrdiff_t const np = dk * nk;
      tm->data.resize(np);
      bytes_allocated += tm->data.size();
      if (verbose) {
        CCTK_VInfo(CCTK_THORNSTRING, "   Current buffer size: %td bytes",
                   bytes_allocated);
      }
      char *const dst = &tm->data[0];

      ptrdiff_t const nis = fd.imax[0] - fd.imin[0] + 1;
      ptrdiff_t const njs = fd.imax[1] - fd.imin[1] + 1;
      ptrdiff_t const nks = fd.imax[2] - fd.imin[2] + 1;
      ptrdiff_t const dis = vartypesize;
      ptrdiff_t const djs = vartypesize * nis;
      ptrdiff_t const dks = vartypesize * nis * njs;
      ptrdiff_t const nps = vartypesize * nis * njs * nks;
      ptrdiff_t const i0s = tm->fragdesc.imin[0] - fd.imin[0];
      ptrdiff_t const j0s = tm->fragdesc.imin[1] - fd.imin[1];
      ptrdiff_t const k0s = tm->fragdesc.imin[2] - fd.imin[2];
      ptrdiff_t const ind0s = i0s * dis + j0s * djs + k0s * dks;
      assert(ind0s < nps);
      char const *const src = &((char const *)data)[ind0s];

#pragma omp parallel for collapse(2)
      for (ptrdiff_t k = 0; k < nk; ++k) {
        for (ptrdiff_t j = 0; j < nj; ++j) {
          ptrdiff_t const ind = j * dj + k * dk;
          ptrdiff_t const inds = j * djs + k * dks;
          assert(ind + ni * vartypesize <= np);
          assert(inds + ni * vartypesize <= nps);
          memcpy(&dst[ind], &src[inds], ni * vartypesize);
        }
      }

      tosend.push_back(tm);
    }
  }
  // Don't enforce this -- mesh refinement boundaries have ghosts
  // that do not overlap with any interior
  // assert(done == mybox);

  return tosend;
}

void scatter_t::write_data(transmission_t *const tm) {
  fragdesc_t const &fd = tm->fragdesc;
  int const groupindex = CCTK_GroupIndexFromVarI(fd.varindex);
  assert(groupindex >= 0);
  int const varoffset = fd.varindex - CCTK_FirstVarIndexI(groupindex);
  assert(varoffset >= 0 and varoffset <= fd.varindex);

  gh const &hh = *Carpet::arrdata.at(groupindex).at(fd.map).hh;
  dh const &dd = *Carpet::arrdata.at(groupindex).at(fd.map).dd;
  ggf &ff = *Carpet::arrdata.at(groupindex).at(fd.map).data.at(varoffset);
  int const lc = hh.get_local_component(fd.reflevel, fd.component);
  gdata &data = *ff.data_pointer(fd.timelevel, fd.reflevel, lc, fd.mglevel);

  ibbox const &baseext = hh.baseextents.AT(fd.mglevel).AT(fd.reflevel);

  ibbox const mybox(baseext.lower() + fd.imin * baseext.stride(),
                    baseext.lower() + fd.imax * baseext.stride(),
                    baseext.stride());

  dh::light_dboxes const &light_box =
      dd.light_boxes.at(fd.mglevel).at(fd.reflevel).at(fd.component);
  ibbox const &intr = light_box.interior;
  assert(mybox.is_contained_in(intr));
  ibbox const &extr = light_box.exterior;

  ptrdiff_t const vartypesize = fd.vartypesize();
  ptrdiff_t const ni = fd.imax[0] - fd.imin[0] + 1;
  ptrdiff_t const nj = fd.imax[1] - fd.imin[1] + 1;
  ptrdiff_t const nk = fd.imax[2] - fd.imin[2] + 1;
  ptrdiff_t const di = vartypesize;
  ptrdiff_t const dj = di * ni;
  ptrdiff_t const dk = dj * nj;
  ptrdiff_t const np = dk * nk;
  char const *const src = (char const *)&tm->data[0];

  ivect const lbnd = (extr.lower() - baseext.lower()) / baseext.stride();
  ivect const lsh = extr.shape() / baseext.stride();
  ptrdiff_t const nid = lsh[0];
  ptrdiff_t const njd = lsh[1];
  ptrdiff_t const nkd = lsh[2];
  ptrdiff_t const did = vartypesize;
  ptrdiff_t const djd = vartypesize * nid;
  ptrdiff_t const dkd = vartypesize * nid * njd;
  ptrdiff_t const npd = vartypesize * nid * njd * nkd;
  ptrdiff_t const i0d = fd.imin[0] - lbnd[0];
  ptrdiff_t const j0d = fd.imin[1] - lbnd[1];
  ptrdiff_t const k0d = fd.imin[2] - lbnd[2];
  ptrdiff_t const ind0d = i0d * did + j0d * djd + k0d * dkd;
  assert(data.has_storage());
  char *const dst = &((char *)data.storage())[ind0d];

#pragma omp parallel for collapse(2)
  for (ptrdiff_t k = 0; k < nk; ++k) {
    for (ptrdiff_t j = 0; j < nj; ++j) {
      ptrdiff_t const indd = j * djd + k * dkd;
      ptrdiff_t const ind = j * dj + k * dk;
      assert(indd + ni * vartypesize <= npd);
      assert(ind + ni * vartypesize <= np);
      memcpy(&dst[indd], &src[ind], ni * vartypesize);
    }
  }
}

} // end namespace CarpetIOF5
