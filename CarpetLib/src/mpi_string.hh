#ifndef MPI_STRING_HH
#define MPI_STRING_HH

#include <cctk.h>

#include <string>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include "defs.hh"
#include "region.hh"

namespace CarpetLib {
using namespace std;

// String communication

vector<string> gather_string(MPI_Comm comm, int root, string const &data);

vector<string> allgather_string(MPI_Comm comm, string const &data);

vector<string> alltoallv_string(MPI_Comm comm, vector<string> const &data);

string broadcast_string(MPI_Comm comm, int root, string const &data);

void send_string(MPI_Comm const comm, int const destination,
                   string const &data);

string recv_string(MPI_Comm const comm, int const source);

// Arbitrary datatypes

template <typename T>
vector<vector<T> > allgatherv(MPI_Comm comm, vector<T> const &data);

template <typename T>
vector<T> allgatherv1(MPI_Comm comm, vector<T> const &data);

template <typename T> vector<T> alltoall(MPI_Comm comm, vector<T> const &data);

template <typename T>
vector<vector<T> > alltoallv(MPI_Comm comm, vector<vector<T> > const &data);

template <typename T>
vector<T> alltoallv1(MPI_Comm comm, vector<vector<T> > const &data);

//////////////////////////////////////////////////////////////////////////////

template <typename T>
vector<vector<T> > allgatherv(MPI_Comm comm, vector<T> const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Exchange the sizes of the data vectors
  int const size_in = data.size();
  assert(size_in >= 0);
  vector<int> sizes_out(num_procs);
  MPI_Allgather(const_cast<int *>(&size_in), 1, MPI_INT, &sizes_out.front(), 1,
                MPI_INT, comm);

  // Allocate space for all data vectors
  vector<int> offsets_out(num_procs + 1);
  offsets_out.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    assert(sizes_out.AT(n) >= 0);
    offsets_out.AT(n + 1) = offsets_out.AT(n) + sizes_out.AT(n);
    assert(offsets_out.AT(n + 1) >= 0);
  }
  int const total_length_out = offsets_out.AT(num_procs);
  vector<T> alldata_buffer_out(total_length_out);

  // Exchange all data vectors
  T dummy;
  MPI_Datatype const type = dist::mpi_datatype(dummy);
  int datatypesize;
  MPI_Type_size(type, &datatypesize);
// datatypesize=" << datatypesize << endl;
#if 0
    MPI_Allgatherv (const_cast <T *> (& data.front()),
                    size_in, type,
                    & alldata_buffer_out.front(),
                    & sizes_out.front(), & offsets_out.front(), type,
                    comm);
#else
  int const typesize = sizeof(T);
  for (int n = 0; n < num_procs; ++n) {
    sizes_out.AT(n) *= typesize;
    offsets_out.AT(n) *= typesize;
  }
  MPI_Allgatherv(const_cast<T *>(&data.front()), size_in * typesize, MPI_CHAR,
                 &alldata_buffer_out.front(), &sizes_out.front(),
                 &offsets_out.front(), MPI_CHAR, comm);
  for (int n = 0; n < num_procs; ++n) {
    sizes_out.AT(n) /= typesize;
    offsets_out.AT(n) /= typesize;
  }
#endif

  // Convert data buffer to vectors
  vector<vector<T> > alldata_out(num_procs);
  {
    typename vector<T>::const_iterator p = alldata_buffer_out.begin();
    for (int n = 0; n < num_procs; ++n) {
      typename vector<T>::const_iterator const pold = p;
      advance(p, sizes_out.AT(n));
      alldata_out.AT(n).assign(pold, p);
    }
    assert(p == alldata_buffer_out.end());
  }

  return alldata_out;
}

template <typename T>
vector<T> allgatherv1(MPI_Comm comm, vector<T> const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Exchange the sizes of the data vectors
  int const size_in = data.size();
  assert(size_in >= 0);
  vector<int> sizes_out(num_procs);
  MPI_Allgather(const_cast<int *>(&size_in), 1, MPI_INT, &sizes_out.front(), 1,
                MPI_INT, comm);

  // Allocate space for all data vectors
  vector<int> offsets_out(num_procs + 1);
  offsets_out.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    assert(sizes_out.AT(n) >= 0);
    offsets_out.AT(n + 1) = offsets_out.AT(n) + sizes_out.AT(n);
    assert(offsets_out.AT(n + 1) >= 0);
  }
  int const total_length_out = offsets_out.AT(num_procs);
  vector<T> alldata_buffer_out(total_length_out);

  // Exchange all data vectors
  T dummy;
  MPI_Datatype const type = dist::mpi_datatype(dummy);
  int datatypesize;
  MPI_Type_size(type, &datatypesize);
#if 0
    MPI_Allgatherv (const_cast <T *> (& data.front()),
                    size_in, type,
                    & alldata_buffer_out.front(),
                    & sizes_out.front(), & offsets_out.front(), type,
                    comm);
#else
  int const typesize = sizeof(T);
  for (int n = 0; n < num_procs; ++n) {
    sizes_out.AT(n) *= typesize;
    offsets_out.AT(n) *= typesize;
  }
  MPI_Allgatherv(const_cast<T *>(&data.front()), size_in * typesize, MPI_CHAR,
                 &alldata_buffer_out.front(), &sizes_out.front(),
                 &offsets_out.front(), MPI_CHAR, comm);
  for (int n = 0; n < num_procs; ++n) {
    sizes_out.AT(n) /= typesize;
    offsets_out.AT(n) /= typesize;
  }
#endif

  return alldata_buffer_out;
}

template <typename T>
vector<T> alltoall(MPI_Comm const comm, vector<T> const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Allocate space for all data
  vector<T> alldata(num_procs);

  // Exchange all data vectors
  T const dummy;
  MPI_Datatype const type = mpi_datatype(dummy);
  MPI_Alltoall(&data.front(), 1, type, &alldata.front(), 1, type, comm);

  return alldata;
}

template <typename T>
vector<vector<T> > alltoallv(MPI_Comm const comm,
                             vector<vector<T> > const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Exchange the sizes of the data vectors
  vector<int> sizes_in(num_procs);
  for (int n = 0; n < num_procs; ++n) {
    sizes_in.AT(n) = data.AT(n).size();
  }
  vector<int> sizes_out(num_procs);
  MPI_Alltoall(&sizes_in.front(), 1, MPI_INT, &sizes_out.front(), 1, MPI_INT,
               comm);

  // Copy vectors to data buffer
  vector<int> offsets_in(num_procs + 1);
  offsets_in.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    offsets_in.AT(n + 1) = offsets_in.AT(n) + sizes_in.AT(n);
  }
  int const total_length_in = offsets_in.AT(num_procs);
  vector<T> alldata_buffer_in;
  alldata_buffer_in.reserve(total_length_in);
  for (int n = 0; n < num_procs; ++n) {
    alldata_buffer_in.insert(alldata_buffer_in.end(), data.AT(n).begin(),
                             data.AT(n).end());
  }

  // Allocate space for all data vectors
  vector<int> offsets_out(num_procs + 1);
  offsets_out.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    offsets_out.AT(n + 1) = offsets_out.AT(n) + sizes_out.AT(n);
  }
  int const total_length_out = offsets_out.AT(num_procs);
  vector<T> alldata_buffer_out(total_length_out);

  // Exchange all data vectors
  T const dummy;
  MPI_Datatype const type = mpi_datatype(dummy);
  MPI_Alltoallv(&alldata_buffer_in.front(), &sizes_in.front(),
                &offsets_in.front(), type, &alldata_buffer_out.front(),
                &sizes_out.front(), &offsets_out.front(), type, comm);

  // Convert data buffer to vectors
  vector<vector<T> > alldata_out(num_procs);
  {
    typename vector<T>::const_iterator p = alldata_buffer_out.begin();
    for (int n = 0; n < num_procs; ++n) {
      typename vector<T>::const_iterator const pold = p;
      advance(p, sizes_out.AT(n));
      alldata_out.AT(n).assign(pold, p);
    }
  }

  return alldata_out;
}

template <typename T>
vector<T> alltoallv1(MPI_Comm const comm, vector<vector<T> > const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Exchange the sizes of the data vectors
  vector<int> sizes_in(num_procs);
  for (int n = 0; n < num_procs; ++n) {
    sizes_in.AT(n) = data.AT(n).size();
  }
  vector<int> sizes_out(num_procs);
  MPI_Alltoall(&sizes_in.front(), 1, MPI_INT, &sizes_out.front(), 1, MPI_INT,
               comm);

#if 0
    // Copy vectors to data buffer
    vector <int> offsets_in (num_procs + 1);
    offsets_in.AT(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets_in.AT(n + 1) = offsets_in.AT(n) + sizes_in.AT(n);
    }
    int const total_length_in = offsets_in.AT(num_procs);
    vector <T> alldata_buffer_in;
    alldata_buffer_in.reserve (total_length_in);
    for (int n = 0; n < num_procs; ++ n)
    {
      alldata_buffer_in.insert (alldata_buffer_in.end(),
                                data.AT(n).begin(), data.AT(n).end());
    }
    
    // Allocate space for all data vectors
    vector <int> offsets_out (num_procs + 1);
    offsets_out.AT(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets_out.AT(n + 1) = offsets_out.AT(n) + sizes_out.AT(n);
    }
    int const total_length_out = offsets_out.AT(num_procs);
    vector <T> alldata_buffer_out (total_length_out);
    
    // Exchange all data vectors
    T const dummy;
    MPI_Datatype const type = mpi_datatype (dummy);
    MPI_Alltoallv (& alldata_buffer_in.front(),
                   & sizes_in.front(), & offsets_in.front(), type,
                   & alldata_buffer_out.front(),
                   & sizes_out.front(), & offsets_out.front(), type,
                   comm);
#endif

  // Allocate space for all data vectors
  vector<int> offsets_out(num_procs + 1);
  offsets_out.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    offsets_out.AT(n + 1) = offsets_out.AT(n) + sizes_out.AT(n);
  }
  int const total_length_out = offsets_out.AT(num_procs);
  vector<T> alldata_buffer_out(total_length_out);

  // Exchange all data vectors
  T dummy;
  MPI_Datatype const type = dist::mpi_datatype(dummy);
  int const tag = 4711;
  vector<MPI_Request> reqs(2 * num_procs);
  int nreqs = 0;
  for (int n = 0; n < num_procs; ++n) {
    if (sizes_out.AT(n) > 0) {
      MPI_Irecv(&alldata_buffer_out.AT(offsets_out.AT(n)), sizes_out.AT(n),
                type, n, tag, comm, &reqs.AT(nreqs));
      ++nreqs;
    }
  }
  for (int n = 0; n < num_procs; ++n) {
    if (sizes_in.AT(n) > 0) {
      MPI_Isend(const_cast<T *>(&data.AT(n).front()), sizes_in.AT(n), type, n,
                tag, comm, &reqs.AT(nreqs));
      ++nreqs;
    }
  }
  MPI_Waitall(nreqs, &reqs.front(), MPI_STATUSES_IGNORE);

  return alldata_buffer_out;
}

} // namespace CarpetLib

#endif // #ifndef MPI_STRING_HH
