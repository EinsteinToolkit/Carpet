#include <algorithm>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "dh.hh"
#include "mpi_string.hh"
#include "region.hh"



namespace CarpetLib
{
  
  using namespace std;
  
  
  
  vector <string>
  gather_string (MPI_Comm const comm,
                 int const root,
                 string const & data)
  {
    // Get my rank
    int rank;
    MPI_Comm_rank (comm, & rank);
    
    if (rank == root) {
      
      // Get the total number of processors
      int num_procs;
      MPI_Comm_size (comm, & num_procs);
      
      // Gather the lengths of the data strings
      int const length = data.length();
      vector <int> lengths (num_procs);
      
      MPI_Gather (const_cast <int *> (& length), 1, MPI_INT,
                  & lengths.front(), 1, MPI_INT,
                  root, comm);
      
      // Allocate space for all data strings
      vector <int> offsets (num_procs + 1);
      offsets.AT(0) = 0;
      for (int n = 0; n < num_procs; ++ n)
      {
        offsets.AT(n + 1) = offsets.AT(n) + lengths.AT(n);
      }
      int const total_length = offsets.AT(num_procs);
      vector <char> alldata_buffer (total_length);
      
      // Gather all data strings
      MPI_Gatherv (const_cast <char *> (data.c_str()), length, MPI_CHAR,
                   & alldata_buffer.front(),
                   const_cast <int *> (& lengths.front()),
                   const_cast <int *> (& offsets.front()),
                   MPI_CHAR,
                   root, comm);
      
      // Convert data buffer with C strings to C++ strings
      vector <string> alldata  (num_procs);
      for (int n = 0; n < num_procs; ++ n)
      {
        alldata.AT(n) =
          string (& alldata_buffer.AT (offsets.AT(n)), lengths.AT(n));
      }
      
      return alldata;
      
    } else {
      
      // Gather the lengths of the data strings
      int const length = data.length();
      
      MPI_Gather (const_cast <int *> (& length), 1, MPI_INT,
                  NULL, 1, MPI_INT,
                  root, comm);
      
      // Gather all data strings
      MPI_Gatherv (const_cast <char *> (data.c_str()), length, MPI_CHAR,
                   NULL, NULL, NULL, MPI_CHAR,
                   root, comm);
      
      // Convert data buffer with C strings to C++ strings
      vector <string> alldata;
      
      return alldata;
      
    }
  }
  
  
  
  vector <string>
  allgather_string (MPI_Comm const comm,
                    string const & data)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (comm, & num_procs);
    
    // Gather the lengths of the data strings
    int const length = data.length();
    vector <int> lengths (num_procs);
    
    MPI_Allgather (const_cast <int *> (& length), 1, MPI_INT,
                   & lengths.front(), 1, MPI_INT,
                   comm);
    
    // Allocate space for all data strings
    vector <int> offsets (num_procs + 1);
    offsets.AT(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets.AT(n + 1) = offsets.AT(n) + lengths.AT(n);
    }
    int const total_length = offsets.AT(num_procs);
    vector <char> alldata_buffer (total_length);
    
    // Gather all data strings
    MPI_Allgatherv (const_cast <char *> (data.c_str()), length, MPI_CHAR,
                    & alldata_buffer.front(),
                    const_cast <int *> (& lengths.front()),
                    const_cast <int *> (& offsets.front()),
                    MPI_CHAR,
                    comm);
    
    // Convert data buffer with C strings to C++ strings
    vector <string> alldata  (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      alldata.AT(n) =
        string (& alldata_buffer.AT (offsets.AT(n)), lengths.AT(n));
    }
    
    return alldata;
  }
  
  
  
  vector <string>
  alltoallv_string (MPI_Comm const comm,
                    vector<string> const & data)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (comm, & num_procs);
    
    // Exchange the lengths of the data strings
    vector <int> lengths_in (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      lengths_in.AT(n) = data.AT(n).length();
    }
    vector <int> lengths (num_procs);
    MPI_Alltoall (& lengths_in.front(), 1, MPI_INT,
                  & lengths.front(), 1, MPI_INT,
                  comm);
    
    // Allocate space for all data strings
    vector <int> offsets_in (num_procs + 1);
    offsets_in.AT(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets_in.AT(n + 1) = offsets_in.AT(n) + lengths_in.AT(n);
    }
    int const total_length_in = offsets_in.AT(num_procs);
    vector <char> alldata_buffer_in (total_length_in);
    
    vector <int> offsets (num_procs + 1);
    offsets.AT(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets.AT(n + 1) = offsets.AT(n) + lengths.AT(n);
    }
    int const total_length = offsets.AT(num_procs);
    vector <char> alldata_buffer (total_length);
    
    // Convert C++ strings to data buffer with C strings
    for (int n = 0; n < num_procs; ++ n)
    {
      memcpy (& alldata_buffer_in.AT (offsets_in.AT(n)),
              data.AT(n).c_str(),
              lengths_in.AT(n));
    }
    
    // Exchange all data strings
    MPI_Alltoallv (& alldata_buffer_in.front(),
                   & lengths_in.front(), & offsets_in.front(), MPI_CHAR,
                   & alldata_buffer.front(),
                   & lengths.front(), & offsets.front(), MPI_CHAR,
                   comm);
    
    // Convert data buffer with C strings to C++ strings
    vector <string> alldata  (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      alldata.AT(n) =
        string (& alldata_buffer.AT (offsets.AT(n)), lengths.AT(n));
    }
    
    return alldata;
  }



  string
  broadcast_string (MPI_Comm const comm,
                    int const root,
                    string const & data)
  {
    // Get my rank
    int rank;
    MPI_Comm_rank (comm, & rank);
    
    if (rank == root) {
      
      // Broadcast the length of the data string
      int const length = data.length();
      MPI_Bcast (const_cast <int *> (& length), 1, MPI_INT, root, comm);
      
      // Broadcast data string
      char const * const buf = data.c_str();
      MPI_Bcast (const_cast <char *> (buf), length, MPI_CHAR, root, comm);
      
      // Return original string
      return data;
      
    } else {
      
      // Broadcast the length of the data string
      int length;
      MPI_Bcast (& length, 1, MPI_INT, root, comm);
      
      // Allocate space for data string
      vector <char> data_buffer (length);
      
      // Broadcast data string
      char * const buf = & data_buffer.front();
      MPI_Bcast (buf, length, MPI_CHAR, root, comm);
      
      // Convert data buffer with C strings to C++ strings
      string const result = string (& data_buffer.front(), length);
      
      return result;
      
    }
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  template <typename T>
  vector <vector <T> >
  allgatherv (MPI_Comm comm,
              vector <T> const & data)
  {
    // cerr << "QQQ: allgatherv[0]" << endl;
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (comm, & num_procs);
    
    // Exchange the sizes of the data vectors
    int const size_in = data.size();
    assert (size_in >= 0);
    vector <int> sizes_out (num_procs);
    // cerr << "QQQ: allgatherv[1] size_in=" << size_in << endl;
    MPI_Allgather (const_cast <int *> (& size_in), 1, MPI_INT,
                   & sizes_out.front(), 1, MPI_INT,
                   comm);
    // cerr << "QQQ: allgatherv[2]" << endl;
    
    // Allocate space for all data vectors
    vector <int> offsets_out (num_procs + 1);
    offsets_out.AT(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      assert (sizes_out.AT(n) >= 0);
      offsets_out.AT(n + 1) = offsets_out.AT(n) + sizes_out.AT(n);
      assert (offsets_out.AT(n + 1) >= 0);
    }
    int const total_length_out = offsets_out.AT(num_procs);
    vector <T> alldata_buffer_out (total_length_out);
    
    // Exchange all data vectors
    T const dummy;
    MPI_Datatype const type = mpi_datatype (dummy);
    int datatypesize;
    MPI_Type_size (type, &datatypesize);
    // cerr << "QQQ: allgatherv[3] total_length_out=" << total_length_out << " datatypesize=" << datatypesize << endl;
#if 0
    MPI_Allgatherv (const_cast <T *> (& data.front()),
                    size_in, type,
                    & alldata_buffer_out.front(),
                    & sizes_out.front(), & offsets_out.front(), type,
                    comm);
#else
    int const typesize = sizeof(T);
    for (int n = 0; n < num_procs; ++ n)
    {
      sizes_out.AT(n) *= typesize;
      offsets_out.AT(n) *= typesize;
    }
    MPI_Allgatherv (const_cast <T *> (& data.front()),
                    size_in * typesize, MPI_CHAR,
                    & alldata_buffer_out.front(),
                    & sizes_out.front(), & offsets_out.front(), MPI_CHAR,
                    comm);
    for (int n = 0; n < num_procs; ++ n)
    {
      sizes_out.AT(n) /= typesize;
      offsets_out.AT(n) /= typesize;
    }
#endif
    // cerr << "QQQ: allgatherv[4]" << endl;
    
    // Convert data buffer to vectors
    vector <vector <T> > alldata_out (num_procs);
    {
      typename vector <T>::const_iterator p = alldata_buffer_out.begin();
      for (int n = 0; n < num_procs; ++ n)
      {
        typename vector <T>::const_iterator const pold = p;
        advance (p, sizes_out.AT(n));
        alldata_out.AT(n).assign (pold, p);
      }
      assert (p == alldata_buffer_out.end());
    }
    
    // cerr << "QQQ: allgatherv[5]" << endl;
    return alldata_out;
  }
  
  
  
  template <typename T>
  vector <T>
  alltoall (MPI_Comm const comm,
            vector <T> const & data)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (comm, & num_procs);
    
    // Allocate space for all data
    vector <T> alldata (num_procs);
    
    // Exchange all data vectors
    T const dummy;
    MPI_Datatype const type = mpi_datatype (dummy);
    MPI_Alltoall (& data.front(), 1, type,
                  & alldata.front(), 1, type,
                  comm);
    
    return alldata;
  }
  
  
  
  template <typename T>
  vector <vector <T> >
  alltoallv (MPI_Comm const comm,
             vector <vector <T> > const & data)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (comm, & num_procs);
    
    // Exchange the sizes of the data vectors
    vector <int> sizes_in (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      sizes_in.AT(n) = data.AT(n).size();
    }
    vector <int> sizes_out (num_procs);
    MPI_Alltoall (& sizes_in.front(), 1, MPI_INT,
                  & sizes_out.front(), 1, MPI_INT,
                  comm);
    
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
    
    // Convert data buffer to vectors
    vector <vector <T> > alldata_out (num_procs);
    {
      typename vector <T>::const_iterator p = alldata_buffer_out.begin();
      for (int n = 0; n < num_procs; ++ n)
      {
        typename vector <T>::const_iterator const pold = p;
        advance (p, sizes_out.AT(n));
        alldata_out.AT(n).assign (pold, p);
      }
    }
    
    return alldata_out;
  }
  
  
  
  template <typename T>
  vector <T>
  alltoallv1 (MPI_Comm const comm,
              vector <vector <T> > const & data)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (comm, & num_procs);
    
    // Exchange the sizes of the data vectors
    vector <int> sizes_in (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      sizes_in.AT(n) = data.AT(n).size();
    }
    vector <int> sizes_out (num_procs);
    // cerr << "QQQ: alltoallv1[1]" << endl;
    MPI_Alltoall (& sizes_in.front(), 1, MPI_INT,
                  & sizes_out.front(), 1, MPI_INT,
                  comm);
    // cerr << "QQQ: alltoallv1[2]" << endl;
    
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
    // cerr << "QQQ: alltoallv1[3]" << endl;
    MPI_Alltoallv (& alldata_buffer_in.front(),
                   & sizes_in.front(), & offsets_in.front(), type,
                   & alldata_buffer_out.front(),
                   & sizes_out.front(), & offsets_out.front(), type,
                   comm);
    // cerr << "QQQ: alltoallv1[4]" << endl;
#endif
    
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
    int const tag = 4711;
    vector <MPI_Request> reqs (2 * num_procs);
    int nreqs = 0;
    // cerr << "QQQ: alltoallv1[5]" << endl;
    for (int n = 0; n < num_procs; ++ n)
    {
      if (sizes_out.AT(n) > 0) {
        MPI_Irecv (& alldata_buffer_out.AT(offsets_out.AT(n)),
                   sizes_out.AT(n),
                   type,
                   n, tag, comm, & reqs.AT(nreqs));
        ++ nreqs;
      }
    }
    // cerr << "QQQ: alltoallv1[6]" << endl;
    for (int n = 0; n < num_procs; ++ n)
    {
      if (sizes_in.AT(n) > 0) {
        MPI_Isend (const_cast <T *> (& data.AT(n).front()),
                   sizes_in.AT(n),
                   type,
                   n, tag, comm, & reqs.AT(nreqs));
        ++ nreqs;
      }
    }
    // cerr << "QQQ: alltoallv1[7]" << endl;
    MPI_Waitall (nreqs, & reqs.front(), MPI_STATUSES_IGNORE);
    // cerr << "QQQ: alltoallv1[8]" << endl;
    
    return alldata_buffer_out;
  }
  
  
  
  template
  vector <vector <dh::light_dboxes> >
  allgatherv (MPI_Comm comm,
              vector <dh::light_dboxes> const & data);
  
  template
  vector <sendrecv_pseudoregion_t>
  alltoallv1 (MPI_Comm comm,
              vector <vector <sendrecv_pseudoregion_t> > const & data);
  
} // namespace CarpetLib
