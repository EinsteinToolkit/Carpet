#include <cassert>
#include <cstring>
//#include <iostream>
//#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "functions.hh"



namespace Carpet
{
  
  using namespace std;
  
  
  
  vector <string>
  AllGatherString (MPI_Comm const world,
                   string const & data)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (world, & num_procs);
    
    // Gather the lengths of the data strings
    int const length = data.length();
    vector <int> lengths (num_procs);
    
    MPI_Allgather (const_cast <int *> (& length), 1, MPI_INT,
                   & lengths.front(), 1, MPI_INT,
                   world);
    
    // Allocate space for all data strings
    vector <int> offsets (num_procs + 1);
    offsets.at(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets.at(n + 1) = offsets.at(n) + lengths.at(n);
    }
    int const total_length = offsets.at(num_procs);
    
    // Gather all data strings
    vector <char> alldata_buffer (total_length);
    
    MPI_Allgatherv (const_cast <char *> (data.c_str()), length, MPI_CHAR,
                    & alldata_buffer.front(),
                    const_cast <int *> (& lengths.front()),
                    const_cast <int *> (& offsets.front()),
                    MPI_CHAR,
                    world);
    
    // Convert data buffer with C strings to C++ strings
    vector <string> alldata  (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      alldata.at(n) =
        string (& alldata_buffer.at (offsets.at(n)), lengths.at(n));
    }
    
    return alldata;
  }
  
} // namespace Carpet
