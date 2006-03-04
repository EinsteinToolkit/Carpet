#include <cassert>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include "cctk.h"

#include "functions.hh"



namespace Carpet
{
  
  using namespace std;
  
  
  
  vector <string> models;             // Model id to model name
  std::map <string, int> model_map;   // Model name to model id
  vector <int> model_ids;             // Processor to model id
  vector <vector <int> > model_procs; // Model id to processors
  
  
  
  vector <string> Models () { return models; }
  std::map <string, int> ModelMap () { return model_map; }
  vector <int> ModelIds () { return model_ids; }
  vector <vector <int> > ModelProcs () { return model_procs; }
  
  string Models (int const id)
  {
    return models.at (id);
  }

  int ModelMap (string const name)
  {
    if (model_map.find (name) != model_map.end())
    {
      return model_map[name];
    }
    else
    {
      return -1;
    }
  }
  
  int ModelId (int const proc)
  {
    return model_ids.at (proc);
  }

  vector <int> ModelProcs (int const proc)
  {
    return model_procs.at (proc);
  }
  
  
  
  void
  SplitUniverse (MPI_Comm const world, string const model, MPI_Comm & comm,
                 bool const verbose)
  {
    // Get the total number of processors
    int num_procs;
    MPI_Comm_size (world, & num_procs);
    int my_proc;
    MPI_Comm_rank (world, & my_proc);
    
    // Gather the lengths of the model strings
    int const length = model.length();
    vector <int> lengths (num_procs);
    
    MPI_Allgather (const_cast <int *> (& length), 1, MPI_INT,
                   & lengths.front(), 1, MPI_INT,
                   world);
    
    // Allocate space for all model strings
    vector <int> offsets (num_procs + 1);
    offsets.at(0) = 0;
    for (int n = 0; n < num_procs; ++ n)
    {
      offsets.at(n + 1) = offsets.at(n) + lengths.at(n);
    }
    int const total_length = offsets.at(num_procs);
    
    // Gather all model strings
    vector <char> models_buffer (total_length);
    
    MPI_Allgatherv (const_cast <char *> (model.c_str()), length, MPI_CHAR,
                    & models_buffer.front(),
                    const_cast <int *> (& lengths.front()),
                    const_cast <int *> (& offsets.front()),
                    MPI_CHAR,
                    world);
    
    // Convert model name buffer with C strings to C++ strings
    models.resize (num_procs);
    for (int n = 0; n < num_procs; ++ n)
    {
      models.at(n)
        = string (& models_buffer.at (offsets.at(n)), lengths.at(n));
    }
    
    // Map model strings to small integers
    int num_models = 0;
    model_ids.resize (num_procs);
    model_map.clear ();
    for (int n = 0; n < num_procs; ++ n)
    {
      if (model_map.find (models.at(n)) != model_map.end())
      {
        model_ids.at(n) = model_map[models.at(n)];
      }
      else
      {
        model_map[models.at(n)] = num_models;
        model_ids.at(n) = num_models;
        ++ num_models;
      }
    }
    
    // Determine processors per model
    vector <int> num_model_procs (num_models, 0);
    for (int n = 0; n < num_procs; ++ n)
    {
      ++ num_model_procs.at (model_ids.at(n));
    }
    
    model_procs.resize (num_models);
    for (int m = 0; m < num_models; ++ m)
    {
      model_procs.at(m).reserve (num_model_procs.at(m));
    }
    for (int n = 0; n < num_procs; ++ n)
    {
      model_procs.at (model_ids.at(n)).push_back (n);
    }
    for (int m = 0; m < num_models; ++ m)
    {
      assert (static_cast<int> (model_procs.at(m).size())
              == num_model_procs.at(m));
    }
    
    // Create a new communicator for each model
    MPI_Comm_split (world, model_ids.at(my_proc), my_proc, & comm);
    
    if (verbose)
    {
      CCTK_INFO ("Multi-Model listing:");
      for (int m = 0; m < num_models; ++ m)
      {
        cout << "   model " << m << ": \"" << models.at(m) << "\"" << endl;
      }
      CCTK_INFO ("Multi-Model processor distribution:");
      for (int n = 0; n < num_procs; ++ n)
      {
        int const m = model_ids.at(n);
        cout << "   processor " << n << ": "
             << "model " << m << " \"" << models.at(m) << "\"" << endl;
      }
      int const my_model = model_ids.at(my_proc);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Multi-Model: This is processor %d, model %d \"%s\"",
                  my_proc, my_model, model.c_str());
    }
  }
  
} // namespace Carpet
