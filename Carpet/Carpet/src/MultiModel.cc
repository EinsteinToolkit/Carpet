#include <cassert>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include <cctk.h>

#include <functions.hh>
#include <mpi_string.hh>



namespace Carpet
{
  
  using namespace std;
  
  
  
  vector <string> models;             // Model id to model name
  std::map <string, int> model_map;   // Model name to model id
  vector <int> model_ids;             // Process to model id
  vector <vector <int> > model_procs; // Model id to processes
  
  
  
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
    // Get the total number of processes
    int num_procs;
    MPI_Comm_size (world, & num_procs);
    int my_proc;
    MPI_Comm_rank (world, & my_proc);
    
    // Gather all model names
    models = allgather_string (world, model);
    
    // Map model strings to small integers
    int num_models = 0;
    model_ids.resize (num_procs);
    model_map.clear ();
    for (int n = 0; n < num_procs; ++ n)
    {
      if (model_map.find (models.AT(n)) != model_map.end())
      {
        model_ids.AT(n) = model_map[models.AT(n)];
      }
      else
      {
        model_map[models.AT(n)] = num_models;
        model_ids.AT(n) = num_models;
        ++ num_models;
      }
    }
    
    // Determine processes per model
    vector <int> num_model_procs (num_models, 0);
    for (int n = 0; n < num_procs; ++ n)
    {
      ++ num_model_procs.at (model_ids.AT(n));
    }
    
    model_procs.resize (num_models);
    for (int m = 0; m < num_models; ++ m)
    {
      model_procs.AT(m).reserve (num_model_procs.AT(m));
    }
    for (int n = 0; n < num_procs; ++ n)
    {
      model_procs.at (model_ids.AT(n)).push_back (n);
    }
    for (int m = 0; m < num_models; ++ m)
    {
      assert (static_cast<int> (model_procs.AT(m).size())
              == num_model_procs.AT(m));
    }
    
    // Create a new communicator for each model
    MPI_Comm_split (world, model_ids.AT(my_proc), my_proc, & comm);
    
    if (verbose)
    {
      CCTK_INFO ("Multi-Model listing:");
      for (int m = 0; m < num_models; ++ m)
      {
        cout << "   model " << m << ": \"" << models.AT(m) << "\"" << endl;
      }
      CCTK_INFO ("Multi-Model process distribution:");
      for (int n = 0; n < num_procs; ++ n)
      {
        int const m = model_ids.AT(n);
        bool const same_model_as_prev =
          n-1 >= 0 and model_ids.AT(n-1) == m;
        bool const same_model_as_next =
          n+1 < num_procs and model_ids.AT(n+1) == m;
        if (same_model_as_next) {
          if (same_model_as_prev) {
            // Output nothing
          } else {
            // This process has the same model as the next one:
            // output only a partial line
            cout << "   processes " << n << "-";
          }
        } else {
          if (same_model_as_prev) {
            // This process has the same model as the previous one:
            // finish a partial line
            cout << n << ": "
                 << "model " << m << " \"" << models.AT(m) << "\"" << endl;
          } else {
            cout << "   process " << n << ": "
                 << "model " << m << " \"" << models.AT(m) << "\"" << endl;
          }
        }
      }
      int const my_model = model_ids.AT(my_proc);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Multi-Model: This is process %d, model %d \"%s\"",
                  my_proc, my_model, model.c_str());
    }
  }
  
} // namespace Carpet
