#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Version.h"

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "CarpetIOHDF5.hh"

/* some macros for HDF5 group names */
#define METADATA_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"


namespace CarpetIOHDF5
{

using namespace std;
using namespace Carpet;

// structure describing a single dataset of an HDF5 file to read from
typedef struct
{
  char *datasetname;

  int vindex;
  int timelevel;
  int mglevel;
  int reflevel;
  int rank;
  int *shape;    // [rank]
  int *iorigin;  // [rank]
} dataset_t;

// structure describing the contents of an HDF5 file to read from
typedef struct
{
  int num_mglevels;
  int num_reflevels;
  int parameter_len;
  int cctk_iteration;
  int main_loop_index;
  CCTK_REAL global_time;
  CCTK_REAL delta_time;
  CCTK_REAL *mgleveltimes;  // [num_mglevels*num_reflevels]

  char *filename;
  hid_t file;
  int num_datasets;
  list<dataset_t> datasets;  // [num_datasets]
  int num_ints;              // total number of integers in datasets[]
} file_t;

static file_t infile = {0, 0, 0, 0, 0, 0, 0, NULL, NULL, -1, -1,
                        list<dataset_t> (), 0};


static int OpenFile (const char *basefilename, file_t *file, int called_from);
static int RecoverVariables (cGH* cctkGH, file_t *file);
static herr_t ReadMetadata (hid_t group, const char *objectname, void *arg);

// callback for I/O parameter parsing routine
static void SetFlag (int vindex, const char* optstring, void* arg);


// Register with the Cactus Recovery Interface
int CarpetIOHDF5_RecoverParameters (void)
{
  int retval = IOUtil_RecoverParameters (Recover, ".h5", "HDF5");

  return (retval);
}


// close a checkpoint/filereader file after recovering grid variables
int CarpetIOHDF5_CloseFile (void)
{
  DECLARE_CCTK_PARAMETERS


  if (infile.num_datasets <= 0)
  {
    return (0);
  }
  infile.num_datasets = -1;

  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "closing file '%s' after recovery",
                infile.filename);
  }

  if (infile.file >= 0)
  {
    HDF5_ERROR (H5Fclose (infile.file));
    infile.file = -1;
  }
  free (infile.filename);
  delete[] infile.mgleveltimes;

  for (list<dataset_t>::iterator dataset = infile.datasets.begin ();
       dataset != infile.datasets.end ();
       dataset++)
  {
    free (dataset->datasetname);
    delete[] dataset->shape;
    delete[] dataset->iorigin;
  }
  infile.datasets.clear ();

  return (0);
}

static int OpenFile (const char *basefilename, file_t *file, int called_from)
{
  hid_t dset = -1;
  const int myproc = CCTK_MyProc (NULL);
  DECLARE_CCTK_PARAMETERS


  // generate filename for an unchunked checkpoint file */
  file->filename = IOUtil_AssembleFilename (NULL, basefilename, "", ".h5",
                                            called_from, 0, 1);
  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "opening %s file '%s'",
                called_from == CP_RECOVER_PARAMETERS ? "checkpoint" : "input",
                file->filename);
  }

  if (myproc == 0)
  {
    // try to open the file (prevent HDF5 error messages if it fails)
    H5E_BEGIN_TRY
    {
      file->file = H5Fopen (file->filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    } H5E_END_TRY;

    if (file->file >= 0)
    {
      if (called_from == CP_RECOVER_PARAMETERS)
      {
        HDF5_ERROR (dset = H5Dopen (file->file,
                                    METADATA_GROUP"/"ALL_PARAMETERS));
        ReadAttribute (dset, "carpet_reflevels", file->num_reflevels);
        ReadAttribute (dset, "numberofmgtimes", file->num_mglevels);
        ReadAttribute (dset, "GH$iteration", file->cctk_iteration);
        ReadAttribute (dset, "main loop index", file->main_loop_index);
        ReadAttribute (dset, "carpet_global_time", file->global_time);
        ReadAttribute (dset, "carpet_delta_time", file->delta_time);
        file->parameter_len = H5Dget_storage_size (dset) + 1;
        assert (file->parameter_len > 1);
      }

      file->num_ints = 0;
      HDF5_ERROR (H5Giterate (file->file, "/", NULL, ReadMetadata, file));
    }
    file->num_datasets = file->datasets.size ();
    if (file->num_datasets <= 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "No valid HDF5 file '%s' found", file->filename);
      free (file->filename);
    }
  }

  // broadcast integer variables
  int *intbuffer = new int[7];
  intbuffer[0] = file->num_datasets;
  intbuffer[1] = file->num_mglevels;
  intbuffer[2] = file->num_reflevels;
  intbuffer[3] = file->cctk_iteration;
  intbuffer[4] = file->main_loop_index;
  intbuffer[5] = file->parameter_len;
  intbuffer[6] = file->num_ints;
  MPI_Bcast (intbuffer, 7, MPI_INT, 0, MPI_COMM_WORLD);
  file->num_datasets    = intbuffer[0];
  file->num_mglevels    = intbuffer[1];
  file->num_reflevels   = intbuffer[2];
  file->cctk_iteration  = intbuffer[3];
  file->main_loop_index = intbuffer[4];
  file->parameter_len   = intbuffer[5];
  file->num_ints        = intbuffer[6];
  delete[] intbuffer;

  // return if no valid checkpoint could be found
  if (file->num_datasets <= 0)
  {
    return (-1);
  }

  // serialize the dataset list metadata into a single MPI_INT buffer
  intbuffer = new int[file->num_ints];
  if (myproc == 0)
  {
    for (list<dataset_t>::iterator dataset = file->datasets.begin ();
         dataset != file->datasets.end ();
         dataset++)
    {
      *intbuffer++ = dataset->vindex;
      *intbuffer++ = dataset->timelevel;
      *intbuffer++ = dataset->mglevel;
      *intbuffer++ = dataset->reflevel;
      *intbuffer++ = dataset->rank;

      for (int i = 0; i < dataset->rank; i++)
      {
        *intbuffer++ = dataset->shape[i];
        *intbuffer++ = dataset->iorigin[i];
      }
    }
    intbuffer -= file->num_ints;
  }

  // broadcast the serialized dataset list metadata
  MPI_Bcast (intbuffer, file->num_ints, MPI_INT, 0, MPI_COMM_WORLD);

  // build the dataset list on non-I/O processors
  if (myproc != 0)
  {
    for (int i = 0; i < file->num_datasets; i++)
    {
      dataset_t dataset;


      dataset.vindex    = *intbuffer++;
      dataset.timelevel = *intbuffer++;
      dataset.mglevel   = *intbuffer++;
      dataset.reflevel  = *intbuffer++;
      dataset.rank      = *intbuffer++;

      dataset.shape     = new int[dataset.rank];
      dataset.iorigin   = new int[dataset.rank];
      for (int j = 0; j < dataset.rank; j++)
      {
        dataset.shape[j]   = *intbuffer++;
        dataset.iorigin[j] = *intbuffer++;
      }
      dataset.datasetname = NULL;

      file->datasets.push_back (dataset);
    }
    intbuffer -= file->num_ints;
  }
  delete[] intbuffer;

  if (called_from == FILEREADER_DATA)
  {
    return (0);
  }

  // leave space at the end for global_time and delta_time
  // so that all double variables can be broadcasted in one go
  int num_times = file->num_mglevels*file->num_reflevels + 2;
  file->mgleveltimes = new CCTK_REAL[num_times];
  if (myproc == 0)
  {
    // FIXME: should store all mgleveltimes in a single contiguous array
    //        to get rid of this loop and save some attributes
    for (int i = 0; i < file->num_mglevels; i++)
    {
      char buffer[32];

      snprintf (buffer, sizeof (buffer), "mgleveltimes %d", i);
      ReadAttribute (dset, buffer, file->mgleveltimes + i*file->num_reflevels,
                     file->num_reflevels);
    }
  }

  // broadcast double variables
  file->mgleveltimes[num_times - 2] = file->global_time;
  file->mgleveltimes[num_times - 1] = file->delta_time;
  MPI_Bcast (file->mgleveltimes, num_times, CARPET_MPI_REAL, 0, MPI_COMM_WORLD);
  file->global_time = file->mgleveltimes[num_times - 2];
  file->delta_time = file->mgleveltimes[num_times - 1];

  char *parameters = new char[file->parameter_len];
  if (myproc == 0)
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Recovering parameters from checkpoint");

    HDF5_ERROR (H5Dread (dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         parameters));
  }
  if (dset >= 0)
  {
    HDF5_ERROR (H5Dclose (dset));
  }

  // broadcast char variables
  MPI_Bcast (parameters, file->parameter_len, MPI_CHAR, 0, MPI_COMM_WORLD);

  // recover parameters
  IOUtil_SetAllParameters (parameters);
  delete[] parameters;

  return (0);
}


int Recover (cGH* cctkGH, const char *basefilename, int called_from)
{
  DECLARE_CCTK_PARAMETERS
  int retval = 0;


  assert (called_from == CP_RECOVER_PARAMETERS ||
          called_from == CP_RECOVER_DATA ||
          called_from == FILEREADER_DATA);

  if (called_from == CP_RECOVER_PARAMETERS ||
      called_from == FILEREADER_DATA)
  {
    // open the file, read and broadcast its metadata information
    // for CP_RECOVER_PARAMETERS: also recover all parameters
    retval = OpenFile (basefilename, &infile, called_from);

    if (called_from == CP_RECOVER_PARAMETERS || retval)
    {
      return (retval == 0 ? 1 : -1);
    }
  }

  // can only proceed with a valid checkpoint file from here on
  assert (infile.num_datasets > 0);

  // set global Cactus/Carpet variables
  if (called_from == CP_RECOVER_DATA)
  {
    global_time = infile.global_time;
    delta_time = infile.delta_time;
    CCTK_SetMainLoopIndex (infile.main_loop_index);

    cctkGH->cctk_iteration = infile.cctk_iteration;
    cctkGH->cctk_time = infile.mgleveltimes[mglevel*infile.num_reflevels +
                                            reflevel];
  }

  // now recover all grid variables on current mglevel and reflevel
  retval = RecoverVariables (cctkGH, &infile);

  if (called_from == CP_RECOVER_DATA)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "restarting simulation at iteration %d (physical time %g)",
                cctkGH->cctk_iteration, (double) cctkGH->cctk_time);
  }
  else
  {
    // FIXME: keep filereader input file open for all mglevels, reflevels
    //        this doesn't work now because we don't know the number of
    //        mglevels and reflevels
    //        also: there may be multiple files to be read in, and they
    //              must not share the same data structures
    CarpetIOHDF5_CloseFile ();
  }

  // now synchronize all variables, in sets of groups of the same vartype
  vector<group_set> groups;

  for (int group = 0; group < CCTK_NumGroups(); group++) {
    if (CCTK_NumVarsInGroupI (group) > 0
        && CCTK_QueryGroupStorageI (cctkGH, group)) {

      group_set newset;
      const int firstvar = CCTK_FirstVarIndexI (group);
      newset.vartype = CCTK_VarTypeI (firstvar);
      assert (newset.vartype >= 0);
      int c;
      for (c = 0; c < groups.size(); c++) {
        if (newset.vartype == groups[c].vartype) {
          break;
        }
      }
      if (c == groups.size()) {
        groups.push_back (newset);
      }
      groups[c].members.push_back (group);
    }
  }

  return (retval);
}


int ReadVar (const cGH* const cctkGH, const int vindex,
             const hid_t dataset, vector<ibset> &regions_read,
             const int called_from_recovery)
{
  DECLARE_CCTK_PARAMETERS;

  const int group = CCTK_GroupIndexFromVarI (vindex);
  assert (group>=0 && group<(int)Carpet::arrdata.size());
  char *fullname = CCTK_FullName (vindex);
  const int n0 = CCTK_FirstVarIndexI(group);
  assert (n0>=0 && n0<CCTK_NumVars());
  const int var = vindex - n0;
  assert (var>=0 && var<CCTK_NumVars());
  int tl = 0;

  bool did_read_something = false;

  // Stuff needed for Recovery

  void *h5data = NULL;

  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "  reading '%s'", fullname);
  }

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, group))
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable \"%s\" because it has no storage",
                fullname);
    free (fullname);
    return 0;
  }

  const int grouptype = CCTK_GroupTypeI(group);
  if ((grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY) && reflevel > 0)
  {
    free (fullname);
    return 0;
  }

  const int gpdim = CCTK_GroupDimI(group);
  const int vartype = CCTK_VarTypeI(vindex);

  int intbuffer[2 + 2*dim];
  int &group_timelevel = intbuffer[0],
      &amr_level       = intbuffer[1];
  int *amr_origin      = &intbuffer[2],
      *amr_dims        = &intbuffer[2+dim];

  if (CCTK_MyProc(cctkGH)==0)
  {
    // get dataset dimensions
    hid_t dataspace;
    HDF5_ERROR (dataspace = H5Dget_space (dataset));
    int rank = (int) H5Sget_simple_extent_ndims (dataspace);
    assert (0 < rank && rank <= dim);
    assert ((grouptype == CCTK_SCALAR ? gpdim+1 : gpdim) == rank);
    vector<hsize_t> shape(rank);
    HDF5_ERROR (H5Sget_simple_extent_dims (dataspace, &shape[0], NULL));
    HDF5_ERROR (H5Sclose (dataspace));

    for (int i = 0; i < dim; i++)
    {
      amr_dims[i]   = 1;
      amr_origin[i] = 0;
    }
    int datalength = 1;
    for (int i = 0; i < rank; i++)
    {
      datalength *= shape[i];
      amr_dims[i] = shape[rank-i-1];
    }

    const hid_t datatype = h5DataType (cctkGH, vartype, 0);

    //cout << "datalength: " << datalength << " rank: " << rank << "\n";
    //cout << shape[0] << " " << shape[1] << " " << shape[2] << "\n";

    // to do: read in an allocate with correct datatype

    h5data = malloc (CCTK_VarTypeSize (vartype) * datalength);
    HDF5_ERROR (H5Dread (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         h5data));

    ReadAttribute (dataset, "level", amr_level);
    ReadAttribute (dataset, "iorigin", amr_origin, rank);

    if(called_from_recovery)
    {
      ReadAttribute(dataset,"group_timelevel", group_timelevel);
      // in old days (before February 2005) timelevels used to be stored
      // as negative numbers
      group_timelevel = abs (group_timelevel);
    }
  } // MyProc == 0

  MPI_Bcast (intbuffer, sizeof (intbuffer) / sizeof (*intbuffer), MPI_INT, 0, dist::comm);

#ifdef CARPETIOHDF5_DEBUG
  cout << "amr_level: " << amr_level << " reflevel: " << reflevel << endl;
#endif

  if (amr_level == reflevel)
  {
    // Traverse all components on all levels
    BEGIN_MAP_LOOP (cctkGH, grouptype)
    {
      BEGIN_COMPONENT_LOOP (cctkGH, grouptype)
      {
        did_read_something = true;

        ggf* ff = 0;

        assert (var < (int)arrdata.at(group).at(Carpet::map).data.size());
        ff = (ggf*)arrdata.at(group).at(Carpet::map).data.at(var);

        if(called_from_recovery) tl = group_timelevel;

        gdata* const data = (*ff) (tl, reflevel, component, mglevel);

        // Create temporary data storage on processor 0
        vect<int,dim> str = vect<int,dim>(maxreflevelfact/reflevelfact);

        if(grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY)
          str = vect<int,dim> (1);

        vect<int,dim> lb = vect<int,dim>::ref(amr_origin) * str;
        vect<int,dim> ub
          = lb + (vect<int,dim>::ref(amr_dims) - 1) * str;

        gdata* const tmp = data->make_typed (vindex);


        cGroup cgdata;
        int ierr = CCTK_GroupData(group,&cgdata);
        assert(ierr==0);
        //cout << "lb_before: " << lb << endl;
        //cout << "ub_before: " << ub << endl;
        if (cgdata.disttype == CCTK_DISTRIB_CONSTANT) {
#ifdef CARPETIOHDF5_DEBUG
          cout << "CCTK_DISTRIB_CONSTANT: " << fullname << endl;
#endif
          assert(grouptype == CCTK_ARRAY || grouptype == CCTK_SCALAR);
          if (grouptype == CCTK_SCALAR)
          {
            lb[0] = arrdata.at(group).at(Carpet::map).hh->processors().at(reflevel).at(component);
            ub[0] = arrdata.at(group).at(Carpet::map).hh->processors().at(reflevel).at(component);
            for(int i=1;i<dim;i++)
            {
              lb[i]=0;
              ub[i]=0;
            }
          }
          else
          {
            const int newlb = lb[gpdim-1] +
              (ub[gpdim-1]-lb[gpdim-1]+1)*
              (arrdata.at(group).at(Carpet::map).hh->processors().at(reflevel).at(component));
            const int newub = ub[gpdim-1] +
              (ub[gpdim-1]-lb[gpdim-1]+1)*
              (arrdata.at(group).at(Carpet::map).hh->processors().at(reflevel).at(component));
            lb[gpdim-1] = newlb;
            ub[gpdim-1] = newub;
          }
#ifdef CARPETIOHDF5_DEBUG
          cout << "lb: " << lb << endl;
          cout << "ub: " << ub << endl;
#endif
        }
        const bbox<int,dim> ext(lb,ub,str);

#ifdef CARPETIOHDF5_DEBUG
        cout << "ext: " << ext << endl;
#endif

        if (CCTK_MyProc(cctkGH)==0)
        {
          tmp->allocate (ext, 0, h5data);
        }
        else
        {
          tmp->allocate (ext, 0);
        }

        // Initialise with what is found in the file -- this does
        // not guarantee that everything is initialised.
        const bbox<int,dim> overlap = tmp->extent() & data->extent();
        regions_read.at(Carpet::map) |= overlap;

#ifdef CARPETIOHDF5_DEBUG
        cout << "working on component: " << component << endl;
        cout << "tmp->extent " << tmp->extent() << endl;
        cout << "data->extent " << data->extent() << endl;
        cout << "overlap " << overlap << endl;
        cout << "-----------------------------------------------------" << endl;
#endif

        // FIXME: is this barrier really necessary ??
        MPI_Barrier(MPI_COMM_WORLD);

        // Copy into grid function
        for (comm_state state(vartype); !state.done(); state.step())
        {
          data->copy_from (state, tmp, overlap);
        }


        // Delete temporary copy
        delete tmp;


      } END_COMPONENT_LOOP;

      if (called_from_recovery)
      {
        arrdata.at(group).at(Carpet::map).tt->set_time(reflevel,mglevel,
                (CCTK_REAL) ((cctkGH->cctk_time - cctk_initial_time)
                             / (delta_time * mglevelfact)) );
      }


    } END_MAP_LOOP;

  } // if amr_level == reflevel

  free (h5data);
  free (fullname);

  return did_read_something;
}


static int InputVarAs (const cGH* const cctkGH, const int vindex,
                       const char* const alias)
{
  DECLARE_CCTK_PARAMETERS;

  char *fullname = CCTK_FullName (vindex);
  const int group = CCTK_GroupIndexFromVarI (vindex);
  assert (group>=0 && group<(int)Carpet::arrdata.size());

  int want_dataset = 0;
  bool did_read_something = false;
  int ndatasets = 0;
  hid_t dataset = 0;

  char datasetname[1024];

  // Check for storage
  if (! CCTK_QueryGroupStorageI(cctkGH, group))
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot input variable \"%s\" because it has no storage",
                fullname);
    free (fullname);
    return 0;
  }

  const int grouptype = CCTK_GroupTypeI(group);
  const int rl = grouptype==CCTK_GF ? reflevel : 0;
  //cout << "want level " << rl << " reflevel " << reflevel << endl;

  // Find the input directory
  const char* const myindir = in_dir;

  // Invent a file name
  ostringstream filenamebuf;
  filenamebuf << myindir << "/" << alias << in_extension;
  string filenamestr = filenamebuf.str();
  const char * const filename = filenamestr.c_str();

  hid_t reader = -1;

  // Read the file only on the root processor
  if (CCTK_MyProc(cctkGH)==0)
  {
    // Open the file
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Opening file \"%s\"", filename);
    }
    reader = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (reader<0)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open file \"%s\" for reading", filename);
    }
    assert (reader>=0);
    // get the number of datasets in the file
    ndatasets=GetnDatasets(reader);
  }

  vector<ibset> regions_read(Carpet::maps);

  // Broadcast number of datasets
  MPI_Bcast (&ndatasets, 1, MPI_INT, 0, dist::comm);
  assert (ndatasets>=0);


  for (int datasetid=0; datasetid<ndatasets; ++datasetid)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Handling dataset #%d", datasetid);
    }

    // Read data
    if (CCTK_MyProc(cctkGH)==0)
    {
      GetDatasetName(reader,datasetid,datasetname);
      //         cout << datasetname << "\n";

      HDF5_ERROR (dataset = H5Dopen (reader, datasetname));
    }

    if (CCTK_MyProc(cctkGH)==0)
    {
      // Read data
      char * name;
      ReadAttribute (dataset, "name", name);
      //        cout << "dataset name is " << name << endl;
      if (CCTK_Equals (verbose, "full") && name)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Dataset name is \"%s\"", name);
      }
      want_dataset = name && CCTK_EQUALS(name, fullname);
      free (name);
    } // myproc == 0

    MPI_Bcast (&want_dataset, 1, MPI_INT, 0, dist::comm);

    if(want_dataset)
    {
      did_read_something = ReadVar(cctkGH,vindex,dataset,regions_read,0);
    } // want_dataset

  } // loop over datasets

  // Close the file
  if (CCTK_MyProc(cctkGH)==0)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Closing file");
    }
    HDF5_ERROR (H5Fclose(reader));
    reader=-1;
  }

  // Was everything initialised?
  if (did_read_something)
  {
    for (int m=0; m<Carpet::maps; ++m)
    {
      dh& thedd = *arrdata.at(group).at(m).dd;
      ibset all_exterior;
      for (size_t c=0; c<thedd.boxes.at(rl).size(); ++c)
      {
        all_exterior |= thedd.boxes.at(mglevel).at(rl).at(c).exterior;
      }
      if (regions_read.at(m) != all_exterior)
      {
#ifdef CARPETIOHDF5_DEBUG
        cout << "read: " << regions_read.at(m) << endl
             << "want: " << all_exterior << endl;
#endif
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Variable \"%s\" could not be initialised from file -- the file may be missing data",
                    fullname);
      }
    }
  } // if did_read_something
  //        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,"stop!");

  return did_read_something ? 0 : -1;
}


int CarpetIOHDF5_ReadData (const cGH* const cctkGH)
{
  int retval = 0;
  DECLARE_CCTK_PARAMETERS


  int numvars = CCTK_NumVars ();
  vector<bool> flags (numvars);

  if (CCTK_TraverseString (in_vars, SetFlag, &flags, CCTK_GROUP_OR_VAR) < 0)
  {
    CCTK_VWarn (strict_io_parameter_check ? 0 : 1,
                __LINE__, __FILE__, CCTK_THORNSTRING,
                "error while parsing parameter 'IOHDF5::in_vars'");
  }

  for (int vindex = 0; vindex < numvars && retval == 0; vindex++)
  {
    if (flags.at (vindex))
    {
      retval = InputVarAs (cctkGH, vindex, CCTK_VarName (vindex));
    }
  }

  return (retval);
}


static herr_t ReadMetadata (hid_t group, const char *objectname, void *arg)
{
  file_t *file = (file_t *) arg;
  dataset_t dataset;
  hid_t dset, dspace;
  H5G_stat_t object_info;


  // we are interested only in datasets
  HDF5_ERROR (H5Gget_objinfo (group, objectname, 0, &object_info));
  if (object_info.type != H5G_DATASET)
  {
    return (0);
  }

  dataset.datasetname = strdup (objectname);
  assert (dataset.datasetname);

  HDF5_ERROR (dset = H5Dopen (group, objectname));
  char *varname = NULL;
  ReadAttribute (dset, "name", varname);
  dataset.vindex = CCTK_VarIndex (varname);
  free (varname);
  ReadAttribute (dset, "level", dataset.timelevel);
  ReadAttribute (dset, "carpet_mglevel", dataset.mglevel);
  ReadAttribute (dset, "carpet_reflevel", dataset.reflevel);
  HDF5_ERROR (dspace = H5Dget_space (dset));
  HDF5_ERROR (dataset.rank = H5Sget_simple_extent_ndims (dspace));
  dataset.iorigin = new int[dataset.rank];
  ReadAttribute (dset, "iorigin", dataset.iorigin, dataset.rank);
  hsize_t *shape = new hsize_t[dataset.rank];
  HDF5_ERROR (H5Sget_simple_extent_dims (dspace, shape, NULL));
  dataset.shape = new int[dataset.rank];
  for (int i = 0; i < dataset.rank; i++)
  {
    dataset.shape[i] = shape[i];
  }
  delete[] shape;
  HDF5_ERROR (H5Sclose (dspace));
  HDF5_ERROR (H5Dclose (dset));

  // add this dataset to our list and count the number of int elements
  file->datasets.push_back (dataset);
  file->num_ints += 5 + 2*dataset.rank;

  return (0);
}


static int RecoverVariables (cGH* cctkGH, file_t *file)
{
  DECLARE_CCTK_PARAMETERS
  int myproc = CCTK_MyProc (cctkGH);
  hid_t dset = -1;


  // Use refinement levels parameter from checkpointing file ?
  if (use_reflevels_from_checkpoint)
  {
    char buffer[32];

    snprintf (buffer, sizeof (buffer), "%d", file->num_reflevels);
    CCTK_ParameterSet ("refinement_levels", "CarpetRegrid", buffer);

    CCTK_VInfo (CCTK_THORNSTRING, "Using %i reflevels read from checkpoint "
                "file. Ignoring value in parameter file.", file->num_reflevels);
  }

#if 0
  double leveltime = MPI_Wtime();
  double comparetime = MPI_Wtime();
  static double totaltime;
  if (reflevel == 0) totaltime = 0;

  comparetime = MPI_Wtime() - comparetime;
  // cout << "Time for string comparison: " << comparetime << endl;
  // cout << "I have for this reflevel " << refleveldatasetnamelist.size() << endl;
#endif

  CCTK_VInfo (CCTK_THORNSTRING,
              "reading grid variables on mglevel %d reflevel %d",
              mglevel, reflevel);

  int num_vars = CCTK_NumVars ();
  for (list<dataset_t>::iterator dataset = file->datasets.begin ();
       dataset != file->datasets.end ();
       dataset++)
  {
    // only recover grid variables for the current mglevel/reflevel
    if (dataset->mglevel != mglevel || dataset->reflevel != reflevel)
    {
      continue;
    }

    if (dataset->vindex <  0 || dataset->vindex >= num_vars)
    {
      if (myproc == 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Ignoring dataset '%s' (invalid variable name)",
                    dataset->datasetname);
      }
      continue;
    }

    if (myproc == 0)
    {
      HDF5_ERROR (dset = H5Dopen (file->file, dataset->datasetname));
      assert (dset >= 0);
    }

    vector<ibset> regions_read (Carpet::maps);
    int did_read_something = ReadVar (cctkGH, dataset->vindex, dset,
                                      regions_read, 1);
    MPI_Bcast (&did_read_something, 1, MPI_INT, 0, dist::comm);

    if (dset >= 0)
    {
      HDF5_ERROR (H5Dclose (dset));
    }

  }

#if 0
  leveltime = MPI_Wtime() - leveltime;
  totaltime += leveltime;

  if (CCTK_Equals (verbose, "full"))
  {
    cout << "Timers: leveltime: " << leveltime << " totaltime: " << totaltime << endl;
  }
#endif

  return (0);
}


static void SetFlag (int vindex, const char* optstring, void* arg)
{
  if (optstring)
  {
    char *fullname = CCTK_FullName (vindex);
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Option string '%s' will be ignored for HDF5 input of "
                "variable '%s'", optstring, fullname);
    free (fullname);
  }
  vector<bool>& flags = *(vector<bool>*)arg;
  flags.at(vindex) = true;
}


} // namespace CarpetIOHDF5
