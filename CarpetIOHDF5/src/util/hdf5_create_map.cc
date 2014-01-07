 /*@@
   @file      hdf5_create_map.cc
   @date      Mon Jan  6 22:58:27 PST 2014
   @author    Roland Haas, Thomas Radke
   @desc
              Utility program to create a map file describing the location of
              each dataset in a set of files.
   @enddesc
 @@*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>

#include <hdf5.h>

using namespace std;

/*****************************************************************************
 *************************     Type Definitions    ***************************
 *****************************************************************************/
struct mapentry {
  int map, mglevel, reflevel, timestep, timelevel, component, rank;
  int iorigin[3], ioffset[3], ioffsetdenom[3], shape[3];
  double origin[3], delta[3], time;
  int nghosts[3];
  int vnameindex, filenameindex;
  char* objectname;
};

/*****************************************************************************
 *************************     Macro Definitions   ***************************
 *****************************************************************************/

// macro to check the return code of calls to the HDF5 library
#define CHECK_HDF5(fn_call) CHECK_HDF5_CALL(fn_call, #fn_call, __LINE__)
static int CHECK_HDF5_CALL(int _retval, const char *fn_call, int line) { \
  if (_retval < 0) {                                                     \
    cerr << "HDF5 call " << fn_call                                      \
         << " in file " << __FILE__ << " line " << line                  \
         << " returned error code " << _retval << endl;                  \
  }                                                                      \
  return retval;
}

/*****************************************************************************
 *************************       Global Data         *************************
 *****************************************************************************/

// output progress information
static bool verbose = false;

// output file id
static hid_t outfile = -1;

// variable names encountered
map<string, int> varnames;
int varnamesmaxvalue = -1;

// file names encountered
map<string, int> filenames;
int filenamesmaxvalue = -1;

// HDF5 type ids for the map entries
hid_t mapentry_t;

/*****************************************************************************
 *************************     Function Prototypes   *************************
 *****************************************************************************/
static herr_t ProcessDataset (hid_t group, const char *name, void *_file);
static void InitializeMappings ();


 /*@@
   @routine    main
   @date       Fri 17 October 2008
   @author     Thomas Radke
   @desc
               Evaluates command line options and opens the input files.
   @enddesc

   @var        argc
   @vdesc      number of command line parameters
   @vtype      int
   @vio        in
   @endvar
   @var        argv
   @vdesc      array of command line parameters
   @vtype      char* const []
   @vio        in
   @endvar
 @@*/
int main (int argc, char *const argv[])
{
  int i;
  bool help = false;
  int num_slab_options = 0;


  // evaluate command line parameters
  for (i = 1; i < argc; i++) {
    if (strcmp (argv[i], "--help") == 0) {
      help = true; break;
    } else if (strcmp (argv[i], "--verbose") == 0) {
      verbose = true;
    } else {
      break;
    }
  }

  /* give some help if called with incorrect number of parameters */
  if (help or i >= argc-1) {
    const string indent (strlen (argv[0]) + 1, ' ');
    cerr << endl << "   -----------------------"
         << endl << "   Carpet HDF5 Map Creator"
         << endl << "   -----------------------" << endl
         << endl
         << "Usage: " << endl
         << argv[0] << " [--help]" << endl
         << indent << "[--verbose]" << endl
         << indent << "<hdf5_infiles> <hdf5_mapfile>" << endl << endl
         << "  where" << endl
         << "    [--help]                         prints this help" << endl
         << "    [--verbose]                      output progress information" << endl
    return (help ? 0 : 1);
  }


  // open or create the output file
  const char *const outfilename = argv[argc-1];
  cout << endl << "  creating output file '" << outfilename << "'" << endl;
  H5E_BEGIN_TRY {
    outfile = H5Fcreate (outfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose (outfile);
  } H5E_END_TRY;
  outfile = H5Fopen (outfilename, H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0) {
    fprintf (stderr, "Could not open output map file \"%s\" for writing",
             outfilename);
    return 1;
  }

  // read existing string<->enum mapping from file
  InitializeMappings();

  // browse though input file(s)
  for (; i < argc-1; i++) {
    hid_t file;

    H5E_BEGIN_TRY {
      file = H5Fopen (argv[i], H5F_ACC_RDONLY, H5P_DEFAULT);
    } H5E_END_TRY;

    if (file < 0) {
      cerr << "Could not open input file '" << argv[i] << "'" << endl << endl;
      break;
    }

    cout << "  iterating through input file '" << argv[i] << "'..." << endl;
    CHECK_HDF5 (H5Giterate (file, "/", NULL, ProcessDataset, &file));

    // close file
    if (file >= 0) {
      CHECK_HDF5 (H5Fclose (file));
    }
  }

  // close output file
  CHECK_HDF5 (H5Fclose (outfile));

  cout << endl << "Done." << endl << endl;

  return (i == argc-1 ? 0 : 1);
}


 /*@@
   @routine    InitializeMappings
   @date       Mon Jan  6 23:09:32 PST 2014
   @author     Roland Haas
   @desc
               Read existing enum types from file and populate string to
               integer mappings from them.
   @enddesc

@@*/
static void InitializeMappings ()
{
  hid_t datatype;

  H5_BEGIN_TRY {
    datatype = H5Topen(outfile, "VARNAMES_ENUM", H5P_DEFAULT);
    if (datatype >= 0) {
      const int nmembers = CHECK_HDF5(H5Tget_nmembers(datatype));
      for(int i = 0 ; i < nmembers ; ++i) {
        int value;
        char* name;
        CHECK_HDF5(datatype, i, &value);
        name = CHECK_HDF5(datatype, i);
        assert (name);
        varnames[name] = value;
        varnamesmaxvalue = max(varnamesmaxvalue, value);
        free(name);
      }
      CHECK_HDF5(H5Tclose(datatype));
    }
  } H5_END_TRY;

  H5_BEGIN_TRY {
    datatype = H5Topen(outfile, "FILENAMES_ENUM", H5P_DEFAULT);
    if (datatype >= 0) {
      const int nmembers = CHECK_HDF5(H5Tget_nmembers(datatype));
      for(int i = 0 ; i < nmembers ; ++i) {
        int value;
        char* name;
        CHECK_HDF5(datatype, i, &value);
        name = CHECK_HDF5(datatype, i);
        assert (name);
        filenames[name] = value;
        filenamesmaxvalue = max(filenamesmaxvalue, value);
        free(name);
      }
      CHECK_HDF5(H5Tclose(datatype));
    }
  } H5_END_TRY;

  // compound type for map entry structure
  hid_t ivec_t = CHECK_HDF5(H5Tarray_create (H5T_NATIVE_INT, 1, dims, NULL));
  hid_t rvec_t = CHECK_HDF5(H5Tarray_create (H5T_NATIVE_DOUBLE, 1, dims, NULL));

  hid_t objectname_t = CHECK_HDF5(H5Tcopy(H5T_C_S1));
  CHECK_HDF5(H5Tset_size(objectname_t. H5T_VARIABLE));

  CHECK_HDF5(mapentry_t = H5Tcreate(H5T_COMPOUND, sizeof(mapentry)));
  CHECK_HDF5(H5Tinsert(mapentry_t, "map", H5OFFSET(mapentry, map),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "mglevel", H5OFFSET(mapentry, mglevel),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "reflevel", H5OFFSET(mapentry, reflevel),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "timestep", H5OFFSET(mapentry, timestep),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "timelevel", H5OFFSET(mapentry, timelevel),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "component", H5OFFSET(mapentry, component),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "rank", H5OFFSET(mapentry, rank),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "iorigin", H5OFFSET(mapentry, iorigin),
                       ivec_t));
  CHECK_HDF5(H5Tinsert(mapentry_t, "ioffset", H5OFFSET(mapentry, ioffset),
                       ivec_t));
  CHECK_HDF5(H5Tinsert(mapentry_t, "ioffsetdenom", H5OFFSET(mapentry,
                                                            ioffsetdenom),
                       ivec_t));
  CHECK_HDF5(H5Tinsert(mapentry_t, "shape", H5OFFSET(mapentry, shape),
                       ivec_t));
  CHECK_HDF5(H5Tinsert(mapentry_t, "origin", H5OFFSET(mapentry, origin),
                       rvec_t));
  CHECK_HDF5(H5Tinsert(mapentry_t, "delta", H5OFFSET(mapentry, delta),
                       rvec_t));
  CHECK_HDF5(H5Tinsert(mapentry_t, "time", H5OFFSET(mapentry, time),
                       H5T_NATIVE_DOUBLE));
  CHECK_HDF5(H5Tinsert(mapentry_t, "varnameindex", H5OFFSET(mapentry,
                                                            varnameindex),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "filenameindex", H5OFFSET(mapentry,
                                                             filenameindex),
                       H5T_NATIVE_INT));
  CHECK_HDF5(H5Tinsert(mapentry_t, "objectname", H5OFFSET(mapentry, objectname),
                       objectname_t));
  // HDF5 keeps the types open since something still refers to them
  CHECK_HDF5(H5Tclose(objectname_t));
  CHECK_HDF5(H5Tclose(rvec_t));
  CHECK_HDF5(H5Tclose(ivec_t));
}

 /*@@
   @routine    ProcessDataset
   @date       Fri 17 October 2008
   @author     Thomas Radke
   @desc
               The worker routine which is called by H5Giterate().
               It checks whether the current HDF5 object is a dataset matching
               the user's slab criteria.
   @enddesc

   @var        group
   @vdesc      HDF5 object to start the iteration
   @vtype      hid_t
   @vio        in
   @endvar
   @var        datasetname
   @vdesc      name of the object at the current iteration
   @vtype      const char *
   @vio        in
   @endvar
   @var        _file
   @vdesc      pointer to the descriptor of the currently opened file
   @vtype      void *
   @vio        in
   @endvar
@@*/
static herr_t ProcessDataset (hid_t group, const char *datasetname, void *_file)
{
  // we are interested in datasets only - skip anything else
  H5G_stat_t object_info;
  CHECK_HDF5 (H5Gget_objinfo (group, datasetname, 0, &object_info));
  if (object_info.type != H5G_DATASET) {
    return (0);
  }

  CHECK_HDF5 (from = H5Dopen (from, objectname));
  CHECK_HDF5 (dataspace = H5Dget_space (from));

  hsize_t shape[3] = {1,1,1};
  hid_t attr, attrtype;
  char *s;

  char vname[1024];
  int map;
  int mglevel;
  int reflevel;
  int timestep;
  int timelevel;
  int component;
  int rank;
  int iorigin[3];
  int ioffset[3];
  int ioffsetdenom[3];   
  double origin[3];
  double delta[3];
  double time;
  int nghosts[3];

  rank = H5Sget_simple_extent_ndims (dataspace);
  assert(rank <= 3);
  // Ian's original writer sets the dimensions to (1,1,1) so we have to try
  // and read this information out of an attribute first
  H5E_BEGIN_TRY {
    attr = H5Aopen_name (from, "h5shape");
  } H5E_END_TRY
  if (attr >= 0) {
    CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_HSIZE, shape));
    CHECK_HDF5 (H5Aclose (attr));
  } else {
    CHECK_HDF5 (H5Sget_simple_extent_dims (dataspace, &shape[0], NULL));
  }

  if ((s = strstr(objectname, "m="))) sscanf(s, "m=%d", &map); else map=0;
  sscanf(strstr(objectname, "it="), "it=%d", &timestep); 
  sscanf(strstr(objectname, "tl="), "tl=%d", &timelevel); 
  if ((s = strstr(objectname, "c="))) sscanf(s, "c=%d", &component); else component = 0;
  if ((s = strstr(objectname, "rl="))) sscanf(s, "rl=%d", &reflevel); else reflevel = 0;

  CHECK_HDF5 (attr = H5Aopen_name (from, "name"));
  CHECK_HDF5 (attrtype = H5Aget_type (attr));
  size_t length = H5Tget_size (attrtype);
  assert (length > 0);
  CHECK_HDF5 (H5Aread (attr, attrtype, &vname));
  vname[length] = '\0';
  CHECK_HDF5 (H5Tclose (attrtype));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "carpet_mglevel"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, &mglevel));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "iorigin"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, iorigin));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "ioffset"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, ioffset));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "ioffsetdenom"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, ioffsetdenom));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "origin"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_DOUBLE, origin));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "delta"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_DOUBLE, delta));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "time"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_DOUBLE, &time));
  CHECK_HDF5 (H5Aclose (attr));
  CHECK_HDF5 (attr = H5Aopen_name (from, "cctk_nghostzones"));
  CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, nghosts));
  CHECK_HDF5 (H5Aclose (attr));

  int varnamesvalue;
  if (varnames.count(vname)) {
    varnamesvalue = varnames[vname];
  } else {
    varnamesvalue = ++varnamesmaxvalue;
    varnames[vname] = varnamesvalue;
  }

  int filenamessvalue;
  if (filenames.count(vname)) {
    filenamesvalue = filenames[vname];
  } else {
    filenamesvalue = ++filenamesmaxvalue;
    filenames[vname] = filenamesvalue;
  }

  #define DOUBLE_TO_INT(d) ((int*)&(d))[0],((int*)&(d))[1]
  assert(2*sizeof(int) == sizeof(double));
  int len;
  int data[] = {
         0x4343, filenum,
         map, mglevel, reflevel, timestep, timelevel, component,
         rank, iorigin[0], iorigin[1], iorigin[2], ioffset[0], ioffset[1],
         ioffset[2], ioffsetdenom[0], ioffsetdenom[1], ioffsetdenom[2],
         (int)shape[0], (int)shape[1], (int)shape[2],
         DOUBLE_TO_INT(origin[0]), DOUBLE_TO_INT(origin[1]), DOUBLE_TO_INT(origin[2]),
         DOUBLE_TO_INT(delta[0]), DOUBLE_TO_INT(delta[1]), DOUBLE_TO_INT(delta[2]),
         DOUBLE_TO_INT(time),
         (int)nghosts[0], (int)nghosts[1], (int)nghosts[2],
         (int)strlen(vname), (int)strlen(objectname)
  };
  fwrite(data, sizeof(data), 1, outfile);
  len = (strlen(vname) + 1 + sizeof(int)-1) / sizeof(int);
  fwrite(vname, len, sizeof(int), outfile);
  len = (strlen(objectname) + 1 + sizeof(int)-1) / sizeof(int);
  fwrite(objectname, len, sizeof(int), outfile);

  CHECK_HDF5 (H5Dclose (from));
  CHECK_HDF5 (H5Sclose (dataspace));

  return 0;
}
