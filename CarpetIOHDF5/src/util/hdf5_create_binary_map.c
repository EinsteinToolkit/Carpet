 /*@@
   @file      hdf5_extract.c
   @date      Thu 19 Feb 2002
   @author    Thomas Radke
   @desc
              This utility program extracts objects from an HDF5 file and writes
              them into a new one.
   @enddesc
   @version   $Id: hdf5_extract.c 37 2009-09-29 14:38:14Z schnetter $
 @@*/

//#include "cctk.h"

#define H5_USE_16_API 1
#include <hdf5.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
//CCTK_FILEVERSION(CactusExternal_HDF5_util_hdf5_extract_c)


/*****************************************************************************/
/*                           macro definitions                               */
/*****************************************************************************/
/* macro to do an HDF5 call, check its return code, and print a warning
   in case of an error */
#define CHECK_ERROR(hdf5_call)                                                \
          do                                                                  \
          {                                                                   \
            int _error_code = hdf5_call;                                      \
                                                                              \
                                                                              \
            if (_error_code < 0)                                              \
            {                                                                 \
              fprintf (stderr, "WARNING: line %d: HDF5 call '%s' returned "   \
                               "error code %d\n",                             \
                                __LINE__, #hdf5_call, _error_code);           \
              nerrors++;                                                      \
            }                                                                 \
          } while (0)


/*****************************************************************************/
/*                           global variables                                */
/*****************************************************************************/
/* NOTE: although it isn't good programming practice
         we make these variables global for convenience
         since they are accessed from recursively or
         indirectly called routines which only get passed
         a single user-supplied argument */
static char *pathname = NULL;            /* pathname of the current object */
static int filenum = -1;                 /* Carpet file number of input file */
static unsigned int nerrors = 0;         /* global error counter */

/*****************************************************************************/
/*                           local function prototypes                       */
/*****************************************************************************/
static herr_t CopyObject (hid_t from, const char *objectname, void *_to);


 /*@@
   @routine    main
   @date       Thu 19 Feb 2002
   @author     Thomas Radke
   @desc
               Main routine of the HDF5 file extractor
   @enddesc

   @calls      CopyObject

   @var        argc
   @vdesc      number of command line arguments
   @vtype      int
   @vio        in
   @endvar
   @var        argv
   @vdesc      command line arguments
   @vtype      char *[]
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success, negative return values indicate an error
   @endreturndesc
@@*/
int main (int argc, char *argv[])
{
  char objectname[256];
  FILE *listfile, *outfile;
  hid_t infile;


  /* give some help if called with incorrect number of parameters */
  if (argc != 4)
  {
    fprintf (stderr, "Usage: %s <listfile> <infile> <outfile>\n", argv[0]);
    fprintf (stderr, "       where <listfile> can be created from the output of the h5ls command\n\n");
    fprintf (stderr, "   eg, h5ls alp.h5 > listfile\n"
                     "       # now edit listfile to select all datasets to be extracted\n"
                     "       # (don't forget to remove backslashes and the trailing 'Dataset {...}')\n"
                     "       %s extract_list alp.h5 alp_extract.h5\n\n", argv[0]);

    return (0);
  }

  /* open the list file */
  listfile = fopen (argv[1], "r");
  if (listfile == NULL)
  {
    fprintf (stderr, "ERROR: Cannot open list file '%s' !\n\n", argv[1]);
    return (-1);
  }

  H5E_BEGIN_TRY
  {
    /* open the input file */
    infile = H5Fopen (argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
    if (infile < 0)
    {
      fprintf (stderr, "ERROR: Cannot open HDF5 input file '%s' !\n\n",argv[2]);
      return (-1);
    }
  } H5E_END_TRY
  outfile = fopen(argv[3], "a");
  if (!outfile)
  {
      fprintf (stderr, "ERROR: Cannot open output file '%s' !\n\n",argv[3]);
      return (-1);
  }

  {
    const char marker[] = ".file_";
    const char *itstring = strstr(argv[2], marker);
    char *dot;
    assert(itstring);
    filenum = strtol(itstring+strlen(marker), &dot, 10);
    assert(strcmp(dot,".h5")==0);
  }

  /* do the copying by iterating over all objects */
  while (fgets (objectname, sizeof (objectname), listfile))
  {
    /* get rid of the trailing newline */
    objectname[strlen (objectname) - 1] = 0;
    pathname = "";

    CopyObject (infile, objectname, outfile);
  }

  /* finally, close all open files */
  CHECK_ERROR (H5Fclose (infile));
  //CHECK_ERROR (H5Fclose (outfile));
  fclose (listfile);

  /* report status */
  if (nerrors != 0)
  {
    fprintf (stderr, "\n\n   *** WARNING: %u errors occured during "
                     "data extraction. ***\n\n", nerrors);
  }

  return (0);
}


/*****************************************************************************/
/*                           local routines                                  */
/*****************************************************************************/
 /*@@
   @routine    CopyObject
   @date       Thu 19 Feb 2002
   @author     Thomas Radke
   @desc
               Copies an object with given name from the input file into the
               output file.
   @enddesc

   @calls      CopyAttribute

   @var        from
   @vdesc      input group handle
   @vtype      hid_t
   @vio        in
   @endvar
   @var        objectname
   @vdesc      name of the object to extract
   @vtype      const char *
   @vio        in
   @endvar
   @var        _to
   @vdesc      output group handle
   @vtype      void *
   @vio        in
   @endvar

   @returntype herr_t
   @returndesc
               0 - continue the iteration for following group objects
               1 - short-curcuit, no further iteration of this group
   @endreturndesc
@@*/
static herr_t CopyObject (hid_t from,
                          const char *objectname,
                          void *_outfile)
{
  hid_t dataspace;
  H5G_stat_t objectinfo;
  char *current_pathname;
  herr_t result;
  FILE *outfile = (FILE*)_outfile;


  /* check whether the requested object exists */
  H5E_BEGIN_TRY
  {
    result = H5Gget_objinfo (from, objectname, 0, &objectinfo);
  } H5E_END_TRY
  if (result < 0)
  {
    fprintf (stderr, "WARNING: object '%s' will not be extracted (object "
                     "doesn't exist)\n", objectname);
    nerrors++;
    return (0);
  }

  /* build the full pathname for the current to object to process */
  current_pathname = pathname;
  pathname = (char *) malloc (strlen (current_pathname) +
                              strlen (objectname) + 2);
  sprintf (pathname, "%s/%s", current_pathname, objectname);

  /* check the type of the current object */
  if (objectinfo.type == H5G_DATASET)
  {
    CHECK_ERROR (from = H5Dopen (from, objectname));
    CHECK_ERROR (dataspace = H5Dget_space (from));

    {
    hsize_t shape[3] = {1,1,1};
    hid_t attr, attrtype;
    char *s;

    char vname[1024];
    int map;
    int mglevel = 0;
    int reflevel;
    int timestep;
    int timelevel;
    int component;
    int rank;
    int iorigin[3];
    int ioffset[3] = {0,0,0};
    int ioffsetdenom[3] = {1,1,1};
    double origin[3] = {0.,0.,0.};
    double delta[3] = {0.,0.,0.};
    double time;
    int nghosts[3] = {0,0,0};
    int bbox[6] = {0,0,0,0,0,0};

    rank = H5Sget_simple_extent_ndims (dataspace);
    assert(rank <= 3);
    // Ian's original writer sets the dimensions to (1,1,1) so we have to try
    // and read this information out of an attribute first
    H5E_BEGIN_TRY {
      attr = H5Aopen_name (from, "h5shape");
    } H5E_END_TRY
    if (attr >= 0) {
      CHECK_ERROR (H5Aread (attr, H5T_NATIVE_HSIZE, shape));
      CHECK_ERROR (H5Aclose (attr));
    } else {
      CHECK_ERROR (H5Sget_simple_extent_dims (dataspace, &shape[0], NULL));
    }

    if ((s = strstr(objectname, "m="))) sscanf(s, "m=%d", &map); else map=0;
    sscanf(strstr(objectname, "it="), "it=%d", &timestep); 
    sscanf(strstr(objectname, "tl="), "tl=%d", &timelevel); 
    if ((s = strstr(objectname, "c="))) sscanf(s, "c=%d", &component); else component = 0;
    if ((s = strstr(objectname, "rl="))) sscanf(s, "rl=%d", &reflevel); else reflevel = 0;

    CHECK_ERROR (attr = H5Aopen_name (from, "name"));
    CHECK_ERROR (attrtype = H5Aget_type (attr));
    size_t length = H5Tget_size (attrtype);
    assert (length > 0);
    CHECK_ERROR (H5Aread (attr, attrtype, &vname));
    vname[length] = '\0';
    CHECK_ERROR (H5Tclose (attrtype));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "carpet_mglevel"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &mglevel));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "iorigin"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, iorigin));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "ioffset"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, ioffset));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "ioffsetdenom"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, ioffsetdenom));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "origin"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, origin));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "delta"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, delta));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "time"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, &time));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "cctk_nghostzones"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, nghosts));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (from, "cctk_bbox"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, bbox));
    CHECK_ERROR (H5Aclose (attr));

    #define DOUBLE_TO_INT(d) ((int*)&(d))[0],((int*)&(d))[1]
    assert(2*sizeof(int) == sizeof(double));
    int len;
    int data[] = {
           0x4444, filenum,
           map, mglevel, reflevel, timestep, timelevel, component,
           rank, iorigin[0], iorigin[1], iorigin[2], ioffset[0], ioffset[1],
           ioffset[2], ioffsetdenom[0], ioffsetdenom[1], ioffsetdenom[2],
           (int)shape[0], (int)shape[1], (int)shape[2],
           DOUBLE_TO_INT(origin[0]), DOUBLE_TO_INT(origin[1]), DOUBLE_TO_INT(origin[2]),
           DOUBLE_TO_INT(delta[0]), DOUBLE_TO_INT(delta[1]), DOUBLE_TO_INT(delta[2]),
           DOUBLE_TO_INT(time),
           (int)nghosts[0], (int)nghosts[1], (int)nghosts[2],
           bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
           (int)strlen(vname), (int)strlen(objectname)
    };
    fwrite(data, sizeof(data), 1, outfile);
    len = (strlen(vname) + 1 + sizeof(int)-1) / sizeof(int);
    fwrite(vname, len, sizeof(int), outfile);
    len = (strlen(objectname) + 1 + sizeof(int)-1) / sizeof(int);
    fwrite(objectname, len, sizeof(int), outfile);

    //printf("filenum=%d,"
    //       "vname=%s,map=%d,mglevel=%d,reflevel=%d,timestep=%d,"
    //       "timelevel=%d,component=%d,rank=%d,iorigin=(%d,%d,%d),"
    //       "ioffset=(%d,%d,%d),ioffsetdenom=(%d,%d,%d),shape=(%d,%d,%d),"
    //       "patchname=%s,\n",
    //       filenum,
    //       vname, map, mglevel, reflevel, timestep, timelevel, component,
    //       rank, iorigin[0], iorigin[1], iorigin[2], ioffset[0], ioffset[1],
    //       ioffset[2], ioffsetdenom[0], ioffsetdenom[1], ioffsetdenom[2],
    //       (int)shape[0], (int)shape[1], (int)shape[2], objectname);
    }

    CHECK_ERROR (H5Dclose (from));
    CHECK_ERROR (H5Sclose (dataspace));
  }

  /* reset the pathname */
  free (pathname);
  pathname = current_pathname;

  return (0);
}
