// This was adopted from Thomas Radke's IOHDF5 thorn.
// Thanks, Thomas!

#ifndef _CARPETIOHDF5_IOHDF5GH_H_
#define _CARPETIOHDF5_IOHDF5GH_H_ 1

#include "StoreNamedData.h"

/* I took basically everything in this file from Thomas' IOHDF5; much of
   below is still unused.. */

/* CARPET IOHDF5 GH extension structure */
typedef struct
{
  /* default number of times to output */
  int out_every_default;

  /* number of times to output for each variable */
  CCTK_INT *out_every;

  /* the last iteration output for each variable */
  int *out_last;

  /* list of variables to output */
  char *out_vars;

  /* I/O request description list (for all variables) */
  ioRequest **requests;

  /* directory in which to output */
  char *out_dir;

  /* filename database for opened files */
  pNamedData *open_output_files;

  /* timer array for checkpointing/recovery */
  // int timers[IOHDF5_NUM_TIMERS];

  /* flag to indicate request for timer output */
  // int print_timing_info;

  /* ring buffer for list of successfully created cp files */
  int    cp_filename_index;
  char **cp_filename_list;

  /* iteration number of the last checkpoint */
  int last_checkpoint_iteration;

  /* hdf5 datatype for stupid complex variables; to be set at run time */
  hid_t HDF5_COMPLEX, HDF5_COMPLEX8, HDF5_COMPLEX16, HDF5_COMPLEX32;
  

} CarpetIOHDF5GH;


#endif  
