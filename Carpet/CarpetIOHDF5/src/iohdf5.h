/* $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetIOHDF5/src/iohdf5.h,v 1.1 2004/03/03 09:44:26 schnetter Exp $ */

#ifndef CARPETIOHDF5_H
#define CARPETIOHDF5_H

#include "cctk_Arguments.h"

#ifdef __cplusplus
namespace CarpetIOHDF5 {
  extern "C" {
#endif
    
    /* Scheduled functions */
    int CarpetIOHDF5Startup (void);
    int CarpetIOHDF5ReadData (CCTK_ARGUMENTS);
    
#ifdef __cplusplus
  } /* extern "C" */
} /* namespace CarpetIOHDF5 */
#endif

#endif /* !defined(CARPETIOHDF5_H) */
