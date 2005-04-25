#include <cassert>

#include "cctk.h"

#include "extending.hh"



namespace CarpetIOF5 {
  
  char const * const extending::
  extension_name = "CarpetIOF5";
  
  void extending::
  create (cGH * cctkGH)
  {
    assert (cctkGH);
    int const handle = CCTK_RegisterGHExtension (extension_name);
    assert (handle >= 0);
    int const ierr = CCTK_RegisterGHExtensionSetupGH (handle, setup);
    assert (! ierr);
  }
  
  void * extending::
  setup (tFleshConfig * const config,
         int const convlevel,
         cGH * const cctkGH)
  {
    assert (config);
    assert (cctkGH);
    return new extension_t;
  }
  
  extending::
  extending (cGH const * cctkGH)
  {
    assert (cctkGH);
    void * const ext = CCTK_GHExtension (cctkGH, extension_name);
    assert (ext != 0);
    m_extension = static_cast<extension_t *> (ext);
  }
  
  set<char const *> extending::
  get_did_truncate()
    const
  {
    return m_extension->did_truncate;
  }

} // namespace CarpetIOF5
