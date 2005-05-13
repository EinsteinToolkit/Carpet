#include <cassert>
#include <cstring>

#include "cctk.h"

#include "extending.hh"



namespace CarpetIOF5 {
  
  char const * const extending_t::
  extension_name = "CarpetIOF5";
  
  void extending_t::
  create ()
  {
    int const handle = CCTK_RegisterGHExtension (extension_name);
    assert (handle >= 0);
    int const ierr = CCTK_RegisterGHExtensionSetupGH (handle, setup);
    assert (! ierr);
  }
  
  void * extending_t::
  setup (tFleshConfig * const fleshconfig,
         int const convlevel,
         cGH * const cctkGH)
  {
    assert (fleshconfig != 0);
    assert (cctkGH != 0);
    return new extension_t;
  }
  
  extending_t::
  extending_t (cGH const * cctkGH)
  {
    assert (cctkGH);
    void * const ext = CCTK_GHExtension (cctkGH, extension_name);
    assert (ext != 0);
    m_extension = static_cast<extension_t *> (ext);
  }
  
  bool extending_t::
  get_did_truncate (char const * const name)
    const
  {
    assert (name != 0);
    return (m_extension->did_truncate.find (name)
            != m_extension->did_truncate.end());
  }
  
  void extending_t::
  set_did_truncate (char const * const name)
  {
    assert (name != 0);
    m_extension->did_truncate.insert (strdup (name));
  }
  
  int extending_t::
  get_last_output_iteration (int const ml, int const rl, int const vi)
    const
  {
    resize_last_output (ml, rl, vi, m_extension->last_output_iteration);
    return m_extension->last_output_iteration.at(ml).at(rl).at(vi);
  }
  
  void extending_t::
  set_last_output_iteration (int const ml, int const rl, int const vi,
                             int const iteration)
  {
    resize_last_output (ml, rl, vi, m_extension->last_output_iteration);
    m_extension->last_output_iteration.at(ml).at(rl).at(vi) = iteration;
  }
  
  CCTK_REAL extending_t::
  get_last_output_time (int const ml, int const rl, int const vi)
    const
  {
    resize_last_output (ml, rl, vi, m_extension->last_output_time);
    return m_extension->last_output_time.at(ml).at(rl).at(vi);
  }
  
  void extending_t::
  set_last_output_time (int const ml, int const rl, int const vi,
                        CCTK_REAL const time)
  {
    resize_last_output (ml, rl, vi, m_extension->last_output_time);
    m_extension->last_output_time.at(ml).at(rl).at(vi) = time;
  }
  
  template<typename T>
  void extending_t::
  resize_last_output (int ml, int rl, int vi,
                      vector<vector<vector<T> > > & array)
  {
    assert (ml >= 0);
    if (ml >= array.size())
    {
      array.resize (ml+1);
    }
    assert (rl >= 0);
    if (rl >= array.at(ml).size())
    {
      array.at(ml).resize (rl+1);
    }
    if (vi >= array.at(ml).at(rl).size())
    {
      array.at(ml).at(rl).resize (vi+1, -1);
    }
  }
  
  template
  void extending_t::
  resize_last_output (int ml, int rl, int vi,
                      vector<vector<vector<int> > > & array);
  template
  void extending_t::
  resize_last_output (int ml, int rl, int vi,
                      vector<vector<vector<CCTK_REAL> > > & array);
  
} // namespace CarpetIOF5
