// $Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/modes.hh,v 1.2 2004/06/14 09:42:53 schnetter Exp $

#ifndef MODES_HH
#define MODES_HH

#include "cctk.h"
  


namespace Carpet {
  
  using namespace std;
  
  
  
  //
  // These are the modes:
  //
  // meta mode:
  // global mode:    mglevel ("convergence level") is defined
  // level mode:     reflevel is defined
  // singlemap mode: map ("map index") is defined
  // local mode:     component ("patch index")
  //
  // maybe missing:
  // convtest mode:  [rl, map, c]
  //
  
  
  
  // Mode indicators
  bool is_meta_mode ();
  bool is_global_mode ();
  bool is_level_mode ();
  bool is_singlemap_mode ();
  bool is_local_mode ();
  
  
  
  // Mode setting
  
  void enter_global_mode (cGH * const cgh, int const ml);
  void leave_global_mode (cGH * const cgh);
  
  void enter_level_mode (cGH * const cgh, int const rl);
  void leave_level_mode (cGH * const cgh);
  
  void enter_singlemap_mode (cGH * const cgh, int const m);
  void leave_singlemap_mode (cGH * const cgh);
  
  void enter_local_mode (cGH * const cgh, int const c);
  void leave_local_mode (cGH * const cgh);
  
  
  
  // Mode iterators
  
  class mglevel_iterator {
    cGH * cgh;
    int ml;
  public:
    mglevel_iterator (cGH const * const cgh);
    ~mglevel_iterator ();
    bool done () const;
    void step ();
  };
  
  class reverse_mglevel_iterator {
    cGH * cgh;
    int ml;
  public:
    reverse_mglevel_iterator (cGH const * const cgh);
    ~reverse_mglevel_iterator ();
    bool done () const;
    void step ();
  };
  
  class reflevel_iterator {
    cGH * cgh;
    int rl;
  public:
    reflevel_iterator (cGH const * const cgh);
    ~reflevel_iterator ();
    bool done () const;
    void step ();
  };
  
  class reverse_reflevel_iterator {
    cGH * cgh;
    int rl;
  public:
    reverse_reflevel_iterator (cGH const * const cgh);
    ~reverse_reflevel_iterator ();
    bool done () const;
    void step ();
  };
  
  // Loop over all maps.  If grouptype is CCTK_GF, then loop over grid
  // function maps.  If grouptype is CCTK_ARRAY (or CCTK_SCALAR), then
  // loop over grid array (or grid scalar) maps.  In the latter case,
  // map denotes the current grid array map, i.e. it cannot be used to
  // access grid functions.
  class map_iterator {
    cGH * cgh;
    int grouptype;
    int m;
  public:
    map_iterator (cGH const * const cgh, int const grouptype);
    ~map_iterator ();
    bool done () const;
    void step ();
  };
    
  // Loop over all components.  If grouptype is CCTK_GF, then loop
  // over grid function components.  If grouptype is CCTK_ARRAY (or
  // CCTK_SCALAR), then loop over grid array (or grid scalar)
  // components.  In the latter case, component denotes the current
  // grid array component, i.e. it cannot be used to index grid
  // functions.
  class component_iterator {
    cGH * cgh;
    int grouptype;
    int c;
  public:
    component_iterator (cGH const * const cgh, int const grouptype);
    ~component_iterator ();
    bool done () const;
    void step ();
  };
    
  class local_component_iterator {
    cGH * cgh;
    int grouptype;
    int c;
  public:
    local_component_iterator (cGH const * const cgh, int const grouptype);
    ~local_component_iterator ();
    bool done () const;
    void step ();
  };
    
    
  
  // Compatibility defines for the mode iterators
  
#define BEGIN_MGLEVEL_LOOP(cgh)                 \
  do {                                          \
    bool mglevel_loop_ = true;                  \
    for (mglevel_iterator mg_iter_(cgh);        \
         !mg_iter_.done();                      \
         mg_iter_.step()) {
#define END_MGLEVEL_LOOP                        \
    }                                           \
    assert (mglevel_loop_);                     \
    mglevel_loop_ = false;                      \
  } while (false)
  
#define BEGIN_REVERSE_MGLEVEL_LOOP(cgh)                 \
  do {                                                  \
    bool reverse_mglevel_loop_ = true;                  \
    for (reverse_mglevel_iterator mg_iter_(cgh);        \
         !mg_iter_.done();                              \
         mg_iter_.step()) {
#define END_REVERSE_MGLEVEL_LOOP                \
    }                                           \
    assert (reverse_mglevel_loop_);             \
    reverse_mglevel_loop_ = false;              \
  } while (false)

#define BEGIN_REFLEVEL_LOOP(cgh)                \
  do {                                          \
    bool reflevel_loop_ = true;                 \
    for (reflevel_iterator ref_iter_(cgh);      \
         !ref_iter_.done();                     \
         ref_iter_.step()) {
#define END_REFLEVEL_LOOP                       \
    }                                           \
    assert (reflevel_loop_);                    \
    reflevel_loop_ = false;                     \
  } while (false)

#define BEGIN_REVERSE_REFLEVEL_LOOP(cgh)                \
  do {                                                  \
    bool reverse_reflevel_loop_ = true;                 \
    for (reverse_reflevel_iterator ref_iter_(cgh);      \
         !ref_iter_.done();                             \
         ref_iter_.step()) {
#define END_REVERSE_REFLEVEL_LOOP               \
    }                                           \
    assert (reverse_reflevel_loop_);            \
    reverse_reflevel_loop_ = false;             \
  } while (false)

#define BEGIN_MAP_LOOP(cgh, grouptype)                  \
  do {                                                  \
    bool map_loop_ = true;                              \
    for (map_iterator map_iter_(cgh, grouptype);        \
         !map_iter_.done();                             \
         map_iter_.step()) {
#define END_MAP_LOOP                            \
    }                                           \
    assert (map_loop_);                         \
    map_loop_ = false;                          \
  } while (false)

#define BEGIN_COMPONENT_LOOP(cgh, grouptype)            \
  do {                                                  \
    bool component_loop_ = true;                        \
    for (component_iterator comp_iter_(cgh, grouptype); \
         !comp_iter_.done();                            \
         comp_iter_.step()) {
#define END_COMPONENT_LOOP                      \
    }                                           \
    assert (component_loop_);                   \
    component_loop_ = false;                    \
  } while (false)

#define BEGIN_LOCAL_COMPONENT_LOOP(cgh, grouptype)              \
  do {                                                          \
    bool local_component_loop_ = true;                          \
    for (local_component_iterator comp_iter_(cgh, grouptype);   \
         !comp_iter_.done();                                    \
         comp_iter_.step()) {
#define END_LOCAL_COMPONENT_LOOP                \
    }                                           \
    assert (local_component_loop_);             \
    local_component_loop_ = false;              \
  } while (false)
  
  
  
  // Mode escapes
  
  class singlemap_escape {
    cGH * cgh;
    int c;
  public:
    singlemap_escape (cGH const * const cgh);
    ~singlemap_escape ();
  };
  
  class level_escape {
    cGH * cgh;
    int m;
    int c;
  public:
    level_escape (cGH const * const cgh);
    ~level_escape ();
  };
  
  class global_escape {
    cGH * cgh;
    int rl;
    int m;
    int c;
  public:
    global_escape (cGH const * const cgh);
    ~global_escape ();
  };
  
  class meta_escape {
    cGH * cgh;
    int ml;
    int rl;
    int m;
    int c;
  public:
    meta_escape (cGH const * const cgh);
    ~meta_escape ();
  };
  
  
  
  // Compatibility defines for the mode escapes
  
#define BEGIN_SINGLEMAP_MODE(cgh)               \
  do {                                          \
    bool singlemap_mode_ = true;                \
    singlemap_escape esc_(cgh);                 \
    {
#define END_SINGLEMAP_MODE                      \
    }                                           \
    assert (singlemap_mode_);                   \
    singlemap_mode_ = false;                    \
  } while (false)
  
#define BEGIN_LEVEL_MODE(cgh)                   \
  do {                                          \
    bool level_mode_ = true;                    \
    level_escape esc_(cgh);                     \
    {
#define END_LEVEL_MODE                          \
    }                                           \
    assert (level_mode_);                       \
    level_mode_ = false;                        \
  } while (false)
  
#define BEGIN_GLOBAL_MODE(cgh)                  \
  do {                                          \
    bool global_mode_ = true;                   \
    global_escape esc_(cgh);                    \
    {
#define END_GLOBAL_MODE                         \
    }                                           \
    assert (global_mode_);                      \
    global_mode_ = false;                       \
  } while (false)
  
#define BEGIN_META_MODE(cgh)                    \
  do {                                          \
    bool meta_mode_ = true;                     \
    meta_escape esc_(cgh);                      \
    {
#define END_META_MODE                           \
    }                                           \
    assert (meta_mode_);                        \
    meta_mode_ = false;                         \
  } while (false)
  
} // namespace Carpet

#endif // !defined(MODES_HH)
