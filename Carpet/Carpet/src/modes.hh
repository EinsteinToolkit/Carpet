#ifndef MODES_HH
#define MODES_HH

#include <cctk.h>
  


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
  bool is_meta_mode () CCTK_ATTRIBUTE_PURE;
  bool is_global_mode () CCTK_ATTRIBUTE_PURE;
  bool is_level_mode () CCTK_ATTRIBUTE_PURE;
  bool is_singlemap_mode () CCTK_ATTRIBUTE_PURE;
  bool is_local_mode () CCTK_ATTRIBUTE_PURE;
  
  
  
  // Mode setting
  
  void enter_global_mode (cGH * cctkGH, int ml);
  void leave_global_mode (cGH * cctkGH);
  
  void enter_level_mode (cGH * cctkGH, int rl);
  void leave_level_mode (cGH * cctkGH);
  
  void enter_singlemap_mode (cGH * cctkGH, int m, int grouptype);
  void leave_singlemap_mode (cGH * cctkGH);
  
  void enter_local_mode (cGH * cctkGH, int c, int lc, int grouptype);
  void leave_local_mode (cGH * cctkGH);
  
  
  
  // Mode iterators
  
  class mglevel_iterator {
    cGH * cctkGH;
    int ml;
  public:
    mglevel_iterator (cGH const * cctkGH);
    ~mglevel_iterator ();
    bool done () const;
    void step ();
  };
  
  class reverse_mglevel_iterator {
    cGH * cctkGH;
    int ml;
  public:
    reverse_mglevel_iterator (cGH const * cctkGH);
    ~reverse_mglevel_iterator ();
    bool done () const;
    void step ();
  };
  
  class reflevel_iterator {
    cGH * cctkGH;
    int rl;
  public:
    reflevel_iterator (cGH const * cctkGH);
    ~reflevel_iterator ();
    bool done () const;
    void step ();
  };
  
  class reverse_reflevel_iterator {
    cGH * cctkGH;
    int rl;
  public:
    reverse_reflevel_iterator (cGH const * cctkGH);
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
    cGH * cctkGH;
    int grouptype;
    int m;
  public:
    map_iterator (cGH const * cctkGH, int grouptype);
    ~map_iterator ();
    bool done () const;
    void step ();
  };
    
  class local_map_iterator {
    cGH * cctkGH;
    int grouptype;
    int m;
  public:
    local_map_iterator (cGH const * cctkGH, int grouptype);
    ~local_map_iterator ();
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
    cGH * cctkGH;
    int grouptype;
    int c;
  public:
    component_iterator (cGH const * cctkGH, int grouptype);
    ~component_iterator ();
    bool done () const;
    void step ();
  };
    
  class local_component_iterator {
    cGH * cctkGH;
    int grouptype;
    int lc;
  public:
    local_component_iterator (cGH const * cctkGH, int grouptype);
    ~local_component_iterator ();
    bool done () const;
    void step ();
  };
  
  
  
  // Compatibility defines for the mode iterators
  
#define BEGIN_MGLEVEL_LOOP(cctkGH)                      \
  do {                                                  \
    bool mglevel_loop_ = true;                          \
    for (Carpet::mglevel_iterator mg_iter_(cctkGH);     \
         not mg_iter_.done();                           \
         mg_iter_.step()) {
#define END_MGLEVEL_LOOP                        \
    }                                           \
    assert (mglevel_loop_);                     \
    mglevel_loop_ = false;                      \
  } while (false)
  
#define BEGIN_REVERSE_MGLEVEL_LOOP(cctkGH)                      \
  do {                                                          \
    bool reverse_mglevel_loop_ = true;                          \
    for (Carpet::reverse_mglevel_iterator mg_iter_(cctkGH);     \
         not mg_iter_.done();                                   \
         mg_iter_.step()) {
#define END_REVERSE_MGLEVEL_LOOP                \
    }                                           \
    assert (reverse_mglevel_loop_);             \
    reverse_mglevel_loop_ = false;              \
  } while (false)

#define BEGIN_REFLEVEL_LOOP(cctkGH)                     \
  do {                                                  \
    bool reflevel_loop_ = true;                         \
    for (Carpet::reflevel_iterator ref_iter_(cctkGH);   \
         not ref_iter_.done();                          \
         ref_iter_.step()) {
#define END_REFLEVEL_LOOP                       \
    }                                           \
    assert (reflevel_loop_);                    \
    reflevel_loop_ = false;                     \
  } while (false)

#define BEGIN_REVERSE_REFLEVEL_LOOP(cctkGH)                     \
  do {                                                          \
    bool reverse_reflevel_loop_ = true;                         \
    for (Carpet::reverse_reflevel_iterator ref_iter_(cctkGH);   \
         not ref_iter_.done();                                  \
         ref_iter_.step()) {
#define END_REVERSE_REFLEVEL_LOOP               \
    }                                           \
    assert (reverse_reflevel_loop_);            \
    reverse_reflevel_loop_ = false;             \
  } while (false)

#define BEGIN_MAP_LOOP(cctkGH, grouptype)                       \
  do {                                                          \
    bool map_loop_ = true;                                      \
    for (Carpet::map_iterator map_iter_(cctkGH, grouptype);     \
         not map_iter_.done();                                  \
         map_iter_.step()) {
#define END_MAP_LOOP                            \
    }                                           \
    assert (map_loop_);                         \
    map_loop_ = false;                          \
  } while (false)

#define BEGIN_LOCAL_MAP_LOOP(cctkGH, grouptype)                         \
  do {                                                                  \
    bool local_map_loop_ = true;                                        \
    for (Carpet::local_map_iterator local_map_iter_(cctkGH, grouptype); \
         not local_map_iter_.done();                                    \
         local_map_iter_.step()) {
#define END_LOCAL_MAP_LOOP                      \
    }                                           \
    assert (local_map_loop_);                   \
    local_map_loop_ = false;                    \
  } while (false)

#define BEGIN_COMPONENT_LOOP(cctkGH, grouptype)                         \
  do {                                                                  \
    bool component_loop_ = true;                                        \
    for (Carpet::component_iterator comp_iter_(cctkGH, grouptype);      \
         not comp_iter_.done();                                         \
         comp_iter_.step()) {
#define END_COMPONENT_LOOP                      \
    }                                           \
    assert (component_loop_);                   \
    component_loop_ = false;                    \
  } while (false)

#define BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, grouptype)                   \
  do {                                                                  \
    bool local_component_loop_ = true;                                  \
    for (Carpet::local_component_iterator comp_iter_(cctkGH, grouptype); \
         not comp_iter_.done();                                         \
         comp_iter_.step()) {
#define END_LOCAL_COMPONENT_LOOP                \
    }                                           \
    assert (local_component_loop_);             \
    local_component_loop_ = false;              \
  } while (false)

#define BEGIN_TIMELEVEL_LOOP(cctkGH)                                    \
  do {                                                                  \
    bool timelevel_loop_ = true;                                        \
    assert (do_allow_past_timelevels);                                  \
    do_allow_past_timelevels = false;                                   \
    assert (timelevel == 0);                                            \
    for (timelevel = 0; timelevel < timelevels; ++ timelevel) {         \
      cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);  \
      {
#define END_TIMELEVEL_LOOP                                              \
      }                                                                 \
    }                                                                   \
    assert (timelevel_loop_);                                           \
    timelevel_loop_ = false;                                            \
    timelevel = 0;                                                      \
    cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);    \
    do_allow_past_timelevels = true;                                    \
  } while (false)
  
  
  
  // Mode escapes
  
  class singlemap_escape {
    cGH * cctkGH;
    int c;
    int lc;
  public:
    singlemap_escape (cGH const * cctkGH);
    ~singlemap_escape ();
  };
  
  class level_escape {
    cGH * cctkGH;
    int grouptype;
    int m;
    int c;
    int lc;
  public:
    level_escape (cGH const * cctkGH);
    ~level_escape ();
  };
  
  class global_escape {
    cGH * cctkGH;
    int rl;
    int grouptype;
    int m;
    int c;
    int lc;
  public:
    global_escape (cGH const * cctkGH);
    ~global_escape ();
  };
  
  class meta_escape {
    cGH * cctkGH;
    int ml;
    int rl;
    int grouptype;
    int m;
    int c;
    int lc;
  public:
    meta_escape (cGH const * cctkGH);
    ~meta_escape ();
  };
  
  
  
  // Compatibility defines for the mode escapes
  
#define BEGIN_SINGLEMAP_MODE(cctkGH)            \
  do {                                          \
    bool singlemap_mode_ = true;                \
    Carpet::singlemap_escape esc_(cctkGH);      \
    {
#define END_SINGLEMAP_MODE                      \
    }                                           \
    assert (singlemap_mode_);                   \
    singlemap_mode_ = false;                    \
  } while (false)
  
#define BEGIN_LEVEL_MODE(cctkGH)                \
  do {                                          \
    bool level_mode_ = true;                    \
    Carpet::level_escape esc_(cctkGH);          \
    {
#define END_LEVEL_MODE                          \
    }                                           \
    assert (level_mode_);                       \
    level_mode_ = false;                        \
  } while (false)
  
#define BEGIN_GLOBAL_MODE(cctkGH)               \
  do {                                          \
    bool global_mode_ = true;                   \
    Carpet::global_escape esc_(cctkGH);         \
    {
#define END_GLOBAL_MODE                         \
    }                                           \
    assert (global_mode_);                      \
    global_mode_ = false;                       \
  } while (false)
  
#define BEGIN_META_MODE(cctkGH)                 \
  do {                                          \
    bool meta_mode_ = true;                     \
    Carpet::meta_escape esc_(cctkGH);           \
    {
#define END_META_MODE                           \
    }                                           \
    assert (meta_mode_);                        \
    meta_mode_ = false;                         \
  } while (false)
  
  
  
  // Mode setters
  
  class mglevel_setter {
    cGH * cctkGH;
  public:
    mglevel_setter (cGH const * cctkGH, int ml);
    ~mglevel_setter ();
  };
  
  class reflevel_setter {
    cGH * cctkGH;
  public:
    reflevel_setter (cGH const * cctkGH, int rl);
    ~reflevel_setter ();
  };
  
  class map_setter {
    cGH * cctkGH;
  public:
    map_setter (cGH const * cctkGH, int m, int grouptype);
    ~map_setter ();
  };
  
  class component_setter {
    cGH * cctkGH;
  public:
    component_setter (cGH const * cctkGH, int c, int grouptype);
    ~component_setter ();
  };
  
  
  
  // Compatibility defines for the mode setters
  
#define ENTER_GLOBAL_MODE(cctkGH, ml)                   \
  do {                                                  \
    Carpet::mglevel_setter mg_setter_(cctkGH, ml);      \
    {
#define LEAVE_GLOBAL_MODE                       \
    }                                           \
  } while (false)
  
#define ENTER_LEVEL_MODE(cctkGH, rl)                    \
  do {                                                  \
    Carpet::reflevel_setter ref_setter_(cctkGH, rl);    \
    {
#define LEAVE_LEVEL_MODE                        \
    }                                           \
  } while (false)
  
#define ENTER_SINGLEMAP_MODE(cctkGH, m, grouptype)      \
  do {                                                  \
    Carpet::map_setter m_setter_(cctkGH, m, grouptype); \
    {
#define LEAVE_SINGLEMAP_MODE                    \
    }                                           \
  } while (false)
  
#define ENTER_LOCAL_MODE(cctkGH, c, grouptype)                  \
  do {                                                          \
    Carpet::component_setter c_setter_(cctkGH, c, grouptype);   \
    {
#define LEAVE_LOCAL_MODE                        \
    }                                           \
  } while (false)
  
  
  
} // namespace Carpet

#endif // #ifndef MODES_HH
