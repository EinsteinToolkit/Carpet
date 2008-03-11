#include <algorithm>
#include <set>
#include <string>

#include <cctk.h>
#include <cctk_Schedule.h>

#include <cctki_GHExtensions.h>
#include <cctki_Schedule.h>

#include "carpet.hh"

using namespace std;



#ifdef CACTUS_HAS_REQUIRES_CLAUSES

// Illegally copied from ScheduleInterface.c:

typedef enum {sched_none, sched_group, sched_function} iSchedType;
typedef enum {schedpoint_misc, schedpoint_analysis} iSchedPoint;

typedef struct
{
  /* Static data */
  char *description;

  /*char *thorn; MOVED TO FunctionData */
  char *implementation;

  iSchedType type;

  cFunctionData FunctionData;

  int n_mem_groups;
  int *mem_groups;
  int *timelevels;

  int n_comm_groups;
  int *comm_groups;

  /* Timer data */

  int timer_handle;

  /* Dynamic data */
  int *CommOnEntry;
  int *StorageOnEntry;

  int done_entry;
  int synchronised;

} t_attribute;



namespace Carpet {
  namespace Requirements {
    
    
    
    // Rules:
    //
    // 1. Everything that is required by a routine must be provided by
    //    another routine which is scheduled earlier.
    //
    // 2. Things can be provided only once, not multiple times.
    //    Except when they are also provided.
    
    
    
    int CheckEntry (void * attribute, void * data);
    int CheckExit  (void * attribute, void * data);
    int CheckWhile (int n_whiles, char ** whiles, void * attribute, void * data, int first);
    int CheckIf    (int n_ifs, char ** ifs, void * attribute, void * data);
    int CheckCall  (void * function, void * attribute, void * data);
    
    void CheckOneGroup (cGH const * cctkGH, char const * where);
    
    
    
    int
    CheckEntry (void * const attribute,
                void * const data)
    {
      if (not attribute) {
        // Nothing to check
        return 1;
      }
      
      // Convert argument types
      cFunctionData & function_data =
        (static_cast <t_attribute *> (attribute))->FunctionData;
      set <string> & active_provisions = * static_cast <set <string> *> (data);
      
      // Gather all required items
      set <string> requires;
      for (int n = 0; n < function_data.n_RequiresClauses; ++ n) {
        requires.insert (string (function_data.RequiresClauses[n]));
      }
      
      // Check whether all required items have already been provided
      set <string> required_but_not_provided;
      set_difference
        (requires.begin(), requires.end(),
         active_provisions.begin(), active_provisions.end(),
         insert_iterator <set <string> >
         (required_but_not_provided, required_but_not_provided.begin()));
      
      // Are there unmet requirements?
      if (not required_but_not_provided.empty()) {
        for (set<string>::const_iterator ri = required_but_not_provided.begin();
             ri != required_but_not_provided.end(); ++ ri)
        {
          string const req = * ri;
          CCTK_VParamWarn (CCTK_THORNSTRING,
                           "Requirement inconsistency:\n"
                           "   Group %s, function %s::%s requires \"%s\" which has not been provided",
                           function_data.where,
                           function_data.thorn, function_data.routine,
                           req.c_str());
        }
      }
      
      // Do traverse this schedule item
      return 1;
    }
    
    
    
    int
    CheckExit (void * const attribute,
               void * const data)
    {
      if (not attribute) {
        // Nothing to check
        return 1;
      }
      
      // Convert argument types
      cFunctionData & function_data =
        (static_cast <t_attribute *> (attribute))->FunctionData;
      set <string> & active_provisions = * static_cast <set <string> *> (data);
      
      // Gather all required and provided items
      set <string> requires;
      for (int n = 0; n < function_data.n_RequiresClauses; ++ n) {
        requires.insert (string (function_data.RequiresClauses[n]));
      }
      set <string> provides;
      for (int n = 0; n < function_data.n_ProvidesClauses; ++ n) {
        provides.insert (string (function_data.ProvidesClauses[n]));
      }
      
      // Check whether any of the providions have already been
      // provided.  (We disallow this as well, so that a routine
      // cannot overwrite what another routine has already set up.)
      set <string> provided_twice;
      set_intersection
        (provides.begin(), provides.end(),
         active_provisions.begin(), active_provisions.end(),
         insert_iterator <set <string> >
         (provided_twice, provided_twice.begin()));
      // But we do allow to provide things which are also required
      set <string> provided_too_often;
      set_difference
        (provided_twice.begin(), provided_twice.end(),
         requires.begin(), requires.end(),
         insert_iterator <set <string> >
         (provided_too_often, provided_too_often.begin()));
      
      // Are there things provided twice?
      if (not provided_too_often.empty()) {
        for (set<string>::const_iterator pi = provided_too_often.begin();
             pi != provided_too_often.end(); ++ pi)
        {
          string const prov = * pi;
          CCTK_VParamWarn (CCTK_THORNSTRING,
                           "Requirement inconsistency:\n"
                           "   Group %s, function %s::%s provides (and does not require) \"%s\" which has already been provided",
                           function_data.where,
                           function_data.thorn, function_data.routine,
                           prov.c_str());
        }
      }
      
      // Add the new provisions
      for (set<string>::const_iterator pi =
             provides.begin(); pi != provides.end(); ++ pi)
      {
        string const prov = * pi;
        active_provisions.insert (prov);
      }
      
      // ???
      return 1;
    }
    
    
    
    int
    CheckWhile (int const n_whiles, char ** const whiles,
                void * const attribute,
                void * const data,
                int const first)
    {
      // Execute item once
      return first;
    }
    
    
    
    int
    CheckIf (int const n_ifs, char ** const ifs,
             void * const attribute,
             void * const data)
    {
      // Execute item
      return 1;
    }
    
    
    
    int CheckCall (void * const function,
                   void * const attribute,
                   void * const data)
    {
      // Do nothing
      return 0;
    }
    
    
    
    // Check one schedule bin
    void
    CheckOneGroup (cGH const * const cctkGH,
                   char const * const where)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Checking requirements of schedule bin %s", where);
      
      // Set up initial provision (none at the moment)
      set <string> active_provisions;
      
      // Output initial provisions
      CCTK_VInfo (CCTK_THORNSTRING,
                  "   Initial provisions:");
      for (set<string>::const_iterator pi =
             active_provisions.begin(); pi != active_provisions.end(); ++ pi)
      {
        string const prov = * pi;
        CCTK_VInfo (CCTK_THORNSTRING,
                    "      %s", prov.c_str());
      }
      
      // Check the schedule bin
      CCTKi_DoScheduleTraverse (where,
                                CheckEntry, CheckExit,
                                CheckWhile, CheckIf,
                                CheckCall,
                                & active_provisions);
      
      // Output the final provisions
      CCTK_VInfo (CCTK_THORNSTRING,
                  "   Final provisions:");
      for (set<string>::const_iterator pi =
             active_provisions.begin(); pi != active_provisions.end(); ++ pi)
      {
        string const prov = * pi;
        CCTK_VInfo (CCTK_THORNSTRING,
                    "      %s", prov.c_str());
      }
    }
    
    
    
    // Check everything
    void
    CheckRequirements (cGH const * const cctkGH)
    {
      Checkpoint ("Checking schedule requirements");
      
      // Check some bins
      CheckOneGroup (cctkGH, "CCTK_WRAGH");
      CheckOneGroup (cctkGH, "CCTK_BASEGRID");
      
      CheckOneGroup (cctkGH, "CCTK_RECOVER_VARIABLES");
      CheckOneGroup (cctkGH, "CCTK_POST_RECOVER_VARIABLES");
      
      CheckOneGroup (cctkGH, "CCTK_PREREGRIDINITIAL");
      CheckOneGroup (cctkGH, "CCTK_POSTREGRIDINITIAL");
      CheckOneGroup (cctkGH, "CCTK_INITIAL");
      CheckOneGroup (cctkGH, "CCTK_POSTRESTRICTINITIAL");
      CheckOneGroup (cctkGH, "CCTK_POSTINITIAL");
      CheckOneGroup (cctkGH, "CCTK_CPINITIAL");
      
      CheckOneGroup (cctkGH, "CCTK_PREREGRID");
      CheckOneGroup (cctkGH, "CCTK_POSTREGRID");
      CheckOneGroup (cctkGH, "CCTK_PRESTEP");
      CheckOneGroup (cctkGH, "CCTK_EVOL");
      CheckOneGroup (cctkGH, "CCTK_POSTSTEP");
      CheckOneGroup (cctkGH, "CCTK_CHECKPOINT");
      CheckOneGroup (cctkGH, "CCTK_ANALYSIS");
      
      CheckOneGroup (cctkGH, "CCTK_TERMINATE");
    }
    
  };                            // namespace Carpet
};                              // namespace Requirements



#else // #ifndef CACTUS_HAS_REQUIRES_CLAUSES



namespace Carpet {
  namespace Requirements {
    // Check one schedule bin
    void
    CheckRequirements (cGH const * const cctkGH)
    {
      Checkpoint ("Skipping check of schedule requirements (no flesh support)");
      // do nothing
    }
  };                            // namespace Carpet
};                              // namespace Requirements



#endif // #ifdef CACTUS_HAS_REQUIRES_CLAUSES
