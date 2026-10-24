#ifndef CARPET_HERMITE_FILL_HH
#define CARPET_HERMITE_FILL_HH

#include <cctk.h>

#include <vector>

namespace CarpetLib {
class ggf;
}

namespace Carpet {

// Borrows an unused time level of each evolved variable, fills it with the
// cubic Hermite value of the parent step, and puts the original storage back
// on destruction.
class HermiteFill {
public:
  HermiteFill(cGH const *cctkGH, std::vector<int> const &groups,
              CCTK_REAL time);
  ~HermiteFill();
  bool source_timelevel(CarpetLib::ggf const *evolved, int &src_tl) const;

private:
  HermiteFill(HermiteFill const &);
  HermiteFill &operator=(HermiteFill const &);
  struct State;
  State *state_;
};

} // namespace Carpet

#endif
