#ifndef CASM_monte2_OccEventProposal
#define CASM_monte2_OccEventProposal

#include <vector>
class MTRand;

namespace CASM {
namespace Monte2 {

struct OccEvent;
class OccLocation;
class OccSwap;

/// Choose a swap type from a list of allowed canonical swap types
OccSwap const &choose_canonical_swap(OccLocation const &occ_location,
                                     std::vector<OccSwap> const &canonical_swap,
                                     MTRand &mtrand);

/// Propose canonical OccEvent, given chosen OccSwap
OccEvent &propose_canonical_event(OccEvent &e, OccLocation const &occ_location,
                                  OccSwap const &swap, MTRand &mtrand);

/// Propose canonical OccEvent
OccEvent &propose_canonical_event(OccEvent &e, OccLocation const &occ_location,
                                  std::vector<OccSwap> const &canonical_swap,
                                  MTRand &mtrand);

/// Choose a swap type from a list of allowed grand canonical swap types
OccSwap const &choose_grand_canonical_swap(
    OccLocation const &occ_location,
    std::vector<OccSwap> const &grand_canonical_swap, MTRand &mtrand);

/// Propose grand canonical OccEvent of particular swap type
OccEvent &propose_grand_canonical_event(OccEvent &e,
                                        OccLocation const &occ_location,
                                        OccSwap const &swap, MTRand &mtrand);

/// Propose grand canonical OccEvent
OccEvent &propose_grand_canonical_event(
    OccEvent &e, OccLocation const &occ_location,
    std::vector<OccSwap> const &grand_canonical_swap, MTRand &mtrand);

}  // namespace Monte2
}  // namespace CASM

#endif
