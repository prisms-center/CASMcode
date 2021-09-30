#include "casm/monte2/events/OccEventProposal.hh"

#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte2/Conversions.hh"
#include "casm/monte2/events/OccCandidate.hh"
#include "casm/monte2/events/OccLocation.hh"

namespace CASM {
namespace Monte2 {

/// Choose a swap type from a list of allowed canonical swap types
///
/// \param occ_location Contains lookup table with occupant locations for
/// efficient choosing in dilute systems.
/// \param canonical_swap List of allowed swap types (OccSwap). Swap types
/// consists of two pairs (asymmetric unit a, species_index a) and (asymmetric
/// unit b, species_index b) defining what will change, but not which sites.
/// For canonical swaps, the species indices must be the different.
/// \param mtrand The random number generator used to stochastically choose the
/// swap.
///
/// Method:
/// Stochastically chooose which swap type will be chosen. To do so, calculate
/// the cumulative number of swaps of each swap type and use a random number to
/// choose which type occurs.
///
OccSwap const &choose_canonical_swap(OccLocation const &occ_location,
                                     std::vector<OccSwap> const &canonical_swap,
                                     MTRand &mtrand) {
  // Calculate m_tsum[i]:
  // - The total number of possible events of swap types [0, i)`
  // - The total number of each swap type is
  //   `cand_size(canonical_swap[i].cand_a) *
  //    cand_size(canonical_swap[i].cand_b)`
  // - For example: System with two asym unit sites (asym1, asym2), and two
  //   species (A, B), allowed on both asym unit sites. The number of swap type
  //   (A on asym1 <-> B on asym1) is equal to the number of species A on asym1
  //   times the number of species B on asym1.
  //
  // Notes:
  // - m_tsum[0] is 0.0
  // - m_tsum[canonical_swap.size()] is total number of possible events
  //
  Index tsize = canonical_swap.size();
  static std::vector<double> m_tsum;
  m_tsum.resize(tsize + 1);

  m_tsum[0] = 0.;
  for (Index i = 0; i < tsize; ++i) {
    m_tsum[i + 1] =
        m_tsum[i] +
        ((double)occ_location.cand_size(canonical_swap[i].cand_a)) *
            ((double)occ_location.cand_size(canonical_swap[i].cand_b));
  }

  if (m_tsum.back() == 0.0) {
    throw std::runtime_error(
        "Error in choose_canonical_swap: No events possible.");
  }

  // Choose a random number on [0, m_tsum[canonical_swap.size()])
  // Swap type canonical_swap[i] occurs if the random number is between
  // m_tsum[i] and m_tsum[i+1]
  double rand = mtrand.randExc(m_tsum.back());

  for (Index i = 0; i < tsize; ++i) {
    if (rand < m_tsum[i + 1]) {
      return canonical_swap[i];
    }
  }

  throw std::runtime_error("Error in choose_canonical_swap");
}

/// Propose canonical OccEvent, given choice of OccSwap
///
/// \param e [out] The OccEvent that will be populated with the proposed event.
/// \param occ_location Contains lookup table with occupant locations for
/// efficient choosing in dilute systems.
/// \param swap: Type of swap, defining what will change, but not which sites.
/// For canonical swaps, the species indices must be the different.
/// \param mtrand The random number generator used to stochastically choose the
/// event.
///
/// Method:
/// Stochastically choose which two sites (consistent with the given swap type)
/// will be swapped by using the list of sites of each candidate type.
OccEvent &propose_canonical_event(OccEvent &e, OccLocation const &occ_location,
                                  OccSwap const &swap, MTRand &mtrand) {
  e.occ_transform.resize(2);
  e.species_traj.resize(0);

  OccTransform &transform_a = e.occ_transform[0];
  Mol const &mol_a = occ_location.choose_mol(swap.cand_a, mtrand);
  transform_a.mol_id = mol_a.id;
  transform_a.l = mol_a.l;
  transform_a.asym = swap.cand_a.asym;
  transform_a.from_species = swap.cand_a.species_index;
  transform_a.to_species = swap.cand_b.species_index;

  OccTransform &transform_b = e.occ_transform[1];
  Mol const &mol_b = occ_location.choose_mol(swap.cand_b, mtrand);
  transform_b.mol_id = mol_b.id;
  transform_b.l = mol_b.l;
  transform_b.asym = swap.cand_b.asym;
  transform_b.from_species = swap.cand_b.species_index;
  transform_b.to_species = swap.cand_a.species_index;

  return e;
}

/// Propose canonical OccEvent
///
/// \param e [out] The OccEvent that will be populated with the proposed event.
/// \param occ_location Contains lookup table with occupant locations for
/// efficient choosing in dilute systems.
/// \param canonical_swap: List of allowed swap types (OccSwap). Swap types
/// consists of two pairs (asymmetric unit a, species_index a) and (asymmetric
/// unit b, species_index b) defining what will change, but not which sites.
/// For canonical swaps, the species indices must be the different. Do not
/// include reverse swaps (i.e. a->b and b->a).
/// \param mtrand The random number generator used to stochastically choose the
/// event.
///
/// Method:
/// - First, stochastically chooose which swap type will be chosen:
///   - To do so, calculate the cumulative number of swaps of each
///     swap type and use a random number to choose which type occurs.
/// - Second, stochastically choose which two sites (consistent with the chosen
///   swap type) will be swapped by using the list of sites of each candidate
///   type.
OccEvent &propose_canonical_event(OccEvent &e, OccLocation const &occ_location,
                                  const std::vector<OccSwap> &canonical_swap,
                                  MTRand &mtrand) {
  auto const &swap =
      choose_canonical_swap(occ_location, canonical_swap, mtrand);
  return propose_canonical_event(e, occ_location, swap, mtrand);
}

/// Choose a swap type from a list of allowed grand canonical swap types
///
/// \param occ_location Contains lookup table with occupant locations for
/// efficient choosing in dilute systems.
/// \param grand_canonical_swap List of allowed swap types (OccSwap). Swap types
/// consists of two pairs (asymmetric unit a, species_index a) and (asymmetric
/// unit b, species_index b) defining what will change, but not which sites.
/// For grand canonical swaps, the species indices must be the different and
/// the asymmetric unit index must be the same. Do include reverse swaps (i.e.
/// a->b and b->a).
/// \param mtrand The random number generator used to stochastically choose the
/// swap.
///
/// Method:
/// Stochastically chooose which swap type will be chosen. To do so, calculate
/// the cumulative number of swaps of each swap type and use a random number to
/// choose which type occurs.
///
OccSwap const &choose_grand_canonical_swap(
    OccLocation const &occ_location,
    std::vector<OccSwap> const &grand_canonical_swap, MTRand &mtrand) {
  // Calculate m_tsum[i]:
  // - The total number of possible events of swap types [0, i)`
  // - The total number of each swap type is
  //   `cand_size(canonical_swap[i].cand_a)`
  // - For example: System with two asym unit sites (asym1, asym2), and two
  //   species (A, B), allowed on both asym unit sites. The number of swap type
  //   (A on asym1 <-> B on asym1) is equal to the number of species A on asym1
  //   times the number of species B on asym1.
  //
  // Notes:
  // - m_tsum[0] is 0.0
  // - m_tsum[canonical_swap.size()] is total number of possible events
  //
  Index tsize = grand_canonical_swap.size();
  static std::vector<double> m_tsum;
  m_tsum.resize(tsize + 1);

  m_tsum[0] = 0.;
  for (Index i = 0; i < tsize; ++i) {
    m_tsum[i + 1] =
        m_tsum[i] +
        ((double)occ_location.cand_size(grand_canonical_swap[i].cand_a));
  }

  if (m_tsum.back() == 0.0) {
    throw std::runtime_error(
        "Error in choose_grand_canonical_swap: No events possible.");
  }

  // Choose a random number on [0, m_tsum[grand_canonical_swap.size()])
  // Swap type grand_canonical_swap[i] occurs if the random number is between
  // m_tsum[i] and m_tsum[i+1]
  double rand = mtrand.randExc(m_tsum.back());

  for (Index i = 0; i < tsize; ++i) {
    if (rand < m_tsum[i + 1]) {
      return grand_canonical_swap[i];
    }
  }

  throw std::runtime_error("Error in choose_grand_canonical_swap");
}

/// Propose grand canonical OccEvent
///
/// Given the choice of OccSwap, this method stochastically chooses which site
/// the swap occurs on.
OccEvent &propose_grand_canonical_event(OccEvent &e,
                                        OccLocation const &occ_location,
                                        OccSwap const &swap, MTRand &mtrand) {
  e.occ_transform.resize(1);
  e.species_traj.resize(0);

  OccTransform &transform = e.occ_transform[0];
  Mol const &mol = occ_location.choose_mol(swap.cand_a, mtrand);
  transform.mol_id = mol.id;
  transform.l = mol.l;
  transform.asym = swap.cand_a.asym;
  transform.from_species = swap.cand_a.species_index;
  transform.to_species = swap.cand_b.species_index;

  return e;
}

/// Propose grand canonical OccEvent
OccEvent &propose_grand_canonical_event(
    OccEvent &e, OccLocation const &occ_location,
    std::vector<OccSwap> const &grand_canonical_swap, MTRand &mtrand) {
  auto const &swap =
      choose_grand_canonical_swap(occ_location, grand_canonical_swap, mtrand);
  return propose_grand_canonical_event(e, occ_location, swap, mtrand);
}

}  // namespace Monte2
}  // namespace CASM
