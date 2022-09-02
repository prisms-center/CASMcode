#ifndef CASM_Monte_OccLocation_HH
#define CASM_Monte_OccLocation_HH

#include <vector>

#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"

class MTRand;

namespace CASM {

class Configuration;
class ConfigDoF;

}  // namespace CASM

namespace CASM {
namespace Monte {

class Conversions;
struct OccCandidate;
class OccSwap;
class OccCandidateList;

/// \brief Represents an indivisible molecule component
struct Species {
  Species(xtal::UnitCellCoord _bijk_begin) : bijk_begin(_bijk_begin) {}

  Index species_index;             ///< Species type index
  Index id;                        ///< Location in OccLocation.m_species
  xtal::UnitCellCoord bijk_begin;  ///< Saves initial position
  Index mol_comp_begin;            ///< Saves initial Mol.component index
};

/// \brief Represents the occupant on a site
///
/// - May be divisible Indexo components or indivisible
struct Mol {
  Index id;             ///< Location in OccLocation.m_mol
  Index l;              ///< Location in config
  Index asym;           ///< Asym unit index (must be consistent with l)
  Index species_index;  ///< Species type index (must be consistent with
                        ///< config.occ(l))
  std::vector<Index>
      component;  ///< Location of component Specie in OccLocation.m_species
  Index loc;      ///< Location in OccLocation.m_loc
};

struct OccTransform {
  Index l;             ///< Config occupant that is being transformed
  Index mol_id;        ///< Location in OccLocation.m_mol
  Index asym;          ///< Asym index
  Index from_species;  ///< Species index before transformation
  Index to_species;    ///< Species index after transformation
};

struct SpeciesLocation {
  Index l;         ///< Config occupant that is being transformed
  Index mol_id;    ///< Location in OccLocation.m_mol
  Index mol_comp;  ///< Location in mol.components
};

struct SpeciesTraj {
  SpeciesLocation from;
  SpeciesLocation to;
};

struct OccEvent {
  /// \brief Linear site indices, indicating on which sites the occupation will
  ///     be modified
  std::vector<Index> linear_site_index;

  /// \brief Occupant indices, indicating the new occupation index on the sites
  ///     being modified
  std::vector<int> new_occ;

  std::vector<OccTransform> occ_transform;
  std::vector<SpeciesTraj> species_traj;
};

/// \brief Stores data to enable efficient proposal and update of occupation
/// mutation
class OccLocation {
 public:
  typedef Index size_type;

  OccLocation(const Conversions &_convert, const OccCandidateList &_cand);

  /// Fill tables with occupation info
  void initialize(const Configuration &config);

  /// Total number of mutating sites
  size_type size() const;

  Mol &mol(Index mol_id);

  const Mol &mol(Index mol_id) const;

  /// Total number of mutating sites, of OccCandidate type, specified by index
  size_type cand_size(Index cand_index) const;

  /// Total number of mutating sites, of OccCandidate type
  size_type cand_size(const OccCandidate &cand) const;

  /// Mol.id of a particular OccCandidate type
  Index mol_id(Index cand_index, Index loc) const;

  /// Mol.id of a particular OccCandidate type
  Index mol_id(const OccCandidate &cand, Index loc) const;

  /// Convert from config index to variable site index
  Index l_to_mol_id(Index l) const;

  /// Propose canonical OccEvent
  OccEvent &propose_canonical(OccEvent &e,
                              const std::vector<OccSwap> &canonical_swap,
                              MTRand &mtrand) const;

  /// Propose grand canonical OccEvent
  OccEvent &propose_grand_canonical(OccEvent &e, const OccSwap &swap,
                                    MTRand &mtrand) const;

  /// Update configdof and this to reflect that event 'e' occurred
  void apply(const OccEvent &e, ConfigDoF &configdof);

 private:
  /// Canonical propose
  OccEvent &_propose(OccEvent &e, const OccSwap &swap, MTRand &mtrand,
                     Index cand_a, Index cand_b, Index size_a,
                     Index size_b) const;

  /// Canonical propose
  OccEvent &_propose(OccEvent &e, const OccSwap &swap, MTRand &mtrand) const;

  const Conversions &m_convert;

  const OccCandidateList &m_cand;

  /// Gives a list of all Mol of the same {asym, species}-type allowed to mutate
  ///   m_loc[cand_index][i] -> m_mol index
  std::vector<std::vector<Index> > m_loc;

  /// Holds Monte::Species objects
  std::vector<Species> m_species;

  /// Holds Mol objects, one for each mutating site in the configuration
  std::vector<Mol> m_mol;

  /// l_to_mol[l] -> Mol.id, m_mol.size() otherwise
  std::vector<Index> m_l_to_mol;

  /// If true, update Species location during apply
  bool m_kmc;

  /// Data structure used store temporaries during apply
  std::vector<Mol> m_tmol;

  /// Data used by propose_canonical
  mutable std::vector<double> m_tsum;
};

}  // namespace Monte
}  // namespace CASM

#endif
