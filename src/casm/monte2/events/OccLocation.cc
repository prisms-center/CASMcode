#include "casm/monte2/events/OccLocation.hh"

#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte2/Conversions.hh"
#include "casm/monte2/events/OccCandidate.hh"

namespace CASM {
namespace Monte2 {

OccLocation::OccLocation(const Conversions &_convert,
                         const OccCandidateList &_candidate_list)
    : m_convert(_convert),
      m_candidate_list(_candidate_list),
      m_loc(_candidate_list.size()),
      m_update_species(false) {}

/// Fill tables with occupation info
void OccLocation::initialize(Eigen::VectorXi const &occupation) {
  m_mol.clear();
  m_species.clear();
  m_l_to_mol.clear();
  for (auto &vec : m_loc) {
    vec.clear();
  }

  Index Nmut = 0;
  for (Index l = 0; l < occupation.size(); ++l) {
    Index asym = m_convert.l_to_asym(l);
    if (m_convert.occ_size(asym) > 1) {
      Nmut++;
    }
  }

  m_mol.resize(Nmut);
  m_l_to_mol.reserve(occupation.size());
  Index mol_id = 0;
  for (Index l = 0; l < occupation.size(); ++l) {
    Index asym = m_convert.l_to_asym(l);
    if (m_convert.occ_size(asym) > 1) {
      Index species_index = m_convert.species_index(asym, occupation[l]);
      Index cand_index = m_candidate_list.index(asym, species_index);

      Mol &mol = m_mol[mol_id];
      mol.id = mol_id;
      mol.l = l;
      mol.asym = asym;
      mol.species_index = species_index;
      mol.loc = m_loc[cand_index].size();

      if (m_update_species) {
        // only atoms now
        for (Index j = 0; j < 1; j++) {
          Species spec(m_convert.l_to_bijk(l));
          spec.species_index = species_index;
          spec.id = m_species.size();
          mol.component.push_back(spec.id);

          m_species.push_back(spec);
        }
      }

      m_loc[cand_index].push_back(mol_id);
      m_l_to_mol.push_back(mol_id);
      mol_id++;
    } else {
      m_l_to_mol.push_back(Nmut);
    }
  }
  if (m_update_species) {
    m_tmol = m_mol;
  }
}

/// Stochastically choose an occupant of a particular OccCandidate type
Mol const &OccLocation::choose_mol(Index cand_index, MTRand &mtrand) const {
  return mol(m_loc[cand_index][mtrand.randInt(m_loc[cand_index].size() - 1)]);
}

/// Stochastically choose an occupant of a particular OccCandidate type
Mol const &OccLocation::choose_mol(OccCandidate const &cand,
                                   MTRand &mtrand) const {
  return choose_mol(m_candidate_list.index(cand), mtrand);
}

/// Update occupation vector and this to reflect that event 'e' occurred
void OccLocation::apply(const OccEvent &e, Eigen::VectorXi &occupation) {
  // copy original Mol.component
  if (m_update_species) {
    for (const auto &occ : e.occ_transform) {
      m_tmol[occ.mol_id].component = m_mol[occ.mol_id].component;
    }
  }

  // update Mol and config occupation
  for (const auto &occ : e.occ_transform) {
    auto &mol = m_mol[occ.mol_id];

    // set config occupation
    occupation[mol.l] = m_convert.occ_index(mol.asym, occ.to_species);

    // remove from m_loc
    Index cand_index = m_candidate_list.index(mol.asym, mol.species_index);
    Index back = m_loc[cand_index].back();
    m_loc[cand_index][mol.loc] = back;
    m_mol[back].loc = mol.loc;
    m_loc[cand_index].pop_back();

    // set Mol.species index
    mol.species_index = occ.to_species;

    if (m_update_species) {
      mol.component.resize(m_convert.components_size(mol.species_index));
    }

    // add to m_loc
    cand_index = m_candidate_list.index(mol.asym, mol.species_index);
    mol.loc = m_loc[cand_index].size();
    m_loc[cand_index].push_back(mol.id);
  }

  if (m_update_species) {
    // update Mol.component
    for (const auto &traj : e.species_traj) {
      m_mol[traj.to.mol_id].component[traj.to.mol_comp] =
          m_tmol[traj.from.mol_id].component[traj.from.mol_comp];
    }
  }
}

/// Total number of mutating sites
OccLocation::size_type OccLocation::mol_size() const { return m_mol.size(); }

Mol &OccLocation::mol(Index mol_id) { return m_mol[mol_id]; }

const Mol &OccLocation::mol(Index mol_id) const { return m_mol[mol_id]; }

/// Total number of mutating sites, of OccCandidate type, specified by index
OccLocation::size_type OccLocation::cand_size(Index cand_index) const {
  return m_loc[cand_index].size();
}

/// Total number of mutating sites, of OccCandidate type
OccLocation::size_type OccLocation::cand_size(const OccCandidate &cand) const {
  return cand_size(m_candidate_list.index(cand));
}

/// The index into the configuration of a particular mutating site
Index OccLocation::mol_id(Index cand_index, Index loc) const {
  return m_loc[cand_index][loc];
}

/// The index into the configuration of a particular mutating site
Index OccLocation::mol_id(const OccCandidate &cand, Index loc) const {
  return mol_id(m_candidate_list.index(cand), loc);
}

/// Convert from config index to variable site index
Index OccLocation::l_to_mol_id(Index l) const { return m_l_to_mol[l]; }

}  // namespace Monte2
}  // namespace CASM
