#include "casm/monte_carlo/OccLocation.hh"
#include "casm/monte_carlo/Conversions.hh"
#include "casm/monte_carlo/OccCandidate.hh"
#include "casm/clex/Configuration.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"


namespace CASM {
  namespace Monte {

    OccLocation::OccLocation(const Conversions &_convert, const OccCandidateList &_cand) :
      m_convert(_convert),
      m_cand(_cand),
      m_loc(_cand.size()),
      m_kmc(false) {}

    /// Fill tables with occupation info
    void OccLocation::initialize(const Configuration &config) {

      m_mol.clear();
      m_species.clear();
      m_l_to_mol.clear();
      for(auto &vec : m_loc) {
        vec.clear();
      }

      Index Nmut = 0;
      for(Index l = 0; l < config.size(); ++l) {
        Index asym = m_convert.l_to_asym(l);
        if(m_convert.occ_size(asym) > 1) {
          Nmut++;
        }
      }

      m_mol.resize(Nmut);
      m_l_to_mol.reserve(config.size());
      Index mol_id = 0;
      for(Index l = 0; l < config.size(); ++l) {
        Index asym = m_convert.l_to_asym(l);
        if(m_convert.occ_size(asym) > 1) {
          Index species_index = m_convert.species_index(asym, config.occ(l));
          Index cand_index = m_cand.index(asym, species_index);

          Mol &mol = m_mol[mol_id];
          mol.id = mol_id;
          mol.l = l;
          mol.asym = asym;
          mol.species_index = species_index;
          mol.loc = m_loc[cand_index].size();

          if(m_kmc) {
            // only atoms now
            for(Index j = 0; j < 1; j++) {
              Species spec;
              spec.species_index = species_index;
              spec.id = m_species.size();
              spec.bijk_begin = m_convert.l_to_bijk(l);
              mol.component.push_back(spec.id);

              m_species.push_back(spec);
            }
          }

          m_loc[cand_index].push_back(mol_id);
          m_l_to_mol.push_back(mol_id);
          mol_id++;
        }
        else {
          m_l_to_mol.push_back(Nmut);
        }
      }
      if(m_kmc) {
        m_tmol = m_mol;
      }
    }

    /// Total number of mutating sites
    OccLocation::size_type OccLocation::size() const {
      return m_mol.size();
    }

    Mol &OccLocation::mol(Index mol_id) {
      return m_mol[mol_id];
    }

    const Mol &OccLocation::mol(Index mol_id) const {
      return m_mol[mol_id];
    }

    /// Total number of mutating sites, of OccCandidate type, specified by index
    OccLocation::size_type OccLocation::cand_size(Index cand_index) const {
      return m_loc[cand_index].size();
    }

    /// Total number of mutating sites, of OccCandidate type
    OccLocation::size_type OccLocation::cand_size(const OccCandidate &cand) const {
      return cand_size(m_cand.index(cand));
    }

    /// The index into the configuration of a particular mutating site
    Index OccLocation::mol_id(Index cand_index, Index loc) const {
      return m_loc[cand_index][loc];
    }

    /// The index into the configuration of a particular mutating site
    Index OccLocation::mol_id(const OccCandidate &cand, Index loc) const {
      return mol_id(m_cand.index(cand), loc);
    }

    /// Convert from config index to variable site index
    Index OccLocation::l_to_mol_id(Index l) const {
      return m_l_to_mol[l];
    }

    /// Propose canonical OccEvent
    OccEvent &OccLocation::propose_canonical(
      OccEvent &e,
      const std::vector<OccSwap> &canonical_swap,
      MTRand &mtrand) const {

      Index tsize = canonical_swap.size();
      m_tsum.resize(tsize + 1);

      m_tsum[0] = 0.;
      for(Index i = 0; i < tsize; ++i) {
        m_tsum[i + 1] = m_tsum[i] +
                        ((double) cand_size(canonical_swap[i].cand_a)) *
                        ((double) cand_size(canonical_swap[i].cand_b));
      }

      double rand = mtrand.randExc(m_tsum.back());

      for(Index i = 0; i < tsize; ++i) {
        if(rand < m_tsum[i + 1]) {
          return _propose(e, canonical_swap[i], mtrand);
        }
      }

      throw std::runtime_error("OccLocation::propose_canonical error");
    }

    /// Propose grand canonical OccEvent
    OccEvent &OccLocation::propose_grand_canonical(OccEvent &e, const OccSwap &swap, MTRand &mtrand) const {
      e.occ_transform.resize(1);
      e.species_traj.resize(0);

      Index index_cand_a = m_cand.index(swap.cand_a);
      Index index_cand_b = m_cand.index(swap.cand_b);

      OccTransform &f_a = e.occ_transform[0];
      f_a.mol_id = m_loc[index_cand_a][mtrand.randInt(cand_size(swap.cand_a) - 1)];
      f_a.l = m_mol[f_a.mol_id].l;
      f_a.asym = m_cand[index_cand_a].asym;
      f_a.from_species = m_cand[index_cand_a].species_index;
      f_a.to_species = m_cand[index_cand_b].species_index;

      return e;
    }

    /// Update configdof and this to reflect that event 'e' occurred
    void OccLocation::apply(const OccEvent &e, ConfigDoF &configdof) {

      // copy original Mol.component
      if(m_kmc) {
        for(const auto &occ : e.occ_transform) {
          m_tmol[occ.mol_id].component = m_mol[occ.mol_id].component;
        }
      }

      // update Mol and config occupation
      for(const auto &occ : e.occ_transform) {
        auto &mol = m_mol[occ.mol_id];

        // set config occupation
        configdof.occ(mol.l) = m_convert.occ_index(mol.asym, occ.to_species);

        // remove from m_loc
        Index cand_index = m_cand.index(mol.asym, mol.species_index);
        Index back = m_loc[cand_index].back();
        m_loc[cand_index][mol.loc] = back;
        m_mol[back].loc = mol.loc;
        m_loc[cand_index].pop_back();

        // set Mol.species index
        mol.species_index = occ.to_species;

        if(m_kmc) {
          mol.component.resize(m_convert.components_size(mol.species_index));
        }

        // add to m_loc
        cand_index = m_cand.index(mol.asym, mol.species_index);
        mol.loc = m_loc[cand_index].size();
        m_loc[cand_index].push_back(mol.id);

      }

      if(m_kmc) {
        // update Mol.component
        for(const auto &traj : e.species_traj) {
          m_mol[traj.to.mol_id].component[traj.to.mol_comp] = m_tmol[traj.from.mol_id].component[traj.from.mol_comp];
        }
      }
    }

    /// Canonical propose
    OccEvent &OccLocation::_propose(OccEvent &e, const OccSwap &swap, MTRand &mtrand, Index cand_a, Index cand_b, Index size_a, Index size_b) const {
      e.occ_transform.resize(2);
      e.species_traj.resize(0);

      OccTransform &f_a = e.occ_transform[0];
      f_a.mol_id = m_loc[cand_a][mtrand.randInt(size_a - 1)];
      f_a.l = m_mol[f_a.mol_id].l;
      f_a.asym = m_cand[cand_a].asym;
      f_a.from_species = m_cand[cand_a].species_index;
      f_a.to_species = m_cand[cand_b].species_index;
      //std::cout << "size_a: " << size_a << "  loc: " << m_mol[f_a.mol_id].loc << std::endl;

      OccTransform &f_b = e.occ_transform[1];
      f_b.mol_id = m_loc[cand_b][mtrand.randInt(size_b - 1)];
      f_b.l = m_mol[f_b.mol_id].l;
      f_b.asym = m_cand[cand_b].asym;
      f_b.from_species = m_cand[cand_b].species_index;
      f_b.to_species = m_cand[cand_a].species_index;
      //std::cout << "size_b: " << size_b << "  loc: " << m_mol[f_b.mol_id].loc << std::endl;

      return e;
    }

    /// Canonical propose
    OccEvent &OccLocation::_propose(OccEvent &e, const OccSwap &swap, MTRand &mtrand) const {
      Index cand_a = m_cand.index(swap.cand_a);
      Index cand_b = m_cand.index(swap.cand_b);
      Index size_a = m_loc[cand_a].size();
      Index size_b = m_loc[cand_b].size();
      return _propose(e, swap, mtrand, cand_a, cand_b, size_a, size_b);
    }
  }
}
