#ifndef CASM_CanonicalEvent_HH
#define CASM_CanonicalEvent_HH

#include <vector>
#include "casm/external/Eigen/Dense"
#include "casm/CASM_global_definitions.hh"
#include "casm/monte_carlo/DoFMod.hh"

namespace CASM {

  namespace {
    Configuration default_prim_config(PrimClex &primclex) {
      auto &scel = primclex.get_supercell(primclex.add_supercell(primclex.get_prim().lattice()));
      Configuration res(scel);
      res.init_occupation();
      return res;
    }
  }

  /// l: linear index into mc_scel
  /// b: prim basis site index
  /// ijk: prim unit cell indices
  /// bijk: prim basis site index + unit cell indices
  /// unitl: ref config basis site index
  /// asym: asymmetric unit index
  class MCConversions {

  public:

    MCConversions(const Supercell &mc_scel) :
      MCConversions(default_prim_config(mc_scel.get_primclex()), mc_scel) {}


    MCConversions(const Configuration &unit_config, const Supercell &mc_scel) :
      m_unit_scel(&unit_config.supercell()),
      m_mc_scel(&mc_scel),
      m_struc_mol(m_mc_scel.get_prim().get_struc_molecule()),
      m_struc_molname(m_mc_scel.get_prim().get_struc_molecule_name()) {

      // make m_unitl_to_asym, m_Nasym
      int asym = 0;
      Index unit_Nsites = m_unit_scel->num_sites();
      m_unitl_to_asym.resize(unit_Nsites, -1);
      std::vector<PermuteIterator> fg = unit_config.factor_group();
      for(int unitl = 0; unitl < unit_Nsites; ++unitl) {
        if(m_unitl_to_asym[unitl] == -1) {
          for(auto it = fg.begin(); it != fg.end(); ++it) {
            // permutation defined by: after[i] = before[it->permute_ind(i)]
            m_unitl_to_asym[it->permute_ind(unitl)] = asym;
          }
          ++asym;
        }
      }
      m_Nasym = asym;

      // make m_asym_to_unitl & m_asym_to_b
      m_asym_to_unitl.resize(Nasym);
      m_asym_to_b.resize(Nasym);
      for(int unitl = 0; unitl < unit_Nsites; ++unitl) {
        int asym = m_unit_to_asym[unitl];
        m_asym_to_unitl[asym].insert(unitl);
        m_asym_to_b[asym].insert(unitl_to_b(unitl));
      }

      // make m_occ_to_species and m_species_to_occ

      // [b][occ] -> species
      auto index_converter = get_index_converter(m_mc_scel.get_prim(), m_struc_molname);

      // [b][species] -> occ, index_converter[b].size() if not allowed
      auto index_converter_inv = get_index_converter_inverse(m_mc_scel.get_prim(), m_struc_molname);

      m_occ_to_species.resize(Nasym);
      m_species_to_occ.resize(Nasym);
      for(int asym = 0; asym < Nasym; ++asym) {
        int b = *m_asym_to_b.begin();
        m_occ_to_species[asym] = index_converter[b];
        m_species_to_occ[asym] = index_converter_inv[b];
      }
    }


    int l_to_b(Index l) const {
      return m_mc_scel->get_b(l);
    }
    UnitCell l_to_ijk(Index l) const {
      return m_mc_scel->prim_grid().unitcell(l % m_mc_scel->volume());
    }
    UnitCellCoord l_to_bijk(Index l) const {
      return m_mc_scel->uccoord(l);
    }
    int l_to_unitl(Index l) const {
      return bijk_to_unitl(l_to_bijk(l));
    }
    int l_to_asym(Index l) const {
      return m_unitl_to_asym[l_to_unitl(l)];
    }

    Index bijk_to_l(const UnitCellCoord &bijk) const {
      return m_mc_scel->find(bijk);
    }
    int bijk_to_unitl(const UnitCellCoord &bijk) const {
      return m_unit_scel->find(bijk);
    }
    int bijk_to_asym(const UnitCellCoord &bijk) const {
      return l_to_asym(bijk_to_l(bijk));
    }

    int unitl_to_b(int unitl) const {
      return m_unit_scel->get_b(unitl);
    }
    UnitCellCoord unitl_to_bijk(int unitl) const {
      return m_unit_scel->uccoord(unitl);
    }
    int unitl_to_asym(int unitl) const {
      return m_unitl_to_asym[unitl];
    }

    int asym_size() const {
      return m_Nasym;
    }
    const std::set<int> &asym_to_b(int asym) const {
      return m_asym_to_b[asym];
    }
    const std::set<int> &asym_to_unitl(int asym) const {
      return m_asym_to_unitl[asym];
    }

    const Supercell &unit_scel() const {
      return *m_unit_scel;
    }
    const Supercell &mc_scel() const {
      return *m_mc_scel;
    }

    int occ_size(int asym) const {
      return m_occ_to_species[asym].size();
    }
    int species_index(int asym, int occ_index) const {
      return m_occ_to_species[asym][occ_index];
    }
    int occ_index(int asym, int species_index) const {
      // returns occ_size(asym) if species not allowed
      return m_species_to_occ[asym][species_index];
    }
    bool species_allowed(int asym, int species_index) const {
      return occ_index(asym, species_index) != occ_size(asym);
    }

    int species_size() const {
      return m_struc_mol.size();
    }
    const Molecule &species_to_mol(int species_index) const {
      return m_struc_mol[species_index];
    }
    const std::string &species_name(int species_index) const {
      return m_struc_mol_name[species_index];
    }

  private:

    const Supercell *m_unit_scel;
    const Supercell *m_mc_scel;
    std::vector<Molecule> m_struc_mol;
    std::vector<std::string> m_struc_molname;

    int m_Nasym;
    std::vector<int> m_unitl_to_asym;
    std::vector<std::set<int> > m_asym_to_unitl;
    std::vector<std::set<int> > m_asym_to_b;

    /// m_occ_to_species[asym][occ_index] -> species_index
    std::vector<std::vector<int> > m_occ_to_species;

    /// m_species_to_occ[asym][species_index] -> occ_index
    std::vector<std::vector<int> > m_species_to_occ;
  };


  /// \brief Represents an indivisible molecule component
  struct MCSpecies {
    int species_index;           ///< Species type index
    Index id;                    ///< Location in MCOccLocation.m_species
    UnitCellCoord bijk_begin;    ///< Saves initial position
    int mol_comp_begin;          ///< Saves initial MCMol.component index
  };

  /// \brief Represents the occupant on a site
  ///
  /// - May be divisible into components or indivisible
  struct MCMol {
    Index id;                             ///< Location in MCOccLocation.m_mol
    Index l;                              ///< Location in config
    int asym;                             ///< Asym unit index (must be consistent with l)
    int species_index;                    ///< Species type index (must be consistent with config.occ(l))
    std::vector<Index> component;         ///< Location of component MCSpecie in MCOccLocation.m_species
    Index loc;                            ///< Location in MCOccLocation.m_loc
  };

  struct MCOccTransform {
    Index l;              ///<Config occupant that is being transformed
    Index mol_id;         ///<Location in MCOccLocation.m_mol
    int asym;             ///<Asym index
    int from_species;     ///<Species index before transformation
    int to_species;       ///<Species index after transformation
  };

  struct MCSpecieLocation {
    Index l;        ///<Config occupant that is being transformed
    Index mol_id;   ///<Location in MCOccLocation.m_mol
    int mol_comp;   ///<Location in MCmol.components
  };

  struct MCSpecieTraj {
    MCSpecieLocation from;
    MCSpecieLocation to;
  };

  struct MCOccEvent {
    std::vector<MCOccTransformation> occ_transform;
    std::vector<MCSpecieTraj> specie_traj;
  };


  /// \brief Stores data to enable efficient proposal and update of occupation mutation
  class MCOccLocation {

  public:

    MCOccLocation(const MCConversions &_convert, const CandidateList &_cand) :
      m_convert(_occ_converter),
      m_cand(_cand),
      m_loc(_cand.size()) {}

    /// Fill tables with occupation info
    void initialize(const Configuration &config) {

      m_mol.clear();
      m_species.clear();
      m_l_to_mol.clear();
      for(auto &vec : m_loc) {
        vec.clear();
      }

      Index Nmut = 0;
      for(Index l = 0; l < config.size(); ++l) {
        int asym = m_convert.l_to_asym(l);
        if(m_convert.occ_size(asym) > 1) {
          Nmut++;
        }
      }

      m_mol.resize(Nmut);
      m_l_to_mol.resize(config.size());
      Index mol_id = 0;
      for(Index l = 0; l < config.size(); ++l) {
        int asym = m_convert.l_to_asym(l);
        if(m_convert.occ_size(asym) > 1) {
          int species_index = m_convert.species_index(asym, config.occ(l));
          int cand_index = m_cand.index(asym, species_index);

          MCMol &mol = m_mol[mol_id];
          mol.id = mol_id;
          mol.l = l;
          mol.asym = asym;
          mol.species_index = species_index;
          mol.cand_id = m_loc[cand_index].size();

          // only atoms now
          for(int j = 0; j < 1; j++) {
            MCSpecies spec;
            spec.species_index = species_index;
            spec.id = m_species.size();
            spec.bijk_begin = m_convert.l_to_bijk(l);
            spec.mol_comp = j;
            mol.component.push_back(spec.id);

            m_species.push_back(spec);
          }

          m_loc[cand_index].push_back(mol_id);
          m_l_to_mol.push_back(mol_id);
          mol_id++;
        }
        else {
          m_l_to_mol.push_back(Nmut);
        }
      }
    }

    /// Total number of mutating sites
    size_type size() const {
      return m_mol.size();
    }

    MCMol &mol(Index mol_id) {
      return m_mol[mol_id];
    }

    const MCMol &mol(Index mol_id) const {
      return m_mol[mol_id];
    }

    /// Total number of mutating sites, of OccCandidate type, specified by index
    size_type cand_size(int cand_index) const {
      return m_loc[cand_index].size();
    }

    /// Total number of mutating sites, of OccCandidate type
    size_type cand_size(const OccCandidate &cand) const {
      return cand_size(m_cand.index(cand));
    }

    /// The index into the configuration of a particular mutating site
    Index mol_id(int cand_index, Index loc) const {
      return m_loc[cand_index][loc];
    }

    /// The index into the configuration of a particular mutating site
    Index mol_id(const OccCandidate &cand, Index index) const {
      return mol_id(m_cand.index(cand), index);
    }

    /// Propose canonical MCOccEvent
    MCOccEvent propose_canonical(const std::vector<MCOccSwap> &canonical_swap, MTRand &mtrand) const {
      double sum = 0.;
      for(const auto &swap : canonical_swap) {
        Index index_cand_a = m_cand.index(swap.cand_a);
        Index index_cand_b = m_cand.index(swap.cand_b);
        Index size_cand_a = m_loc[index_cand_a].size();
        Index size_cand_b = m_loc[index_cand_b].size();

        sum += ((double) size_cand_a) * ((double) size_cand_b);
      }

      double rand = mtrand.randExc(sum);

      sum = 0.0;
      for(const auto &swap : canonical_swap) {
        Index index_cand_a = m_cand.index(swap.cand_a);
        Index index_cand_b = m_cand.index(swap.cand_b);
        Index size_cand_a = m_loc[index_cand_a].size();
        Index size_cand_b = m_loc[index_cand_b].size();

        sum += ((double) size_cand_a) * ((double) size_cand_b);
        if(rand < sum) {
          return _propose(swap, mtrand, index_cand_a, index_cand_b, size_cand_a, size_cand_b);
        }
      }
      throw std::runtime_error("MCOccLocation::propose_canonical error");
    }

    /// Propose grand canonical MCOccEvent
    MCOccEvent propose_grand_canonical(const MCOccSwap &swap, MTRand &mtrand) const {
      MCOccEvent e;
      e.occ_transform.resize(1);

      int index_cand_a = m_cand.index(swap.cand_a);
      int index_cand_b = m_cand.index(swap.cand_b);

      MCOccTransform &f_a = e.occ_transform[0];
      f_a.mol_id = m_loc[swap.cand_a][mtrand.randInt(cand_size(swap.cand_a) - 1)];
      f_a.l = m_mol[f_a.mol_id].l;
      f_a.asym = m_cand[index_cand_a].asym;
      f_a.from_species = m_cand[index_cand_a].species_index;
      f_a.to_species = m_cand[index_cand_b].species_index;

      return e;
    }

    /// Update config and this to reflect that event 'e' occurred
    void apply(const MCOccEvent &e, Configuration &config) {

      // copy original MCMol.components
      std::map<Index, std::vector<Index> > tmol;
      for(const auto &occ : occ_transform) {
        tmol.insert({ occ.mol_id, m_mol[occ.mol_id].components});
      }

      // update MCMol and config occupation
      for(const auto &occ : e.occ_transform) {
        auto &mol = m_mol[occ.mol_id];

        // set config occupation
        config.set_occ(mol.l, m_convert.occ_index(mol.asym, mol.species_index));

        // remove from m_loc
        int cand_index = m_cand.index(mol.asym, mol.species_index);
        Index loc_before = mol.loc;
        m_loc[cand_index][loc_before] = m_loc[cand_index].back();
        m_loc[cand_index].pop_back();

        // set MCMol.species index && resize MCMol.components
        mol.species_index = occ.to_species;
        mol.components.resize(m_convert.components_size(mol.species_index));

        // add to m_loc
        cand_index = m_cand.index(mol.asym, mol.species_index);
        mol.loc = m_loc[cand_index].size();
        m_loc[cand_index].push_back(mol.id);
      }

      // update MCMol.components
      for(const auto &traj : e.specie_traj) {
        m_mol[traj.to.mol_id].component[traj.to.mol_comp] = tmol[traj.from.mol_id][traj.from.mol_comp];
      }
    }

  private:

    /// Canonical propose
    MCOccEvent &_propose(MCOccSwap &swap, MTRand &mtrand, Index cand_a, Index cand_b, Index size_a, Index size_b) const {
      MCOccEvent e;
      e.occ_transform.resize(2);

      MCOccTransform &f_a = e.occ_transform[0];
      f_a.mol_id = m_loc[cand_a][mtrand.randInt(size_a - 1)];
      f_a.l = m_mol[f_a.mol_id].l;
      f_a.asym = m_cand[cand_a].asym;
      f_a.from_species = m_cand[cand_a].species_index;
      f_a.to_species = m_cand[cand_b].species_index;

      MCOccTransform &f_b = e.occ_transform[1];
      f_b.mol_id = m_loc[cand_b][mtrand.randInt(size_b - 1)];
      f_b.l = m_mol[f_b.mol_id].l;
      f_b.asym = m_cand[cand_b].asym;
      f_b.from_species = m_cand[cand_b].species_index;
      f_b.to_species = m_cand[cand_a].species_index;

      return e;
    }

    const MCConversions &m_converter;

    const CandidateList &m_cand;

    /// Gives a list of all MCMol of the same {asym, species}-type allowed to mutate
    ///   m_loc[asym][i] -> m_mol index
    std::vector<std::vector<Index> > m_loc;

    /// Holds MCSpecies objects
    std::vector<MCSpecies> m_species;

    /// Holds MCMol objects, one for each mutating site in the configuration
    std::vector<MCMol> m_mol;

    /// l_to_mol[l] -> MCMol.id, m_mol.size() otherwise
    std::vector<Index> m_l_to_mol;

  };


  struct OccCandidate : public Comparsisons<OccCandidate> {

    OccCandidate(int _asym, int _species_index) :
      asym(_asym),
      species_index(_species_index) {}

    int asym;
    int species_index;

    bool operator<(OccCandidate A, OccCandidate B) const {
      if(A.asym != B.asym) {
        return A.asym < B.asym;
      }
      return A.species_index < B.species_index;
    }
  };

  /// List of asym / species_index pairs indicating allowed variable occupation dof
  class CandidateList {

  public:

    typedef const_iterator std::vector<OccCandidate>::const_iterator;

    CandidateList() {}

    CandidateList(const MCConversions &convert) {

      // create set of 'candidate' asym / species pairs
      m_candidate.clear();
      for(int asym = 0; asym < convert.asym_size(); ++asym) {

        // hard code allowed sublattices: >1 allowed occupant
        if(convert.occ_size(asym) < 2) {
          continue;
        }

        // add candidates
        for(int i = 0; i < convert.occ_size(asym); ++i) {
          m_candidate.push_back(OccCandidate(asym, convert.species_index(asym, i)));
        }
      }

      // create lookup table of asym, species_index -> candidate index,
      //   will return {Nasym, Nspecies} if {asym, species_index} not allowed
      int Nspecies = convert.species_size();
      int Nasym = convert.asym_size();
      m_end = _candidate.size();
      std::vector<int> unallowed(Nspecies, end);
      m_species_to_cand_index = std::vector<std::vector<int> >(Nasym, unallowed);

      int index = 0;
      for(const auto &cand : _candidate) {
        m_species_to_cand_index[cand.asym][cand.species_index] = index;
        ++index;
      }

      // make canonical and grand canonical swaps
      _make_possible_swaps(convert);
    }

    /// Return index into std::vector<OccCandidate>, or _candidate.size() if not allowed
    int index(const OccCandidate &cand) const {
      return m_species_to_cand_index[cand.asym][cand.species_index];
    }

    /// Return index into std::vector<OccCandidate>, or _candidate.size() if not allowed
    int index(int asym, int species_index) const {
      return m_species_to_cand_index[asym][species_index];
    }

    const OccCandidate &operator[](int candidate_index) const {
      return m_candidate[candidate_index];
    }

    const_iterator begin() const {
      return m_candidate.begin();
    }

    const_iterator end() const {
      return m_candidate.end();
    }

    Index size() const {
      return m_end;
    }

    const std::vector<MCOccSwap> &canonical_swap() const {
      return m_canonical_swap;
    }

    const std::vector<MCOccSwap> &grand_canonical_swap() const {
      return m_grand_canonical_swap;
    }

  private:

    /// \brief Construct m_canonical_swaps, m_grand_canonical_swaps
    ///
    /// - Currently settings is not used, but we could add restrictions
    void Canonical::_make_possible_swaps(const MCConversions &convert) {

      // construct canonical and grand canonical swaps
      m_canonical_swaps.clear();
      m_grand_canonical_swaps.clear();

      // check that species are different and allowed on both sites
      auto allowed_canonical_swap = [&](OccCandidate cand_a, OccCandidate cand_b) {
        return cand_a.species_index != cand_b.species_index &&
               convert.species_allowed(cand_a.asym_unit, cand_b.species_index) &&
               convert.species_allowed(cand_b.asym_unit, cand_a.species_index);
      };

      // check that asym_unit is the same and species_index is different
      auto allowed_grand_canonical_swap = [&](OccCandidate cand_a, OccCandidate cand_b) {
        return cand_a.asym_unit == cand_b.asym_unit &&
               cand_a.species_index != cand_b.species_index;
      };

      // for each pair of candidates, check if they are allowed to swap
      for(const auto &cand_a : m_candidate) {
        for(const auto &cand_b : m_candidate) {

          // don't repeat a->b, b->a
          // and check that cand_b's species is allowed on cand_a's sublat && vice versa
          if(cand_a < cand_b && allowed_canonical_swap(cand_a, cand_b)) {
            m_canonical_swaps.push_back(MCOccSwap(cand_a, cand_b));
          }

          // allow a->b, b->a
          // check that asym_unit is the same and species_index is different
          if(allowed_grand_canonical_swap(cand_a, cand_b)) {
            m_grand_canonical_swaps.push_back(MCOccSwap(cand_a, cand_b));
          }
        }
      }
    }

    /// m_converter[asym][species_index] -> candidate_index
    std::vector<std::vector<int> > m_species_to_cand_index;

    std::vector<OccCandidate> m_candidate;

    /// Number of allowed candidates, what is returned if a candidate is not allowed
    int m_end;

    /// vector of allowed canonical swaps
    std::vector<MCOccSwap> m_canonical_swap;

    /// vector of allowed grand canonical swaps
    std::vector<MCOccSwap> m_grand_canonical_swap;

  };


  /// \brief Hold lookup tables of OccCandidate -> config index
  ///
  /// - This avoids proposing events that are not allowed or not possible
  class OccLocation {

  public:

    void prepare(const std::vector<CanonicalSwap> &possible) {
      m_canonical_swap = possible;
    }

    /// Canonical propose
    CanonicalSwap &propose(MTRand &mtrand) const {
      double sum = 0.;
      for(const auto &swap : m_canonical_swap) {
        Index index_cand_a = m_cand.index(swap.cand_a);
        Index index_cand_b = m_cand.index(swap.cand_b);
        Index size_cand_a = m_loc[index_cand_a].size();
        Index size_cand_b = m_loc[index_cand_b].size();

        sum += ((double) size_cand_a) * ((double) size_cand_b);
      }

      double rand = mtrand.randExc(sum);

      sum = 0.0;
      for(const auto &swap : m_canonical_swap) {
        Index index_cand_a = m_cand.index(swap.cand_a);
        Index index_cand_b = m_cand.index(swap.cand_b);
        Index size_cand_a = m_loc[index_cand_a].size();
        Index size_cand_b = m_loc[index_cand_b].size();

        sum += ((double) size_cand_a) * ((double) size_cand_b);
        if(rand < sum) {
          return _propose(swap, mtrand, index_cand_a, index_cand_b, size_cand_a, size_cand_b);
        }
      }
      throw std::runtime_error("Propose canonical swap error");
    }

    /// Canonical propose
    CanonicalSwap &propose(GrandCanonicalSwap &swap, MTRand &mtrand) const {
      Index index_cand_a = m_cand.index(swap.cand_a);
      Index index_cand_b = m_cand.index(swap.cand_b);
      Index size_cand_a = m_loc[index_cand_a].size();
      Index size_cand_b = m_loc[index_cand_b].size();

      return _propose(swap, mtrand, index_cand_a, index_cand_b, size_cand_a, size_cand_b);
    }

    /// Update this after a canonical swap took place
    void update(const CanonicalSwap &swap) {

      // m_loc[index_a] gets moved from asym_a, species_a,  (cand_a) -> asym_a, species_b (cand_a2)
      // m_loc[index_a] gets moved from asym_b, species_b,  (cand_b) -> asym_b, species_a (cand_b2)

      Index tmp;
      Index index_a = swap.mutating_site_index_a;
      Index index_b = swap.mutating_site_index_b;

      int cand_index_a = m_cand.index(swap.cand_a);
      int cand_index_a2 = m_cand.index(swap.cand_a.asym, swap.cand_b.species_index);
      int cand_index_b = m_cand.index(swap.cand_b);
      int cand_index_b2 = m_cand.index(swap.cand_b.asym, swap.cand_a.species_index);

      // store (cand_a)[index_a] -> tmp (mutating_site_a)
      tmp = m_loc[cand_index_a][index_a];

      // (cand_a) remove mutating_site_a
      m_loc[cand_index_a][index_a] = m_loc[cand_index_a].back();
      m_loc[cand_index_a].pop_back();

      // (cand_a2) add remove mutating_site_a
      m_loc[cand_index_a2].push_back(tmp);


      // store (cand_b)[index_b] -> tmp (mutating_site_b)
      tmp = m_loc[cand_index_b][index_b];

      // (cand_b) remove mutating_site_b
      m_loc[cand_index_b][index_b] = m_loc[cand_index_b].back();
      m_loc[cand_index_b].pop_back();

      // (cand_b2) add mutating_site_b
      m_loc[cand_index_ba2].push_back(tmp);

    }

    void prepare(const std::vector<GrandCanonicalSwap> &possible) {
      m_grand_canonical_swap = possible;

      // create lookup table of cand_a index -> possible cand_b
      m_gc_cand_b.clear();
      m_gc_cand_b.resize(m_cand.size());
      for(const auto &swap : m_grand_canonical_swap) {
        m_gc_cand_b[m_cand.index(swap.cand_a)].push_back(swap.cand_b);
      }
    }

    /// Grand canonical propose
    GrandCanonicalSwap &propose(MTRand &mtrand) {
      Index rand = mtrand.randInt(m_size - 1);

      for(int cand_index = 0; cand_index < m_loc.size(); ++cand_index) {
        Index cand_size = m_loc[cand_index].size();
        if(rand < cand_size) {
          m_grand_canonical_swap[cand_index].mutating_site_index = rand;
          m_grand_canonical_swap[cand_index].mutating_site = mutating_site(cand_index, rand);
          int rand2 = mtrand.randInt(m_gc_cand_b[cand_index].size() - 1);
          m_grand_canonical_swap[cand_index].cand_b = m_gc_cand_b[cand_index][rand];
        }
        else {
          rand -= cand_size;
        }
      }

      throw std::runtime_error("Propose grand canonical swap error");
    }

    /// Given a GrandCanonicalSwap, pick the mutating_site
    GrandCanonicalSwap &propose(GrandCanonicalSwap &swap, MTRand &mtrand) const {
      Index index_cand_a = m_cand.index(swap.cand_a);
      Index size_cand_a = m_loc[index_cand_a].size();
      return _propose(swap, mtrand, index_cand_a, size_cand_a);
    }

    /// Update this after a grand canonical swap took place
    void update(const GrandCanonicalSwap &swap) {

      Index index = swap.mutating_site_index;
      int cand_index_a = m_cand.index(swap.cand_a);
      int cand_index_b = m_cand.index(swap.cand_b);

      // the mutating_site
      Index tmp = m_loc[cand_index_a][index];

      // (cand_a) remove mutating_site
      m_loc[cand_index_a][index] = m_loc[cand_index].back();
      m_loc[cand_index].pop_back();

      // (cand_b2) add mutating_site
      m_loc[cand_index_b].push_back(tmp_mutating_site);
    }


  private:

    /// Canonical propose
    CanonicalSwap &_propose(CanonicalSwap &swap, MTRand &mtrand, Index cand_a, Index cand_b, Index size_a, Index size_b) const {
      swap.mutating_site_index_a = mtrand.randInt(size_a - 1);
      swap.mutating_site_a = m_loc[cand_a][swap.mutating_site_index_a];
      swap.mutating_site_index_b = mtrand.randInt(size_b - 1);
      swap.mutating_site_b = m_loc[cand_b][swap.mutating_site_index_b];
      return swap;
    }

    /// Given a GrandCanonicalSwap, pick the mutating_site
    GrandCanonicalSwap &_propose(GrandCanonicalSwap &swap, MTRand &mtrand, Index index_cand_a, Index size_cand_a) const {
      swap.mutating_site_index = mtrand.randInt(size_cand_a - 1);
      swap.mutating_site = m_loc[index_cand_a][swap.mutating_site_index];
      return swap;
    }

    /// candidate_index -> std::vector<mutating_site>
    std::vector<std::vector<Index> > m_loc;

    /// total number of variable sites
    Index m_size;

    /// same as m_size;
    double m_total;
  };


  /// \brief Store swap type, mutating sites, and info for keeping OccLocation up-to-date
  class MCOccSwap : public Comparisons<MCOccSwap> {

  public:

    MCOccSwap(const OccCandidate &_cand_a, const OccCandidate &_cand_b) :
      cand_a(_cand_a)
      cand_b(_cand_b) {}

    OccCandidate cand_a;
    OccCandidate cand_b;

    void reverse() {
      using std::swap;
      std::swap(cand_a, cand_b);
    }

    CanonicalSwap &sort() {
      CanonicalSwap B(*this);
      B.reverse();

      if(B._lt(*this)) {
        *this = B;
      }
      return *this;
    }

    CanonicalSwap sorted() const {
      CanonicalSwap res(*this);
      res.sort();
      return res;
    }

    bool operator<(const CanonicalSwap &B) const {
      return this->sorted()._lt(B.sorted());
    }


  private:

    bool _lt(const CanonicalSwap &B) const {
      return this->tuple() < B.tuple();
    }

    typedef std::tuple<OccCandidate, OccCandidate> tuple_type;

    tuple_type tuple() const {
      return std::make_tuple(cand_a, cand_b);
    }

  };


  /// \brief Data structure for storing information regarding a proposed grand canonical Monte Carlo event
  class CanonicalEvent {

  public:

    typedef Index size_type;

    /// \brief Default constructor
    CanonicalEvent() {}

    /// \brief Constructor
    ///
    /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
    /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
    ///
    CanonicalEvent(size_type Nspecies, size_type Ncorr);


    /// \brief Set the change in (extensive) formation energy associated with this event
    void set_dEf(double dE);

    /// \brief Return change in (extensive) formation energy associated with this event
    double dEf() const;


    /// \brief const Access change in number of species per supercell. Zeros, size of CompositionConverter::components().
    const Eigen::VectorXl &dN() const;

    /// \brief Return change in number of species in supercell. Zeros, size of CompositionConverter::components().
    long int dN(size_type species_type_index) const;


    /// \brief Return change in (extensive) potential energy, dEpot = dEf
    double dEpot() const;

    /// \brief Access the changes in (extensive) correlations associated with this event
    Eigen::VectorXd &dCorr();

    /// \brief const Access the changes in (extensive) correlations associated with this event
    const Eigen::VectorXd &dCorr() const;


    /// \brief Access the data describing this event
    MCOccEvent &event();

    /// \brief const Access the data describing this event
    const MCOccEvent &event() const;


  private:

    /// \brief Change in (extensive) correlations due to this event
    Eigen::VectorXd m_dCorr;

    /// \brief Change in (extensive) formation energy due to this event
    double m_dEf;

    /// \brief Change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    double m_dEpot;

    /// \brief Change in number of each species in supercell due to this event.
    ///        Zeros, size of primclex.get_param_comp().get_components()
    Eigen::VectorXl m_dN;

    /// \brief The modifications performed by this event
    MCOccEvent &m_event;

  };


  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
  /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
  ///
  inline CanonicalEvent::CanonicalEvent(size_type Nspecies, size_type Ncorr) :
    m_dN(Eigen::VectorXl::Zero(Nspecies)),
    m_dCorr(Eigen::VectorXd(Ncorr)) { }


  /// \brief Set the change in total (formation) energy associated with this event
  inline void CanonicalEvent::set_dEf(double dEf) {
    m_dEf = dEf;
  }

  /// \brief Return change in total (formation) energy associated with this event
  inline double CanonicalEvent::dEf() const {
    return m_dEf;
  }


  /// \brief const Access change in number of all species (extensive). Order as in CompositionConverter::components().
  inline const Eigen::VectorXl &CanonicalEvent::dN() const {
    return m_dN;
  }

  /// \brief Return change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
  inline long int CanonicalEvent::dN(size_type species_type_index) const {
    return m_dN(species_type_index);
  }


  /// \brief Return change in potential energy: dEpot = dEf
  inline double CanonicalEvent::dEpot() const {
    return m_dEf;
  }

  /// \brief Access the changes in correlations associated with this event
  inline Eigen::VectorXd &CanonicalEvent::dCorr() {
    return m_dCorr;
  }

  /// \brief const Access the changes in correlations associated with this event
  inline const Eigen::VectorXd &CanonicalEvent::dCorr() const {
    return m_dCorr;
  }

  /// \brief Access the swap for this event
  inline CanonicalSwap &CanonicalEvent::canonical_swap() {
    return m_swap;
  }

  /// \brief const Access the swap for this event
  inline const CanonicalSwap &CanonicalEvent::canonical_swap() const {
    return m_swap;
  }

}

#endif
