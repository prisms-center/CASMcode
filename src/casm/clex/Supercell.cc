#include "casm/clex/Supercell.hh"

#include <math.h>
#include <map>
#include <vector>
#include <stdlib.h>

#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"

namespace CASM {

  /*****************************************************************/

  ReturnArray<int> Supercell::max_allowed_occupation() const {
    Array<int> max_allowed;

    // Figures out the maximum number of occupants in each basis site, to initialize counter with
    for(Index i = 0; i < prim().basis.size(); i++) {
      max_allowed.append(Array<int>(volume(), prim().basis[i].site_occupant().size() - 1));
    }
    //std::cout << "max_allowed_occupation is:  " << max_allowed << "\n\n";
    return max_allowed;
  }

  /*****************************************************************/

  const Structure &Supercell::prim() const {
    return m_primclex->prim();
  }

  /// \brief Returns the SuperNeighborList
  const SuperNeighborList &Supercell::nlist() const {

    // if any additions to the prim nlist, must update the super nlist
    if(primclex().nlist().size() != m_nlist_size_at_construction) {
      m_nlist.unique().reset();
    }

    // lazy construction of neighbor list
    if(!m_nlist) {
      m_nlist_size_at_construction = primclex().nlist().size();
      m_nlist = notstd::make_cloneable<SuperNeighborList>(
                  m_prim_grid,
                  primclex().nlist()
                );
    }
    return *m_nlist;
  };

  /*****************************************************************/

  // begin and end iterators for iterating over configurations
  Supercell::config_iterator Supercell::config_begin() {
    return config_iterator(m_primclex, m_id, 0);
  }

  Supercell::config_iterator Supercell::config_end() {
    return ++config_iterator(m_primclex, m_id, m_config_list.size() - 1);
  }

  // begin and end const_iterators for iterating over configurations
  Supercell::config_const_iterator Supercell::config_cbegin() const {
    return config_const_iterator(m_primclex, m_id, 0);
  }

  Supercell::config_const_iterator Supercell::config_cend() const {
    return ++config_const_iterator(m_primclex, m_id, m_config_list.size() - 1);
  }

  /*****************************************************************/

  const SymGroup &Supercell::factor_group() const {
    if(!m_factor_group.size())
      generate_factor_group();
    return m_factor_group;
  }

  /*****************************************************************/

  // permutation_symrep() populates permutation symrep if needed
  const Permutation &Supercell::factor_group_permute(Index i) const {
    return *(permutation_symrep().get_permutation(factor_group()[i]));
  }
  /*****************************************************************/

  // PrimGrid populates translation permutations if needed
  const Permutation &Supercell::translation_permute(Index i) const {
    return m_prim_grid.translation_permutation(i);
  }

  /*****************************************************************/

  // PrimGrid populates translation permutations if needed
  const Array<Permutation> &Supercell::translation_permute() const {
    return m_prim_grid.translation_permutations();
  }

  /*****************************************************************/
  /* //Example usage case:
   *  Supercell my_supercell;
   *  Configuration my_config(my_supercell, configuration_info);
   *  ConfigDoF my_dof=my_config.configdof();
   *  my_dof.is_canonical(my_supercell.permute_begin(),my_supercell.permute_end());
   */
  Supercell::permute_const_iterator Supercell::permute_begin() const {
    return permute_const_iterator(SymGroupRep::RemoteHandle(this->factor_group(), this->permutation_symrep_ID()),
                                  m_prim_grid,
                                  0, 0); // starting indices
  }

  /*****************************************************************/

  Supercell::permute_const_iterator Supercell::permute_end() const {
    return permute_const_iterator(SymGroupRep::RemoteHandle(factor_group(), permutation_symrep_ID()),
                                  m_prim_grid,
                                  factor_group().size(), 0); // one past final indices
  }


  /*****************************************************************/

  //Printing config_index_to_bijk
  void Supercell::print_bijk(std::ostream &stream) {
    for(Index i = 0; i < num_sites(); i++) {
      stream << uccoord(i);
    }
  }

  //*******************************************************************************
  /**
   *   Checks if the Configuration 'config' is contained in Supercell::m_config_list.
   *     Only checks Configuration::occupation for equivalence.
   *     Does not check for symmetrically equivalent Configurations, so put your
   *     'config' in canonical form first.
   */
  //*******************************************************************************
  bool Supercell::contains_config(const Configuration &config) const {
    Index index;
    return contains_config(config, index);
  };

  //*******************************************************************************
  /**
   *   Checks if the Configuration 'config' is contained in Supercell::m_config_list.
   *     Only checks Configuration::configdof for equivalence.
   *     Does not check for symmetrically equivalent Configurations, so put your
   *     'config' in canonical form first.
   *
   *   If equivalent found, 'index' contains it's index into m_config_list, else
   *     'index' = m_config_list.size().
   */
  //*******************************************************************************
  bool Supercell::contains_config(const Configuration &config, Index &index) const {
    for(Index i = 0; i < m_config_list.size(); i++)
      if(config.configdof() == m_config_list[i].configdof()) {
        index = i;
        return true;
      }

    index = m_config_list.size();
    return false;
  };

  //*******************************************************************************
  /**
   *   Converts 'config' to canonical form, then adds to m_config_list if not already
   *     present. Location in m_config_list is stored in 'index'.
   *     Permutation that resulted in canonical form is stored in 'permute_it'.
   *     Return 'true' if new config, 'false' otherwise.
   *
   *   Might want to rewrite without using new canon_config for memory/speed issues.
   */
  //*******************************************************************************

  bool Supercell::add_config(const Configuration &config) {
    Index index;
    Supercell::permute_const_iterator permute_it;
    return add_config(config, index, permute_it);
  }

  bool Supercell::add_config(const Configuration &config, Index &index, Supercell::permute_const_iterator &permute_it) {
    // 'canon_config' is 'config' permuted to canonical form
    //    std::cout << "get canon_config" << std::endl;
    Configuration canon_config = config.canonical_form(permute_begin(), permute_end(), permute_it);

    // std::cout << "    config: " << config.occupation() << std::endl;
    // std::cout << "     canon: " << canon_config.occupation() << std::endl;

    return add_canon_config(canon_config, index);
  }

  //*******************************************************************************
  /**
   *   Assumes 'canon_config' is in canonical form, adds to m_config_list if not already there.
   *     Location in m_config_list is stored in 'index'.
   */
  //*******************************************************************************
  bool Supercell::add_canon_config(const Configuration &canon_config, Index &index) {

    // Add 'canon_config' to 'm_config_list' if it doesn't already exist
    //   store it's index into 'm_config_list' in 'm_config_list_index'
    //std::cout << "check if canon_config is in m_config_list" << std::endl;
    if(!contains_config(canon_config, index)) {
      //std::cout << "new config" << std::endl;
      m_config_list.push_back(canon_config);
      m_config_list.back().set_id(m_config_list.size() - 1);
      m_config_list.back().set_selected(false);
      return true;
      //std::cout << "    added" << std::endl;
    }
    else {
      m_config_list[index].push_back_source(canon_config.source());
    }
    return false;
  }

  //*******************************************************************************

  void Supercell::read_config_list(const jsonParser &json) {

    // Provide an error check
    if(m_config_list.size() != 0) {
      std::cerr << "Error in Supercell::read_configuration." << std::endl;
      std::cerr << "  config_list().size() != 0, only use this once" << std::endl;
      exit(1);
    }

    if(!json.contains("supercells")) {
      return;
    }

    if(!json["supercells"].contains(name())) {
      return;
    }

    // Read all configurations for this supercell. They should be numbered sequentially, so read until not found.
    Index configid = 0;
    while(true) {
      std::stringstream ss;
      ss << configid;

      if(json["supercells"][name()].contains(ss.str())) {
        m_config_list.push_back(Configuration(json, *this, configid));
      }
      else {
        return;
      }
      configid++;
    }
  }


  //*******************************************************************************

  //Copy constructor is needed for proper initialization of m_prim_grid
  Supercell::Supercell(const Supercell &RHS) :
    m_primclex(RHS.m_primclex),
    m_real_super_lattice(RHS.m_real_super_lattice),
    m_prim_grid((*m_primclex).prim().lattice(), m_real_super_lattice, (*m_primclex).prim().basis.size()),
    m_name(RHS.m_name),
    m_nlist(RHS.m_nlist),
    m_config_list(RHS.m_config_list),
    m_transf_mat(RHS.m_transf_mat) {
  }

  //*******************************************************************************

  Supercell::Supercell(PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &transf_mat_init) :
    m_primclex(_prim),
    m_real_super_lattice((*m_primclex).prim().lattice().lat_column_mat() * transf_mat_init.cast<double>()),
    m_prim_grid((*m_primclex).prim().lattice(), m_real_super_lattice, (*m_primclex).prim().basis.size()),
    m_transf_mat(transf_mat_init) {
    generate_name();
    //    fill_reciprocal_supercell();
  }

  //*******************************************************************************

  Supercell::Supercell(PrimClex *_prim, const Lattice &superlattice) :
    m_primclex(_prim),
    m_real_super_lattice(superlattice),
    m_prim_grid((*m_primclex).prim().lattice(), m_real_super_lattice, (*m_primclex).prim().basis.size()) {

    auto res = is_supercell(superlattice, prim().lattice(), primclex().settings().lin_alg_tol());
    if(!res.first) {
      std::cerr << "Error in Supercell(PrimClex *_prim, const Lattice &superlattice)" << std::endl
                << "  Bad supercell, the transformation matrix is not integer." << std::endl;
      throw std::invalid_argument("Error constructing Supercell: the transformation matrix is not integer");
    }
    m_transf_mat = res.second;

    generate_name();

  }

  //*******************************************************************************
  /**
   * Run through every selected Configuration in *this and call write() on it. This will
   * update all the JSON files and also rewrite POS, DoF etc. Meant for when
   * you calculated some properties (e.g. formation energies or correlations) and
   * want it outputted, but didn't generate any new configurations.
   */

  jsonParser &Supercell::write_config_list(jsonParser &json) {
    for(Index c = 0; c < m_config_list.size(); c++) {
      m_config_list[c].write(json);
    }
    return json;
  }

  //***********************************************************

  void Supercell::generate_factor_group()const {
    m_real_super_lattice.find_invariant_subgroup(prim().factor_group(), m_factor_group);
    m_factor_group.set_lattice(m_real_super_lattice);
    return;
  }

  //***********************************************************

  void Supercell::generate_permutations()const {
    if(!m_perm_symrep_ID.empty()) {
      std::cerr << "WARNING: In Supercell::generate_permutations(), but permutations data already exists.\n"
                << "         It will be overwritten.\n";
    }
    m_perm_symrep_ID = m_prim_grid.make_permutation_representation(factor_group(), prim().basis_permutation_symrep_ID());
    //m_trans_permute = m_prim_grid.make_translation_permutations(basis_size()); <--moved to PrimGrid

    /*
      std::cerr << "For SCEL " << " -- " << name() << " Translation Permutations are:\n";
      for(int i = 0; i < m_trans_permute.size(); i++)
      std::cerr << i << ":   " << m_trans_permute[i].perm_array() << "\n";

      std::cerr << "For SCEL " << " -- " << name() << " factor_group Permutations are:\n";
      for(int i = 0; i < m_factor_group.size(); i++){
    std::cerr << "Operation " << i << ":\n";
    m_factor_group[i].print(std::cerr,FRAC);
    std::cerr << '\n';
    std::cerr << i << ":   " << m_factor_group[i].get_permutation_rep(m_perm_symrep_ID)->perm_array() << '\n';

    }
    std:: cerr << "End permutations for SCEL " << name() << '\n';
    */

    return;
  }

  //***********************************************************

  void Supercell::generate_name() {
    m_name = CASM::generate_name(m_transf_mat);
    return;
  }

  //***********************************************************

  fs::path Supercell::path() const {
    return primclex().dir().supercell_dir(m_name);
  }

  /*
   * Run through the configuration list and count how many of them
   * have been selected then return value.
   */

  Index Supercell::amount_selected() const {
    Index amount_selected = 0;
    for(Index c = 0; c < m_config_list.size(); c++) {
      if(m_config_list[c].selected()) {
        amount_selected++;
      }
    }
    return amount_selected;
  }

  //***********************************************************
  /**  Check if a Structure fits in this Supercell
   *  - Checks that 'structure'.lattice is supercell of 'real_super_lattice'
   *  - Does *NOT* check basis sites
   */
  //***********************************************************
  bool Supercell::is_supercell_of(const Structure &structure) const {
    Eigen::Matrix3d mat;
    return is_supercell_of(structure, mat);
  };

  //***********************************************************
  /**  Check if a Structure fits in this Supercell
   *  - Checks that 'structure'.lattice is supercell of 'real_super_lattice'
   *  - Does *NOT* check basis sites
   */
  //***********************************************************
  bool Supercell::is_supercell_of(const Structure &structure, Eigen::Matrix3d &mat) const {
    Structure tstruct = structure;
    SymGroup point_group;
    tstruct.lattice().generate_point_group(point_group);
    return m_real_super_lattice.is_supercell_of(tstruct.lattice(), point_group, mat);
  };

  //***********************************************************
  /**  Generate a Configuration from a Structure
   *  - Generally expected the user will first call
   *      Supercell::is_supercell_of(const Structure &structure, Matrix3<double> multimat)
   *  - tested OK for perfect prim coordinates, not yet tested with relaxed coordinates using 'tol'
   */
  //***********************************************************
  Configuration Supercell::configuration(const BasicStructure<Site> &structure_to_config, double tol) {
    //Because the user is a fool and the supercell may not be a supercell (This still doesn't check the basis!)
    Eigen::Matrix3d transmat;
    if(!structure_to_config.lattice().is_supercell_of(prim().lattice(), prim().factor_group(), transmat)) {
      std::cerr << "ERROR in Supercell::configuration" << std::endl;
      std::cerr << "The provided structure is not a supercell of the PRIM. Tranformation matrix was:" << std::endl;
      std::cerr << transmat << std::endl;
      exit(881);
    }

    std::cerr << "WARNING in Supercell::config(): This routine has not been tested on relaxed structures using 'tol'" << std::endl;
    //std::cout << "begin config()" << std::endl;
    //std::cout << "  mat:\n" << mat << std::endl;

    const Structure &prim = (*m_primclex).prim();

    // create a 'superstruc' that fills '*this'
    BasicStructure<Site> superstruc = structure_to_config.create_superstruc(m_real_super_lattice);

    //std::cout << "superstruc:\n";
    //superstruc.print(std::cout);
    //std::cout << " " << std::endl;

    // Set the occuation state of a Configuration from superstruc
    //   Allow Va on sites where Va are allowed
    //   Do not allow interstitials, print an error message and exit
    Configuration config(*this);

    // Initially set occupation to -1 (for unknown) on every site
    config.set_occupation(Array<int>(num_sites(), -1));

    Index _linear_index, b;
    int val;

    // For each site in superstruc, set occ index
    for(Index i = 0; i < superstruc.basis.size(); i++) {
      //std::cout << "i: " << i << "  basis: " << superstruc.basis[i] << std::endl;
      _linear_index = linear_index(Coordinate(superstruc.basis[i]), tol);
      b = sublat(_linear_index);

      // check that we're not over-writing something already set
      if(config.occ(_linear_index) != -1) {
        std::cerr << "Error in Supercell::config." << std::endl;
        std::cerr << "  Adding a second atom on site: linear index: " << _linear_index << " bijk: " << uccoord(_linear_index) << std::endl;
        exit(1);
      }

      // check that the Molecule in superstruc is allowed on the site in 'prim'
      if(!prim.basis[b].contains(superstruc.basis[i].occ_name(), val)) {
        std::cerr << "Error in Supercell::config." << std::endl;
        std::cerr << "  The molecule: " << superstruc.basis[i].occ_name() << " is not allowed on basis site " << b << " of the Supercell prim." << std::endl;
        exit(1);
      }
      config.set_occ(_linear_index, val);
    }

    // Check that vacant sites are allowed
    for(Index i = 0; i < config.size(); i++) {
      if(config.occ(i) == -1) {
        b = sublat(i);

        if(prim.basis[b].contains("Va", val)) {
          config.set_occ(i, val);
        }
        else {
          std::cerr << "Error in Supercell::config." << std::endl;
          std::cerr << "  Missing atom.  Vacancies are not allowed on the site: " << uccoord(i) << std::endl;
          exit(1);
        }
      }
    }

    return config;

  };

  //***********************************************************
  /**  Returns a Structure equivalent to the Supercell
   *  - basis sites are ordered to agree with Supercell::config_index_to_bijk
   *  - occupation set to prim default, not curr_state
   */
  //***********************************************************
  Structure Supercell::superstructure() const {
    // create a 'superstruc' that fills '*this'
    Structure superstruc = (*m_primclex).prim().create_superstruc(m_real_super_lattice);

    // sort basis sites so that they agree with config_index_to_bijk
    //   This sorting may not be necessary,
    //   but it depends on how we construct the config_index_to_bijk,
    //   so I'll leave it in for now just to be safe
    for(Index i = 0; i < superstruc.basis.size(); i++) {
      superstruc.basis.swap_elem(i, linear_index(superstruc.basis[i]));
    }

    //superstruc.reset();

    //set_site_internals() is better than Structure::reset(), because
    //it doesn't destroy all the info that
    //Structure::create_superstruc makes efficiently
    superstruc.set_site_internals();
    return superstruc;

  }

  //***********************************************************
  /**  Returns a Structure equivalent to the Supercell
   *  - basis sites are ordered to agree with Supercell::config_index_to_bijk
   *  - occupation set to config
   *  - prim set to (*m_primclex).prim
   */
  //***********************************************************
  Structure Supercell::superstructure(const Configuration &config) const {
    // create a 'superstruc' that fills '*this'
    Structure superstruc = superstructure();

    // set basis site occupants
    for(Index i = 0; i < superstruc.basis.size(); i++) {
      superstruc.basis[i].set_occ_value(config.occ(i));
    }

    // setting the occupation changes symmetry properties, so must reset
    superstruc.reset();

    return superstruc;

  }

  /**
   * This is a safer version that takes an Index instead of an actual Configuration.
   * It might be better to have the version that takes a Configuration private,
   * that way you can't pass it anything that's incompatible.
   */

  Structure Supercell::superstructure(Index config_index) const {
    if(config_index >= m_config_list.size()) {
      std::cerr << "ERROR in Supercell::superstructure" << std::endl;
      std::cerr << "Requested superstructure of configuration with index " << config_index << " but there are only " << m_config_list.size() << " configurations" << std::endl;
      exit(185);
    }
    return superstructure(m_config_list[config_index]);
  }

  //***********************************************************
  /**  Returns an Array<int> consistent with
   *     Configuration::occupation that is all vacancies.
   *     A site which can not contain a vacancy is set to -1.
   */
  //***********************************************************
  ReturnArray<int> Supercell::vacant() const {
    Array<int> occupation = Array<int>(num_sites(), -1);
    int b, index;
    for(Index i = 0; i < num_sites(); i++) {
      b = sublat(i);
      if(prim().basis[b].contains("Va", index)) {
        occupation[i] = index;
      }
    }
    return occupation;
  }


}

