#include "casm/clex/Supercell.hh"

//#include <math.h>
#include <vector>
//#include <stdlib.h>

#include "casm/app/ProjectSettings.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/database/ScelDatabase.hh"

namespace CASM {

  bool ConfigMapCompare::operator()(const Configuration *A, const Configuration *B) const {
    return *A < *B;
  }

  //Copy constructor is needed for proper initialization of m_prim_grid
  Supercell::Supercell(const Supercell &RHS) :
    Named(RHS.primclex()),
    m_lattice(RHS.m_lattice),
    m_prim_grid(primclex().prim().lattice(), m_lattice, primclex().prim().basis.size()),
    m_nlist(RHS.m_nlist),
    m_canonical(nullptr),
    m_transf_mat(RHS.m_transf_mat) {
  }

  Supercell::Supercell(const PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &transf_mat_init) :
    Named(*_prim),
    m_lattice(primclex().prim().lattice().lat_column_mat() * transf_mat_init.cast<double>()),
    m_prim_grid(primclex().prim().lattice(), m_lattice, primclex().prim().basis.size()),
    m_canonical(nullptr),
    m_transf_mat(transf_mat_init) {
    //    fill_reciprocal_supercell();
  }

  Supercell::Supercell(const PrimClex *_prim, const Lattice &superlattice) :
    Named(*_prim),
    m_lattice(superlattice),
    m_canonical(nullptr),
    m_prim_grid(primclex().prim().lattice(), m_lattice, primclex().prim().basis.size()) {

    auto res = is_supercell(superlattice, prim().lattice(), primclex().settings().crystallography_tol());
    if(!res.first) {
      _prim->err_log() << "Error in Supercell(PrimClex *_prim, const Lattice &superlattice)" << std::endl
                       << "  Bad supercell, the transformation matrix is not integer." << std::endl;
      _prim->err_log() << "superlattice: \n" << superlattice.lat_column_mat() << std::endl;
      _prim->err_log() << "prim lattice: \n" << prim().lattice().lat_column_mat() << std::endl;
      _prim->err_log() << "lin_alg_tol: " << primclex().settings().lin_alg_tol() << std::endl;
      _prim->err_log() << "transformation matrix: \n" << prim().lattice().lat_column_mat().inverse() * superlattice.lat_column_mat() << std::endl;
      throw std::invalid_argument("Error constructing Supercell: the transformation matrix is not integer");
    }
    m_transf_mat = res.second;
  }

  Supercell::~Supercell() {}

  /// \brief Return the sublattice index for a linear index
  ///
  /// Linear indices are grouped by sublattice, then ordered as determined by
  /// PrimGrid. This function is equivalent to:
  /// \code
  /// linear_index / volume();
  /// \endcode
  Index Supercell::sublat(Index linear_index) const {
    return linear_index / volume();
  }

  /// \brief Given a Coordinate and tolerance, return linear index into Configuration
  ///
  ///   This may be slow, first converts Coordinate -> UnitCellCoord,
  ///   then gets linear_index from UnitCellCoord
  ///
  /// Implementation:
  /// \code
  /// Coordinate tcoord(coord);
  /// tcoord.within();
  /// return linear_index(UnitCellCoord(prim(), coord, tol));
  /// \endcode
  Index Supercell::linear_index(const Coordinate &coord, double tol) const {
    Coordinate tcoord(coord);
    tcoord.within();
    return linear_index(UnitCellCoord(prim(), coord, tol));
  };

  /// \brief Return the linear index corresponding to integral coordinates
  ///
  /// Linear indices are grouped by sublattice, then ordered as determined by
  /// PrimGrid. This function is equivalent to:
  /// \code
  /// bijk[0] * volume() + m_prim_grid.find(bijk.unitcell());
  /// \endcode
  Index Supercell::linear_index(const UnitCellCoord &bijk) const {
    return bijk[0] * volume() + m_prim_grid.find(bijk.unitcell());
  }

  /// \brief Return the linear index corresponding to integral coordinates
  ///
  /// Equivalent to:
  /// \code
  /// uccoord(linear_index).coordinate()
  /// \endcode
  Coordinate Supercell::coord(Index linear_index) const {
    return uccoord(linear_index).coordinate();
  }

  /// \brief Return the integral coordinates corresponding to a linear index
  ///
  /// Linear indices are grouped by sublattice, then ordered as determined by
  /// PrimGrid. This function is equivalent to:
  /// \code
  /// UnitCellCoord(prim(), sublat(linear_index), m_prim_grid.unitcell(linear_index % volume()))
  /// \endcode
  UnitCellCoord Supercell::uccoord(Index linear_index) const {
    return UnitCellCoord(prim(), sublat(linear_index), m_prim_grid.unitcell(linear_index % volume()));
  };

  std::vector<int> Supercell::max_allowed_occupation() const {
    std::vector<int> max_allowed;

    // Figures out the maximum number of occupants in each basis site, to initialize counter with
    for(Index i = 0; i < prim().basis.size(); i++) {
      std::vector<int> tmp(volume(), prim().basis[i].site_occupant().size() - 1);
      max_allowed.insert(max_allowed.end(), tmp.begin(), tmp.end());
    }
    //std::cout << "max_allowed_occupation is:  " << max_allowed << "\n\n";
    return max_allowed;
  }

  bool Supercell::is_canonical() const {
    return lattice().is_canonical(
             prim().point_group(),
             primclex().crystallography_tol());
  }

  SymOp Supercell::to_canonical() const {
    return lattice().to_canonical(
             prim().point_group(),
             primclex().crystallography_tol());
  }

  SymOp Supercell::from_canonical() const {
    return lattice().from_canonical(
             prim().point_group(),
             primclex().crystallography_tol());
  }

  const Supercell &Supercell::canonical_form() const {
    if(!m_canonical) {
      m_canonical = &*insert().first;
    }
    return *m_canonical;
  }

  //***********************************************************
  /**  Generate a Configuration from a Structure
   *  - Generally expected the user will first call
   *      Supercell::is_supercell_of(const Structure &structure, Matrix3<double> multimat)
   *  - tested OK for perfect prim coordinates, not yet tested with relaxed coordinates using 'tol'
   */
  //***********************************************************
  Configuration Supercell::configuration(const BasicStructure<Site> &structure_to_config, double tol) const {
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

    const Structure &prim = primclex().prim();

    // create a 'superstruc' that fills '*this'
    BasicStructure<Site> superstruc = structure_to_config.create_superstruc(m_lattice);

    //std::cout << "superstruc:\n";
    //superstruc.print(std::cout);
    //std::cout << " " << std::endl;

    // Set the occuation state of a Configuration from superstruc
    //   Allow Va on sites where Va are allowed
    //   Do not allow interstitials, print an error message and exit
    Configuration config(*this);

    // Initially set occupation to -1 (for unknown) on every site
    config.set_occupation(std::vector<int>(num_sites(), -1));

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

  ///  Returns a Structure equivalent to the Supercell
  ///  - basis sites are ordered to agree with Supercell::config_index_to_bijk
  ///  - occupation set to prim default, not curr_state
  ///
  Structure Supercell::superstructure() const {
    // create a 'superstruc' that fills '*this'
    Structure superstruc = primclex().prim().create_superstruc(m_lattice);

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

  ///  Returns a Structure equivalent to the Supercell
  ///  - basis sites are ordered to agree with Supercell::config_index_to_bijk
  ///  - occupation set to config
  ///  - prim set to (*m_primclex).prim
  ///
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

  /// \brief Get the PrimClex crystallography_tol
  double Supercell::crystallography_tol() const {
    return primclex().crystallography_tol();
  }

  const PrimGrid &Supercell::prim_grid() const {
    return m_prim_grid;
  }

  const Structure &Supercell::prim() const {
    return primclex().prim();
  }

  ///Return number of primitive cells that fit inside of *this
  Index Supercell::volume() const {
    return m_prim_grid.size();
  };

  Index Supercell::basis_size() const {
    return prim().basis.size();
  }

  Index Supercell::num_sites() const {
    return volume() * basis_size();
  };

  // the permutation_symrep is the SymGroupRep of prim().factor_group() that describes how
  // operations of m_factor_group permute sites of the Supercell.
  // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
  //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
  //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
  //       SymGroupRep::get_representation(m_factor_group[i]) or SymGroupRep::get_permutation(m_factor_group[i]),
  //       so that you don't encounter the gaps (i.e., the representation can be indexed using the
  //       SymOps of m_factor_group
  SymGroupRepID Supercell::permutation_symrep_ID() const {
    if(m_perm_symrep_ID.empty()) {
      _generate_permutations();
    }
    return m_perm_symrep_ID;
  }

  SymGroupRep const &Supercell::permutation_symrep() const {
    return prim().factor_group().representation(permutation_symrep_ID());
  }

  const Eigen::Matrix3i &Supercell::transf_mat() const {
    return m_transf_mat;
  };

  const Lattice &Supercell::lattice() const {
    return m_lattice;
  };

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

  const SymGroup &Supercell::factor_group() const {
    if(!m_factor_group.size()) {
      _generate_factor_group();
    }
    return m_factor_group;
  }

  // permutation_symrep() populates permutation symrep if needed
  const Permutation &Supercell::factor_group_permute(Index i) const {
    return *(permutation_symrep().get_permutation(factor_group()[i]));
  }

  // PrimGrid populates translation permutations if needed
  const Permutation &Supercell::translation_permute(Index i) const {
    return m_prim_grid.translation_permutation(i);
  }

  // PrimGrid populates translation permutations if needed
  const std::vector<Permutation> &Supercell::translation_permute() const {
    return m_prim_grid.translation_permutations();
  }

  /// \brief Begin iterator over translation permutations
  Supercell::permute_const_iterator Supercell::translate_begin() const {
    return permute_begin();
  }

  /// \brief End iterator over translation permutations
  Supercell::permute_const_iterator Supercell::translate_end() const {
    return permute_begin().begin_next_fg_op();
  }

  /// Example usage case:
  ///  Supercell my_supercell;
  ///  Configuration my_config(my_supercell, configuration_info);
  ///  ConfigDoF my_dof=my_config.configdof();
  ///  my_dof.is_canonical(my_supercell.permute_begin(),my_supercell.permute_end());
  ///
  Supercell::permute_const_iterator Supercell::permute_begin() const {
    return permute_it(0, 0); // starting indices
  }

  Supercell::permute_const_iterator Supercell::permute_end() const {
    return permute_it(factor_group().size(), 0); // one past final indices
  }

  Supercell::permute_const_iterator Supercell::permute_it(Index fg_index, Index trans_index) const {
    return permute_const_iterator(SymGroupRep::RemoteHandle(factor_group(), permutation_symrep_ID()),
                                  m_prim_grid,
                                  fg_index, trans_index); // one past final indices
  }

  bool Supercell::operator<(const Supercell &B) const {
    if(&primclex() != &B.primclex()) {
      throw std::runtime_error(
        "Error using Supercell::operator<(const Supercell& B): "
        "Only Supercell with the same PrimClex may be compared this way.");
    }
    if(volume() != B.volume()) {
      return volume() < B.volume();
    }
    return lattice() < B.lattice();
  }

  /// \brief Insert the canonical form of this into the database
  std::pair<DB::DatabaseIterator<Supercell>, bool> Supercell::insert() const {
    return primclex().db<Supercell>().emplace(
             & primclex(),
             canonical_equivalent_lattice(
               lattice(),
               prim().point_group(),
               crystallography_tol()));
  }

  ///  Check if a Structure fits in this Supercell
  ///  - Checks that 'structure'.lattice is supercell of 'lattice'
  ///  - Does *NOT* check basis sites
  ///
  bool Supercell::is_supercell_of(const Structure &structure) const {
    Eigen::Matrix3d mat;
    return is_supercell_of(structure, mat);
  };

  ///  Check if a Structure fits in this Supercell
  ///  - Checks that 'structure'.lattice is supercell of 'lattice'
  ///  - Does *NOT* check basis sites
  ///
  bool Supercell::is_supercell_of(const Structure &structure, Eigen::Matrix3d &mat) const {
    Structure tstruct = structure;
    SymGroup point_group;
    tstruct.lattice().generate_point_group(point_group);
    return m_lattice.is_supercell_of(tstruct.lattice(), point_group, mat);
  };

  ///  Returns an std::vector<int> consistent with
  ///    Configuration::occupation that is all vacancies.
  ///    A site which can not contain a vacancy is set to -1.
  ///
  std::vector<int> Supercell::vacant() const {
    std::vector<int> occupation(num_sites(), -1);
    int b, index;
    for(Index i = 0; i < num_sites(); i++) {
      b = sublat(i);
      if(prim().basis[b].contains("Va", index)) {
        occupation[i] = index;
      }
    }
    return occupation;
  }

  bool Supercell::_eq(const Supercell &B) const {
    if(&primclex() != &B.primclex()) {
      throw std::runtime_error(
        "Error using Supercell::operator==(const Supercell& B): "
        "Only Supercell with the same PrimClex may be compared this way.");
    }
    return transf_mat() == B.transf_mat();
  }

  void Supercell::_generate_factor_group()const {
    m_lattice.find_invariant_subgroup(prim().factor_group(), m_factor_group);
    m_factor_group.set_lattice(m_lattice);
    return;
  }

  void Supercell::_generate_permutations()const {
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

  /// \brief Return supercell name
  ///
  /// - If lattice is the canonical equivalent, then return 'SCELV_A_B_C_D_E_F'
  /// - Else, return 'SCELV_A_B_C_D_E_F.$FG_INDEX', where $FG_INDEX is the index of the first
  ///   symmetry operation in the primitive structure's factor group such that the lattice
  ///   is equivalent to `apply(fg_op, canonical equivalent)`
  std::string Supercell::_generate_name() const {
    return CASM::generate_name(m_transf_mat);
  }

  Supercell &apply(const SymOp &op, Supercell &scel) {
    return scel = copy_apply(op, scel);
  }

  Supercell copy_apply(const SymOp &op, const Supercell &scel) {
    return Supercell(&scel.primclex(), copy_apply(op, scel.lattice()));
  }

  std::string generate_name(const Eigen::Matrix3i &transf_mat) {
    std::string name_str;

    Eigen::Matrix3i H = hermite_normal_form(transf_mat).first;
    name_str = "SCEL";
    std::stringstream tname;
    //Consider using a for loop with HermiteCounter_impl::_canonical_unroll here
    tname << H(0, 0)*H(1, 1)*H(2, 2)
          << "_" << H(0, 0) << "_" << H(1, 1) << "_" << H(2, 2)
          << "_" << H(1, 2) << "_" << H(0, 2) << "_" << H(0, 1);
    name_str.append(tname.str());

    return name_str;
  }

}

