#include "casm/clex/Supercell_impl.hh"

//#include <math.h>
#include <vector>
//#include <stdlib.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include "casm/clex/ChemicalReference.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/casm_io/stream_io/container.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/database/Named_impl.hh"
#include "casm/database/ScelDatabase.hh"


namespace CASM {

  template class SupercellCanonicalForm<CRTPBase<Supercell> >;
  template class HasPrimClex<DB::Named<Comparisons<SupercellCanonicalForm<CRTPBase<Supercell> > > > >;

  namespace DB {
    template class DB::Named<Comparisons<SupercellCanonicalForm<CRTPBase<Supercell> > > >;
  }

  bool ConfigMapCompare::operator()(const Configuration *A, const Configuration *B) const {
    return *A < *B;
  }

  //Copy constructor is needed for proper initialization of m_prim_grid
  Supercell::Supercell(const Supercell &RHS) :
    m_primclex(&RHS.primclex()),
    m_lattice(RHS.m_lattice),
    m_prim_grid(prim().lattice(), m_lattice, prim().basis().size()),
    m_nlist(RHS.m_nlist),
    m_transf_mat(RHS.m_transf_mat) {
  }

  Supercell::Supercell(const PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &transf_mat_init) :
    m_primclex(_prim),
    m_lattice(prim().lattice().lat_column_mat() * transf_mat_init.cast<double>()),
    m_prim_grid(prim().lattice(), m_lattice, prim().basis().size()),
    m_transf_mat(transf_mat_init) {
    //    fill_reciprocal_supercell();
  }

  Supercell::Supercell(const PrimClex *_prim, const Lattice &superlattice) :
    m_primclex(_prim),
    m_lattice(superlattice),
    m_prim_grid(prim().lattice(), m_lattice, prim().basis().size()) {

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

  const PrimClex &Supercell::primclex() const {
    return *m_primclex;
  }

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

  /// \brief Return the coordinate corresponding to linear index in the supercell
  ///
  Coordinate Supercell::coord(Index linear_index) const {
    Coordinate tcoord(m_prim_grid.coord(linear_index % volume(), SCEL));
    tcoord.cart() += prim().basis()[linear_index / volume()].cart();
    return tcoord;
    // return uccoord(linear_index).coordinate();
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
    for(Index i = 0; i < prim().basis().size(); i++) {
      std::vector<int> tmp(volume(), prim().basis()[i].site_occupant().size() - 1);
      max_allowed.insert(max_allowed.end(), tmp.begin(), tmp.end());
    }
    //std::cout << "max_allowed_occupation is:  " << max_allowed << "\n\n";
    return max_allowed;
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
      default_err_log() << "ERROR in Supercell::configuration" << std::endl;
      default_err_log() << "The provided structure is not a supercell of the PRIM. Tranformation matrix was:" << std::endl;
      default_err_log() << transmat << std::endl;
      exit(881);
    }

    default_err_log() << "WARNING in Supercell::config(): This routine has not been tested on relaxed structures using 'tol'" << std::endl;

    // create a 'superstruc' that fills '*this'
    BasicStructure<Site> superstruc = structure_to_config.create_superstruc(m_lattice);


    // Set the occuation state of a Configuration from superstruc
    //   Allow Va on sites where Va are allowed
    //   Do not allow interstitials, print an error message and exit
    Configuration config(*this);

    // Initially set occupation to -1 (for unknown) on every site
    config.set_occupation(std::vector<int>(num_sites(), -1));

    Index _linear_index, b;
    int val;

    // For each site in superstruc, set occ index
    for(Index i = 0; i < superstruc.basis().size(); i++) {
      _linear_index = linear_index(Coordinate(superstruc.basis()[i]), tol);
      b = sublat(_linear_index);

      // check that we're not over-writing something already set
      if(config.occ(_linear_index) != -1) {
        default_err_log() << "Error in Supercell::config." << std::endl;
        default_err_log() << "  Adding a second atom on site: linear index: " << _linear_index << " bijk: " << uccoord(_linear_index) << std::endl;
        throw std::runtime_error("Error in Supercell::configuration: multiple molecule map to same site");
      }

      // check that the Molecule in superstruc is allowed on the site in 'prim'
      if(!prim().basis()[b].contains(superstruc.basis()[i].occ_name(), val)) {
        default_err_log() << "Error in Supercell::config." << std::endl;
        default_err_log() << "  The molecule: " << superstruc.basis()[i].occ_name() << " is not allowed on basis site " << b << " of the Supercell prim." << std::endl;
        throw std::runtime_error("Error in Supercell::configuration: molecule site mapping not allowed");
      }
      config.set_occ(_linear_index, val);
    }

    // Check that vacant sites are allowed
    for(Index i = 0; i < config.size(); i++) {
      if(config.occ(i) == -1) {
        b = sublat(i);

        if(prim().basis()[b].contains("Va", val)) {
          config.set_occ(i, val);
        }
        else {
          default_err_log() << "Error in Supercell::config." << std::endl;
          default_err_log() << "  Missing atom.  Vacancies are not allowed on the site: " << uccoord(i) << std::endl;
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
    Structure superstruc = prim().create_superstruc(m_lattice);

    // sort basis sites so that they agree with config_index_to_bijk
    //   This sorting may not be necessary,
    //   but it depends on how we construct the config_index_to_bijk,
    //   so I'll leave it in for now just to be safe
    //for(Index i = 0; i < superstruc.basis().size(); i++) {
    //superstruc.basis.swap_elem(i, linear_index(superstruc.basis()[i]));
    //}

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
  ///  - prim set to prim()
  ///
  Structure Supercell::superstructure(const Configuration &config) const {
    // create a 'superstruc' that fills '*this'
    Structure superstruc = superstructure();

    // set basis site occupants
    for(Index i = 0; i < superstruc.basis().size(); i++) {
      superstruc.set_occ(i, config.occ(i));
    }

    // setting the occupation changes symmetry properties, so must reset
    superstruc.reset();

    return superstruc;

  }

  const PrimGrid &Supercell::prim_grid() const {
    return m_prim_grid;
  }

  ///Return number of primitive cells that fit inside of *this
  Index Supercell::volume() const {
    return m_prim_grid.size();
  };

  Index Supercell::basis_size() const {
    return prim().basis().size();
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
                                  fg_index, trans_index);
  }

  Supercell::permute_const_iterator Supercell::permute_it(Index fg_index, UnitCell trans) const {
    return permute_it(fg_index, prim_grid().find(trans));
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
  ///
  /// Note: does not commit the change in the database
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
      if(prim().basis()[b].contains("Va", index)) {
        occupation[i] = index;
      }
    }
    return occupation;
  }

  bool Supercell::eq_impl(const Supercell &B) const {
    if(this == &B) {
      return true;
    }
    if(&primclex() != &B.primclex()) {
      throw std::runtime_error(
        "Error using Supercell::operator==(const Supercell& B): "
        "Only Supercell with the same PrimClex may be compared this way.");
    }
    return transf_mat() == B.transf_mat();
  }

  void Supercell::_generate_factor_group()const {
    m_factor_group = m_lattice.invariant_subgroup(prim().factor_group());
  }

  void Supercell::_generate_permutations()const {
    if(!m_perm_symrep_ID.empty()) {
      default_err_log() << "WARNING: In Supercell::generate_permutations(), but permutations data already exists.\n"
                        << "         It will be overwritten.\n";
    }
    m_perm_symrep_ID = m_prim_grid.make_permutation_representation(factor_group(), prim().basis_permutation_symrep_ID());
    //m_trans_permute = m_prim_grid.make_translation_permutations(basis_size()); <--moved to PrimGrid

    /*
      default_err_log() << "For SCEL " << " -- " << name() << " Translation Permutations are:\n";
      for(int i = 0; i < m_trans_permute.size(); i++)
      default_err_log() << i << ":   " << m_trans_permute[i].perm_array() << "\n";

      default_err_log() << "For SCEL " << " -- " << name() << " factor_group Permutations are:\n";
      for(int i = 0; i < m_factor_group.size(); i++){
    default_err_log() << "Operation " << i << ":\n";
    m_factor_group[i].print(default_err_log(),FRAC);
    default_err_log() << '\n';
    default_err_log() << i << ":   " << m_factor_group[i].get_permutation_rep(m_perm_symrep_ID)->perm_array() << '\n';

    }
    std:: cerr << "End permutations for SCEL " << name() << '\n';
    */

    return;
  }

  std::ostream &Supercell::write_pos(std::ostream &sout) const {
    sout << lattice().lat_column_mat() << std::endl;
    return sout;
  }


  void Supercell::write_pos() const {
    const auto &dir = primclex().dir();
    try {
      fs::create_directories(dir.configuration_dir(name()));
    }
    catch(const fs::filesystem_error &ex) {
      default_err_log() << "Error in Supercell::write_pos()." << std::endl;
      default_err_log() << ex.what() << std::endl;
    }

    fs::ofstream file(dir.LAT(name()));
    write_pos(file);
    return;
  }


  /// \brief Return supercell name
  ///
  /// For supercells that are equivalent to the canonical supercell:
  /// - EQUIV_SCEL_NAME = `$CANON_SCELNAME` = `SCELV_A_B_C_D_E_F`
  /// - where 'V' is supercell volume (number of unit cells), and
  ///   'A-F' are the six non-zero elements of the hermite normal form of the
  ///   supercell transformation matrix (T00*T11*T22, T00, T11, T22, T12, T02, T01)
  /// - CANON_SCEL is found in the supercell database (or constructed using the HNF
  ///   for the tranformation matrix and then making the lattice canonical)
  /// For supercells that are not equivalent to the canonical supercell:
  /// - NONEQUIV_SCEL_NAME = `$CANON_SCELNAME.$FG_INDEX`
  /// - The CANON_SCEL is constructed,
  ///   then the FG_INDEX-th prim factor_group operation is applied
  ///
  std::string Supercell::generate_name_impl() const {
    if(lattice().is_equivalent(canonical_form().lattice())) {
      return CASM::generate_name(m_transf_mat);
    }
    else {
      return canonical_form().name() + "." + std::to_string(from_canonical().index());
    }
  }

  /// \brief Get canonical supercell from name. If not yet in database, construct and insert.
  ///
  /// Note: does not commit the change in the database
  const Supercell &make_supercell(const PrimClex &primclex, std::string name) {

    // check if scel is in database
    const auto &db = primclex.db<Supercell>();
    auto it = db.find(name);

    // if already in database, return ref
    if(it != db.end()) {
      return *it;
    }

    // else construct transf_mat from name (make sure to remove any empty tokens)
    std::vector<std::string> tmp, tokens;
    boost::split(tmp, name, boost::is_any_of("SCEL_"), boost::token_compress_on);
    std::copy_if(tmp.begin(), tmp.end(), std::back_inserter(tokens),
    [](const std::string & val) {
      return !val.empty();
    });
    if(tokens.size() != 7) {
      std::string format = "SCELV_T00_T11_T22_T12_T02_T01";
      primclex.err_log().error("In make_supercell");
      primclex.err_log() << "expected format: " << format << "\n";
      primclex.err_log() << "name: |" << name << "|" << std::endl;
      primclex.err_log() << "tokens: " << tokens << std::endl;
      for(const auto &val : tokens) {
        std::cout << "|" << val << "|" << std::endl;
      }
      primclex.err_log() << "tokens.size(): " << tokens.size() << std::endl;
      throw std::invalid_argument("Error in make_supercell: supercell name format error");
    }
    Eigen::Matrix3i T;
    auto cast = [](std::string val) {
      return boost::lexical_cast<Index>(val);
    };
    T << cast(tokens[1]), cast(tokens[6]), cast(tokens[5]),
    0, cast(tokens[2]), cast(tokens[4]),
    0, 0, cast(tokens[3]);

    // construct supercell, insert into database, and return result
    Supercell scel(&primclex, T);
    return *(scel.insert().first);
  }

  /// \brief Construct non-canonical supercell from name. Uses equivalent niggli lattice.
  std::shared_ptr<Supercell> make_shared_supercell(const PrimClex &primclex, std::string name) {

    // tokenize name
    std::vector<std::string> tokens;
    boost::split(tokens, name, boost::is_any_of("."), boost::token_compress_on);

    // validate name
    if(tokens.size() != 2) {
      std::string format = "$CANON_SCEL_NAME.$PRIM_FG_OP";
      primclex.err_log().error("In make_shared_supercell");
      primclex.err_log() << "expected format: " << format << "\n";
      primclex.err_log() << "name: " << name << std::endl;
      primclex.err_log() << "tokens: " << tokens << std::endl;
      throw std::invalid_argument("Error in make_shared_supercell: supercell name format error");
    }

    // generate scel lattice, and put in niggli form
    Index fg_op_index = boost::lexical_cast<Index>(tokens[1]);
    Lattice hnf_lat = copy_apply(
                        primclex.prim().factor_group()[fg_op_index],
                        make_supercell(primclex, tokens[0]).lattice());
    Lattice niggli_lat = niggli(hnf_lat, primclex.crystallography_tol());

    // construct Supercell
    return std::make_shared<Supercell>(&primclex, niggli_lat);
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

