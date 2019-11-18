#include "casm/clex/Supercell_impl.hh"

//#include <math.h>
#include <vector>
//#include <stdlib.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include "casm/clex/ChemicalReference.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/crystallography/BasicStructure_impl.hh"
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
    m_sym_info(make_supercell_sym_info(prim(), RHS.lattice())),
    //m_nlist(RHS.m_nlist),
    m_nlist_size_at_construction(-1) {

  }

  Supercell::Supercell(const PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &transf_mat_init) :
    m_primclex(_prim),
    m_sym_info(make_supercell_sym_info(prim(), Lattice(prim().lattice().lat_column_mat() * transf_mat_init.cast<double>(), _prim->crystallography_tol()))),
    m_nlist_size_at_construction(-1) {

  }

  Supercell::Supercell(const PrimClex *_prim, const Lattice &superlattice) :
    m_primclex(_prim),
    m_sym_info(make_supercell_sym_info(prim(), superlattice)),
    m_nlist_size_at_construction(-1) {

    auto res = xtal::is_superlattice(superlattice, prim().lattice(), primclex().settings().crystallography_tol());
    if(!res.first) {
      _prim->err_log() << "Error in Supercell(PrimClex *_prim, const Lattice &superlattice)" << std::endl
                       << "  Bad supercell, the transformation matrix is not integer." << std::endl;
      _prim->err_log() << "superlattice: \n" << superlattice.lat_column_mat() << std::endl;
      _prim->err_log() << "prim lattice: \n" << prim().lattice().lat_column_mat() << std::endl;
      _prim->err_log() << "lin_alg_tol: " << primclex().settings().lin_alg_tol() << std::endl;
      _prim->err_log() << "transformation matrix: \n" << prim().lattice().lat_column_mat().inverse() * superlattice.lat_column_mat() << std::endl;
      throw std::invalid_argument("Error constructing Supercell: the transformation matrix is not integer");
    }
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
    return prim_grid().sublat(linear_index);
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
    return linear_index(UnitCellCoord::from_coordinate(this->primclex().prim(), coord, tol));
  }

  /// \brief Return the linear index corresponding to integral coordinates
  ///
  /// Linear indices are grouped by sublattice, then ordered as determined by
  /// PrimGrid. This function is equivalent to:
  /// \code
  /// bijk[0] * volume() + prim_grid().find(bijk.unitcell());
  /// \endcode
  Index Supercell::linear_index(const UnitCellCoord &bijk) const {
    return bijk[0] * volume() + prim_grid().find(bijk.unitcell());
  }

  /// \brief Return the coordinate corresponding to linear index in the supercell
  ///
  Coordinate Supercell::coord(Index linear_index) const {
    Coordinate tcoord(prim_grid().scel_coord(linear_index % volume()));
    tcoord.cart() += prim().basis()[sublat(linear_index)].cart();
    return tcoord;
    // return uccoord(linear_index).coordinate();
  }

  /// \brief Return the integral coordinates corresponding to a linear index
  ///
  /// Linear indices are grouped by sublattice, then ordered as determined by
  /// PrimGrid. This function is equivalent to:
  /// \code
  /// UnitCellCoord(prim(), sublat(linear_index), prim_grid().unitcell(linear_index % volume()))
  /// \endcode
  UnitCellCoord Supercell::uccoord(Index linear_index) const {
    return UnitCellCoord(sublat(linear_index), prim_grid().unitcell(linear_index % volume()));
  }

  /// \brief returns Supercell-compatible configdof with zeroed DoF values and user-specified tolerance
  ConfigDoF Supercell::zero_configdof(double tol) const {
    return ConfigDoF(basis_size(),
                     volume(),
                     global_dof_info(prim()),
                     local_dof_info(prim()),
                     occ_symrep_IDs(prim()),
                     tol);
  }

  std::vector<int> Supercell::max_allowed_occupation() const {
    std::vector<int> max_allowed;

    // Figures out the maximum number of occupants in each basis site, to initialize counter with
    for(Index i = 0; i < prim().basis().size(); i++) {
      std::vector<int> tmp(volume(), prim().basis()[i].occupant_dof().size() - 1);
      max_allowed.insert(max_allowed.end(), tmp.begin(), tmp.end());
    }
    //std::cout << "max_allowed_occupation is:  " << max_allowed << "\n\n";
    return max_allowed;
  }


  const PrimGrid &Supercell::prim_grid() const {
    return sym_info().prim_grid();
  }

  ///Return number of primitive cells that fit inside of *this
  Index Supercell::volume() const {
    return prim_grid().size();
  }

  Index Supercell::basis_size() const {
    return prim().basis().size();
  }

  Index Supercell::num_sites() const {
    return volume() * basis_size();
  }

  Eigen::Matrix3i Supercell::transf_mat() const {
    return iround(prim_grid().trans_mat());
    //return CASM::transf_mat(primclex().prim().lattice(), lattice(), primclex().crystallography_tol());
  }

  const Lattice &Supercell::lattice() const {
    return prim_grid().scel_lattice();
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
                  prim_grid(),
                  primclex().nlist()
                );
    }
    return *m_nlist;
  }

  // Factor group of this supercell
  const SymGroup &Supercell::factor_group() const {
    return sym_info().factor_group();
  }

  // SymInfo object of this supercell
  const SupercellSymInfo &Supercell::sym_info() const {
    return m_sym_info;
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
             xtal::canonical::equivalent(
               lattice(),
               prim().point_group(),
               crystallography_tol()));
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



  std::string pos_string(Supercell const &_scel) {
    std::stringstream ss;
    ss << _scel.lattice().lat_column_mat() << std::endl;
    return ss.str();
  }


  void write_pos(Supercell const &_scel) {
    const auto &dir = _scel.primclex().dir();
    try {
      fs::create_directories(dir.configuration_dir(_scel.name()));
    }
    catch(const fs::filesystem_error &ex) {
      default_err_log() << "Error in Supercell::write_pos()." << std::endl;
      default_err_log() << ex.what() << std::endl;
    }

    fs::ofstream file(dir.LAT(_scel.name()));
    file << pos_string(_scel);
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
    return scelname(prim(), lattice());
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

    std::vector<std::string> tmp, tokens;
    try {
      // else construct transf_mat from name (make sure to remove any empty tokens)
      boost::split(tmp, name, boost::is_any_of("SCEL_"), boost::token_compress_on);
      std::copy_if(tmp.begin(), tmp.end(), std::back_inserter(tokens),
      [](const std::string & val) {
        return !val.empty();
      });
      if(tokens.size() != 7) {
        throw std::invalid_argument("Error in make_supercell: supercell name format error");
      }
    }
    catch(std::exception &e) {
      std::string format = "SCELV_T00_T11_T22_T12_T02_T01";
      primclex.err_log().error("In make_supercell");
      primclex.err_log() << "expected format: " << format << "\n";
      primclex.err_log() << "name: |" << name << "|" << std::endl;
      primclex.err_log() << "tokens: " << tokens << std::endl;
      primclex.err_log() << "tokens.size(): " << tokens.size() << std::endl;
      primclex.err_log() << e.what() << std::endl;
      throw e;
    }

    Eigen::Matrix3i T;
    try {
      auto cast = [](std::string val) {
        return boost::lexical_cast<Index>(val);
      };
      T << cast(tokens[1]), cast(tokens[6]), cast(tokens[5]),
      0, cast(tokens[2]), cast(tokens[4]),
      0, 0, cast(tokens[3]);
    }
    catch(std::exception &e) {
      primclex.err_log().error("In make_supercell");
      primclex.err_log() << "Could not construct transformation matrix from supercell name" << std::endl;
      primclex.err_log() << "  name: " << name << std::endl;
      primclex.err_log() << "  tokens: " << tokens << std::endl;
      primclex.err_log() << e.what() << std::endl;
      throw e;
    }

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

  SupercellSymInfo make_supercell_sym_info(Structure const &_prim, Lattice const &_slat) {
    std::map<DoFKey, SymGroupRepID> global_dof_symrep_IDs;
    for(auto const &key : global_dof_types(_prim))
      global_dof_symrep_IDs.emplace(std::make_pair(key, _prim.global_dof(key).symrep_ID()));

    std::map<DoFKey, std::vector<SymGroupRepID> > local_dof_symrep_IDs;
    for(auto const &key : continuous_local_dof_types(_prim)) {
      std::vector<SymGroupRepID> treps(_prim.basis().size());
      for(Index b = 0; b < _prim.basis().size(); ++b) {
        if(_prim.basis()[b].has_dof(key))
          treps[b] = _prim.basis()[b].dof(key).symrep_ID();
      }
      local_dof_symrep_IDs.emplace(std::make_pair(key, std::move(treps)));
    }


    return SupercellSymInfo(_prim.lattice(),
                            _slat,
                            _prim.basis().size(),
                            _prim.factor_group(),
                            _prim.basis_permutation_symrep_ID(),
                            global_dof_symrep_IDs,
                            occ_symrep_IDs(_prim),
                            local_dof_symrep_IDs);


  }



  Supercell &apply(const SymOp &op, Supercell &scel) {
    return scel = copy_apply(op, scel);
  }

  Supercell copy_apply(const SymOp &op, const Supercell &scel) {
    return Supercell(&scel.primclex(), copy_apply(op, scel.lattice()));
  }

  Eigen::Matrix3i transf_mat(const Lattice &prim_lat, const Lattice &super_lat, double tol) {
    auto res = xtal::is_superlattice(super_lat, prim_lat, tol);
    if(!res.first) {
      std::stringstream err_msg;
      err_msg <<
              "Error finding supercell transformation matrix:\n" <<
              "  Bad supercell, the transformation matrix is not integer.\n\n" <<
              "superlattice: \n" << super_lat.lat_column_mat() << "\n" <<
              "prim lattice: \n" << prim_lat.lat_column_mat() << "\n" <<
              "tolerance: " << tol << "\n" <<
              "transformation matrix: \n" << prim_lat.lat_column_mat().inverse() * super_lat.lat_column_mat() << "\n";
      throw std::invalid_argument(err_msg.str());
    }
    return iround(res.second);
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

  std::string scelname(const Structure &prim, const Lattice &superlat) {
    const SymGroup &pg = prim.point_group();
    Lattice canon_lat = xtal::canonical::equivalent(superlat, pg);
    std::string result = CASM::generate_name(transf_mat(prim.lattice(), canon_lat, prim.lattice().tol()));
    if(!xtal::is_equivalent(superlat, canon_lat)) {
      auto to_canonical_ix = xtal::canonical::operation_index(superlat, pg);
      result += ("." + std::to_string(pg[to_canonical_ix].inverse().index()));
    }
    return result;
  }

  std::string canonical_scelname(const Structure &prim, const Lattice &superlat) {
    const SymGroup &pg = prim.point_group();
    return CASM::generate_name(transf_mat(prim.lattice(), xtal::canonical::equivalent(superlat, pg), prim.lattice().tol()));
  }


}
