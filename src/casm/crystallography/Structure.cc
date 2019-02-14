#include "casm/crystallography/Structure.hh"

#include <sstream>
#include <string>
#include <exception>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/misc/algorithm.hh"
#include "casm/container/algorithm.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/BasicStructure_impl.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/DoFIsEquivalent_impl.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/casm_io/Log.hh"


namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Structure::Structure(const fs::path &filepath) : BasicStructure<Site>() {
    if(!fs::exists(filepath)) {
      default_err_log() << "Error in Structure::Structure(const fs::path &filepath)." << std::endl;
      default_err_log() << "  File does not exist at: " << filepath << std::endl;
      exit(1);
    }
    fs::ifstream infile(filepath);

    read(infile);
  }

  //***********************************************************

  Structure::Structure(const Structure &RHS) :
    BasicStructure<Site>(RHS) {

    copy_attributes_from(RHS);

  };

  //***********************************************************

  Structure::~Structure() {}

  //***********************************************************

  Structure &Structure::operator=(const Structure &RHS) {
    BasicStructure<Site>::operator=(RHS);

    //Following gets done by base class
    //lattice = RHS.lattice;
    //basis = RHS.basis;
    //title = RHS.title;

    /*
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].set_lattice(lattice);
    }
    */

    copy_attributes_from(RHS);

    return *this;
  }


  //***********************************************************

  void Structure::copy_attributes_from(const Structure &RHS) {

    BasicStructure<Site>::copy_attributes_from(RHS);

    m_basis_perm_rep_ID = RHS.m_basis_perm_rep_ID; //this *should* work
    //assert(0);
    m_factor_group = RHS.m_factor_group;
    m_factor_group.set_lattice(lattice());
  }

  //***********************************************************


  void Structure::generate_factor_group_slow() const {
    m_factor_group.clear();
    m_factor_group.set_lattice(lattice());
    BasicStructure<Site>::generate_factor_group_slow(m_factor_group);
    return;
  }

  //************************************************************
  void Structure::generate_factor_group() const {
    m_factor_group.clear();
    m_factor_group.set_lattice(lattice());
    BasicStructure<Site>::generate_factor_group(m_factor_group);
    _generate_basis_symreps();
    _generate_global_symreps();
    return;
  }

  //************************************************************
  const MasterSymGroup &Structure::factor_group() const {
    if(!m_factor_group.size())
      generate_factor_group();
    return m_factor_group;
  }

  //************************************************************
  const SymGroup &Structure::point_group() const {
    return factor_group().point_group();
  }

  //***********************************************************

  SymGroupRep const *Structure::basis_permutation_symrep() const {
    return &(factor_group().representation(basis_permutation_symrep_ID()));
  }

  //***********************************************************

  SymGroupRepID Structure::basis_permutation_symrep_ID() const {
    if(m_basis_perm_rep_ID.empty())
      _generate_basis_symreps();

    return m_basis_perm_rep_ID;
  }

  //************************************************************
  std::vector<AtomSpecies> Structure::struc_species() const {
    return CASM::struc_species(*this);
  }
  //************************************************************
  std::vector<Molecule> Structure::struc_molecule() const {
    return CASM::struc_molecule(*this);
  }
  //************************************************************
  std::vector<std::string> Structure::struc_species_name() const {
    return CASM::struc_species_name(*this);
  }

  //************************************************************
  std::vector<std::string> Structure::struc_molecule_name() const {
    return CASM::struc_molecule_name(*this);
  }
  //************************************************************
  Eigen::VectorXi Structure::num_each_species() const {
    return CASM::num_each_species(*this);
  }
  //************************************************************
  Eigen::VectorXi Structure::num_each_molecule() const {
    return CASM::num_each_molecule(*this);
  }
  //************************************************************

  void Structure::fg_converge(double small_tol, double large_tol, double increment) {
    BasicStructure<Site>::fg_converge(m_factor_group, small_tol, large_tol, increment);
    return;
  }

  //************************************************************
  void Structure::fg_converge(double large_tol) {
    BasicStructure<Site>::fg_converge(m_factor_group, lattice().tol(), large_tol, (large_tol - lattice().tol()) / 10.0);
    return;
  }

  //************************************************************
  /*void Structure::print_factor_group(std::ostream &stream) const {
    stream << "Factor Group of " << title << ", containing "
           << factor_group().size() << " symmetry operations:\n";

    for(Index i = 0; i < m_factor_group.size(); i++) {
      m_factor_group[i].print(stream);
    }

    return;
  }
  */
  //***********************************************************
  /**
   * It is NOT wise to use this function unless you have already
   * initialized a superstructure with lattice vectors.
   *
   * It is more wise to use the two methods that call this method:
   * Either the overloaded * operator which does:
   *  SCEL_Lattice * Prim_Structrue = New_Superstructure
   *       --- or ---
   *  New_Superstructure=Prim_Structure.create_superstruc(SCEL_Lattice);
   *
   *  Both of these will return NEW superstructures.
   */
  //***********************************************************

  void Structure::fill_supercell(const Structure &prim) {
    Index i, j;

    SymGroup latvec_pg;
    m_lattice.generate_point_group(latvec_pg);

    m_SD_flag = prim.m_SD_flag;
    PrimGrid prim_grid(prim.lattice(), lattice());

    m_basis.clear();
    Coordinate tcoord(lattice());

    //loop over basis sites of prim
    for(j = 0; j < prim.basis().size(); j++) {

      //loop over prim_grid points
      for(i = 0; i < prim_grid.size(); i++) {

        //push back translated basis site of prim onto superstructure basis
        push_back(prim.basis()[j] + prim_grid.coord(i, PRIM));

        m_basis.back().within();
        for(Index k = 0; k < basis().size() - 1; k++) {
          if(basis()[k].compare(basis().back())) {
            m_basis.pop_back();
            break;
          }
        }
      }
    }
    //trans_and_expand primitive factor_group
    for(i = 0; i < prim.factor_group().size(); i++) {
      if(latvec_pg.find_no_trans(prim.factor_group()[i]) == latvec_pg.size()) {
        continue;
      }
      else {
        for(Index j = 0; j < prim_grid.size(); j++) {
          m_factor_group.push_back(within_cell(SymOp::translation(prim_grid.coord(j, SCEL).const_cart())*prim.factor_group()[i],
                                               lattice(),
                                               PERIODIC));
        }
      }
    }
    if(m_factor_group.size() > 200) {// how big is too big? this is at least big enough for FCC conventional cell
#ifndef NDEBUG
      default_err_log() << "WARNING: You have a very large factor group of a non-primitive structure. Certain symmetry features will be unavailable.\n";
#endif
      m_factor_group.invalidate_multi_tables();
    }
    update();

    return;
  }

  //***********************************************************
  /**
   * Operates on the primitive structure and takes as an argument
   * the supercell lattice.  It then returns a new superstructure.
   *
   * This is similar to the Lattice*Primitive routine which returns a
   * new superstructure.  Unlike the fill_supercell routine which takes
   * the primitive structure, this WILL fill the sites.
   */
  //***********************************************************

  Structure Structure::create_superstruc(const Lattice &scel_lat) const {
    Structure tsuper(scel_lat);
    tsuper.fill_supercell(*this);
    return tsuper;
  }

  //***********************************************************
  /*shorttag allows the user to decide whether or not to print the
   * entire symmetry matrix.  If shorttag=0, then only the name of the
   * operation and the eigenvector will be printed.  If it equals anything else,
   * then the symmetry operation matrix will be printed out.
   */
  /*
    void Structure::print_site_symmetry(std::ostream &stream, COORD_TYPE mode, int shorttag = 1) {
      GenericOrbitBranch<SiteCluster> asym_unit(lattice());
      asym_unit.generate_asymmetric_unit(basis, factor_group());

      stream <<  " Printing symmetry operations that leave each site unchanged:\n \n";
      for(Index i = 0; i < asym_unit.size(); i++) {
        for(Index j = 0; j < asym_unit[i].size(); j++) {
          stream << "Site: ";
          asym_unit.at(i).at(j).at(0).print(stream); //Print the site (not cluster) for the asymmetric unit
          stream << " Total Symmetry Operations: " << asym_unit[i][j].clust_group().size() << "\n";

          for(Index nc = 0; nc < asym_unit[i][j].clust_group().size(); nc++) {
            //print_short is a new fxn I added in SymOp to not print out the whole symmetry matrix
            if(shorttag == 0) {
              stream <<  std::setw(4)  << nc + 1 << ": ";
              //asym_unit[i][j].clust_group()[nc].print_short(stream); // I turned this off because the "print_short" function does not appear to exist...
            }


            else {
              stream.flags(std::ios::left);
              stream << "\n" <<  std::setw(3)  << nc + 1 << ": ";
              stream.unsetf(std::ios::left);
              if(mode == CART)
                asym_unit[i][j].clust_group()[nc].print(stream, Eigen::Matrix3d::Identity());
              else
                asym_unit[i][j].clust_group()[nc].print(stream, lattice().inv_lat_column_mat());
            }
          }
          stream << "\n  ---------------------------------------------------------------------   \n\n";
        }
      }
      return;
    }
  */
  //***********************************************************

  void Structure::reset() {
    //std::cout << "begin reset() " << this << std::endl;

    for(Index nb = 0; nb < basis().size(); nb++) {
      m_basis[nb].set_basis_ind(nb);
    }
    within();
    m_factor_group.clear();

    /** Should we also invalidate the occupants?
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].site_occupant.set_value(-1);
    }
    */
    return;
  }

  //*********************************************************

  void Structure::map_superstruc_to_prim(Structure &prim) {

    Index prim_to_scel;
    Site shifted_site(prim.lattice());


    //Check that (*this) is actually a supercell of the prim
    if(!lattice().is_supercell_of(prim.lattice(), prim.point_group())) {
      default_err_log() << "*******************************************\n"
                        << "ERROR in Structure::map_superstruc_to_prim:\n"
                        << "The structure \n";
      default_err_log() << jsonParser(*this) << std::endl;
      default_err_log() << "is not a supercell of the given prim!\n";
      default_err_log() << jsonParser(prim) << std::endl;
      default_err_log() << "*******************************************\n";
      exit(1);
    }

    //Get prim grid of supercell to get the lattice translations
    //necessary to stamp the prim in the superstructure
    PrimGrid prim_grid(prim.lattice(), lattice());

    // Translate each of the prim atoms by prim_grid translation
    // vectors, and map that translated atom in the supercell.
    for(Index pg = 0; pg < prim_grid.size(); pg++) {
      for(Index pb = 0; pb < prim.basis().size(); pb++) {
        shifted_site = prim.basis()[pb];
        shifted_site += prim_grid.coord(pg, PRIM);
        shifted_site.set_lattice(lattice(), CART);
        shifted_site.within();

        // invalidate asym_ind and basis_ind because when we use
        // Structure::find, we don't want a comparison using the
        // basis_ind and asym_ind; we want a comparison using the
        // cartesian and Specie type.

        shifted_site.set_basis_ind(-1);
        prim_to_scel = find(shifted_site);

        if(prim_to_scel == basis().size()) {
          default_err_log() << "*******************************************\n"
                            << "ERROR in Structure::map_superstruc_to_prim:\n"
                            << "Cannot associate site \n"
                            << shifted_site << "\n"
                            << "with a site in the supercell basis. \n"
                            << "*******************************************\n";
          default_err_log() << "The basis_ind is "
                            << shifted_site.basis_ind() << "\n ";
          exit(2);
        }

        // Set ind_to_prim of the basis site
        //basis[prim_to_scel].ind_to_prim = pb;
      }
    }
  };

  //***********************************************************
  //
  void Structure::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {
    bool is_equiv(lattice() == new_lat);

    m_lattice = new_lat;

    for(Index nb = 0; nb < basis().size(); nb++) {
      m_basis[nb].set_lattice(lattice(), mode);
    }

    if(is_equiv)
      m_factor_group.set_lattice(lattice());
    else
      reset();
  }

  //***********************************************************
  /**
   * Loop through basis and rearrange atoms by type. Uses bubble
   * sort algorithm by comparing integer values of the strings
   * assigned to the basis occupants.
   *
   * The basis gets sorted in a sort of alphabetical way, so be
   * mindful of any POTCARs you might have if you run this.
   */
  //***********************************************************

  void Structure::sort_basis() {

    for(Index i = 0; i < basis().size(); i++) {
      for(Index j = 0; j < basis().size() - 1; j++) {

        if(basis()[j].occ_name() > basis()[j + 1].occ_name()) {
          std::iter_swap(m_basis.begin() + j, m_basis.begin() + j + 1);
        }
      }
    }
    set_site_internals();
    return;
  }

  //***********************************************************
  /**
   * Given a symmetry group, the basis of the structure will have
   * each operation applied to it. The resulting set of basis
   * from performing these operations will be averaged out,
   * yielding a new average basis that will replace the current one.
   */
  //***********************************************************

  void Structure::symmetrize(const SymGroup &relaxed_factors) {
    //First make a copy of your current basis
    //This copy will eventually become the new average basis.
    reset();
    std::vector<Site> avg_basis = basis();

    //Loop through given symmetry group an fill a temporary "operated basis"
    std::vector<Site> operbasis;
    for(Index rf = 0; rf < relaxed_factors.size(); rf++) {
      operbasis.clear();
      for(Index b = 0; b < basis().size(); b++) {
        operbasis.push_back(relaxed_factors[rf]*basis()[b]);
      }

      //Now that you have a transformed basis, find the closest mapping of atoms
      //Then average the distance and add it to the average basis
      for(Index b = 0; b < basis().size(); b++) {
        double smallest = 1000000;
        Coordinate bshift(lattice()), tshift(lattice());
        for(Index ob = 0; ob < operbasis.size(); ob++) {
          double dist = operbasis[ob].min_dist(basis()[b], tshift);
          if(dist < smallest) {
            bshift = tshift;
            smallest = dist;
          }
        }
        bshift.cart() *= (1.0 / relaxed_factors.size());
        avg_basis[b] += bshift;
      }

    }
    set_basis(avg_basis);
    //generate_factor_group();
    update();
    return;
  }

  //***********************************************************
  /**
   * Same as the other symmetrize routine, except this one assumes
   * that the symmetry group you mean to use is the factor group
   * of your structure within a certain tolerance.
   *
   * Notice that the tolerance is also used on your point group!!
   */
  //***********************************************************

  void Structure::symmetrize(const double &tolerance) {
    double orig_tol = lattice().tol();
    m_lattice.set_tol(tolerance);
    generate_factor_group();
    SymGroup g = factor_group();
    symmetrize(g);
    m_lattice.set_tol(orig_tol);
    return;
  }

  //This function gets the permutation representation of the
  // factor group operations of the structure. It first applies
  // the factor group operation to the structure, and then tries
  // to map the new position of the basis atom to the various positions
  // before symmetry was applied. It only checks the positions after
  // it brings the basis within the crystal.
  void Structure::_generate_basis_symreps(bool verbose) const {
    std::string clr(100, ' ');
    if(factor_group().size() <= 0) {
      default_err_log() << "ERROR in generate_basis_permutation_representation" << std::endl;
      default_err_log() << "Factor group is empty." << std::endl;;
      exit(1);
    }

    std::vector<UnitCellCoord> sitemap;

    m_basis_perm_rep_ID = m_factor_group.allocate_representation();

    for(std::string const &dof : continuous_local_dof_types(*this)) {
      for(Site const &site : basis()) {
        if(site.has_dof(dof))
          site.dof(dof).allocate_symrep(m_factor_group);
      }
    }

    // The sitemap specifies that op*basis(b) -> sitemap[b] (which is a UnitCellCoord)
    // The new dofs of site specified by UCC sitemap[b] will be transformations of the dofs
    // that previously resided at basis(b). As such, for dofs, we use the inverse permutation and
    //   basis()[b].symrep(doftype.name()) = basis()[b].dof(doftype.name()).basis().transpose()
    //                                       * doftype.symop_to_matrix(op)
    //                                       * basis()[sitemap[b].sublat()].dof(doftype.name().basis())
    Eigen::MatrixXd trep, trepblock;
    for(Index s = 0; s < m_factor_group.size(); ++s) {
      SymOp const &op = m_factor_group[s];
      if(verbose) {
        if(op.index() % 100 == 0)
          std::cout << '\r' << clr.c_str() << '\r' << "Find permute rep for symOp " << op.index() << "/" << m_factor_group.size() << std::flush;
      }

      sitemap = symop_site_map(op, *this);
      op.set_rep(m_basis_perm_rep_ID, SymBasisPermute(op, lattice(), sitemap));

      for(Index b = 0; b < basis().size(); ++b) {
        // copy_aply(symop,dofref_from) = P.permute(dofref_to);
        auto const &dofref_to = basis()[sitemap[b].sublat()].site_occupant();
        auto const &dofref_from = basis()[b].site_occupant();
        OccupantDoFIsEquivalent<Molecule> eq(dofref_from);
        if(eq(op, dofref_to)) {
          if(dofref_from.symrep_ID().is_identity()) {
            if(!eq.perm().is_identity()) {
              dofref_from.allocate_symrep(m_factor_group);
              Index s2;
              for(s2 = 0; s2 < s; ++s2) {
                m_factor_group[s2].set_rep(dofref_from.symrep_ID(), SymPermutation(sequence<Index>(0, dofref_from.size())));
              }
              m_factor_group[s2].set_rep(dofref_from.symrep_ID(), SymPermutation(eq.perm().inverse()));
            }
          }
          else {
            op.set_rep(dofref_from.symrep_ID(), SymPermutation(eq.perm().inverse()));
          }
        }
        else throw std::runtime_error("In Structure::_generate_basis_symreps(), Sites originally identified as equivalent cannot be mapped by symmetry.");
      }

      for(auto const &dof_dim : local_dof_dims(*this)) {
        for(Index b = 0; b < basis().size(); ++b) {
          if(!basis()[b].has_dof(dof_dim.first))
            continue;

          // copy_aply(symop,dofref_from) = dofref_to * U
          DoFSet const &dofref_to = basis()[sitemap[b].sublat()].dof(dof_dim.first);
          DoFSet const &dofref_from = basis()[b].dof(dof_dim.first);
          DoFIsEquivalent eq(dofref_from);
          if(!eq(op, dofref_to)) {
            throw std::runtime_error("While generating symmetry representation for local DoF \""
                                     + dof_dim.first
                                     + "\", a symmetry operation was identified that invalidates the degree of freedom. "
                                     + "Degrees of freedom must be fully specified before performing symmetry analyses.");
          }

          trep.setIdentity(dof_dim.second, dof_dim.second);
          trep.topLeftCorner(dofref_from.size(), dofref_from.size()) = eq.U();
          op.set_rep(dofref_from.symrep_ID(), SymMatrixXd(trep));
        }
      }
    }

    if(verbose) std::cout << '\r' << clr.c_str() << '\r' << std::flush;
    return;
  }

  //***********************************************************

  void Structure::_generate_global_symreps(bool verbose) const {
    if(factor_group().size() <= 0) {
      default_err_log() << "ERROR in generate_global_dof_representations" << std::endl;
      default_err_log() << "Factor group is empty." << std::endl;
      exit(1);
    }
    //std::cout << "INSIDE _generate_global_symreps\n";
    for(auto const &dof : m_dof_map) {
      dof.second.allocate_symrep(m_factor_group);
      for(SymOp const &op : m_factor_group) {
        DoFIsEquivalent eq(dof.second);
        if(!eq(op)) {
          throw std::runtime_error("While generating symmetry representation for global DoF \""
                                   + dof.first
                                   + "\", a symmetry operation was identified that invalidates the degree of freedom. "
                                   + "Degrees of freedom must be fully specified before performing symmetry analyses.");
        }
        op.set_rep(dof.second.symrep_ID(), SymMatrixXd(eq.U().transpose()));
      }

    }
  }
  //***********************************************************

  //Sets the occupants in the basis sites to those specified by occ_index
  void Structure::set_occs(std::vector <int> occ_index) {
    if(occ_index.size() != basis().size()) {
      default_err_log() << "The size of the occ index and basis index do not match!\nEXITING\n";
      exit(1);
    }
    for(Index i = 0; i < basis().size(); i++) {
      m_basis[i].set_occ_value(occ_index[i]);
    }
  }

  //***********************************************************
  Structure &Structure::operator+=(const Coordinate &shift) {

    for(Index i = 0; i < basis().size(); i++) {
      m_basis[i] += shift;
    }

    m_factor_group += shift.cart();
    return (*this);
  }


  //***********************************************************
  Structure &Structure::operator-=(const Coordinate &shift) {

    for(Index i = 0; i < basis().size(); i++) {
      m_basis[i] -= shift;
    }
    m_factor_group -= shift.cart();
    return (*this);
  }

  //****************************************************

  jsonParser &Structure::to_json(jsonParser &json) const {

    // class Structure : public BasicStructure<Site>
    BasicStructure<Site>::to_json(json);

    // mutable MasterSymGroup m_factor_group;
    json["factor_group"] = m_factor_group;

    return json;
  }

  //****************************************************

  // Assumes constructor CoordType::CoordType(Lattice) exists
  void Structure::from_json(const jsonParser &json) {


    // class Structure : public BasicStructure<Site>
    BasicStructure<Site> &basic = *this;
    basic.from_json(json);

    // mutable MasterSymGroup m_factor_group;
    m_factor_group.clear();
    m_factor_group.from_json(json["factor_group"]);

    update();
  }

  //****************************************************

  jsonParser &to_json(const Structure &structure, jsonParser &json) {
    return structure.to_json(json);
  }

  void from_json(Structure &structure, const jsonParser &json) {
    structure.from_json(json);
  }

  //***********************************************************
  /*
  Structure operator*(const SymOp &LHS, const Structure &RHS) { //AAB

    return Structure(RHS).apply_sym(LHS);
  }
  */
  //***********************************************************
  Structure operator*(const Lattice &LHS, const Structure &RHS) {
    Structure tsuper(LHS);
    tsuper.fill_supercell(RHS);
    return tsuper;
  }






  /// Helper Functions


  /// Returns 'converter' which converts Site::site_occupant indices to 'mol_list' indices:
  ///   mol_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > make_index_converter(const Structure &struc, std::vector<Molecule> mol_list) {

    std::vector< std::vector<Index> > converter(struc.basis().size());

    for(Index i = 0; i < struc.basis().size(); i++) {
      converter[i].resize(struc.basis()[i].site_occupant().size());

      for(Index j = 0; j < struc.basis()[i].site_occupant().size(); j++) {
        converter[i][j] = find_index(mol_list, struc.basis()[i].site_occupant()[j]);
      }
    }

    return converter;

  }

  /// Returns 'converter' which converts Site::site_occupant indices to 'mol_name_list' indices:
  ///   mol_name_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > make_index_converter(const Structure &struc, std::vector<std::string> mol_name_list) {

    std::vector< std::vector<Index> > converter(struc.basis().size());

    for(Index i = 0; i < struc.basis().size(); i++) {
      converter[i].resize(struc.basis()[i].site_occupant().size());

      for(Index j = 0; j < struc.basis()[i].site_occupant().size(); j++) {
        converter[i][j] = find_index(mol_name_list, struc.basis()[i].site_occupant()[j].name());
      }
    }

    return converter;

  }

  /// Returns 'converter_inverse' which converts 'mol_name_list' indices to Site::site_occupant indices:
  ///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
  ///
  /// If mol is not allowed on basis_site, return struc.basis()[basis_site].site_occupant().size()
  std::vector< std::vector<Index> > make_index_converter_inverse(const Structure &struc, std::vector<std::string> mol_name_list) {

    std::vector< std::vector<Index> > converter_inv(struc.basis().size());

    for(Index i = 0; i < struc.basis().size(); i++) {
      converter_inv[i].resize(mol_name_list.size());

      std::vector<std::string> site_occ_name_list;
      for(Index j = 0; j < struc.basis()[i].site_occupant().size(); j++) {
        site_occ_name_list.push_back(struc.basis()[i].site_occupant()[j].name());
      }

      for(Index j = 0; j < mol_name_list.size(); j++) {
        converter_inv[i][j] = find_index(site_occ_name_list, mol_name_list[j]);
      }
    }

    return converter_inv;

  }

};


