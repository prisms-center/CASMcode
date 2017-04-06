#ifndef SUPERCELL_HH
#define SUPERCELL_HH

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {


  template<typename T, typename U> class ConfigIterator;
  class PermuteIterator;
  class PrimClex;
  class Clexulator;

  struct ConfigMapCompare {

    bool operator()(const Configuration *A, const Configuration *B) const {
      return *A < *B;
    }

  };

  /** \defgroup Supercell
   *  \ingroup Clex
   *  \brief Represents a supercell of the primitive parent crystal structure
   *  @{
   */

  /// \brief Represents a supercell of the primitive parent crystal structure
  ///
  class Supercell : public Comparisons<Supercell> {

  public:

    typedef boost::container::stable_vector<Configuration> ConfigList;
    //typedef boost::container::stable_vector<Transition> TransitionList;

    typedef ConfigIterator<Configuration, PrimClex> config_iterator;
    typedef ConfigIterator<const Configuration, const PrimClex> config_const_iterator;
    //typedef ConfigIterator<Transition, PrimClex> transition_iterator;
    //typedef ConfigIterator<const Transition, const PrimClex> transition_const_iterator;
    typedef PermuteIterator permute_const_iterator;

  private:
    // pointer to Primcell containing all the cluster expansion data
    PrimClex *primclex;

    // lattice of supercell in real space
    Lattice real_super_lattice;

    // reciprocal of real_super_lattice (grid of recip scell lattice)
    Lattice recip_prim_lattice;

    PrimGrid m_prim_grid;
    //Grid in reciprocal space of the supercell that perfectly tiles
    //the prim cell
    PrimGrid recip_grid;

    // m_perm_symrep_ID is the ID of the SymGroupRep of get_prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).get_prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       the the Supercel::factor_group_permute(int) method, so that you don't encounter the gaps
    //       OR, see note for Supercell::permutation_symrep() below.
    mutable SymGroupRepID m_perm_symrep_ID;

    // m_factor_group is factor group of the super cell, found by identifying the subgroup of
    // (*this).get_prim().factor_group() that leaves the supercell lattice vectors unchanged
    // if (*this).get_prim() is actually primitive, then m_factor_group.size() <= 48
    // NOTE: This is different from the SymGroup found by doing (*this).superstruc().factor_group()
    //       if Tprim is the translation group formed by the primitive cell lattice vectors, then
    //       m_factor_group is the group formed by the cosets of Tprim in the supercell space group
    //       if Tsuper is the translation group formed by the supercell lattice vectors, then,
    //       m_occupation(init_config.get_occupation()),
    //       m_displacement(init_config.get_displacement()),
    //       m_strain(init_config.get_supercell().strain
    //       (*this).superstruc().factor_group() is the group formed by the cosets of Tsuper in the supercell space group
    mutable SymGroup m_factor_group;

    /// unique name of the supercell based on hermite normal form (see _generate_name() )
    mutable std::string m_name;


    ///************************************************************************************************
    /// STRUCTURE FACTOR ROUTINES
    /// Fourier matrix is arranged as :
    ///           [exp(-i k1 r1)  exp(- k2 r1) ... exp(-i kn r1)]
    ///           [exp(-i k1 r2)  exp(- k2 r2) ... exp(-i kn r2)]
    ///           ...
    ///           [exp(-i k1 rn)  exp(- k2 rn) ... exp(-i kn rn)]
    /// r_{i} are the real space coordinates of the prim grid points
    /// k_{i} are k-points commensurate with the supercell
    Eigen::MatrixXcd m_fourier_matrix;

    /// The phase factor matrix is arranged as:
    ///           [exp(-i k1 tau1)  exp(-i k2 tau1) ... exp(-i kn tau1)]
    ///           [exp(-i k1 tau2)  exp(-i k2 tau2) ... exp(-i kn tau2)]
    ///           ...
    ///           [exp(-i k1 taun)  exp(-i k2 taun) ... exp(-i kn taun)]
    /// tau_{i} is the shift vector for the ith basis atom
    Eigen::MatrixXcd m_phase_factor;

    /// Calculating the structure factor:
    ///    vectors you need : Eigen::VectorXd intensities //this needs to be as long
    ///                                                   //as the number of sites in the supercell
    ///    Calculations (these are done in Configuration):
    ///       intensities.segment<volume()>(i*volume()) * fourier_matrix = Q.column(i)
    ///       Q is arranged as: [Q1(k1) Q2(k1) ... Qn(k1)]
    ///                         [Q1(k2) Q2(k2) ... Qn(k2)]
    ///                         ...
    ///                         [Q1(kn) Q2(kn) ... Qn(kn)]
    ///       Structure factors are then calculated as S:
    ///       S = Q * m_phase_factor
    ///       S is arranged as: [S(k1) S(k2) ... S(kn)]
    ///       In the code: Q is called sublat_sf

    /// k-point mesh for the Fourier Transform
    /// An Eigen::MatrixXd that should be nx3
    /// each row is a coordinate in reciprocal space
    Eigen::MatrixXd m_k_mesh;

    /// Structure factor calculation routines -> Use only the public versions of these functions
    /// The calculations, and math is explained as comments above m_fourier_matrix and m_phase_factor
    void generate_fourier_matrix(const Eigen::MatrixXd &real_coordinates, const Eigen::MatrixXd &recip_coordinates, const bool &override);
    void generate_phase_factor(const Eigen::MatrixXd &shift_vectors, const Array<bool> &is_commensurate, const bool &override);
    ///************************************************************************************************

    /// SuperNeighborList, mutable for lazy construction
    mutable notstd::cloneable_ptr<SuperNeighborList> m_nlist;

    /// Store size of PrimNeighborList at time of construction of SuperNeighborList
    /// to enable checking if SuperNeighborList should be re-constructed
    mutable Index m_nlist_size_at_construction;

    /// Store a pointer to the canonical equivalent Supercell
    mutable Supercell *m_canonical;

    // Could hold either enumerated configurations or any 'saved' configurations
    ConfigList config_list;

    // Improve performance of 'contains_config' by sorting Configuration
    std::map<const Configuration *, Index, ConfigMapCompare> m_config_map;

    Eigen::Matrix3i transf_mat;

    double scaling;

    /// index into PrimClex::supercell_list
    Index m_id;

  public:

    //The current 'state' of the supercell in real space
    //Configuration curr_state;

    // **** Constructors ****

    //Supercell(PrimClex *_prim);
    Supercell(const Supercell &RHS);
    Supercell(PrimClex *_prim, const Lattice &superlattice);
    Supercell(PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &superlattice_matrix);


    // **** Coordinates ****
    Index get_linear_index(const Site &site, double tol = TOL) const;
    Index get_linear_index(const Coordinate &coord, double tol = TOL) const;
    Index find(const UnitCellCoord &bijk) const;
    Coordinate coord(const UnitCellCoord &bijk) const;
    Coordinate coord(Index linear_ind) const;

    // returns maximum allowed occupation bitstring -- used for initializing enumeration counters
    ReturnArray<int> max_allowed_occupation() const;

    Configuration configuration(const BasicStructure<Site> &structure_to_config, double tol = TOL);

    // return Structure corresponding to this supercell
    //    w/ basis site occupation as per primclex.prim
    Structure superstructure() const;
    //    w/ basis site occupation as per config
    Structure superstructure(const Configuration &config) const;        //This should be private
    ///Returns a structure corresponding to the specified configuration.
    Structure superstructure(Index config_index) const;
    //    w/ basis site occupation as per config
    //    and itself as prim, not primclex->prim

    //Structure structure(const Configuration &config) const;
    //Routines that generate real and reciprocal coordinates
    //for *this supercell
    Eigen::MatrixXd real_coordinates() const;
    Eigen::MatrixXd recip_coordinates() const;

    // **** Accessors ****

    PrimClex &get_primclex() const {
      return *primclex;
    }

    const PrimGrid &prim_grid() const {
      return m_prim_grid;
    }

    const Structure &get_prim() const;

    ///Return number of primitive cells that fit inside of *this
    Index volume()const {
      return m_prim_grid.size();
    };

    Index basis_size() const {
      return get_prim().basis.size();
    }

    Index num_sites()const {
      return volume() * basis_size();
    };

    UnitCellCoord uccoord(Index i) const {
      UnitCellCoord t_bijk = m_prim_grid.uccoord(i % volume());
      t_bijk[0] = get_b(i);
      return t_bijk;
    };
    const Eigen::MatrixXcd &fourier_matrix() const {
      return m_fourier_matrix;
    }

    const Eigen::MatrixXcd &phase_factor() const {
      return m_phase_factor;
    }

    const Eigen::MatrixXd &k_mesh() const {
      return m_k_mesh;
    }

    // the permutation_symrep is the SymGroupRep of get_prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).get_prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       SymGroupRep::get_representation(m_factor_group[i]) or SymGroupRep::get_permutation(m_factor_group[i]),
    //       so that you don't encounter the gaps (i.e., the representation can be indexed using the
    //       SymOps of m_factor_group
    SymGroupRepID permutation_symrep_ID()const {
      if(m_perm_symrep_ID.empty())
        generate_permutations();
      return m_perm_symrep_ID;
    }

    SymGroupRep const &permutation_symrep() const {
      return get_prim().factor_group().representation(permutation_symrep_ID());
    }

    Index get_b(Index i) const {
      return i / volume();
    }

    const Eigen::Matrix3i &get_transf_mat() const {
      return transf_mat;
    };

    const Lattice &get_real_super_lattice() const {
      return real_super_lattice;
    };

    const Lattice &get_recip_prim_lattice() const {
      return recip_prim_lattice;
    };

    /// \brief Returns the SuperNeighborList
    const SuperNeighborList &nlist() const;


    ConfigList &get_config_list() {
      return config_list;
    };

    const ConfigList &get_config_list() const {
      return config_list;
    };

    const Configuration &get_config(Index i) const {
      return config_list[i];
    };

    Configuration &get_config(Index i) {
      return config_list[i];
    }

    // begin and end iterators for iterating over configurations
    config_iterator config_begin();
    config_iterator config_end();

    // begin and end const_iterators for iterating over configurations
    config_const_iterator config_cbegin() const;
    config_const_iterator config_cend() const;

    Index get_id() const {
      return m_id;
    }

    /// \brief Return supercell name
    std::string get_name() const;

    // Populates m_factor_group (if necessary) and returns it.
    const SymGroup &factor_group() const;

    // Returns the permutation representation of the i'th element of m_factor_group
    const Permutation &factor_group_permute(Index i) const;

    // Returns the i'th element of m_trans_permute
    // Populates m_trans_permute if needed
    const Permutation &translation_permute(Index i) const;

    // Const access of m_trans_permute
    // Populates m_trans_permute if needed
    const Array<Permutation> &translation_permute() const;

    /// \brief Begin iterator over pure translational permutations
    permute_const_iterator translate_begin() const;

    /// \brief End iterator over pure translational permutations
    permute_const_iterator translate_end() const;

    // begin and end iterators for iterating over translation and factor group permutations
    permute_const_iterator permute_begin() const;
    permute_const_iterator permute_end() const;
    permute_const_iterator permute_it(Index fg_index, Index trans_index) const;

    ///Return path to supercell directory
    fs::path get_path() const;

    ///Count how many configs are selected in *this
    Index amount_selected() const;

    bool is_canonical() const;

    SymOp to_canonical() const;

    SymOp from_canonical() const;

    Supercell &canonical_form() const;

    bool is_equivalent(const Supercell &B) const;

    bool operator<(const Supercell &B) const;

    // **** Mutators ****

    void set_id(Index id) {
      m_id = id;
    }

    // **** Generating functions ****

    // Populate m_factor_group -- probably should be private
    void generate_factor_group() const;

    // Populate m_trans_permute -- probably should be private
    void generate_permutations() const;

    /// Calculate reference properties for a configuration (must have reference states in appropriate directories)
    void generate_reference_config_props(Index config_index);
    /// Calculate reference properties for each configuration in *this (see above)
    void generate_all_reference_config_props();

    /// Calculate delta properties for a configuration (must have read in calculated and reference properties)
    void generate_delta_config_props(Index config_index);
    /// Calculate delta properties for each configuration in *this (see above)
    void generate_all_delta_config_props();

    /// Structure Factor
    void generate_fourier_matrix();
    void generate_fourier_matrix(const Eigen::MatrixXd &real_coordinates, const Eigen::MatrixXd &recip_coordinates);
    Array< bool > is_commensurate_kpoint(const Eigen::MatrixXd &recip_coordinates, double tol = TOL);
    void populate_structure_factor();
    void populate_structure_factor(const Index &config_index);

    // **** Enumerating functions ****

    //Functions for enumerating configurations that are perturbations of a 'background' structure
    void enumerate_perturb_configurations(const std::string &background, fs::path CSPECS, double tol = TOL, bool verbose = false, bool print = false);
    void enumerate_perturb_configurations(Configuration background_config, fs::path CSPECS, double tol = TOL, bool verbose = false, bool print = false);
    void enumerate_perturb_configurations(const Structure &background, fs::path CSPECS, double tol = TOL, bool verbose = false, bool print = false);

    void enumerate_perturb_configurations(Configuration background_config,
                                          const SiteOrbitree &background_tree,
                                          Array< Array< Array<Index> > > &config_index,
                                          Array< Array< Array<permute_const_iterator> > > &config_symop_index,
                                          jsonParser &jsonsrc,
                                          double tol = TOL);

    bool contains_config(const Configuration &config) const;
    bool contains_config(const Configuration &config, Index &index) const;
    bool add_config(const Configuration &config);
    bool add_config(const Configuration &config, Index &index, Supercell::permute_const_iterator &permute_it);
    bool add_canon_config(const Configuration &config, Index &index);

    std::pair<config_const_iterator, bool> insert_config(const Configuration &config);
    std::pair<config_const_iterator, bool> insert_canon_config(const Configuration &config);
    config_const_iterator find(const Configuration &config) const;

    void read_config_list(const jsonParser &json);

    template<typename ConfigIterType>
    void add_unique_canon_configs(ConfigIterType it_begin, ConfigIterType it_end);

    template<typename ConfigIterType>
    void add_configs(ConfigIterType it_begin, ConfigIterType it_end);

    // **** Other ****
    // Reads a relaxed structure and calculates the strains and stretches using the reference structure
    void read_relaxed_structure(Index configNum, const Lattice &home_lattice);
    void read_relaxed_structure(Index configNum);
    void read_clex_relaxations(const Lattice &home_lattice);

    bool is_supercell_of(const Structure &structure) const;
    bool is_supercell_of(const Structure &structure, Eigen::Matrix3d &multimat) const;
    ReturnArray<int> vacant()const;

    // **** Printing ****

    void print_bijk(std::ostream &stream);
    //   void print_clex_correlations(std::ostream &corrFile);
    ///Old CASM style corr.in output for all the configurations in *this supercell
    //   void print_global_correlations_simple(std::ostream &corrstream) const;
    void print_sublat_to_comp(std::ostream &stream);
    void print_PERTURB_json(std::ofstream &file,
                            const Configuration &background_config,
                            const Array< Array< Array<Index > > > &perturb_config_index,
                            const Array< Array< Array<permute_const_iterator> > > &perturb_config_symop_index,
                            bool print_config_name) const;

    ///Call Configuration::write out every configuration in supercell
    jsonParser &write_config_list(jsonParser &json);

    void printUCC(std::ostream &stream, COORD_TYPE mode, UnitCellCoord ucc, char term = 0, int prec = 7, int pad = 5) const;
    //\Michael 241013

  private:

    friend Comparisons<Supercell>;

    bool _eq(const Supercell &B) const;

    void _add_canon_config(const Configuration &config);

    void _generate_name() const;

  };

  Supercell &apply(const SymOp &op, Supercell &scel);

  Supercell copy_apply(const SymOp &op, const Supercell &scel);


  //*******************************************************************************
  // Warning: Assumes configurations are in canonical form
  template<typename ConfigIterType>
  void Supercell::add_unique_canon_configs(ConfigIterType it_begin, ConfigIterType it_end) {
    // Remember existing configs, to avoid duplicates
    //   Enumerated configurations are added after existing configurations
    Index N_existing = config_list.size();
    Index N_existing_enumerated = 0;
    Index index;
    //std::cout << "ADDING CONFIGS TO SUPERCELL; N_exiting: " << N_existing << " N_enumerated: " << N_existing_enumerated << "\n";
    //std::cout << "beginning iterator: " << it_begin->occupation() << "\n";
    // Loops through all possible configurations
    for(; it_begin != it_end; ++it_begin) {
      //std::cout << "Attempting to add configuration: " << it_begin->occupation() << "\n";
      // Adds the configuration to the list, if not among previously existing configurations

      bool add = true;
      if(N_existing_enumerated != N_existing) {
        if(contains_config(*it_begin, index)) {
          config_list[index].push_back_source(it_begin->source());
          add = false;
          N_existing_enumerated++;
        }
      }
      if(add) {
        _add_canon_config(*it_begin);
      }
    }

  }

  //*******************************************************************************

  template<typename ConfigIterType>
  void Supercell::add_configs(ConfigIterType it_begin, ConfigIterType it_end) {
    for(; it_begin != it_end; ++it_begin) {
      add_config(*it_begin);
    }
  }

  std::string generate_name(const Eigen::Matrix3i &transf_mat);

  /** @} */
}
#endif
