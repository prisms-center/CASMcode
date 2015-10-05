#ifndef SUPERCELL_HH
#define SUPERCELL_HH

#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigEnumIterator.hh"
#include "casm/clex/ConfigDoF.hh"

namespace CASM {


  template<typename T, typename U> class ConfigIterator;
  class PermuteIterator;
  class PrimClex;
  class Clexulator;

  class Supercell {

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
    mutable Index m_perm_symrep_ID;

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

    // m_trans_permute describes how translations by a primitive lattice translation permutes the sites
    // of the Supercell. m_trans_permute[i] is the effect of translating all sites by (*this).uccoord[i]
    // and mapping them back within the Supercell. m_trans_permute.size() == (*this).volume()
    // NOTE: m_trans_permute can be thought of as the translational factor group formed by the the cosets
    //       of the group Tsuper in the group Tprim, as they are defined above
    //mutable Array<Permutation> m_trans_permute;

    //Indices that map the linear index to the bijk point in real space
    //  Generate in void Supercell::fill_supercell() which is called in constructors
    //  This always exists and is populated
    //Array < UnitCellCoord > config_index_to_bijk;

    //Indices that map the linear index to the apqr point in reciprocal space
    //Array < UnitCellCoord > config_index_to_apqr;

    //Array of coordinates for both real and reciprocal space
    //Array<Coordinate> real_coords;  //primitive=primitive
    //Array<Coordinate> recip_coords; //primitive=recip(primitive)

    //The current 'state' of the supercell in reciprocal space
    //  Configuration recip_curr_state;

    /// unique name of the supercell based on hermite normal form (see generate_name() )
    std::string name;


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

    Array< Array<Index> >    nlists;  //[scell site][ nbor_indices]

    // Could hold either enumerated configurations or any 'saved' configurations
    ConfigList config_list;

    Matrix3 < int > transf_mat;

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
    Supercell(PrimClex *_prim, const Matrix3<int> &superlattice_matrix);
    //Supercell(PrimClex *_prim, const Eigen::Matrix3i &superlattice_matrix);   //I wish
    //Supercell(PrimClex *_prim, const Matrix3<int> &superlattice_matrix, const Lattice &superlattice, const int warningFlag = 1);
    //Supercell(PrimClex *_prim, std::string _name);
    //Supercell(PrimClex *_prim, const Lattice &superlattice, std::string _name);
    //Supercell(PrimClex *_prim, const Matrix3<int> &superlattice_matrix, std::string _name);
    //Supercell(PrimClex *_prim, const Matrix3<int> &superlattice_matrix, const Lattice &superlattice, std::string _name, const int warningFlag = 1);

    // **** Coordinates ****
    //UnitCellCoord get_bijk(Vector3<double> cartesian_coord, int b) const;
    Index get_linear_index(const Site &site, double tol = TOL) const;
    Index get_linear_index(const Coordinate &coord, double tol = TOL) const;
    Index find(const UnitCellCoord &bijk) const;
    Coordinate coord(const UnitCellCoord &bijk) const;
    Coordinate coord(Index linear_ind) const;

    // only used for populate_bijk_l_map
    //void location_within(const UnitCellCoord &supercell_point, Array< Array < Array <Array <Index > > > > &linear_index, const UnitCellCoord &centering);

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
    //PrimClex &get_primclex() {
    //return *primclex;
    //}

    const PrimClex &get_primclex() const {
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
    Index permutation_symrep_ID()const {
      if(m_perm_symrep_ID == Index(-1))
        generate_permutations();
      return m_perm_symrep_ID;
    }

    SymGroupRep const *permutation_symrep() const {
      return get_prim().factor_group().representation(permutation_symrep_ID());
    }

    Index get_b(Index i) const {
      return i / volume();
    }

    Matrix3<int> get_transf_mat() const {
      return transf_mat;
    };
    Lattice get_real_super_lattice() const {
      return real_super_lattice;
    };
    Lattice get_recip_prim_lattice() const {
      return recip_prim_lattice;
    };

    // get indices of neighbor sites ('nlist_index') in Configuration to some 'site'
    Index get_nlist_l(Index pivot_l, Index nlist_index) const {
      return nlists[pivot_l][nlist_index];
    };

    const Array<Index> &get_nlist(Index pivot_l) const {
      return nlists[pivot_l];
    };


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

    /*
    TransitionList &transition_list() {
      return m_transition_list;
    };

    const TransitionList &transition_list() const {
      return m_transition_list;
    };

    Transition &transition(int i) {
      return m_transition_list[i];
    }

    const Transition &transition(int i) const {
      return m_transition_list[i];
    };

    // begin and end iterators for iterating over transitions
    trans_iterator trans_begin();
    trans_iterator config_end();

    // begin and end const_iterators for iterating over transitions
    trans_const_iterator trans_cbegin() const;
    trans_const_iterator trans_cend() const;
    */

    Index get_id() const {
      return m_id;
    }

    std::string get_name() const {
      return name;
    };

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

    // begin and end iterators for iterating over translation and factor group permutations
    permute_const_iterator permute_begin() const;
    permute_const_iterator permute_end() const;

    ///Return path to supercell directory
    fs::path get_path() const;

    ///Count how many configs are selected in *this
    Index amount_selected() const;

    // **** Mutators ****

    void set_id(Index id) {
      m_id = id;
    }

    // **** Generating functions ****

    // Populate m_factor_group -- probably should be private
    void generate_factor_group() const;

    // Populate m_trans_permute -- probably should be private
    void generate_permutations() const;

    //void fill_supercell();
    //void populate_bijk_l_map(Array< Array < Array <Array <Index > > > > &linear_index, UnitCellCoord &centering);
    void generate_neighbor_list();

    ///Return true if the Supercell is smaller than the neighborhood of the sites, causing periodic overlap
    bool neighbor_image_overlaps() const;

    //\John G 070713
    void generate_name();

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

    /// Loop over all configurations enumerated by 'enumerator' and add them to the Supercell, if they are not already present
    void add_enumerated_configurations(ConfigEnum<Configuration> &enumerator);

    /// Loop over all configurations enumerated by some (invisible) enumerator from 'it_begin' to 'it_end', adding them to the Supercell, if they are not already present
    void add_enumerated_configurations(ConfigEnumIterator<Configuration> it_begin, ConfigEnumIterator<Configuration> it_end);

    /// Enumerate all possible occupation configurations that are symmetrically equivalent and fit inside this supercell (but cannot be described by a smaller supercell)
    void enumerate_all_occupation_configurations();

    /// Enumerate 'Nstep' configurations that linearly interpolate deformation and displacement from 'initial' configuration to 'final' configuration
    /// 'initial' and 'final' must either have the same occupation or have unspecified occupation
    /// The range can be adjusted using 'being_delta' and 'end_delta' (which can be positive or negative). begin_delta<0 indicates interpolation starts
    /// abs(begin_delta) number of steps *before* 'initial'. positive values indicate interpolation starts the specified number of steps *after* initial.
    /// semantics are similar for 'end_delta', but for the final configuration.
    void enumerate_interpolated_configurations(Supercell::config_const_iterator initial, Supercell::config_const_iterator final,
                                               long Nstep, long begin_delta = 0, long end_delta = 0);


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
    void read_config_list(const jsonParser &json);

    template<typename ConfigIterType>
    void add_configs(ConfigIterType it_begin, ConfigIterType it_end);

    // **** Correlations ****
    //void calc_correlations_curr();
    //void calc_corr_all();
    //void calc_curr_energy();
    //void calc_clex_correlations();

    //Fill up values in the correlations of all configurations
    //    void populate_correlations(Clexulator &clexulator);
    //Fill up values in the correlations of configuration specified by index
    //    void populate_correlations(Clexulator &clexulator, const Index &config_index);


    // **** Other ****
    // Reads a relaxed structure and calculates the strains and stretches using the reference structure
    void read_relaxed_structure(Index configNum, const Lattice &home_lattice);
    void read_relaxed_structure(Index configNum);
    void read_clex_relaxations(const Lattice &home_lattice);

    bool is_supercell_of(const Structure &structure) const;
    bool is_supercell_of(const Structure &structure, Matrix3<double> &multimat) const;
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

    // Va_mode		description
    // 0			print no information about the vacancies
    // 1			print only the coordinates of the vacancies
    // 2			print the number of vacancies and the coordinates of the vacancies
    void print(const Configuration &config, std::ostream &stream, COORD_TYPE mode, int Va_mode = 0, char term = 0, int prec = 7, int pad = 5) const;
    void print(std::ostream &stream, Configuration tconfig) const;

    ///Call Configuration::write out every configuration in supercell
    jsonParser &write_config_list(jsonParser &json);

    void printUCC(std::ostream &stream, COORD_TYPE mode, UnitCellCoord ucc, char term = 0, int prec = 7, int pad = 5) const;
    //\Michael 241013

    // **** ParamComposition Calculators ****
    //    std::map < std::string , double > composition_calculate(int ConfigNum);
    //    std::map < std::string , double > true_composition_calculate(int ConfigNum);
    //    std::map < int , std::map <std::string , double > > sublattice_composition_calculate(int ConfigNum);
    //Array< Array< double > > sublattice_composition_calculate(int ConfigNum);
    //    std::map < std::string , double > calculate_composition(int ConfigName);
    //    std::map < std::string , double > calculate_true_composition(int ConfigName);
    //void populate_sublat_to_comp();

  };

  template<typename ConfigIterType>
  void Supercell::add_configs(ConfigIterType it_begin, ConfigIterType it_end) {
    if(ConfigIterType::is_canonical_iter()) {
      for(; it_begin != it_end; ++it_begin) {
        add_canon_config(*it_begin);
      }
    }
    else {
      for(; it_begin != it_end; ++it_begin) {
        add_config(*it_begin);
      }
    }
  }

}
#endif
