#ifndef CASM_Configuration
#define CASM_Configuration

#include <map>

#include "casm/external/boost.hh"

#include "casm/misc/Comparisons.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/container/Array.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/Properties.hh"
#include "casm/clex/Correlation.hh"
#include "casm/clusterography/Orbitree.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ConfigIterator.hh"

namespace CASM {

  class PrimClex;
  class Supercell;
  class UnitCellCoord;
  class Clexulator;
  class FillSupercell;
  template<typename ConfigType, typename PrimClexType>
  class ConfigIterator;


  namespace ConfigIO {
    class Selected;
  }
  template <bool IsConst> class ConfigSelection;
  typedef ConfigSelection<true> ConstConfigSelection;

  template<typename DataObject> struct QueryTraits;

  template<>
  struct QueryTraits<Configuration> {
    static const std::string name;
    typedef ConfigIO::Selected Selected;
    typedef ConstConfigSelection Selection;
  };

  /// \defgroup Configuration
  ///
  /// \brief A Configuration represents the values of all degrees of freedom in a Supercell
  ///
  /// \ingroup Clex
  ///
  /// @{


  /// \brief Holds results of Configuration::insert
  ///
  /// - 'canonical' refers to the canonical form of the Configuration in it's
  ///   canonical equivalent Supercell.  The canonical form may be primitive or
  ///   non-primitive
  /// - 'primitive' refers to the primitive canonical Configuration.
  ///
  struct ConfigInsertResult {

    typedef ConfigIterator<const Configuration, const PrimClex> config_const_iterator;

    /// True if primitive did not exist before insertion
    bool insert_primitive;

    /// Iterator pointing at primitive
    config_const_iterator primitive_it;

    /// True if canonical configuration did not exist before insertion
    bool insert_canonical;

    /// Iterator pointing at canonical, if existing
    config_const_iterator canonical_it;

  };

  /// \brief A Configuration represents the values of all degrees of freedom in a Supercell
  ///
  class Configuration : public Comparisons<Configuration> {
  private:

    /// Configuration DFT data is expected in:
    ///   casmroot/supercells/SCEL_NAME/CONFIG_ID/CURR_CALCTYPE/properties.calc.json

    /// POS files are written to:
    ///  casmroot/supercells/SCEL_NAME/CONFIG_ID/POS


    /// Identification

    // Configuration id is the index into Supercell::config_list
    std::string id;

    /// const pointer to the (non-const) Supercell for this Configuration
    Supercell *supercell;

    ///
    std::shared_ptr<Supercell> m_supercell_ptr;

    /// a jsonParser object indicating where this Configuration came from
    jsonParser m_source;
    bool source_updated;


    // symmetric multiplicity (i.e., size of configuration's factor_group)
    int multiplicity;


    /// Degrees of Freedom

    // 'occupation' is a list of the indices describing the occupants in each crystal site.
    //   get_prim().basis[ get_b(i) ].site_occupant[ occupation[i]] -> Molecule on site i
    //   This means that for the background structure, 'occupation' is all 0

    // Configuration sites are arranged by basis, and then prim:
    //   occupation: [basis0                |basis1               |basis2          |...] up to prim.basis.size()
    //       basis0: [prim0|prim1|prim2|...] up to supercell.volume()
    //
    ConfigDoF m_configdof;


    /// Properties
    ///Keeps track of whether the Configuration properties change since reading. Be sure to set to true in your routine if it did!
    /// PROPERTIES (AS OF 07/27/15)
    /*  calculated:
     *    calculated["energy"]
     *    calculated["relaxed_energy"]
     *
     *  generated:
     *    generated["is_groundstate"]
     *	  generated["dist_from_hull"]
     *    generated["sublat_struct_fact"]
     *    generated["struct_fact"]
     */
    bool prop_updated;
    Properties calculated;  //Stuff you got directly from your DFT calculations
    Properties generated;   //Everything else you came up with through casm


    bool m_selected;

    /// Remember how to copy into the canonical Supercell
    mutable notstd::cloneable_ptr<FillSupercell> m_fill_canonical;

    /// Remember name
    mutable std::string m_name;

  public:
    typedef ConfigDoF::displacement_matrix_t displacement_matrix_t;
    typedef ConfigDoF::displacement_t displacement_t;
    typedef ConfigDoF::const_displacement_t const_displacement_t;

    //********* CONSTRUCTORS *********

    /// Construct a default Configuration
    Configuration(Supercell &_supercell, const jsonParser &source = jsonParser(), const ConfigDoF &_dof = ConfigDoF());

    /// Construct by reading from main data file (json)
    Configuration(const jsonParser &json, Supercell &_supercell, Index _id);


    /// Construct a Configuration with occupation specified by string 'con_name'
    //Configuration(Supercell &_supercell, std::string con_name, bool select, const jsonParser &source = jsonParser());


    //********** DESTRUCTORS *********

    //********** MUTATORS  ***********

    /// \brief set symmetric multiplicity (i.e., size of configuration's factor_group)
    void set_multiplicity(int m) {
      multiplicity = m;
    }

    void set_id(Index _id);

    void set_source(const jsonParser &source);

    void push_back_source(const jsonParser &source);

    // ** Degrees of Freedom **
    //
    // ** Note: Properties and correlations are not automatically updated when dof are changed, **
    // **       nor are the written records automatically updated                               **

    /// \brief Clear all DoF
    ///
    /// - This will invalidate the Configuration's id
    void clear();


    /// \brief Set occupant variables to background structure
    ///
    /// - This will invalidate the Configuration's id
    void init_occupation();

    /// \brief Set occupant variables
    ///
    /// With one value for each site in the Configuration, this Array describes
    /// which occupant is at each of the 'N' sites of the configuration. The
    /// occupant on site l can be obtained from the occupation variable using:
    /// \code
    /// Molecule on_site_l = config.get_prim().basis[ config.get_b(l) ].site_occupant[ config.occupation()[l]];
    /// \endcode
    /// - For a CASM project, the occupation variables will be ordered according
    /// to the occupant DoF in a "prim.json" file. This means that for the
    /// background structure, 'occupation' is all 0
    /// - This will invalidate the Configuration's id
    ///
    /// \throws If \code newoccupation.size() != this->size() \endcode
    ///
    void set_occupation(const Array<int> &newoccupation);

    /// \brief Set occupant variable on site l
    ///
    /// The occupant on site l can be obtained from the occupation variable using:
    /// \code
    /// Molecule on_site_l = config.get_prim().basis[ config.get_b(l) ].site_occupant[config.occ(l)];
    /// \endcode
    /// - For a CASM project, the occupation variables will be ordered according
    /// to the occupant DoF in a "prim.json" file. This means that for the
    /// background structure, 'occupation' is all 0
    /// - No check is performed as to whether val is in the correct range
    /// - This will invalidate the Configuration's id
    ///
    void set_occ(Index site_l, int val);

    /// \brief Clear occupation
    ///
    /// - This will invalidate the Configuration's id
    void clear_occupation();


    /// \brief Set all occupant displacements to (0.,0.,0.)
    ///
    /// - This will invalidate the Configuration's id
    void init_displacement();

    /// \brief Set occupant displacements
    ///
    /// A displacement_t vector for each site in the Configuration to describe
    /// displacements condensed in matrix form. This a 3xN matrix whose columns
    /// are the displacement of each of the N sites of the configuration.
    /// - Displacements are applied before strain.
    /// - This will invalidate the Configuration's id
    ///
    /// \throws If \code _disp.cols() != this->size() \endcode
    ///
    void set_displacement(const displacement_matrix_t &_disp);

    /// \brief Set occupant displacements
    ///
    /// - A displacement_t vector to describe displacement of the occupant on site l.
    /// - Displacements are applied before strain.
    /// - This will invalidate the Configuration's id
    ///
    void set_disp(Index site_l, const Eigen::VectorXd &_disp);

    /// \brief Clear displacement
    ///
    /// - This will invalidate the Configuration's id
    void clear_displacement();


    /// \brief Set applied strain to Eigen::Matrix3d::Zero()
    ///
    /// - This will invalidate the Configuration's id
    void init_deformation();

    /// \brief Set applied strain
    ///
    /// Set strain applied applied to the Configuration.
    /// This is the matrix that relates the reference lattice vectors to the
    /// deformed lattice vectors via
    /// \code
    /// L_deformed = m_deformation * L_reference
    /// \endcode
    /// where L is a 3x3 matrix whose columns are the lattice vectors.
    /// - Strain is applied after displacement.
    /// - This will invalidate the Configuration's id
    ///
    void set_deformation(const Eigen::Matrix3d &_deformation);

    /// \brief Clear applied strain
    ///
    /// - This will invalidate the Configuration's id
    void clear_deformation();


    /// \brief Check if this is a primitive Configuration
    bool is_primitive() const;

    /// \brief Returns a PermuteIterator corresponding to the first non-zero pure
    /// translation that maps the Configuration onto itself.
    PermuteIterator find_translation() const;

    /// \brief Return the primitive Configuration
    Configuration primitive() const;


    /// \brief Check if Configuration is in the canonical form
    bool is_canonical() const;

    /// \brief Returns the operation that applied to *this returns the canonical form
    PermuteIterator to_canonical() const;

    /// \brief Returns the operation that applied to the the canonical form returns *this
    PermuteIterator from_canonical() const;

    /// \brief Returns the canonical form Configuration in the same Supercell
    Configuration canonical_form() const;

    /// \brief Returns the canonical form Configuration in the canonical Supercell
    Configuration in_canonical_supercell() const;

    /// \brief Insert this in the canonical Supercell
    ConfigInsertResult insert(bool primitive_only) const;

    /// \brief Returns the subgroup of the Supercell factor group that leaves the
    ///        Configuration unchanged
    std::vector<PermuteIterator> factor_group() const;

    /// \brief Returns the point group that leaves the Configuration unchanged
    SymGroup point_group() const;

    /// \brief Fills supercell 'scel' with reoriented configuration, as if by apply(op,*this)
    Configuration fill_supercell(Supercell &scel, const SymOp &op) const;

    /// \brief Fills supercell 'scel' with reoriented configuration, as if by apply(op,*this)
    Configuration fill_supercell(Supercell &scel, const SymGroup &g) const;


    // ** Properties **
    //
    // ** Note: DeltaProperties are automatically updated, but not written upon changes **

    /// Read calculation results into the configuration
    //void read_calculated();

    void set_calc_properties(const jsonParser &json);

    bool read_calc_properties(jsonParser &parsed_props) const;

    void set_selected(bool _selected) {
      m_selected = _selected;
    }


    //********** ACCESSORS ***********

    const Lattice &ideal_lattice()const;

    std::string get_id() const;

    /// \brief Get symmetric multiplicity (i.e., size of configuration's factor_group)
    ///
    /// - Must first be set via Configuration::set_multiplicity
    int get_multiplicity()const {
      return multiplicity;
    }

    /// \brief SCELV_A_B_C_D_E_F/i
    std::string name() const;

    std::string calc_status() const;

    std::string failure_type() const;

    const jsonParser &source() const;

    fs::path get_path() const;

    ///Returns number of sites, NOT the number of primitives that fit in here
    Index size() const;

    /// \brief Get the primitive Structure for this Configuration
    const Structure &get_prim() const;

    /// \brief True if this Configuration is currently selected in the MASTER config list
    bool selected() const {
      return m_selected;
    }

    /// \brief Get the PrimClex for this Configuration
    PrimClex &get_primclex() const;

    /// \brief Get the Supercell for this Configuration
    Supercell &get_supercell() const;

    /// \brief Get the PrimClex crystallography_tol
    double crystallography_tol() const;

    /// \brief Get the UnitCellCoord for a given linear site index
    UnitCellCoord get_uccoord(Index site_l) const;

    /// \brief Get the basis site index for a given linear linear site index
    int get_b(Index site_l) const;

    /// \brief const Access the DoF
    const ConfigDoF &configdof() const {
      return m_configdof;
    }

    /// \brief Access the DoF
    ///
    /// - This will invalidate the Configuration's id
    ConfigDoF &configdof() {
      _invalidate_id();
      return m_configdof;
    }

    /// \brief True if Configuration has occupation DoF
    bool has_occupation() const {
      return configdof().has_occupation();
    }

    /// \brief Occupant variables
    ///
    /// With one value for each site in the Configuration, this Array describes
    /// which occupant is at each of the 'N' sites of the configuration. The
    /// occupant on site l can be obtained from the occupation variable using:
    /// \code
    /// Molecule on_site_l = config.get_prim().basis[ config.get_b(l) ].site_occupant[ config.occupation()[l]];
    /// \endcode
    /// - For a CASM project, the occupation variables will be ordered according
    /// to the occupant DoF in a "prim.json" file. This means that for the
    /// background structure, 'occupation' is all 0
    ///
    const Array<int> &occupation() const {
      return configdof().occupation();
    }

    /// \brief Occupant variable on site l
    ///
    /// The occupant on site l can be obtained from the occupation variable using:
    /// \code
    /// Molecule on_site_l = config.get_prim().basis[ config.get_b(l) ].site_occupant[config.occ(l)];
    /// \endcode
    /// - For a CASM project, the occupation variables will be ordered according
    /// to the occupant DoF in a "prim.json" file. This means that for the
    /// background structure, 'occupation' is all 0
    ///
    const int &occ(Index site_l) const {
      return configdof().occ(site_l);
    }

    /// \brief Molecule on site l
    ///
    /// Equivalent to:
    /// \code
    /// config.get_prim().basis[ config.get_b(l) ].site_occupant[ config.occupation()[l]];
    /// \endcode
    ///
    const Molecule &get_mol(Index site_l) const;


    /// \brief True if Configuration has displacement DoF
    bool has_displacement() const {
      return configdof().has_displacement();
    }

    /// \brief Occupant displacements
    ///
    /// A displacement_t vector for each site in the Configuration to describe
    /// displacements condensed in matrix form. This a 3xN matrix whose columns
    /// are the displacement of each of the N sites of the configuration.
    /// - Displacements are applied before strain.
    ///
    const displacement_matrix_t &displacement() const {
      return configdof().displacement();
    }

    /// \brief Occupant displacement
    ///
    /// - A displacement_t vector describes displacement of the occupant on site l.
    /// - Displacements are applied before strain.
    ///
    const_displacement_t disp(Index site_l) const {
      return configdof().disp(site_l);
    }

    /// \brief Applied strain
    ///
    /// Describes possible strains that may have been applied to the Configuration.
    /// This is the matrix that relates the reference lattice vectors to the
    /// deformed lattice vectors via
    /// \code
    /// L_deformed = m_deformation * L_reference
    /// \endcode
    /// where L is a 3x3 matrix whose columns are the lattice vectors.
    ///
    /// - Strain is applied after displacement.
    ///
    const Eigen::Matrix3d &deformation() const {
      return configdof().deformation();
    }

    /// \brief True if Configuration has strain DoF
    bool has_deformation() const {
      return configdof().has_deformation();
    }


    //fs::path get_reference_state_dir() const;

    //const Properties &ref_properties() const;

    const Properties &calc_properties() const;

    //const DeltaProperties &delta_properties() const;

    const Properties &generated_properties() const;


    //const Correlation &get_correlations() const;


    // Returns composition on each sublattice: sublat_comp[ prim basis site / sublattice][ molecule_type]
    //   molucule_type is ordered as in the Prim structure's site_occupant list for that basis site (includes vacancies)
    ReturnArray< Array < double > > get_sublattice_composition() const;

    // Returns number of each molecule by sublattice:
    //   sublat_num_each_molecule[ prim basis site / sublattice ][ molecule_type]
    //   molucule_type is ordered as in the Prim structure's site_occupant list for that basis site
    ReturnArray< Array<int> > get_sublat_num_each_molecule() const;

    // Returns composition, not counting vacancies
    //    composition[ molecule_type ]: molecule_type ordered as prim structure's get_struc_molecule(), with [Va]=0.0
    ReturnArray<double> get_composition() const;

    // Returns composition, including vacancies
    //    composition[ molecule_type ]: molecule_type ordered as prim structure's get_struc_molecule()
    ReturnArray<double> get_true_composition() const;

    /// Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::get_struc_molecule()
    ReturnArray<int> get_num_each_molecule() const;

    /// Returns parametric composition, as calculated using PrimClex::param_comp
    Eigen::VectorXd get_param_composition() const;

    /// Returns num_each_component[ component_type] per prim cell,
    ///   where 'component_type' is ordered as ParamComposition::get_components
    Eigen::VectorXd get_num_each_component() const;

    //-----------------------------------
    //Structure Factor
    Eigen::VectorXd get_struct_fact_intensities() const;
    Eigen::VectorXd get_struct_fact_intensities(const Eigen::VectorXd &component_intensities) const;

    void calc_sublat_struct_fact();
    void calc_struct_fact();
    void calc_sublat_struct_fact(const Eigen::VectorXd &intensities);
    void calc_struct_fact(const Eigen::VectorXd &intensities);

    Eigen::MatrixXcd sublat_struct_fact();
    Eigen::MatrixXd struct_fact();

    //********* IO ************

    /// Writes the Configuration to the correct casm directory
    ///   Uses PrimClex's current settings to write the appropriate
    ///   Properties, DeltaProperties and Correlations files
    jsonParser &write(jsonParser &json) const;

    /// Write the POS file to get_pos_path
    void write_pos() const;

    // Va_mode		description
    // 0			print no information about the vacancies
    // 1			print only the coordinates of the vacancies
    // 2			print the number of vacancies and the coordinates of the vacancies
    void print(std::ostream &stream, COORD_TYPE mode, int Va_mode = 0, char term = '\n', int prec = 10, int pad = 5) const;

    void print_occupation(std::ostream &stream) const;

    void print_config_list(std::ostream &stream, int composition_flag) const;

    void print_composition(std::ostream &stream) const;

    void print_true_composition(std::ostream &stream) const;

    void print_sublattice_composition(std::ostream &stream) const;


    fs::path calc_dir() const;
    fs::path calc_properties_path() const;
    fs::path calc_status_path() const;
    /// Path to various files
    fs::path get_pos_path() const;

    /// \brief Check if Configuration are equivalent wrt the prim's factor group
    ///
    /// Equivalent to:
    /// \code
    /// this->primitive() == B.primitive()
    /// \endcode
    /// - Must have the same PrimClex
    ///
    bool is_equivalent(const Configuration &B) const;

    /// \brief Compare Configuration, via ConfigCompare
    ///
    /// - Must have the same Supercell
    bool operator<(const Configuration &B) const;

    /// \brief Split configuration name string into scelname and config index
    static std::pair<std::string, Index> split_name(std::string configname);

  private:
    /// Convenience accessors:
    int &_occ(Index site_l) {
      return m_configdof.occ(site_l);
    }

    displacement_t _disp(Index site_l) {
      return m_configdof.disp(site_l);
    }

    void _invalidate_id() {
      id = "none";
      m_name.clear();
    }

    void _generate_name() const;

    friend Comparisons<Configuration>;

    /// \brief Equality comparison of Configuration, via ConfigEqual
    ///
    /// - Must have the same Supercell
    /// - Checks that all DoF are the same, within tolerance
    bool _eq(const Configuration &B) const;

    /// Reads the Configuration from the expected casm directory
    ///   Uses PrimClex's current settings to read in the appropriate
    ///   Properties, DeltaProperties and Correlations files if they exist
    ///
    /// This is private, because it is only called from the constructor:
    ///   Configuration(const Supercell &_supercell, Index _id)
    ///   It's called from the constructor because of the Supercell pointer
    ///
    void read(const jsonParser &json);

    /// Functions used to perform read()
    void read_dof(const jsonParser &json);
    void read_properties(const jsonParser &json);

    /// Functions used to perform write to config_list.json:
    jsonParser &write_dof(jsonParser &json) const;
    jsonParser &write_source(jsonParser &json) const;
    jsonParser &write_pos(jsonParser &json) const;
    jsonParser &write_param_composition(jsonParser &json) const;
    jsonParser &write_properties(jsonParser &json) const;

    //bool reference_states_exist() const;
    //void read_reference_states(Array<Properties> &ref_state_prop, Array<Eigen::VectorXd> &ref_state_comp) const;
    //void generate_reference_scalar(std::string propname, const Array<Properties> &ref_state_prop, const Array<Eigen::VectorXd> &ref_state_comp);

  };

  // Calculate transformed ConfigDoF from PermuteIterator via
  //   apply(permute_iterator, dof)
  Configuration &apply(const PermuteIterator &it, Configuration &config);

  Configuration sub_configuration(Supercell &sub_scel,
                                  const Configuration &super_config,
                                  const UnitCell &origin = UnitCell(0, 0, 0));

  /// \brief Make Configuration from name string
  Configuration make_configuration(PrimClex &primclex, std::string name);

  /// \brief Returns correlations using 'clexulator'.
  Correlation correlations(const Configuration &config, Clexulator &clexulator);

  /// Returns parametric composition, as calculated using PrimClex::param_comp
  Eigen::VectorXd comp(const Configuration &config);

  /// \brief Returns the composition, as number of each species per unit cell
  Eigen::VectorXd comp_n(const Configuration &config);

  /// \brief Returns the vacancy composition, as number per unit cell
  double n_vacancy(const Configuration &config);

  /// \brief Returns the total number species per unit cell
  double n_species(const Configuration &config);

  /// \brief Returns the composition as species fraction, with [Va] = 0.0, in
  ///        the order of Structure::get_struc_molecule
  Eigen::VectorXd species_frac(const Configuration &config);

  /// \brief Returns the composition as site fraction, in the order of Structure::get_struc_molecule
  Eigen::VectorXd site_frac(const Configuration &config);

  /// \brief Returns the relaxed energy, normalized per unit cell
  double relaxed_energy(const Configuration &config);

  /// \brief Returns the relaxed energy, normalized per species
  double relaxed_energy_per_species(const Configuration &config);

  /// \brief Returns the reference energy, normalized per unit cell
  double reference_energy(const Configuration &config);

  /// \brief Returns the reference energy, normalized per species
  double reference_energy_per_species(const Configuration &config);

  /// \brief Returns the formation energy, normalized per unit cell
  double formation_energy(const Configuration &config);

  /// \brief Returns the formation energy, normalized per species
  double formation_energy_per_species(const Configuration &config);

  /// \brief Returns the formation energy, normalized per unit cell
  double clex_formation_energy(const Configuration &config);

  /// \brief Returns the formation energy, normalized per species
  double clex_formation_energy_per_species(const Configuration &config);

  /// \brief Return true if all current properties have been been calculated for the configuration
  bool is_calculated(const Configuration &config);

  /// \brief Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)
  double rms_force(const Configuration &_config);

  /// \brief Cost function that describes the degree to which basis sites have relaxed
  double basis_deformation(const Configuration &_config);

  /// \brief Cost function that describes the degree to which lattice has relaxed
  double lattice_deformation(const Configuration &_config);

  /// \brief Change in volume due to relaxation, expressed as the ratio V/V_0
  double volume_relaxation(const Configuration &_config);

  /// \brief Returns the relaxed magnetic moment, normalized per unit cell
  double relaxed_magmom(const Configuration &_config);

  /// \brief Returns the relaxed magnetic moment, normalized per species
  double relaxed_magmom_per_species(const Configuration &_config);

  /// \brief Returns the relaxed magnetic moment of each basis site
  Eigen::VectorXd relaxed_mag_basis(const Configuration &_config);
  /* std::vector<double> relaxed_mag_basis(const Configuration &_config); */

  /// \brief Returns the relaxed magnetic moment for each molecule
  Eigen::VectorXd relaxed_mag(const Configuration &_config);



  /// \brief returns true if _config describes primitive cell of the configuration it describes
  bool is_primitive(const Configuration &_config);

  /// \brief returns true if _config no symmetry transformation applied to _config will increase its lexicographic order
  bool is_canonical(const Configuration &_config);

  /// \brief Status of calculation
  inline
  std::string calc_status(const Configuration &_config) {
    return _config.calc_status();
  }

  // \brief Reason for calculation failure.
  inline
  std::string failure_type(const Configuration &_config) {
    return _config.failure_type();
  }

  bool has_relaxed_energy(const Configuration &_config);

  bool has_reference_energy(const Configuration &_config);

  bool has_formation_energy(const Configuration &_config);

  bool has_rms_force(const Configuration &_config);

  bool has_basis_deformation(const Configuration &_config);

  bool has_lattice_deformation(const Configuration &_config);

  bool has_volume_relaxation(const Configuration &_config);

  bool has_relaxed_magmom(const Configuration &_config);

  bool has_relaxed_mag_basis(const Configuration &_config);

  inline
  bool has_calc_status(const Configuration &_config) {
    return !_config.calc_status().empty();
  }

  inline
  bool has_failure_type(const Configuration &_config) {
    return !_config.failure_type().empty();
  }

  class FillSupercell {

  public:

    /// \brief Constructor
    FillSupercell(Supercell &_scel, const SymOp &_op);

    /// \brief Find first SymOp in the prim factor group such that apply(op, motif)
    ///        can be used to fill the Supercell
    FillSupercell(Supercell &_scel, const Configuration &motif, double tol);

    Configuration operator()(const Configuration &motif) const;

    /// \brief Find first SymOp in the prim factor group such that apply(op, motif)
    ///        can be used to fill the Supercell
    const SymOp *find_symop(const Configuration &motif, double tol);

    const SymOp &symop() const {
      return *m_op;
    }


  private:

    void _init(Supercell &_motif_scel) const;

    Supercell *m_scel;
    const SymOp *m_op;

    mutable Supercell *m_motif_scel;
    mutable std::vector<std::vector<Index> > m_index_table;

  };

  std::ostream &operator<<(std::ostream &sout, const Configuration &c);


  inline
  void reset_properties(Configuration &_config) {
    _config.set_calc_properties(jsonParser());
  }

  /** @} */

}

#endif
