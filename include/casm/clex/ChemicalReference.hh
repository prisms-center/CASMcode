#ifndef CASM_ChemicalReference
#define CASM_ChemicalReference

#include "casm/clex/Reference.hh"

namespace CASM {

  class PrimClex;
  class Structure;

  /** \ingroup Reference
   *
   *  @{
   */

  /// \brief Stores the composition and energy in a single reference state
  ///
  /// - Should not include Va
  struct ChemicalReferenceState {

    ChemicalReferenceState() {}

    /// \brief Construct using the results of n(config), and e(config)
    ChemicalReferenceState(const Configuration &config,
                           std::function<Eigen::VectorXd(Configuration)> n,
                           std::function<double(Configuration)> e);

    /// Map of Molecule name : number of each species in reference state
    std::map<std::string, double> species_num;

    /// Energy in this reference state
    double energy_per_species;
  };


  class ChemicalReference : public HyperPlaneReferenceBase {

  public:

    typedef std::vector<ChemicalReferenceState> RefStateVec;
    typedef std::map<std::string, RefStateVec> RefStateMap;
    typedef Index size_type;

    static const std::string Name;
    static const std::string Desc;


    /// \brief Constructor
    explicit ChemicalReference(const Structure &prim,
                               const Eigen::VectorXd &_global_ref,
                               SpecializedRef _supercell_ref = SpecializedRef(),
                               SpecializedRef _config_ref = SpecializedRef());

    /// \brief Construct global reference via range ChemicalReferenceState
    ///
    template<typename RefStateIterator>
    explicit ChemicalReference(const Structure &prim,
                               RefStateIterator begin,
                               RefStateIterator end,
                               double tol);


    /// \brief Clone
    std::unique_ptr<ChemicalReference> clone() const;

    /// \brief Get primitive Structure
    const Structure &prim() const;

    // --- Global reference ---

    /// \brief const Access the global reference
    ///
    const Eigen::VectorXd &global() const;

    /// \brief Set global hyperplane reference
    void set_global(const Eigen::VectorXd &ref);

    /// \brief Set global hyperplane reference
    template<typename RefStateIterator>
    void set_global(RefStateIterator begin,
                    RefStateIterator end,
                    double tol);

    /// \brief const Access a map of configname to RefStateVec for Supercell
    ///        specialized references
    const RefStateVec &global_ref_states() const;


    // --- Supercell specialized references ---

    /// \brief const Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    const std::map<std::string, Eigen::VectorXd> &supercell() const;

    /// \brief Set hyperplane reference specialized for a Supercell
    void set_supercell(const std::string &scelname, const Eigen::VectorXd &ref);

    /// \brief Set hyperplane reference specialized for a Supercell
    template<typename RefStateIterator>
    void set_supercell(const std::string &scelname,
                       RefStateIterator begin,
                       RefStateIterator end,
                       double tol);

    /// \brief Erase hyperplane reference specialized for a Supercell
    size_type erase_supercell(const std::string &scelname);

    /// \brief const Access a map of configname to RefStateVec for Supercell
    ///        specialized references
    const RefStateMap &supercell_ref_states() const;


    // --- Configuration specialized references ---


    /// \brief const Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    const std::map<std::string, Eigen::VectorXd> &config() const;

    /// \brief Set hyperplane reference specialized for a Configuration
    void set_config(const std::string &configname, const Eigen::VectorXd &ref);

    /// \brief Set hyperplane reference specialized for a Configuration
    template<typename RefStateIterator>
    void set_config(const std::string &configname,
                    RefStateIterator begin,
                    RefStateIterator end,
                    double tol);

    /// \brief Erase hyperplane reference specialized for a Configuration
    size_type erase_config(const std::string &configname);

    /// \brief const Access a map of configname to RefStateVec for Configuration
    ///        specialized references
    const RefStateMap &config_ref_states() const;


    /// \brief Convert a set of ChemicalReferenceState to a hyperplane, including checks
    template<typename RefStateIterator>
    static Eigen::VectorXd hyperplane(const Structure &prim,
                                      RefStateIterator begin,
                                      RefStateIterator end,
                                      double tol);

  private:

    /// \brief Access the global reference
    ///
    Eigen::VectorXd &_global();

    /// \brief const Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &_supercell();

    /// \brief const Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &_config();


    // --- For use in hyperplane() ------

    /// \brief Convert a set of ChemicalReferenceState to a hyperplane, including checks
    static Eigen::VectorXd _calc_hyperplane(
      const Structure &prim,
      const std::vector<std::string> &struc_mol_name,
      Eigen::MatrixXd N,
      Eigen::VectorXd E,
      double tol);

    /// \brief Clone
    ChemicalReference *_clone() const;

    // \brief non-owning pointer to const primitive Structure
    const Structure *m_prim;


    // --- Store ChemicalReferenceState if known ----

    RefStateVec m_global_ref_vec;
    RefStateMap m_supercell_ref_map;
    RefStateMap m_config_ref_map;
  };

  /// \brief Automatically set ChemicalReference using calculated Configurations
  ///        with 'extreme' compositions
  ChemicalReference auto_chemical_reference(const PrimClex &primclex, double lin_alg_tol);

  /// \brief Structure to help print ChemicalReference
  struct ChemicalReferencePrinter {

    // -- constructor --
    ChemicalReferencePrinter(std::ostream &_stream,
                             const ChemicalReference &_ref,
                             int _indent = 0,
                             int _indent_incr = 2);

    // -- data --
    std::ostream &stream;
    int indent;
    int indent_incr;
    const ChemicalReference &ref;
    std::vector<std::string> struc_mol_name;


    void incr();

    void decr();

    // print regular string
    void print(const std::string &str);

    // print plane as '<indent>mol_i(1): energy_per_species\n', for each molecule
    void print(const Eigen::VectorXd &plane);

    // print plane as:
    // \code
    // <indent>mol_i(N_i)mol_j(N_j)...: energy_per_species //for each ref state
    // ...
    // \endcode
    void print(const std::vector<ChemicalReferenceState> &ref_state_vec);

    // print supercell/config specific plane as:
    // \code
    // <indent>name:
    // <indent>  mol_i: energy_per_species // for each molecule
    // \endcode
    void print(const std::pair<std::string, Eigen::VectorXd> &_pair);

    // print supercell/config specific ref states as:
    // \code
    // <indent>name:
    // <indent>  mol_i(N_i)mol_j(N_j)...: energy_per_species // for each ref state
    // \endcode
    void print(const std::pair<std::string, std::vector<ChemicalReferenceState> > &_pair);

    void print_global();

    void print_supercell();

    void print_supercell(const std::string &name);

    void print_config();

    void print_config(const std::string &name);

    void print_all();
  };

  /** @} */
}

#endif
