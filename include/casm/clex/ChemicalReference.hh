#ifndef CASM_ChemicalReference
#define CASM_ChemicalReference

#include "casm/clex/Reference.hh"
#include "casm/misc/algorithm.hh"
#include "casm/clex/ConfigIO.hh"

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
    ///
    /// \param _global_ref An Eigen::VectorXd giving the intercepts of the
    ///                    hyperplane used for the global reference
    /// \param _supercell_ref A map of scelname to Eigen::VectorXd specializing the
    ///                       the reference value by Supercell
    /// \param _config_ref A map of configname to Eigen::VectorXd specializing the
    ///                       the reference value by Configuration
    ///
    /// A hyperplane reference, R, maps vector species_frac, x, to output
    /// energy_per_species, y:
    /// - \code y = R.dot(x) \endcode
    ///
    /// The global reference, '_global_ref', is required, but may be specialized
    /// to give a different R for a particular Supercell or Configuration via
    /// optional '_supercell_ref' and '_config_ref'.
    ///
    explicit ChemicalReference(const Structure &prim,
                               const Eigen::VectorXd &_global_ref,
                               SpecializedRef _supercell_ref = SpecializedRef(),
                               SpecializedRef _config_ref = SpecializedRef()) :
      HyperPlaneReferenceBase(Name, Desc, _global_ref, ConfigIO::SpeciesFrac(), _supercell_ref, _config_ref),
      m_prim(&prim) {}

    /// \brief Construct global reference via range ChemicalReferenceState
    ///
    template<typename RefStateIterator>
    explicit ChemicalReference(const Structure &prim,
                               RefStateIterator begin,
                               RefStateIterator end,
                               double tol) :
      HyperPlaneReferenceBase(Name, Desc, Eigen::VectorXd::Zero(0), ConfigIO::AtomFrac()),
      m_prim(&prim) {
      set_global(begin, end, tol);
    }


    /// \brief Clone
    std::unique_ptr<ChemicalReference> clone() const {
      return notstd::make_unique<ChemicalReference>(*this->_clone());
    }

    /// \brief Get primitive Structure
    const Structure &prim() const {
      return *m_prim;
    }

    // --- Global reference ---

    /// \brief const Access the global reference
    ///
    const Eigen::VectorXd &global() const {
      return HyperPlaneReferenceBase::global();
    }

    /// \brief Set global hyperplane reference
    ///
    /// - Erases associated RefStateVec
    void set_global(const Eigen::VectorXd &ref) {
      _global() = ref;
      m_global_ref_vec.clear();
    }

    /// \brief Set global hyperplane reference
    ///
    /// \param prim The Structure defining the composition space the reference should span
    /// \param begin,end Iterators over a range of ChemicalReferenceState
    /// \param tol Tolerance for checking that input spans the prim composition
    ///            space and a solution for the hyperplane is found
    ///
    /// sets global refrence to be the Eigen::VectorXd, R, that solves:
    /// \code energy = R.dot(atom_frac) \endcode for each ChemicalReferenceState
    ///
    ///
    /// - Inserts associated RefStateVec
    template<typename RefStateIterator>
    void set_global(RefStateIterator begin,
                    RefStateIterator end,
                    double tol) {
      Eigen::VectorXd ref = hyperplane(*m_prim, begin, end, tol);
      m_global_ref_vec = RefStateVec(begin, end);
      _global() = ref;
    }

    /// \brief const Access a map of configname to RefStateVec for Supercell
    ///        specialized references
    ///
    /// - There may not be global reference states (maybe only the hyperplane is
    ///   known), in which case this is 'empty' / has 'size() == 0'
    const RefStateVec &global_ref_states() const {
      return m_global_ref_vec;
    }


    // --- Supercell specialized references ---

    /// \brief const Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    const std::map<std::string, Eigen::VectorXd> &supercell() const {
      return HyperPlaneReferenceBase::supercell();
    }

    /// \brief Set hyperplane reference specialized for a Supercell
    ///
    /// - Erases associated RefStateVec
    void set_supercell(const std::string &scelname, const Eigen::VectorXd &ref) {
      _supercell()[scelname] = ref;
      m_supercell_ref_map.erase(scelname);
    }

    /// \brief Set hyperplane reference specialized for a Supercell
    ///
    /// - Inserts associated RefStateVec
    template<typename RefStateIterator>
    void set_supercell(const std::string &scelname,
                       RefStateIterator begin,
                       RefStateIterator end,
                       double tol) {
      Eigen::VectorXd ref = hyperplane(*m_prim, begin, end, tol);
      m_supercell_ref_map[scelname] = RefStateVec(begin, end);
      _supercell()[scelname] = ref;
    }

    /// \brief Erase hyperplane reference specialized for a Supercell
    size_type erase_supercell(const std::string &scelname) {
      auto result = _supercell().erase(scelname);
      m_supercell_ref_map.erase(scelname);
      return result;
    }

    /// \brief const Access a map of configname to RefStateVec for Supercell
    ///        specialized references
    ///
    /// - A configuration with a specialized reference need not have an associated
    ///   RefStateVec
    const RefStateMap &supercell_ref_states() const {
      return m_supercell_ref_map;
    }


    // --- Configuration specialized references ---


    /// \brief const Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    const std::map<std::string, Eigen::VectorXd> &config() const {
      return HyperPlaneReferenceBase::config();
    }

    /// \brief Set hyperplane reference specialized for a Configuration
    ///
    /// - Erases associated RefStateVec
    void set_config(const std::string &configname, const Eigen::VectorXd &ref) {
      _config()[configname] = ref;
      m_config_ref_map.erase(configname);
    }

    /// \brief Set hyperplane reference specialized for a Configuration
    ///
    /// - Inserts associated RefStateVec
    template<typename RefStateIterator>
    void set_config(const std::string &configname,
                    RefStateIterator begin,
                    RefStateIterator end,
                    double tol) {
      Eigen::VectorXd ref = hyperplane(*m_prim, begin, end, tol);
      m_config_ref_map[configname] = RefStateVec(begin, end);
      _config()[configname] = ref;
    }

    /// \brief Erase hyperplane reference specialized for a Configuration
    size_type erase_config(const std::string &configname) {
      auto result = _config().erase(configname);
      m_config_ref_map.erase(configname);
      return result;
    }

    /// \brief const Access a map of configname to RefStateVec for Configuration
    ///        specialized references
    ///
    /// - A configuration with a specialized reference need not have an associated
    ///   RefStateVec
    const RefStateMap &config_ref_states() const {
      return m_config_ref_map;
    }


    /// \brief Convert a set of ChemicalReferenceState to a hyperplane, including checks
    ///
    /// \param prim The Structure defining the composition space the reference should span
    /// \param begin,end Iterators over a range of ChemicalReferenceState
    /// \param tol Tolerance for checking that input spans the prim composition
    ///            space and a solution for the hyperplane is found
    ///
    /// \returns Eigen::VectorXd, R, that solves: energy = R.dot(atom_frac) for
    ///          each ChemicalReferenceState
    ///
    template<typename RefStateIterator>
    static Eigen::VectorXd hyperplane(const Structure &prim,
                                      RefStateIterator begin,
                                      RefStateIterator end,
                                      double tol);

  private:

    /// \brief Access the global reference
    ///
    Eigen::VectorXd &_global() {
      return HyperPlaneReferenceBase::global();
    }

    /// \brief const Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &_supercell() {
      return HyperPlaneReferenceBase::supercell();
    }

    /// \brief const Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &_config() {
      return HyperPlaneReferenceBase::config();
    }


    // --- For use in hyperplane() ------

    /// \brief Convert a set of ChemicalReferenceState to a hyperplane, including checks
    static Eigen::VectorXd _calc_hyperplane(
      const Structure &prim,
      const std::vector<std::string> &struc_mol_name,
      Eigen::MatrixXd N,
      Eigen::VectorXd E,
      double tol);

    /// \brief Clone
    ChemicalReference *_clone() const {
      return new ChemicalReference(*this);
    }

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
                             int _indent_incr = 2) :
      stream(_stream),
      indent(_indent),
      indent_incr(_indent_incr),
      ref(_ref),
      struc_mol_name(ref.prim().get_struc_molecule_name()) {}

    // -- data --
    std::ostream &stream;
    int indent;
    int indent_incr;
    const ChemicalReference &ref;
    std::vector<std::string> struc_mol_name;


    void incr() {
      indent += indent_incr;
    }

    void decr() {
      indent -= indent_incr;
    }

    // print regular string
    void print(const std::string &str) {
      stream << std::string(indent, ' ') << str << "\n";
    }

    // print plane as '<indent>mol_i(1): energy_per_species\n', for each molecule
    void print(const Eigen::VectorXd &plane) {
      for(int i = 0; i < plane.size(); ++i) {
        if(!is_vacancy(struc_mol_name[i])) {
          stream << std::string(indent, ' ') << struc_mol_name[i] << "(1): " << plane(i) << "\n";
        }
      }
    }

    // print plane as:
    // \code
    // <indent>mol_i(N_i)mol_j(N_j)...: energy_per_species //for each ref state
    // ...
    // \endcode
    void print(const std::vector<ChemicalReferenceState> &ref_state_vec) {
      for(auto it = ref_state_vec.begin(); it != ref_state_vec.end(); ++it) {

        const ChemicalReferenceState &ref_state = *it;

        stream << std::string(indent, ' ');
        for(auto it = ref_state.species_num.begin(); it != ref_state.species_num.end(); ++it) {
          double num = it->second;
          if(almost_zero(num, 1e-14)) {
            continue;
          }

          stream << it->first << "(";
          // print integer if ~integer
          if(almost_equal(std::round(num), num, 1e-14)) {
            stream << std::lround(num);
          }
          else {
            stream << num;
          }
          stream << ")";

        }
        stream << ": " << ref_state.energy_per_species << "\n";
      }
    }

    // print supercell/config specific plane as:
    // \code
    // <indent>name:
    // <indent>  mol_i: energy_per_species // for each molecule
    // \endcode
    void print(const std::pair<std::string, Eigen::VectorXd> &_pair) {
      stream << std::string(indent, ' ') << _pair.first << ":\n";
      incr();
      print(_pair.second);
      decr();
    }

    // print supercell/config specific ref states as:
    // \code
    // <indent>name:
    // <indent>  mol_i(N_i)mol_j(N_j)...: energy_per_species // for each ref state
    // \endcode
    void print(const std::pair<std::string, std::vector<ChemicalReferenceState> > &_pair) {
      stream << std::string(indent, ' ') << _pair.first << ":\n";
      incr();
      print(_pair.second);
      decr();
    }

    void print_global() {
      print("Global chemical reference:");
      incr();
      (ref.global_ref_states().empty()) ? print(ref.global()) : print(ref.global_ref_states());
      decr();
    }

    void print_supercell() {
      if(ref.supercell().size()) {
        print("Supercell specific chemical references:");
        for(auto it = ref.supercell().begin(); it != ref.supercell().end(); ++it) {
          print_supercell(it->first);
        }
      }
    }

    void print_supercell(const std::string &name) {
      auto it = ref.supercell().find(name);
      auto res = ref.supercell_ref_states().find(name);
      auto ref_state_end = ref.supercell_ref_states().end();
      (res != ref_state_end) ? print(*res) : print(*it);
    }

    void print_config() {
      if(ref.config().size()) {
        print("Configuration specific chemical references:");
        for(auto it = ref.config().begin(); it != ref.config().end(); ++it) {
          print_config(it->first);
        }
      }
    }

    void print_config(const std::string &name) {
      auto it = ref.config().find(name);
      auto res = ref.config_ref_states().find(name);
      auto ref_state_end = ref.config_ref_states().end();
      (res != ref_state_end) ? print(*res) : print(*it);
    }

    void print_all() {
      print_global();
      if(ref.supercell().size()) {
        print("");
        print_supercell();
      }
      if(ref.config().size()) {
        print("");
        print_config();
      }
    }
  };


  /// \brief Convert a set of ChemicalReferenceState to a hyperplane, including checks
  ///
  /// \param prim The Structure defining the composition space the reference should span
  /// \param begin,end Iterators over a range of ChemicalReferenceState
  /// \param tol Tolerance for checking that input spans the prim composition
  ///            space and a solution for the hyperplane is found
  ///
  /// \returns Eigen::VectorXd, R, that solves: energy = R.dot(atom_frac) for
  ///          each ChemicalReferenceState
  ///
  template<typename RefStateIterator>
  Eigen::VectorXd ChemicalReference::hyperplane(const Structure &prim,
                                                RefStateIterator begin,
                                                RefStateIterator end,
                                                double tol) {

    // store Molecule names in vector
    std::vector<std::string> struc_mol_name = prim.get_struc_molecule_name();

    // --- find any Molecule not in the prim, add to end of vector -------------

    // increase struc_mol_name to include all Molecule names in input
    // ensure no vacancies included
    for(auto it = begin; it != end; ++it) {
      for(auto mol_it = it->species_num.begin(); mol_it != it->species_num.end(); ++mol_it) {
        if(is_vacancy(mol_it->first)) {
          throw std::runtime_error("Error in ChemicalReference::hyperplane: Input should not include vacancies");
        }
        if(!contains(struc_mol_name, mol_it->first)) {
          struc_mol_name.push_back(mol_it->first);
        }
      }
    }

    // --- initialize vectors, matrices ----------------------------------------

    // reference 'relaxed_energy_per_species'
    Eigen::VectorXd E = Eigen::VectorXd::Zero(std::distance(begin, end));

    // column vector matrix of number of each Molecule in each reference state
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(struc_mol_name.size(), std::distance(begin, end));


    // --- get input values ---------------------------------------------------

    // populate E, N
    Index index = 0;
    for(auto it = begin; it != end; ++it, ++index) {
      E(index) = it->energy_per_species;
      for(auto mol_it = it->species_num.begin(); mol_it != it->species_num.end(); ++mol_it) {
        N(find_index(struc_mol_name, mol_it->first), index) = mol_it->second;
      }
    }

    // use E, N to calculate hyperplane
    return _calc_hyperplane(prim, struc_mol_name, N, E, tol);
  }

  /** @} */
}

#endif
