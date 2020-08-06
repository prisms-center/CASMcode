#ifndef CASM_ClusterSpecsParser
#define CASM_ClusterSpecsParser

#include <boost/lexical_cast.hpp>

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"

namespace CASM {

  struct OrbitBranchSpecsParser : public KwargsParser {

    int max_branch;

    OrbitBranchSpecsParser(jsonParser &_input, fs::path _path, bool _required);

    jsonParser::const_iterator branch(int branch_i) const;

  protected:
    int branch_to_int(jsonParser::const_iterator it) const;

    std::string branch_to_string(int branch_i) const;

    bool is_integer(jsonParser::const_iterator it) const;

    /// warn for unexpected non-integer 'branch'
    bool warn_non_integer_branch(jsonParser::const_iterator it, const std::set<std::string> &expected = {});

    /// add warning if option is unnecessary
    bool warn_unnecessary(jsonParser::const_iterator it, const std::set<std::string> &expected);

    /// add warning if branch is unnecessary
    bool warn_unnecessary_branch(jsonParser::const_iterator it);

    jsonParser::const_iterator previous(jsonParser::const_iterator it);

    template<typename RequiredType>
    bool require_previous(jsonParser::const_iterator it, std::string option);

    template<typename RequiredType>
    bool require_nonincreasing(jsonParser::const_iterator it, std::string option);

    int max_orbit_branch() const;

  };

  /// Component used by ClusterSpecs parsers
  struct PrimPeriodicOrbitBranchSpecsParser : OrbitBranchSpecsParser {
    PrimPeriodicOrbitBranchSpecsParser(jsonParser &_input, fs::path _path, bool _required);

    double max_length(int branch_i) const;
  };

  /// Component used by ClusterSpecs parsers
  ///
  /// Notes:
  /// - "coordinate_mode" (COORD_TYPE): (optional, default="Integral")
  ///     Specifies how the cluster sites are represented. May be
  ///     "Integral" (default), "Direct" or "Fractional", or "Cartesian"
  /// - "prototype" or "sites" (2d array):
  ///     Array of sites in the prototype cluster.
  /// - "include_subclusters": (optional, default=true)
  ///     Also include orbits of subclusters of the specified prototype cluster
  /// Example:
  /// \code
  /// "orbit_specs": [
  ///   {
  ///     "coordinate_mode" : "Direct",
  ///     "prototype" : [
  ///       [ 0.000000000000, 0.000000000000, 0.000000000000 ],
  ///       [ 1.000000000000, 0.000000000000, 0.000000000000 ],
  ///       [ 2.000000000000, 0.000000000000, 0.000000000000 ],
  ///       [ 3.000000000000, 0.000000000000, 0.000000000000 ]],
  ///     "include_subclusters" : true
  ///   },
  ///   ...
  /// ]
  /// \endcode
  ///
  struct PrimPeriodicOrbitSpecsParser : KwargsParser {

    const PrimClex &primclex;
    typedef PrimPeriodicOrbit<IntegralCluster> OrbitType;
    OrbitGenerators<OrbitType> custom_generators;

    PrimPeriodicOrbitSpecsParser(
      const PrimClex &_primclex,
      const SymGroup &_generating_grp,
      const PrimPeriodicSymCompare<IntegralCluster> &_sym_compare,
      jsonParser &_input,
      fs::path _path,
      bool _required);
  };

  /// Checks:
  /// - no minimum orbit_branch_specs, branch 1 is always assumed
  /// - warn for non-integer orbit_branch_specs
  /// - error if missing any in range [2, max(branch)]
  /// - for branch=1: ignore with warning
  /// - for branch 2+: max_length required,
  ///      max_length must be <= max_length for branch-1
  /// - that 'orbit_specs' are readable
  ///
  /// Usage:
  /// - get options ("max_length", "cutoff_radious"):
  ///   - template<RequiredType> RequiredType get_option(int branch_i, std::string option)
  /// - get custom orbit generator:
  ///   -
  ///
  /// \code
  /// {
  ///   "cluster_specs": {
  ///     "method": "PrimPeriodicClustersByMaxLength",
  ///     "kwargs": {
  ///       "orbit_branch_specs": {
  ///         "2": {
  ///           "max_length": 6.0
  ///         },
  ///         "3": {
  ///           "max_length": 6.0
  ///         }
  ///       },
  ///       "orbit_specs": [
  ///         {
  ///           "coordinate_mode" : "Direct",
  ///           "prototype" : [
  ///             [ 0.000000000000, 0.000000000000, 0.000000000000 ],
  ///             [ 1.000000000000, 0.000000000000, 0.000000000000 ],
  ///             [ 2.000000000000, 0.000000000000, 0.000000000000 ],
  ///             [ 3.000000000000, 0.000000000000, 0.000000000000 ]],
  ///           "include_subclusters" : true
  ///         },
  ///         ...
  ///       ]
  ///     }
  ///   }
  /// }
  /// \endcode
  struct PrimPeriodicClustersByMaxLength : InputParser<std::nullptr_t> {

    typedef PrimPeriodicOrbit<IntegralCluster> OrbitType;

    static std::string cspecs_help();

    PrimPeriodicClustersByMaxLength(
      const PrimClex &_primclex,
      const SymGroup &_generating_grp,
      const PrimPeriodicSymCompare<IntegralCluster> &_sym_compare,
      jsonParser &_input,
      fs::path _path,
      bool _required);

    int max_branch() const;

    double max_length(int branch_i) const;

    const OrbitGenerators<PrimPeriodicOrbit<IntegralCluster>> &custom_generators() const;

    const PrimPeriodicOrbitBranchSpecsParser &orbit_branch_specs() const;

    const PrimPeriodicOrbitSpecsParser &orbit_specs() const;

  };


  /// --- LocalClustersByMaxLength ---

  struct LocalOrbitBranchSpecsParser : OrbitBranchSpecsParser {
    LocalOrbitBranchSpecsParser(jsonParser &_input, fs::path _path, bool _required);

    bool max_length_including_phenomenal() const;
    double max_length(int branch_i) const;
    double cutoff_radius(int branch_i) const;
  };

  struct LocalOrbitSpecsParser : KwargsParser {

    struct Data {
      IntegralCluster cluster;
      bool include_subclusters;
    };

    const PrimClex &primclex;
    std::vector<Data> prototypes;

    LocalOrbitSpecsParser(const PrimClex &_primclex, jsonParser &_input, fs::path _path, bool _required);

    /// For all custom clusters, insert 'op*prototype' into custom_generators
    ///
    /// \note Use ClusterEquivalenceParser to determine if custom clusters apply
    /// to a given cluster, and determine 'op'.
    ///
    template<typename OrbitType, typename CopyApplySymToElementType>
    OrbitGenerators<OrbitType> &insert_custom_generators(
      const SymOp &op,
      OrbitGenerators<OrbitType> &custom_generators,
      const CopyApplySymToElementType &copy_apply) const;
  };

  /// \brief Specify equivalence type for customizing input based on a particular phenomenal cluster
  ///
  /// Parses input of form:
  /// \code
  /// {
  ///   "coordinate_mode": "<>", // coordinate mode for 'sites'
  ///   "sites": [...], // phenomenal cluster (or whatever is expected for PhenomenalType)
  ///   "equivalence_type": "prim", // one of "prim" (default), "scel", "config"
  ///   "configname": "<configname>" // required if equivalence_type is "config" (must have appropriate Supercell for context)
  /// }
  ///
  /// - For equiv_type == EQUIVALENCE_TYPE::PRIM, use prim factor group
  /// - For equiv_type == EQUIVALENCE_TYPE::SCEL, use factor group of relevant Supercell
  /// - For equiv_type == EQUIVALENCE_TYPE::CONFIG, use factor group of Configuration named
  ///   by 'configname' (must have same Supercell as given to constructor)
  ///
  template<typename PhenomenalType>
  struct ClusterEquivalenceParser : KwargsParser {

    const PrimClex &primclex;
    const Supercell &scel;
    std::unique_ptr<PhenomenalType> phenom;
    EQUIVALENCE_TYPE equiv_type;

    // for EQUIVALENCE_TYPE::PRIM
    std::unique_ptr<PrimPeriodicSymCompare<PhenomenalType>> prim_sym_compare;

    // for EQUIVALENCE_TYPE::SCEL
    std::unique_ptr<ScelPeriodicSymCompare<PhenomenalType>> scel_sym_compare;

    // for EQUIVALENCE_TYPE::CONFIG
    std::string configname;
    std::unique_ptr<ScelPeriodicSymCompare<PhenomenalType>> config_sym_compare;
    std::vector<PermuteIterator> config_fg;


    ClusterEquivalenceParser(const Supercell &_scel, jsonParser &_input, fs::path _path, bool _required);

    /// \brief Check if test is equivalent to phenom
    ///
    /// \returns (phenom.apply_sym(op) is_equivalent to test, op)
    ///
    /// For equiv_type == EQUIVALENCE_TYPE::PRIM, use prim factor group
    /// For equiv_type == EQUIVALENCE_TYPE::SCEL, use factor group of Supercell
    /// For equiv_type == EQUIVALENCE_TYPE::CONFIG, use factor group of Configuration named
    ///   by 'configname' (must have same Supercell as given to constructor)
    ///
    std::pair<bool, SymOp> is_equivalent(const PhenomenalType &test) const;

  private:

    void _init_prim_equivalence();

    void _init_scel_equivalence();

    void _init_config_equivalence();

    template<typename SymCompareType>
    std::pair<bool, SymOp> _is_equivalent(
      const PhenomenalType &test,
      const SymGroup &g,
      const SymCompareType &sym_compare) const;

    template<typename SymCompareType>
    std::pair<bool, SymOp> _is_equivalent(
      const PhenomenalType &test,
      const SymCompareType &sym_compare,
      PermuteIterator begin,
      PermuteIterator end) const;

    template<typename SymCompareType, typename PermuteIteratorIt>
    std::pair<bool, SymOp> _is_equivalent(
      const PhenomenalType &test,
      const SymCompareType &sym_compare,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;
  };

  template<typename PhenomenalType>
  struct CustomLocalClustersByMaxLength : InputParser<std::nullptr_t> {

    struct CustomSpecs {
      std::shared_ptr<ClusterEquivalenceParser<PhenomenalType>> phenom;
      std::shared_ptr<LocalOrbitBranchSpecsParser> orbit_branch_specs;
      std::shared_ptr<LocalOrbitSpecsParser> orbit_specs;
    };
    std::vector<CustomSpecs> data;
    typedef typename std::vector<CustomSpecs>::const_iterator custom_specs_iterator;

    CustomLocalClustersByMaxLength(const Supercell &_scel, jsonParser &_input, fs::path _path, bool _required);

    /// Find if phenom is equivalent to one of the custom phenomenal clusters
    std::pair<custom_specs_iterator, SymOp> find(const PhenomenalType &phenom) const;

  };

  /// Specs for generating local clusters
  ///
  /// - Is supercell specific, all local orbits are generated with a subgroup
  ///   of the supercell factor group
  /// - "standard" orbit_branch_specs give the default generating specs
  /// - "custom" specs provide a way to specify orbit branch and orbit specs
  ///   for particular phenomenal clusters
  /// - For any particular test phenomenal cluster, all custom specs are checked
  ///   to see if one of the custom phenomenal clusters is equivalent (equivalence
  ///   may be by prim, scel, or config symmetry, see ClusterEquivalenceParser).
  /// - If a custom phenomenal cluster is equivalent, then the specified custom
  ///   prototype clusters are transformed by the same symmetry op and used as a
  ///   generator for local clusters:
  ///     op*custom_phenomneal = test  -->  op*custom_prototype = generator_cluster
  ///
  /// Checks:
  /// - no minimum orbit_branch_specs
  /// - warn for non-integer orbit_branch_specs
  /// - error if missing any in range [2, max(branch)]
  /// - for branch 1:
  ///   - max_length required if "max_length_including_phenomenal", else ignore w/ warning
  ///     and throw if not given and requested
  ///   - cutoff_radius required
  /// - for branch 2+: max_length & cutoff_radius required,
  ///   - max_length must be <= max_length for branch-1
  ///   - cutoff_radious must be <= cutoff_radius for branch-1
  ///
  /// - ClusterFilter: max_length, include_phenomenal
  ///
  /// \code
  /// {
  ///   "cluster_specs": {
  ///     "method": "LocalClustersByMaxLength",
  ///     "kwargs": {
  ///       "standard": {
  ///         "orbit_branch_specs": {
  ///           "max_length_including_phenomenal": true, // default==false
  ///           "1": {
  ///             "max_length": 12.0, // max length in cluster
  ///             "cutoff_radius": 6.0 // cutoff 'dist_to_path' for neighborhood
  ///           },
  ///           "2": {
  ///             "max_length": 12.0,
  ///             "cutoff_radius": 6.0
  ///           },
  ///           "3": {
  ///             "max_length": 12.0,
  ///             "cutoff_radius": 6.0
  ///           }
  ///         }
  ///       },
  ///       "custom": [
  ///         {
  ///           "phenomenal": {
  ///             "coordinate_mode": "<>",
  ///             "sites": [...],
  ///             "equivalence_type": "prim", // "scel", "config"
  ///             "configname": "<configname>"
  ///           },
  ///           "orbit_branch_specs": {}
  ///           "orbit_specs": [...(coordinates relative to phenom)...]
  ///         },
  ///         ...
  ///       ]
  ///     }
  ///   }
  /// }
  /// \endcode
  ///
  template<typename PhenomenalType>
  struct LocalClustersByMaxLength : InputParser<std::nullptr_t> {

    std::shared_ptr<LocalOrbitBranchSpecsParser> standard;
    typedef CustomLocalClustersByMaxLength<PhenomenalType> CustomParserType;
    std::shared_ptr<CustomParserType> custom;
    typedef typename CustomParserType::custom_specs_iterator custom_specs_iterator;


    LocalClustersByMaxLength(const Supercell &_scel, jsonParser &_input, fs::path _path, bool _required);

    /// Find if phenom is equivalent to one of the custom phenomenal clusters
    std::pair<custom_specs_iterator, SymOp> find(const PhenomenalType &phenom) const;

    bool max_length_including_phenomenal(custom_specs_iterator it) const;

    int max_branch(custom_specs_iterator it) const;

    double max_length(custom_specs_iterator it, int branch_i) const;

    double cutoff_radius(custom_specs_iterator it, int branch_i) const;

    // *INDENT-OFF*

    /// Get custom local cluster generators
    ///
    /// \throws std::invalid_argument if find_res.first == custom.end()
    template<typename OrbitType, typename CopyApplySymToElementType>
    OrbitGenerators<OrbitType> &
    insert_custom_generators(
      std::pair<custom_specs_iterator, SymOp> find_res,
      OrbitGenerators<OrbitType> &custom_generators,
      const CopyApplySymToElementType& copy_apply) const;

    // *INDENT-ON*
  };

}

#endif
