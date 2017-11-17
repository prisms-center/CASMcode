#ifndef CASM_ClusterSpecsParser
#define CASM_ClusterSpecsParser

#include <boost/lexical_cast.hpp>

#include "casm/casm_io/InputParser_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"

namespace CASM {

  struct OrbitBranchSpecsParser : public KwargsParser {

    int max_branch;

    OrbitBranchSpecsParser(jsonParser &_input, fs::path _path, bool _required);

    jsonParser::const_iterator branch(int branch_i) const;

  protected:
    int branch_to_int(jsonParser::const_iterator it) const;

    std::string branch_to_string(int branch_i) const;

    bool is_integer(jsonParser::const_iterator it) const;

    /// warn for non-integer 'branch'
    bool warn_non_integer_branch(jsonParser::const_iterator it);

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
  ///     "method": "ClustersByMaxLength",
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
  struct PrimPeriodicClustersByMaxLength : InputParser {

    typedef PrimPeriodicOrbit<IntegralCluster> OrbitType;
    fs::path path;

    PrimPeriodicClustersByMaxLength(
      const PrimClex &_primclex,
      const SymGroup &_generating_grp,
      const PrimPeriodicSymCompare<IntegralCluster> &_sym_compare,
      const jsonParser &_input,
      fs::path _path,
      bool _required);

    int max_branch() const;

    double max_length(int branch_i) const;

    const OrbitGenerators<PrimPeriodicOrbit<IntegralCluster>> &custom_generators() const;

    const PrimPeriodicOrbitBranchSpecsParser &orbit_branch_specs() const;

    const PrimPeriodicOrbitSpecsParser &orbit_specs() const;

  };
}

#endif
