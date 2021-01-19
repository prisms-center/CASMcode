#ifndef CASM_BasisFunctionSpecs
#define CASM_BasisFunctionSpecs

#include <algorithm>
#include <map>
#include <vector>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

// --- Generic BasisFunctionSpecs ---

/// Specify the Clexulator underlying data structure type.
///
/// Choices:
/// - DEFAULT: For BasisClexParamPack, the default.
/// - DIFF: For DiffClexParamPack, to make use of FADBAD++ automatic
/// differentiation library
enum class PARAM_PACK_TYPE { DEFAULT, DIFF };

/// Provides DoF-particular specifications for constructing basis functions
///
/// - This gets specialized (i.e. OccupationDoFSpecs for OccupationDoFTraits)
/// for DoFTraits types
///   that need additional information for constructing basis functions
/// - Stored in BasisFunctionSpecs::dof_specs by DoF type name, so that it is
/// accessible when
///   ClexBasis attempts to construct basis functions for each DoF type
///
class DoFSpecs : public notstd::Cloneable {
  ABSTRACT_CLONEABLE(DoFSpecs)
 public:
  DoFKey name() const;

 private:
  virtual DoFKey _name() const = 0;
};

/// Specify how to construct basis functions
struct BasisFunctionSpecs {
  typedef Index MaxPolyOrder;
  typedef Index OrbitBranchSize;

  /// Construct BasisFunctionSpecs
  ///
  /// \param _dof_keys Names of DoF types to include in basis functions.
  /// \param _dof_specs Provides DoF-particular specifications for constructing
  /// basis functions. \param _global_max_poly_order  Maximum polynomial order
  /// for cluster basis functions of
  ///     continuous DoF. Use -1 for default behaviour: max_poly_order ==
  ///     cluster size.
  /// \param _orbit_branch_max_poly_order  Allows specifying the max_poly_order
  /// per orbit
  ///      branch. The key is the orbit branch index (equal to number of sites
  ///      in a cluster), and the value is the maximum polynomial order for
  ///      orbits in that branch.
  /// \param _include_functions  Allows specifying that only ClexBasis functions
  /// with these
  ///      indices should be written in the produced Clexulator. Generally, this
  ///      requires generating the full ClexBasis once (without writing the
  ///      Clexulator), inspecting the basis functions, and then specifying this
  ///      value and re-generating the ClexBasis.
  /// \param _exclude_functions  Allows specifying that ClexBasis functions with
  /// these
  ///      indices should be not written in the produced Clexulator. Generally,
  ///      this requires generating the full ClexBasis once (without writing the
  ///      Clexulator), inspecting the basis functions, and then specifying this
  ///      value and re-generating the ClexBasis.
  /// \param _param_pack_type  Specify the Clexulator underlying data structure
  /// type.
  ///
  /// Notes on _dof_specs:
  /// - Not all DoF types require their own DoFSpecs. See documentation for a
  /// particular DoFTraits
  ///   class to determine if it is required (i.e. OccupationDoFSpecs is
  ///   required for OccupationDoFTraits) for the "construct_site_bases" method,
  ///   which should throw an error message if required specs are not found.
  /// - Example: OccupationDoFSpecs allows specifying the type of basis
  /// functions to generate for
  ///   occupation cluster functions from among CHEBYCHEV, OCCUPATION,
  ///   COMPOSITION. It also allows specifying the particular sublat composition
  ///   to expand about for COMPOSITION basis functions.
  ///
  /// Notes on _include_functions and _exclude_functions:
  /// - These parameters should be seen as mutually exclusive. Including both
  /// will result in an
  ///   error.
  ///
  /// Notes on _param_pack_type choices:
  /// - DEFAULT: For BasisClexParamPack, the default.
  /// - DIFF: For DiffClexParamPack, to make use of FADBAD++ automatic
  /// differentiation library
  ///
  BasisFunctionSpecs(
      std::vector<DoFKey> _dof_keys = {},
      std::vector<notstd::cloneable_ptr<DoFSpecs>> _dof_specs = {},
      MaxPolyOrder _global_max_poly_order = -1,
      std::map<OrbitBranchSize, MaxPolyOrder> _orbit_branch_max_poly_order = {},
      std::vector<Index> _include_functions = {},
      std::vector<Index> _exclude_functions = {},
      PARAM_PACK_TYPE _param_pack_type = PARAM_PACK_TYPE::DEFAULT);

  /// Which DoF types to include in the basis functions
  std::vector<DoFKey> dof_keys;

  /// Provides DoF-particular specifications for constructing basis functions.
  std::vector<notstd::cloneable_ptr<DoFSpecs>> dof_specs;

  /// Maximum polynomial order for basis functions of continuous DoF.
  /// Use -1 for context-dependent default behaviour (usually max_poly_order ==
  /// orbit branch cluster size) The global value can be overriden per orbit
  /// branch with orbit_branch_max_poly_order
  MaxPolyOrder global_max_poly_order;

  /// Maximum polynomial order for basis functions of continuous DoF, by orbit
  /// branch.
  std::map<OrbitBranchSize, MaxPolyOrder> orbit_branch_max_poly_order;

  // Note: Could allow more customization of max_poly_order, by orbit
  // std::map<IntegralCluster, Index> orbit_max_poly_order;

  /// If not empty, specifies that only ClexBasis functions with these indices
  /// should be included in the produced Clexulator.
  std::vector<Index> include_functions;

  /// If not empty, specifies that ClexBasis functions with these indices should
  /// be excluded from the produced Clexulator.
  std::vector<Index> exclude_functions;

  /// Specify the Clexulator underlying data structure type.
  PARAM_PACK_TYPE param_pack_type;
};

/// Find DoFSpecs and static_cast to DoFSpecsType, or throw exception if not
/// found
template <typename DoFSpecsType>
DoFSpecsType const &get(DoFKey const &key,
                        BasisFunctionSpecs const &basis_function_specs);

}  // namespace CASM

#include <sstream>

#include "casm/global/errors.hh"

namespace CASM {

inline DoFKey DoFSpecs::name() const { return this->_name(); }

inline BasisFunctionSpecs::BasisFunctionSpecs(
    std::vector<DoFKey> _dof_keys,
    std::vector<notstd::cloneable_ptr<DoFSpecs>> _dof_specs,
    MaxPolyOrder _global_max_poly_order,
    std::map<OrbitBranchSize, MaxPolyOrder> _orbit_branch_max_poly_order,
    std::vector<Index> _include_functions,
    std::vector<Index> _exclude_functions, PARAM_PACK_TYPE _param_pack_type)
    : dof_keys(std::move(_dof_keys)),
      dof_specs(std::move(_dof_specs)),
      global_max_poly_order(_global_max_poly_order),
      orbit_branch_max_poly_order(_orbit_branch_max_poly_order),
      include_functions(_include_functions),
      exclude_functions(_exclude_functions),
      param_pack_type(_param_pack_type) {}

template <typename DoFSpecsType>
DoFSpecsType const &get(DoFKey const &key,
                        BasisFunctionSpecs const &basis_function_specs) {
  auto it = std::find_if(basis_function_specs.dof_specs.begin(),
                         basis_function_specs.dof_specs.end(),
                         [&](notstd::cloneable_ptr<DoFSpecs> const &ptr) {
                           return ptr->name() == key;
                         });
  if (it == basis_function_specs.dof_specs.end()) {
    std::stringstream ss;
    ss << "DoFSpecs not found in BasisFunctionSpecs for " << key;
    throw libcasm_runtime_error{ss.str()};
  }
  return static_cast<DoFSpecsType const &>(**it);
}

}  // namespace CASM
#endif
