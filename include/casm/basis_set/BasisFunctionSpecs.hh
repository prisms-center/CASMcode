#ifndef CASM_BasisFunctionSpecs
#define CASM_BasisFunctionSpecs

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
  /// - DIFF: For DiffClexParamPack, to make use of FADBAD++ automatic differentiation library
  enum class PARAM_PACK_TYPE {
    DEFAULT, DIFF
  };

  /// Provides DoF-particular specifications for constructing basis functions
  ///
  /// - This gets specialized (i.e. OccupationDoFSpecs for OccupationDoFTraits) for DoFTraits types
  ///   that need additional information for constructing basis functions
  /// - Stored in BasisFunctionSpecs::dof_specs by DoF type name, so that it is accessible when
  ///   ClexBasis attempts to construct basis functions for each DoF type
  ///
  class DoFSpecs : public notstd::Cloneable {
    ABSTRACT_CLONEABLE(DoFSpecs)
  public:

    DoFKey key() const;

  private:
    virtual DoFKey _key() const = 0;
  };

  /// Specify how to construct basis functions
  ///
  /// Only the default constructor Factory functions for constructing BasisFunctionSpecs
  struct BasisFunctionSpecs {

    /// Construct BasisFunctionSpecs
    ///
    /// \param _dof_keys Names of DoF types to include in basis functions.
    /// \param _dof_specs Provides DoF-particular specifications for constructing basis functions.
    /// \param _max_poly_order  Maximum polynomial order for basis functions of continuous DoF.
    ///     Use -1 for context-dependent default behaviour.
    /// \param _parampack_type  Specify the Clexulator underlying data structure type.
    ///
    /// Notes on _dof_specs_begin, _dof_specs_end:
    /// - Not all DoF types require their own DoFSpecs. See documentation for a particular DoFTraits
    ///   class to determine if it is required (i.e. OccupationDoFSpecs is required for
    ///   OccupationDoFTraits) for the "construct_site_bases" method, which should throw an error
    ///   message if required specs are not found.
    /// - Example: OccupationDoFSpecs allows specifying the type of basis functions to generate for
    ///   occupation cluster functions from among CHEBYCHEV, OCCUPATION, COMPOSITION. It also allows
    ///   specifying the particular sublat composition to expand about for COMPOSITION basis
    ///   functions.
    ///
    /// Notes on _parampack_type choices:
    /// - DEFAULT: For BasisClexParamPack, the default.
    /// - DIFF: For DiffClexParamPack, to make use of FADBAD++ automatic differentiation library
    ///
    BasisFunctionSpecs(
      std::vector<DoFKey> _dof_keys = {},
      std::vector<std::unique_ptr<DoFSpecs>> _dof_specs = {},
      Index _max_poly_order = -1,
      PARAM_PACK_TYPE _parampack_type = PARAM_PACK_TYPE::DEFAULT);


    /// Which DoF types to include in the basis functions
    std::vector<DoFKey> dof_keys;

    /// Provides DoF-particular specifications for constructing basis functions.
    std::vector<std::unique_ptr<DoFSpecs>> dof_specs;

    /// Maximum polynomial order for basis functions of continuous DoF.
    /// Use -1 for context-dependent default behaviour.
    Index max_poly_order;

    /// Specify the Clexulator underlying data structure type.
    PARAM_PACK_TYPE parampack_type;

  };

  /// Find DoFSpecs and static_cast to DoFSpecsType, or throw exception if not found
  template<typename DoFSpecsType>
  DoFSpecsType const &get(DoFKey const &key, BasisFunctionSpecs const &basis_function_specs);

}


#include "casm/global/errors.hh"

namespace CASM {

  inline BasisFunctionSpecs::BasisFunctionSpecs(
    std::vector<DoFKey> _dof_keys,
    std::vector<std::unique_ptr<DoFSpecs>> _dof_specs,
    Index _max_poly_order,
    PARAM_PACK_TYPE _parampack_type):
    dof_keys(std::move(_dof_keys)),
    dof_specs(std::move(_dof_specs)),
    max_poly_order(_max_poly_order),
    parampack_type(_parampack_type) {}

  template<typename DoFSpecsType>
  DoFSpecsType const &get(DoFKey const &key, BasisFunctionSpecs const &basis_function_specs) {
    auto it = std::find(
                basis_function_specs.dof_specs.begin(),
                basis_function_specs.dof_specs.end(),
    [&](std::unique_ptr<DoFSpecs> const & ptr) {
      return ptr->key() == key;
    });
    if(it == basis_function_specs.dof_specs.end()) {
      std::stringstream ss;
      ss << "DoFSpecs not found in BasisFunctionSpecs for " << key;
      throw libcasm_runtime_error {ss.str()};
    }
    return static_cast<DoFSpecsType const &>(**it);
  }

}
#endif
