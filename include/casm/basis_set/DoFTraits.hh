#ifndef CASM_DoFTraits
#define CASM_DoFTraits

#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/symmetry/OrbitDecl.hh"

namespace CASM {
namespace xtal {
class Site;
class BasicStructure;
class SimpleStructure;
class UnitCellCoord;
struct SymOp;
}  // namespace xtal
using xtal::BasicStructure;
using xtal::SimpleStructure;
using xtal::Site;
using xtal::UnitCellCoord;

namespace clexulator {
class PrimNeighborList;
}
using clexulator::PrimNeighborList;

template <typename T>
class InputParser;
class jsonParser;
struct BasisFunctionSpecs;
class BasisSet;
class Structure;

template <typename T>
class ParsingDictionary;

class ConfigDoF;

namespace DoFType {
class Traits;
struct ParamAllocation;

/// Parsing dictionary for obtaining the correct DoFType::Traits given a name
using TraitsDictionary = ParsingDictionary<Traits>;

/// Access the global dictionary of DoFType::Traits
///
/// Note: This may eventually be managed by ProjectSettings
TraitsDictionary &traits_dict();

/// Insert new DoFType::Traits into the global dictionary
void register_traits(Traits const &_traits);

/// Lookup DoFType::Traits in the global dictionary
Traits const &traits(std::string const &dof_key);

/// Access DoF::BasicTraits for a given DoF type
DoF::BasicTraits const &basic_traits(std::string const &dof_key);

/// \brief Collection of all the traits specific to a DoF type
class Traits {
 public:
  static std::string class_desc() { return "DoFType::Traits"; }

  Traits(AnisoValTraits const &_val_traits, bool _requires_site_basis = false)
      : m_val_traits(_val_traits),
        m_requires_site_basis(_requires_site_basis) {}

  AnisoValTraits const &val_traits() const { return m_val_traits; }

  std::string const &name() const { return val_traits().name(); }

  std::string site_basis_name() const { return name() + "_site_func"; }

  bool requires_site_basis() const { return m_requires_site_basis; }

  /// \brief Allow destruction through base pointer
  virtual ~Traits() {}

  /// \brief Retrieve the standard values for a DoF from dictionary of
  /// properties from properties.calc.json
  ///  Returns matrix with standard values, and names of properties that were
  ///  used to construct the matrix
  virtual std::pair<Eigen::MatrixXd, std::set<std::string> > find_values(
      std::map<std::string, Eigen::MatrixXd> const &values) const;

  /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given
  /// its site
  virtual std::vector<BasisSet> construct_site_bases(
      Structure const &_prim,
      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
      BasisFunctionSpecs const &_basis_function_specs) const = 0;

  /// \brief Populate @param _in from JSON
  virtual void from_json(DoFSet &_in, jsonParser const &_json) const {}

  /// \brief Output @param _in to JSON
  virtual void to_json(DoFSet const &_out, jsonParser &_json) const;

  /// \brief Transforms SimpleSructure @param _struc by applying DoF values
  /// contained in @param _dof in a type-specific way
  virtual void apply_dof(ConfigDoF const &_dof,
                         BasicStructure const &_reference,
                         SimpleStructure &_struc) const;

  /// \brief Serialize type-specific DoF values from ConfigDoF
  virtual jsonParser dof_to_json(ConfigDoF const &_dof,
                                 BasicStructure const &_reference) const;

  // ** The following functionality is utilized for controlling clexulator
  // printing. It only needs to be overridden in special cases **

  /// \brief
  virtual std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors(
      std::string const &nlist_specifier = "%n") const;

  virtual std::vector<std::unique_ptr<FunctionVisitor> >
  clust_function_visitors() const;

  virtual std::string site_basis_description(BasisSet site_bset, Site site,
                                             Index site_ix) const;

  virtual std::vector<ParamAllocation> param_pack_allocation(
      Structure const &_prim, std::vector<BasisSet> const &_bases) const;

  virtual std::string clexulator_constructor_string(
      Structure const &_prim, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_point_prepare_string(
      Structure const &_prim,
      std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
      PrimNeighborList &_nlist, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_global_prepare_string(
      Structure const &_prim,
      std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
      PrimNeighborList &_nlist, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_member_declarations_string(
      Structure const &_prim, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_private_method_declarations_string(
      Structure const &_prim, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_public_method_declarations_string(
      Structure const &_prim, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_private_method_definitions_string(
      Structure const &_prim, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  virtual std::string clexulator_public_method_definitions_string(
      Structure const &_prim, std::vector<BasisSet> const &site_bases,
      std::string const &indent) const;

  /// Parse DoF-specific basis function specs & validate. Default does nothing.
  ///
  /// This function may be overridden in any DoFType::Traits-derived class
  /// (i.e. OccupationDoFTraits) to read in parameters for the corresponding
  /// DoFSpecs-derived data structure (i.e. OccupationDoFSpecs). Expect that
  /// parser.value is already constructed.
  ///
  /// Example JSON format for BasisFunctionSpecs with OccupationDoFSpecs:
  /// \code
  /// "basis_function_specs": {
  ///   "dofs": ["occ", ...],
  ///   ...
  ///   "dof_specs": {
  ///     "occ": <OccupationDoFSpecs JSON>,
  ///     ...
  ///   }
  /// }
  /// \endcode
  ///
  virtual void parse_dof_specs(InputParser<BasisFunctionSpecs> &parser,
                               Structure const &prim) const {}

  /// Output DoF-specific basis function specs to json. Default does nothing.
  virtual void dof_specs_to_json(BasisFunctionSpecs const &basis_function_specs,
                                 jsonParser &json,
                                 Structure const &prim) const {}

  /// \brief non-virtual method to obtain copy through Traits pointer
  std::unique_ptr<Traits> clone() const {
    return std::unique_ptr<Traits>(_clone());
  }

 private:
  virtual Traits *_clone() const = 0;

  AnisoValTraits m_val_traits;
  bool m_requires_site_basis;
};

struct ParamAllocation {
 public:
  ParamAllocation(std::string const &_param_name, Index _param_dim,
                  Index _num_param, bool _independent)
      : param_name(_param_name),
        param_dim(_param_dim),
        num_param(_num_param),
        independent(_independent) {}

  const std::string param_name;
  const Index param_dim;
  const Index num_param;
  const bool independent;
};

}  // namespace DoFType

}  // namespace CASM
#endif
