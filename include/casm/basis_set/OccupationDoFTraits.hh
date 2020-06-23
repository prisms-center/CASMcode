#ifndef CASM_OccupationDoFTraits
#define CASM_OccupationDoFTraits
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/misc/Validator.hh"

#include "casm/casm_io/enum/io_traits.hh"
#include "casm/casm_io/enum/json_io.hh"
#include "casm/casm_io/enum/stream_io.hh"

namespace CASM {
  namespace DoF_impl {

    // --- OccupationDoFSpecs ---

    enum class SITE_BASIS_FUNCTION_TYPE {
      CHEBYCHEV, OCCUPATION, COMPOSITION
    };

    /// Specifies the composition on one or more equivalent sublattices
    ///
    /// Example, specify sublattices 0 and 1 are "A":0.5, "B":0.5 with:
    /// \code
    /// SublatComp sublat_comp { {0, 1}, {{"A", 0.5}, {"B", 0.5}} };
    /// \endcode
    struct SublatComp {
      SublatComp() {}
      SublatComp(
        std::initializer_list<Index> _indices,
        std::initializer_list<std::pair<std::string, double>> _values):
        indices(_indices), values(_values.begin(), _values.end()) {}
      SublatComp(
        std::set<Index> _indices,
        std::map<std::string, double> _values):
        indices(_indices), values(_values) {}

      std::set<Index> indices;
      std::map<std::string, double> values;
    };

    /// Use to specify how to construct site basis functions for occupation DoF
    ///
    /// The OccupationDoFSpecs are added to BasisFunctionSpecs to provide input to the
    /// OccupationDoFTraits::construct_site_bases method
    ///
    /// Examples, inserting OccupationDoFSpecs into BasisFunctionSpecs:
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.insert(DoFType::chebychev_bfuncs());
    /// \endcode
    ///
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.insert(DoFType::chebychev_bfuncs());
    /// \endcode
    ///
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.insert(
    ///   DoFType::composition_bfuncs( {
    ///     { {0, 1}, {{"A", 0.5}, {"B", 0.5}} },
    ///     { {2, 3}, {{"C", 0.25}, {"D", 0.75}} }
    ///   }));
    /// \endcode
    ///
    struct OccupationDoFSpecs : public DoFSpecs {

      /// Constructor for any SITE_BASIS_FUNCTION_TYPE
      OccupationDoFSpecs(SITE_BASIS_FUNCTION_TYPE _site_basis_function_type):
        site_basis_function_type(_site_basis_function_type) {}

      /// Specialized constructor for SITE_BASIS_FUNCTION_TYPE::COMPOSITION
      OccupationDoFSpecs(std::vector<SublatComp> _sublat_composition):
        site_basis_function_type(SITE_BASIS_FUNCTION_TYPE::COMPOSITION),
        sublat_composition(_sublat_composition) {}

      SITE_BASIS_FUNCTION_TYPE site_basis_function_type;
      std::vector<SublatComp> sublat_composition;

      CLONEABLE(OccupationDoFSpecs)

    private:
      std::string _name() const override;
    };

    std::vector<double> chebychev_sublat_prob_vec(Index occupant_dof_size);

    std::vector<double> occupation_sublat_prob_vec(Index occupant_dof_size);

    std::vector<double> composition_sublat_prob_vec(
      const OccupationDoFSpecs &occ_specs,
      Index sublat_index,
      const std::vector<std::string> &allowed_occupants);

    /// Validates OccupationDoFSpecs against a prim Structure
    ///
    /// Validates OccupationDoFSpecs for:
    /// - valid sublattice indices: in range [0, prim.basis.size()), no duplicates, no missing
    /// - check that molecule names are allowed on sublattices
    /// - check that all sublattices have composition
    Validator validate(OccupationDoFSpecs const &occ_specs, const Structure &prim);


    class OccupationDoFTraits : public DoFType::Traits {
    public:
      using BasicTraits = DoF::BasicTraits;
      OccupationDoFTraits():
        DoFType::Traits(BasicTraits("occ",
      {},
      BasicTraits::LOCAL),
      /*_requires_site_basis = */ true) {

      }

      std::string site_basis_description(BasisSet site_bset, Site site, Index site_ix) const override;

      std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors(std::string const &nlist_specifier) const override;

      std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const override;

      std::vector<DoFType::ParamAllocation> param_pack_allocation(Structure const &_prim,
                                                                  std::vector<BasisSet> const &_bases) const override;

      std::string clexulator_constructor_string(Structure const &_prim,
                                                std::vector<BasisSet> const &site_bases,
                                                std::string const &indent) const override;

      std::string clexulator_point_prepare_string(Structure const &_prim,
                                                  std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                  PrimNeighborList &_nlist,
                                                  std::vector<BasisSet> const &site_bases,
                                                  std::string const &indent) const override;

      std::string clexulator_global_prepare_string(Structure const &_prim,
                                                   std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                   PrimNeighborList &_nlist,
                                                   std::vector<BasisSet> const &site_bases,
                                                   std::string const &indent) const override;

      std::string clexulator_member_declarations_string(Structure const &_prim,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const override;

      std::string clexulator_private_method_declarations_string(Structure const &_prim,
                                                                std::vector<BasisSet> const &site_bases,
                                                                std::string const &indent) const override;

      std::string clexulator_public_method_declarations_string(Structure const &_prim,
                                                               std::vector<BasisSet> const &site_bases,
                                                               std::string const &indent) const override {
        //todo
        return std::string();
      }

      std::string clexulator_private_method_definitions_string(Structure const &_prim,
                                                               std::vector<BasisSet> const &site_bases,
                                                               std::string const &indent) const override {
        // todo
        return std::string();
      }

      std::string clexulator_public_method_definitions_string(Structure const &_prim,
                                                              std::vector<BasisSet> const &site_bases,
                                                              std::string const &indent) const override {
        // todo
        return std::string();
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 BasisFunctionSpecs const &_basis_function_specs) const override;

      /// Parse DoF-specific basis function specs & validate.
      void parse_dof_specs(InputParser<BasisFunctionSpecs> &parser, Structure const &prim) const override;

      /// Output DoF-specific basis function specs to json.
      void dof_specs_to_json(BasisFunctionSpecs const &basis_function_specs, jsonParser &json, Structure const &prim) const override;


    protected:
      DoFType::Traits *_clone() const override;
    };
  }

  namespace DoFType {

    DoF_impl::OccupationDoFTraits occupation();

    /// Specify use of Chebychev site basis functions (random alloy, w/ equal compositions)
    ///
    /// Example, inserting OccupationDoFSpecs into BasisFunctionSpecs:
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.push_back(DoFType::chebychev_bfuncs());
    /// \endcode
    std::unique_ptr<DoFSpecs> chebychev_basis_function_specs();

    /// Specify use of "occupation" site basis functions (ordered alloy) for occupation DoF
    ///
    /// Example, inserting OccupationDoFSpecs into BasisFunctionSpecs:
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.push_back(DoFType::occupation_bfuncs());
    /// \endcode
    std::unique_ptr<DoFSpecs> occupation_basis_function_specs();

    /// Specify site basis functions (random alloy @ specified composition) for occupation DoF
    ///
    /// Example, inserting OccupationDoFSpecs into BasisFunctionSpecs:
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.push_back(
    ///   DoFType::composition_bfuncs( {
    ///     { {0, 1}, {{"A", 0.5}, {"B", 0.5}} },
    ///     { {2, 3}, {{"C", 0.25}, {"D", 0.75}} }
    ///   },
    ///   prim));
    /// \endcode
    std::unique_ptr<DoFSpecs> composition_basis_function_specs(std::vector<DoF_impl::SublatComp> sublat_composition);

  }

  // -- OccupationDoFSpecs IO --

  ENUM_TRAITS(DoF_impl::SITE_BASIS_FUNCTION_TYPE)
  ENUM_IO_DECL(DoF_impl::SITE_BASIS_FUNCTION_TYPE)
  ENUM_JSON_IO_DECL(DoF_impl::SITE_BASIS_FUNCTION_TYPE)

  /// Reads site_basis_functions from bspecs JSON and validate
  ///
  /// \code
  /// site_basis_functions: string or array (required if "occ" dof included)
  ///     Must be one of "chebychev" or "occupation" or an array specifying sublat compositions.
  ///     Example sublat composition specification:
  ///        [
  ///          {
  ///            "sublat_indices": [0, 1],
  ///            "composition": {"A": 0.25, "B": 0.75}
  ///          },
  ///          {
  ///            "sublat_indices": [2, 3],
  ///            "composition": {"A": 0.75, "B": 0.25}
  ///          }
  ///       ]
  /// \endcode
  void parse(InputParser<DoF_impl::OccupationDoFSpecs> &parser, const Structure &prim);

  void to_json(const DoF_impl::OccupationDoFSpecs &occupation_dof_specs, jsonParser &json);

}
#endif
