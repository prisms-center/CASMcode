#ifndef CASM_MagSpinDoFTraits
#define CASM_MagSpinDoFTraits
#include "casm/basis_set/DoFTraits.hh"

#include "casm/basis_set/BasisSet.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {
  namespace DoF_impl {

    struct MagSpinDoFSpecs : public DoFSpecs {

      /// Constructor for any SITE_BASIS_FUNCTION_TYPE
      MagSpinDoFSpecs(std::string _flavor_name, Index _max_poly_order = -1):
        max_poly_order(_max_poly_order),
        m_flavor_name(_flavor_name) {}

      Index max_poly_order;

      CLONEABLE(MagSpinDoFSpecs)

    private:

      std::string m_flavor_name;

      std::string _name() const override {
        return m_flavor_name + "magspin";
      }
    };


    class MagSpinDoFTraits : public DoFType::Traits {
    public:
      // Initialize with particular 'flavor' of magnetic spin. Options are:
      //  collinear: "C" and "Cunit"
      //  non-collinear (w/o spin-orbit): "NC" and "NCunit"
      //  non-collinear (with spin-orbit): "SO" and "SOunit"
      MagSpinDoFTraits(std::string const &flavor_name):
        DoFType::Traits(AnisoValTraits(flavor_name + "magspin"), true) {

      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 BasisFunctionSpecs const &_basis_function_specs) const override;

      void parse_dof_specs(InputParser<BasisFunctionSpecs> &parser, Structure const &prim) const override;

      void dof_specs_to_json(BasisFunctionSpecs const &basis_function_specs, jsonParser &json, Structure const &prim) const override;

    protected:
      DoFType::Traits *_clone() const override {
        return new MagSpinDoFTraits(*this);
      }
    };
  }

  namespace DoFType {
    /// Initialize with particular 'flavor' of magnetic spin. Options are:
    ///  collinear: "C" and "Cunit"
    ///  non-collinear (w/o spin-orbit): "NC" and "NCunit"
    ///  non-collinear (with spin-orbit): "SO" and "SOunit"
    DoF_impl::MagSpinDoFTraits magspin(std::string const &flavor_name);

    /// Specify max polynomial order for magspin site basis functions
    ///
    /// Example, inserting MagSpinDoFSpecs into BasisFunctionSpecs:
    /// \code
    /// BasisFunctionSpecs bspecs;
    /// bspecs.dof_specs.push_back(DoFType::magspin_basis_function_specs(flavor_name, max_poly_order));
    /// \endcode
    std::unique_ptr<DoFSpecs> magspin_specs(std::string const &flavor_name, Index max_poly_order);
  }

  // -- MagSpinDoFSpecs IO --

  /// Specify max polynomial order for magspin site basis functions
  ///
  /// Options:
  /// \code
  /// {
  ///   "max_poly_order": <int, optional, default=-1>
  /// \endcode
  ///
  void parse(InputParser<DoF_impl::MagSpinDoFSpecs> &parser, const Structure &prim);

  /// Specify max polynomial order for magspin site basis functions
  ///
  /// Options:
  /// \code
  /// {
  ///   "max_poly_order": <int, optional, default=-1>
  /// \endcode
  ///
  void to_json(const DoF_impl::MagSpinDoFSpecs &occ_specs, jsonParser &json);

}
#endif
