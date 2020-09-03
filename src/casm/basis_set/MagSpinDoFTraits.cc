#include "casm/basis_set/Adapter.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/MagSpinDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/NeighborList.hh"

#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {
  namespace DoF_impl {
    /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
    std::vector<BasisSet> MagSpinDoFTraits::construct_site_bases(Structure const &_prim,
                                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                                 BasisFunctionSpecs const &_basis_function_specs) const {
      std::vector<BasisSet> result(_prim.basis().size());

      auto const &mag_spin_dof_specs = get<MagSpinDoFSpecs>(name(), _basis_function_specs);

      for(Index b = 0; b < _prim.basis().size(); b++) {

        if(!_prim.basis()[b].has_dof(name()))
          continue;
        BasisSet tresult;
        CASM::DoFSet adapted_dofset = adapter::Adapter<CASM::DoFSet, xtal::SiteDoFSet>()(_prim.basis()[b].dof(name()), _prim.site_dof_symrep_IDs()[b].at(name()), b);
        tresult.set_variable_basis(adapted_dofset);
        Array<BasisSet const *> tsubs(1, &tresult);
        result[b].construct_harmonic_polynomials(tsubs, mag_spin_dof_specs.max_poly_order, 1, false);
        result[b].get_symmetry_representation(_prim.factor_group());

        result[b].set_name(name() + "_site_func");
        //std::cout << "+:+:+:+Created variable set for site " << b << ", size " << result[b].size() << "\n";
      }
      return result;
    }

    void MagSpinDoFTraits::parse_dof_specs(InputParser<BasisFunctionSpecs> &parser, Structure const &prim) const {
      fs::path dof_specs_path {"dof_specs"};
      fs::path subpath = dof_specs_path / this->name();
      auto subparser = parser.subparse<MagSpinDoFSpecs>(subpath, prim);
      if(subparser->value) {
        parser.value->dof_specs.push_back(std::move(subparser->value));
      }
    }

    void MagSpinDoFTraits::dof_specs_to_json(const BasisFunctionSpecs &basis_function_specs, jsonParser &json, Structure const &prim) const {
      MagSpinDoFSpecs const &magspin_dof_specs = get<MagSpinDoFSpecs>(
                                                   this->name(), basis_function_specs);
      CASM::to_json(magspin_dof_specs, json["dof_specs"][this->name()]);
    }

  }

  namespace DoFType {

    DoF_impl::MagSpinDoFTraits magspin(std::string const &flavor_name) {
      return DoF_impl::MagSpinDoFTraits(flavor_name);
    }

    std::unique_ptr<DoFSpecs> magspin_bfuncs(std::string const &flavor_name, Index max_poly_order = -1) {
      return notstd::make_unique<DoF_impl::MagSpinDoFSpecs>(flavor_name, max_poly_order);
    }

  }

  // -- MagSpinDoFSpecs IO --

  void parse(InputParser<DoF_impl::MagSpinDoFSpecs> &parser, const Structure &prim) {
    std::string full_name = parser.name();
    std::string flavor_name = full_name.substr(0, full_name.find("magspin"));
    Index max_poly_order;
    parser.optional_else<Index>(max_poly_order, "max_poly_order", -1);
    if(parser.valid()) {
      parser.value = notstd::make_unique<DoF_impl::MagSpinDoFSpecs>(flavor_name, max_poly_order);
    }
  }

  void to_json(const DoF_impl::MagSpinDoFSpecs &magspin_dof_specs, jsonParser &json) {
    json.put_obj();
    json["max_poly_order"] = magspin_dof_specs.max_poly_order;
  }


}
