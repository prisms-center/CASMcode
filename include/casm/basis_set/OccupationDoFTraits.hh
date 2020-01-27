#ifndef CASM_OccupationDoFTraits
#define CASM_OccupationDoFTraits
#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {
  namespace DoF_impl {
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
                                                 jsonParser const &_bspecs) const override;
    protected:
      DoFType::Traits *_clone() const override;
    };
  }

  namespace DoFType {
    DoF_impl::OccupationDoFTraits occupation();
  }

}
#endif
