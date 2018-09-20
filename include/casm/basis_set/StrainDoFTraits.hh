#ifndef CASM_OccupationDoFTraits
#define CASM_OccupationDoFTraits

#include "casm/basis_set/DoFTraits.hh"
namespace CASM {
  namespace DoF_impl {
    class StrainDoFTraits : public DoFType::Traits {
    public:
      StrainDoFTraits(std::string _metric) :
        DoFType::Traits(_metric + "strain",
                        std::vector<std::string>(6, "e"),
                        CONTINUOUS,
                        GLOBAL),
        m_metric(_metric) {
      }


      /// \brief Output @param _in to JSON
      void to_json(DoFSet const &_out, jsonParser &_json) const override {
        throw std::runtime_error("StrainDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for the supporting vector space
      Eigen::MatrixXd symop_to_matrix(SymOp const &op) const override;


      std::string site_basis_description(BasisSet site_bset, Site site) const override;

      std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors() const override;

      std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const override;

      std::vector<std::pair<std::string, Index> > param_pack_allocation(std::vector<BasisSet> const &_bases) const override;

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

      std::string clexulator_private_method_implementations_string(Structure const &_prim,
                                                                   std::vector<BasisSet> const &site_bases,
                                                                   std::string const &indent) const override {
        // todo
        return std::string();
      }

      std::string clexulator_public_method_implementations_string(Structure const &_prim,
                                                                  std::vector<BasisSet> const &site_bases,
                                                                  std::string const &indent) const override {
        // todo
        return std::string();
      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;
    protected:
      DoFType::BasicTraits *_clone() const override {
        new StrainDoFTraits(*this);
      }

      std::string const m_metric;
    };

    namespace DoFType {
      inline
      notstd::cloneable_ptr<typename DoF_impl::BasicTraits> GLstrain() {
        return DoF_impl::StrainDoFTraits("GL").clone();
      }
      inline
      notstd::cloneable_ptr<typename DoF_impl::BasicTraits> EAstrain() {
        return DoF_impl::StrainDoFTraits("EA").clone();
      }
      inline
      notstd::cloneable_ptr<typename DoF_impl::BasicTraits> Hstrain() {
        return DoF_impl::StrainDoFTraits("H").clone();
      }
    }
  }
}
#endif
