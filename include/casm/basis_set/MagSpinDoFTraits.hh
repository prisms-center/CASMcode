#ifndef CASM_MagSpinDoFTraits
#define CASM_MagSpinDoFTraits

#include "casm/basis_set/DoFTraits.hh"
namespace CASM {
  namespace DoF_impl {
    class MagSpinDoFTraits : public DoFType::Traits {
    public:
      MagSpinDoFTraits():
        DoFType::Traits("magspin", {
        "sx", "sy", "sz"
      },
      DoFType::CONTINUOUS,
      DoFType::LOCAL) {
      }

      bool time_reversal_active() const override {
        return true;
      }

      /// \brief Output @param _in to JSON
      void to_json(DoFSet const &_out, jsonParser &_json) const override {
        throw std::runtime_error("MagSpinDoFTraits::to_json not implemented!");
      }

      Eigen::MatrixXd symop_to_matrix(SymOp const &op) const override;

      std::string site_basis_description(BasisSet site_bset, Site site) const override {
        return "";
      }

      std::vector<std::unique_ptr<FunctionVisitor> > site_function_visitors() const override {
        return std::vector<std::unique_ptr<FunctionVisitor> >();
      }

      std::vector<std::unique_ptr<FunctionVisitor> > clust_function_visitors() const override {
        return std::vector<std::unique_ptr<FunctionVisitor> >();
      }

      std::vector<std::pair<std::string, Index> > param_pack_allocation(std::vector<BasisSet> const &_bases) const override {
        return std::vector<std::pair<std::string, Index> >();
      }

      std::string clexulator_constructor_string(Structure const &_prim,
                                                std::vector<BasisSet> const &site_bases,
                                                std::string const &indent) const override {
        return "";
      }

      std::string clexulator_point_prepare_string(Structure const &_prim,
                                                  std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                  PrimNeighborList &_nlist,
                                                  std::vector<BasisSet> const &site_bases,
                                                  std::string const &indent) const override {
        return "";
      }

      std::string clexulator_global_prepare_string(Structure const &_prim,
                                                   std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                   PrimNeighborList &_nlist,
                                                   std::vector<BasisSet> const &site_bases,
                                                   std::string const &indent) const override {
        return "";
      }

      std::string clexulator_member_declarations_string(Structure const &_prim,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const override {
        return "";
      }

      std::string clexulator_private_method_declarations_string(Structure const &_prim,
                                                                std::vector<BasisSet> const &site_bases,
                                                                std::string const &indent) const override {
        return "";
      }

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
                                                 jsonParser const &_bspecs) const override {
        return std::vector<BasisSet>();
      }
    protected:
      DoFType::BasicTraits *_clone() const override {
        return new MagSpinDoFTraits(*this);
      }
    };
  }

  namespace DoFType {
    DoF_impl::MagSpinDoFTraits magspin();
  }

}
#endif
