#ifndef CASM_MagSpinDoFTraits
#define CASM_MagSpinDoFTraits

#include "casm/basis_set/DoFTraits.hh"
namespace CASM{
  class MagSpinDoFTraits : public Traits {
  public:
    MagSpinDoFTraits():
      Traits("magspin",
      {"sx","sy","sz"},
             CONTINUOUS,
             LOCAL) {
    }

    bool time_reversal_active const override{
      return true;
    }
      
    /// \brief Output @param _in to JSON
    void to_json(DoFSet const &_out, jsonParser &_json) const override {
      throw std::runtime_error("MagSpinDoFTraits::to_json not implemented!");
    }

    /// \brief Generate a symmetry representation for the supporting vector space
    Eigen::MatrixXd symop_to_matrix(SymOp const &op) const override{
      return op.matrix().determinant()*op.matrix();
    }
    
    /// \brief Generate a symmetry representation for this DoF
    SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                  Structure const &_prim,
                                  Index _nb) const override {
      throw std::runtime_error("MagSpinDoFTraits::generate_symrep not implemented!");
    }

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
    BasicTraits *_clone() const override;
  };

  namespace DoFType {

    inline
    notstd::cloneable_ptr<typename DoF_impl::BasicTraits> magspin() {
      return DoF_impl::MagSpinDoFTraits().clone();
    }
  }
}
#endif
