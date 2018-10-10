#ifndef CASM_DisplacementDoFTraits
#define CASM_DisplacementDoFTraits

#include "casm/basis_set/DoFTraits.hh"
namespace CASM {
  namespace DoF_impl {
    class DisplacementDoFTraits : public DoFType::Traits {
    public:
      DisplacementDoFTraits():
        DoFType::Traits("disp", {
        "x", "y", "z"
      },
      DoFType::CONTINUOUS,
      DoFType::LOCAL) {
      }


      /// \brief Output @param _in to JSON
      void to_json(DoFSet const &_out, jsonParser &_json) const override {
        throw std::runtime_error("DisplacementDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for the supporting vector space
      Eigen::MatrixXd symop_to_matrix(SymOp const &op) const override;

      /// \brief Generate a symmetry representation for this DoF
      /*
      SymGroupRepID generate_symrep(MasterSymGroup const &_group,
                                    Structure const &_prim,
                                    Index _nb) const override {
        SymGroupRepID resultID=_group.allocate_representation();
        SymGroupRepHandle bprep(_group,_prim.basis_permutation_symrep_ID());
        for(Index i=0; i<_group.size(); i++){
          Index new_b=bprep[i].get_ucc_permutation(_nb);
          _group[i].set_rep(resultID,SymMatrixXd(_prim.basis(_nb).dof(type_name()).var_basis().transpose()*_group[i].matrix()*_prim.basis(_nb).dof(type_name()).var_basis()));
        }
        return resultID;
        }*/

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
      DoFType::BasicTraits *_clone() const override;
    };
  }

  namespace DoFType {
    DoF_impl::DisplacementDoFTraits displacement();
  }

}
#endif
