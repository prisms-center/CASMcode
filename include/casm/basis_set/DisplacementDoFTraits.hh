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
      DoFType::LOCAL,
      /*_requires_site_basis = */ false,
      /*_unit_length = */ false) {
      }

      /// \brief Generate a symmetry representation for the supporting vector space
      Eigen::MatrixXd symop_to_matrix(SymOp const &op) const override;


      /// \brief Apply DoF values for this DoF to _struc
      void apply_dof(ConfigDoF const &_dof, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const override;


      /// \brief Return list of DoFs that *must* be applied after this DoF is applied
      std::set<std::string> after_dof_apply() const override {
        return {"atomize"};
      }

      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
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
