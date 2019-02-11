#ifndef CASM_StrainDoFTraits
#define CASM_StrainDoFTraits
#include "casm/basis_set/DoFTraits.hh"

namespace CASM {
  namespace DoF_impl {
    class StrainDoFTraits : public DoFType::Traits {
    public:
      StrainDoFTraits(std::string _metric) :
        DoFType::Traits(_metric + "strain",
                        std::vector<std::string>({"e_1", "e_2", "e_3", "e_4", "e_5", "e_6"}),
      DoFType::CONTINUOUS,
      DoFType::GLOBAL,
      false),
      m_metric(_metric) {
      }


      /// \brief Output @param _in to JSON
      void to_json(DoFSet const &_out, jsonParser &_json) const override {
        throw std::runtime_error("StrainDoFTraits::to_json not implemented!");
      }

      /// \brief Generate a symmetry representation for the supporting vector space
      Eigen::MatrixXd symop_to_matrix(SymOp const &op) const override;


      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;

      /// \brief Apply DoF values for this DoF to _struc
      void apply_dof(ConfigDoF const &_dof, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const override;

      /// \brief Return list of DoFs that *must* be applied before this DoF is applied
      std::set<std::string> before_dof_apply() const override {
        return {"atomize", "disp"};
      }

    protected:

      DoFType::BasicTraits *_clone() const override {
        return new StrainDoFTraits(*this);
      }

      std::string const m_metric;
    };
  }

  namespace DoFType {
    DoF_impl::StrainDoFTraits GLstrain();

    DoF_impl::StrainDoFTraits EAstrain();

    DoF_impl::StrainDoFTraits Hstrain();

  }
}
#endif
