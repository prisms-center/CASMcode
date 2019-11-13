#ifndef CASM_StrainDoFTraits
#define CASM_StrainDoFTraits
#include "casm/basis_set/DoFTraits.hh"

namespace CASM {
  namespace DoF_impl {
    class StrainDoFTraits : public DoFType::Traits {
    public:

      StrainDoFTraits(std::string const &_metric) :
        DoFType::Traits(AnisoValTraits::strain(_metric)),
        m_metric(_metric) {
      }


      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;

      /// \brief Serialize strain DoF values from ConfigDoF
      jsonParser dof_to_json(ConfigDoF const &_dof, BasicStructure<Site> const &_reference) const override;

      /// \brief Transforms SimpleSructure @param _struc by applying strain DoF values contained in @param _dof
      void apply_dof(ConfigDoF const &_dof, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const override;

    protected:

      DoFType::Traits *_clone() const override {
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
