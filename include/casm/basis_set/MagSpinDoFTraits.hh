#ifndef CASM_MagSpinDoFTraits
#define CASM_MagSpinDoFTraits
#include "casm/basis_set/DoFTraits.hh"

#include "casm/basis_set/BasisSet.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {
  namespace DoF_impl {
    class MagSpinDoFTraits : public DoFType::Traits {
    public:
      MagSpinDoFTraits():
        DoFType::Traits(AnisoValTraits::magspin(), true) {

      }

      /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
      std::vector<BasisSet> construct_site_bases(Structure const &_prim,
                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                 jsonParser const &_bspecs) const override;
    protected:
      DoFType::Traits *_clone() const override {
        return new MagSpinDoFTraits(*this);
      }
    };
  }

  namespace DoFType {
    DoF_impl::MagSpinDoFTraits magspin();
  }

}
#endif
