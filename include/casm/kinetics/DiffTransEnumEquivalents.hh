#ifndef CASM_DiffTransEnumEquivalents
#define CASM_DiffTransEnumEquivalents

#include "casm/misc/cloneable_ptr.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/EnumEquivalents.hh"

namespace CASM {
  namespace Kinetics {
  /// \brief Enumerate unique DiffusionTransformations in DiffTransOrbit.
  ///
  /// - The 'begin' DiffusionTransformation is always the prototype of the DiffTransOrbit
  /// - The background Configuration must is assumed to be in primitive form
  ///
  /// \ingroup EnumEquivalents
  ///
  class DiffTransEnumEquivalents : public EnumEquivalents<DiffusionTransformation, PermuteIterator> {
    /// Need input for bg_config, but not sure if it needs to be Configuration or some sort of factor-group type.
    // -- Required -------------------

  public:

    DiffTransEnumEquivalents(const DiffusionTransformation &diff_trans, PermuteIterator begin, PermuteIterator end, const Configuration &bg_config_prim);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;

  };
  }
}

#endif
