#ifndef CASM_DoFTransformation_impl
#define CASM_DoFTransformation_impl

#include "casm/kinetics/DoFTransformation.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  namespace Kinetics {

    // DoFTransformation

    template<typename Base>
    Configuration &DoFTransformation<Base>::apply_reverse_to(Configuration &config) const {
      return derived().apply_reverse_to_impl(config);
    }

    template<typename Base>
    Configuration &DoFTransformation<Base>::apply_reverse_to_impl(Configuration &config) const {
      MostDerived tmp {derived()};
      tmp.reverse();
      tmp.apply_to(config);
      return config;
    }

  }
}

#endif
