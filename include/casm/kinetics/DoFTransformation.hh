#ifndef CASM_DoFTransformation
#define CASM_DoFTransformation

#include <memory>
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  class Configuration;

  namespace Kinetics {
    template<typename Base>
    class DoFTransformation;
  }

  template<typename Base>
  using DoFTransformation = Kinetics::DoFTransformation<Base>;

  namespace Kinetics {

    /// \brief CRTP base class for transformations
    template<typename Base>
    class DoFTransformation : public Base {

    public:

      typedef typename Base::MostDerived MostDerived;
      using Base::derived;

      // --- required in MostDerived: ---

      //Configuration &apply_to(Configuration &config) const;
      //MostDerived &apply_sym(const SymOp &op);
      //MostDerived &apply_sym(const PermuteIterator &op);
      //void reverse();


      // --- have optional default implementation: ---

      Configuration &apply_reverse_to(Configuration &config) const;



    protected:

      // customizable functions:
      Configuration &apply_reverse_to_impl(Configuration &config) const;

    };

  }
}

#endif
