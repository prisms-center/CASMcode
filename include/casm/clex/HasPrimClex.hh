#ifndef CASM_HasPrimClex
#define CASM_HasPrimClex

namespace CASM {

  namespace xtal {
    class Structure;
  }

  class PrimClex;

  /// \brief Implements PrimClex dependent functions
  ///
  /// - Useful if MostDerived class has const PrimClex& MostDerived::primclex() const;
  template <typename Base>
  class HasPrimClex : public Base {
  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    const xtal::Structure &prim() const;

    double crystallography_tol() const;
  };

} // namespace CASM

#endif
