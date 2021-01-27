#ifndef CASM_Comparisons
#define CASM_Comparisons

#include "casm/misc/CRTPBase.hh"

namespace notstd {

/// \brief Implements other comparisons in terms of '<'
///
/// If
/// \code
/// Derived::operator<(const MostDerived& B) const
/// \endcode
/// is implemented, then Comparisons implements:
/// - '>', '<=', '>=', '==', '!='
/// - '==' and '!=' can be specialized in MostDerived by implementing private
///   methods '_eq' and '_ne'
///
/// The MostDerived class definition needs to include:
/// \code
/// template<typename T> friend class Comparisons;
/// \endcode
///
template <typename Base>
struct Comparisons : public Base {
  using Base::derived;
  typedef typename Base::MostDerived MostDerived;

  bool operator>(const MostDerived &B) const { return B < derived(); };

  bool operator<=(const MostDerived &B) const { return !(B < derived()); };

  bool operator>=(const MostDerived &B) const { return !(derived() < B); };

  bool operator==(const MostDerived &B) const { return derived().eq_impl(B); };

  bool operator!=(const MostDerived &B) const { return derived().ne_impl(B); };

 protected:
  bool eq_impl(const MostDerived &B) const {
    return (!(derived() < B) && !(B < derived()));
  };

  bool ne_impl(const MostDerived &B) const { return !derived().eq_impl(B); };
};

}  // namespace notstd

namespace CASM {
template <typename Base>
using Comparisons = notstd::Comparisons<Base>;
}

#endif
