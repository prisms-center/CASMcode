#ifndef CASM_Comparisons
#define CASM_Comparisons

namespace CASM {

  /// \brief Implements other comparisons in terms of '<'
  ///
  /// If
  /// \code
  /// Derived::operator<(const Derived& B) const
  /// \endcode
  /// is implemented, then Comparisons<Derived> implements:
  /// - '>', '<=', '>=', '==', '!='
  /// - '==' and '!=' can be specialized in Derived by implementing private
  ///   methods '_eq' and '_ne'
  ///
  /// The Derived class definition needs to include:
  /// \code
  /// friend Comparisons<Derived>;
  /// \endcode
  ///
  template<typename Derived>
  struct Comparisons {

    bool operator>(const Derived &B) const {
      return B < derived();
    };

    bool operator<=(const Derived &B) const {
      return !(B < derived());
    };

    bool operator>=(const Derived &B) const {
      return !(derived() < B);
    };

    bool operator==(const Derived &B) const {
      return derived()._eq(B);
    };

    bool operator!=(const Derived &B) const {
      return derived()._ne(B);
    };


  protected:

    const Derived &derived() const {
      return *static_cast<const Derived *>(this);
    }

    bool _eq(const Derived &B) const {
      return (!(derived() < B) && !(B < derived()));
    };

    bool _ne(const Derived &B) const {
      return !derived()._eq(B);
    };

  };

}

#endif
