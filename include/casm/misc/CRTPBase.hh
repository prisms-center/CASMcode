#ifndef CASM_CRTPBase
#define CASM_CRTPBase

namespace notstd {

/// \brief Base class for CRTP pattern
template <typename _MostDerived>
class CRTPBase {
 public:
  typedef _MostDerived MostDerived;

 protected:
  MostDerived &derived() { return *static_cast<MostDerived *>(this); }

  const MostDerived &derived() const {
    return *static_cast<const MostDerived *>(this);
  }
};
}  // namespace notstd

namespace CASM {

template <typename MostDerived>
using CRTPBase = notstd::CRTPBase<MostDerived>;

namespace CASM_TMP {
template <typename MostDerived>
using CRTPBase = notstd::CRTPBase<MostDerived>;
}
}  // namespace CASM

#endif
