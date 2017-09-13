#ifndef CASM_CRTPBase
#define CASM_CRTPBase

namespace CASM {

  namespace CASM_TMP {

    /// \brief Base class for CRTP pattern
    template<typename _MostDerived>
    class CRTPBase {

    public:

      typedef _MostDerived MostDerived;

    protected:

      MostDerived &derived() {
        return *static_cast<MostDerived *>(this);
      }

      const MostDerived &derived() const {
        return *static_cast<const MostDerived *>(this);
      }
    };
  }

  template<typename MostDerived>
  using CRTPBase = CASM_TMP::CRTPBase<MostDerived>;
}

#endif
