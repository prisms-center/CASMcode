#ifndef CASM_DatabaseTypeOrbit
#define CASM_DatabaseTypeOrbit

#include "casm/symmetry/Orbit.hh"
#include "casm/database/Named.hh"
#include "casm/clex/HasPrimClex.hh"

namespace CASM {

  // -- DatabaseTypeOrbit -------------------------------------

  /// \brief Specialize GenericOrbit for Orbit types that will be stored in a database
  ///
  template<typename _SymCompareType>
  class DatabaseTypeOrbit :
    public Orbit<_SymCompareType>,
    public HasPrimClex<DB::Indexed<CRTPBase<DatabaseTypeOrbit<_SymCompareType>>>> {
  public:
    using Element = typename _SymCompareType::Element;
    using SymCompareType = _SymCompareType;

    DatabaseTypeOrbit(Orbit<_SymCompareType> const &_orbit, PrimClex const *_primclex);

    const PrimClex &primclex() const;

  private:

    friend DB::Named<CRTPBase<DatabaseTypeOrbit<_SymCompareType>>>;

    std::string generate_name_impl() const;

    void set_primclex(const PrimClex *_primclex);

    const PrimClex *m_primclex;
  };


  // -- Orbit Helpers --------------------

  template<typename _SymCompareType>
  void write_pos(DatabaseTypeOrbit<_SymCompareType> const &_el);

  template<typename _SymCompareType>
  std::string pos_string(DatabaseTypeOrbit<_SymCompareType> const &_el);

}

#endif
