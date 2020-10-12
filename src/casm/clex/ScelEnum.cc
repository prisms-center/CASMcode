#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {

  /// Always true for ScelEnumByProps
  template<>
  bool is_guaranteed_for_database_insert(ScelEnumByProps const &enumerator) {
    return true;
  }

  /// \brief Construct with shared prim Structure and ScelEnumProps settings
  ///
  /// \param shared_prim A shared prim Structure for which to enumerate Supercells
  /// \param enum_props Specifies which Supercells to enumerate
  ///
  /// Note: This variant does not require a PrimClex, and there for cannot insert Supercells into
  /// a Supercell database automatically.
  ///
  ScelEnumByProps::ScelEnumByProps(std::shared_ptr<const Structure> const &shared_prim,
                                   const xtal::ScelEnumProps &enum_props) :
    m_shared_prim(shared_prim) {

    auto const &pg = m_shared_prim->point_group();
    m_lattice_enum.reset(new xtal::SuperlatticeEnumerator(
                           pg.begin(),
                           pg.end(),
                           m_shared_prim->lattice(),
                           enum_props
                         ));

    m_lat_it = m_lattice_enum->begin();
    m_lat_end = m_lattice_enum->end();

    if(m_lat_it != m_lat_end) {
      double xtal_tol = m_shared_prim->lattice().tol();
      Lattice canonical_lattice = xtal::canonical::equivalent(*m_lat_it, pg, xtal_tol);
      m_current = notstd::make_unique<Supercell>(m_shared_prim, canonical_lattice);
      this->_initialize(&(*m_current));
    }
    else {
      this->_invalidate();
    }
  }

  std::string ScelEnumByProps::name() const {
    return ScelEnumByProps::enumerator_name;
  }

  const std::string ScelEnumByProps::enumerator_name = "ScelEnumByProps";

  /// Implements increment over supercells
  void ScelEnumByProps::increment() {
    ++m_lat_it;

    if(m_lat_it != m_lat_end) {
      double xtal_tol = m_shared_prim->lattice().tol();
      auto const &pg = m_shared_prim->point_group();
      Lattice canonical_lattice = xtal::canonical::equivalent(*m_lat_it, pg, xtal_tol);
      m_current = notstd::make_unique<Supercell>(m_shared_prim, canonical_lattice);
      this->_initialize(&(*m_current));
      this->_increment_step();
    }
    else {
      this->_invalidate();
    }
  }

}
