#ifndef CASM_DatabaseTypeOrbit_impl
#define CASM_DatabaseTypeOrbit_impl

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/symmetry/Orbit_impl.hh"

#include "casm/clex/PrimClex.hh"
#include "casm/database/Named.hh"
#include "casm/database/Database.hh"
#include "casm/app/DirectoryStructure.hh"

namespace CASM {

  template<typename _SymCompareType>
  DatabaseTypeOrbit<_SymCompareType>::DatabaseTypeOrbit(
    Orbit<_SymCompareType> const &_orbit, PrimClex const *_primclex) :
    Orbit<_SymCompareType>(_orbit),
    m_primclex(_primclex) {}

  template<typename _SymCompareType>
  const PrimClex &DatabaseTypeOrbit<_SymCompareType>::primclex() const {
    if(!m_primclex) {
      throw std::runtime_error("DatabaseTypeOrbit primclex pointer was not set");
    }
    return *m_primclex;
  }

  template<typename _SymCompareType>
  std::string DatabaseTypeOrbit<_SymCompareType>::generate_name_impl() const {
    return OrbitTraits<_SymCompareType>::generate_name_impl(*this);
  }

  template<typename _SymCompareType>
  void DatabaseTypeOrbit<_SymCompareType>::set_primclex(const PrimClex *_primclex) {
    m_primclex = _primclex;
  }


  template<typename _SymCompareType>
  void write_pos(DatabaseTypeOrbit<_SymCompareType> const &_el) {
    const auto &dir = _el.primclex().dir();
    try {
      fs::create_directories(dir.configuration_dir(_el.name()));
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in DatabaseTypeOrbit::write_pos(): could not create_directories" << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream file(dir.POS(_el.name()));
    file << pos_string(_el);
  }

  template<typename _SymCompareType>
  std::string pos_string(DatabaseTypeOrbit<_SymCompareType> const &_el) {
    std::stringstream ss;
    OrbitTraits<_SymCompareType>::write_pos(_el, ss);
    return ss.str();
  }


}

#endif
