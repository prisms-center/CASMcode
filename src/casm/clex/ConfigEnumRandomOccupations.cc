#include "casm/clex/ConfigEnumRandomOccupations.hh"

#include "casm/clex/Supercell.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

ConfigEnumRandomOccupationsParams::ConfigEnumRandomOccupationsParams(
    MTRand &_mtrand, Index _n_config)
    : mtrand(_mtrand), n_config(_n_config) {}

ConfigEnumRandomOccupations::ConfigEnumRandomOccupations(
    ConfigEnumInput const &_in_config,
    ConfigEnumRandomOccupationsParams const &params)
    : ConfigEnumRandomOccupations(_in_config, params.n_config, params.mtrand) {}

ConfigEnumRandomOccupations::ConfigEnumRandomOccupations(
    ConfigEnumInput const &_in_config, Index _n_config, MTRand &_mtrand)
    : m_n_config(_n_config),
      m_mtrand(_mtrand),
      m_max_allowed(
          _in_config.configuration().supercell().max_allowed_occupation()),
      m_site_selection(_in_config.sites().begin(), _in_config.sites().end()) {
  if (m_n_config < 0) {
    throw std::runtime_error(
        "Error in ConfigEnumRandomOccupations: n_config < 0");
  }
  if (m_n_config == 0) {
    this->_invalidate();
    return;
  }

  m_current = notstd::make_cloneable<Configuration>(_in_config.configuration());

  reset_properties(*m_current);
  this->_initialize(&(*m_current));

  // Make initial random config
  this->randomize();
  _set_step(0);
  m_current->set_source(this->source(step()));
}

std::string ConfigEnumRandomOccupations::name() const {
  return enumerator_name;
}

const std::string ConfigEnumRandomOccupations::enumerator_name =
    "ConfigEnumRandomOccupations";

/// Set m_current to correct value at specified step and return a reference to
/// it
void ConfigEnumRandomOccupations::increment() {
  this->_increment_step();
  if (step() < m_n_config) {
    this->randomize();
    m_current->set_source(this->source(step()));
  } else {
    this->_invalidate();
  }
}

void ConfigEnumRandomOccupations::randomize() {
  for (Index l : m_site_selection) {
    m_current->set_occ(l, m_mtrand.randInt(m_max_allowed[l]));
  }
}

}  // namespace CASM
