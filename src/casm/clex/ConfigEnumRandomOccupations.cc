#include "casm/clex/ConfigEnumRandomOccupations.hh"
#include "casm/container/Enumerator_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/clex/FilteredConfigIterator.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumRandomOccupations_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumRandomOccupations>();
  }
}

namespace CASM {

  const std::string ConfigEnumRandomOccupations::enumerator_name = "ConfigEnumRandomOccupations";

  const std::string ConfigEnumRandomOccupations::interface_help =
    "ConfigEnumRandomOccupations: \n\n"

    "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
    "    Indicate supercells to enumerate all occupational configurations in. May \n"
    "    be a JSON array of supercell names, or a JSON object specifying          \n"
    "    supercells in terms of size and unit cell. By default, all existing      \n"
    "    supercells are used. See 'ScelEnum' description for details.         \n\n"

    "  n_config: integer (optional, default=100) \n"
    "    How many random configurations to generate. Includes duplicate and pre-\n"
    "    existing configurations.                                                 \n\n"

    "  primitive_only: bool (default=true)\n"
    "    If true, only the primitive form of a configuration is saved in the      \n"
    "    configuration list. Otherwise, both primitive and non-primitive          \n"
    "    configurations are saved. \n\n"

    "  filter: string (optional, default=None)\n"
    "    A query command to use to filter which Configurations are kept.          \n\n"

    "  Examples:\n"
    "    To enumerate 200 random occupations in supercells up to and including size 4:\n"
    "      casm enum --method ConfigEnumRandomOccupations -i \n"
    "        '{\"supercell\":{\"max\":4}, \"n_config\": 200}' \n"
    "\n"
    "    To enumerate 200 random occupations in all existing supercells:\n"
    "      casm enum --method ConfigEnumRandomOccupations -i '{\"n_config\": 200}' \n"
    "\n"
    "    To enumerate 100 random occupations in all existing supercells:\n"
    "      casm enum --method ConfigEnumRandomOccupations' \n"
    "\n"
    "    To enumerate 200 random occupations in particular supercells:\n"
    "      casm enum --method ConfigEnumRandomOccupations -i \n"
    "      '{ \n"
    "        \"supercells\": { \n"
    "          \"name\": [\n"
    "            \"SCEL1_1_1_1_0_0_0\",\n"
    "            \"SCEL2_1_2_1_0_0_0\",\n"
    "            \"SCEL4_1_4_1_0_0_0\"\n"
    "          ]\n"
    "        }, \n"
    "        \"n_config\": 200\n"
    "      }' \n\n";

  int ConfigEnumRandomOccupations::run(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    std::unique_ptr<ScelEnum> scel_enum = make_enumerator_scel_enum(primclex, _kwargs, enum_opt);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);
    MTRand mtrand;

    Index n_config;
    _kwargs.get_else<Index>(n_config, "n_config", 100);

    bool primitive_only = true;
    _kwargs.get_if(primitive_only, "primitive_only");


    auto lambda = [&](Supercell & scel) {
      return notstd::make_unique<ConfigEnumRandomOccupations>(scel, n_config, mtrand);
    };

    int returncode = insert_configs(
                       enumerator_name,
                       primclex,
                       scel_enum->begin(),
                       scel_enum->end(),
                       lambda,
                       filter_expr,
                       primitive_only);

    return returncode;
  }

  /// \brief Constructor
  ///
  /// \param _initial,_final Initial and final configurations to interpolate between
  /// \param _size The total number of configurations to enumerate, including
  ///              the initial and final configurations
  ///
  /// - The `final` configuration is *not* pointed at by the end iterator,
  ///   which points past-the-final element, as is typical
  /// - `_size` will be equal to \code std::distance(this->begin(), this->end()) \endcode
  ConfigEnumRandomOccupations::ConfigEnumRandomOccupations(
    Supercell &_scel,
    Index _n_config,
    MTRand &_mtrand):
    m_n_config(_n_config),
    m_mtrand(&_mtrand),
    m_max_allowed(_scel.max_allowed_occupation()) {

    if(m_n_config < 0) {
      throw std::runtime_error("Error in ConfigEnumRandomOccupations: n_config < 0");
    }
    if(m_n_config == 0) {
      this->_invalidate();
      return;
    }

    m_current = notstd::make_cloneable<Configuration>(_scel, this->source(0));
    m_current->init_occupation();

    reset_properties(*m_current);
    this->_initialize(&(*m_current));

    // Make initial random config
    this->randomize();
    _set_step(0);
    _current().set_source(this->source(step()));
  }

  /// Set m_current to correct value at specified step and return a reference to it
  void ConfigEnumRandomOccupations::increment() {

    this->_increment_step();
    if(step() < m_n_config) {
      this->randomize();
      _current().set_source(this->source(step()));
    }
    else {
      this->_invalidate();
    }

  }

  void ConfigEnumRandomOccupations::randomize() {
    for(Index i = 0; i < m_current->size(); ++i) {
      m_current->set_occ(i, m_mtrand->randInt(m_max_allowed[i]));
    }
  }

}
