#include "casm/clex/ConfigEnumRandomOccupations.hh"
#include "casm/app/enum/EnumInterface_impl.hh"
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

  std::string ConfigEnumRandomOccupations::interface_help() {
    return
      "ConfigEnumRandomOccupations: \n\n"

      "  confignames: Array of strings (optional) \n"
      "    Names of configurations to be used as initial state of occupation enumeration.\n"
      "    All specified sublattices or sites will be enumerated on and all other DoFs will\n"
      "    maintain the values of the initial state.\n"
      "    Ex: \"confignames\" : [\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]\n\n"

      "  scelnames: Array of strings (optional) \n"
      "    Names of supercells used as initial state of occupation enumeration. All\n"
      "    sites not being enumerated will be set to the first listed occupant, and all\n"
      "    other DoFs will be set to zero.\n"
      "    Ex: \"scelnames\" : [\"SCEL1_1_1_1_0_0_0\",\"SCEL2_2_1_1_0_0_0\"]\n\n"

      "  sublats: array of integers (optional, default none) \n"
      "    Restricts enumeration to specified sublattices. Each sublattice index corresponds\n"
      "    to a basis site in prim.json, indexed from 0.\n"
      "    Ex: \"sublats\" : [0,2]\n\n"

      "  sites: array of 4-entry integer arrays (optional, default none) \n"
      "    Restricts normal coordinate determination to specified sites. Sites are specified\n"
      "    in [b,i,j,k] convention, where 'b' is sublattice index and [i,j,k] spedifies line-\n"
      "    ar combinations of primitive-cell lattice vectors.\n"
      "    Ex: \"sites\" : [[0,0,0,0],\n"
      "                   [2,0,0,0]]\n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
      "    Indicate supercells to use as initial states of enumeration in terms of size\n"
      "    and unit cell via a JSON object conforming to the format of 'ScelEnum' JSON\n"
      "    settings. \"scelnames\" will override \"supercells\", but if neither is specified\n"
      "    all existing supercells are used by default. See 'ScelEnum' description for details.\n\n"

      "  n_config: integer (optional, default=100) \n"
      "    How many random configurations to generate. Includes duplicate and pre-\n"
      "    existing configurations.                                                 \n\n"

      "  primitive_only: bool (default=true)\n"
      "    If true, only the primitive form of a configuration is saved in the      \n"
      "    configuration list. Otherwise, both primitive and non-primitive          \n"
      "    configurations are saved. \n\n"

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
  }

  int ConfigEnumRandomOccupations::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt,
    EnumeratorMap const *interface_map) {

    std::unique_ptr<ScelEnum> scel_enum = make_enumerator_scel_enum(primclex, _kwargs, enum_opt);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);
    MTRand mtrand;

    Index n_config;
    _kwargs.get_else<Index>(n_config, "n_config", 100);

    bool primitive_only = true;
    _kwargs.get_if(primitive_only, "primitive_only");


    auto lambda = [&](const ConfigEnumInput & _input) {
      return notstd::make_unique<ConfigEnumRandomOccupations>(_input, n_config, mtrand);
    };

    int returncode = insert_configs(
                       enumerator_name,
                       primclex,
                       scel_enum->begin(),
                       scel_enum->end(),
                       lambda,
                       filter_expr,
                       primitive_only,
                       CASM::dry_run(_kwargs, enum_opt));

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
    ConfigEnumInput const &_in_config,
    Index _n_config,
    MTRand &_mtrand):
    m_n_config(_n_config),
    m_mtrand(&_mtrand),
    m_max_allowed(_in_config.configuration().supercell().max_allowed_occupation()),
    m_site_selection(_in_config.sites().begin(), _in_config.sites().end()) {

    if(m_n_config < 0) {
      throw std::runtime_error("Error in ConfigEnumRandomOccupations: n_config < 0");
    }
    if(m_n_config == 0) {
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

  /// Set m_current to correct value at specified step and return a reference to it
  void ConfigEnumRandomOccupations::increment() {

    this->_increment_step();
    if(step() < m_n_config) {
      this->randomize();
      m_current->set_source(this->source(step()));
    }
    else {
      this->_invalidate();
    }

  }

  void ConfigEnumRandomOccupations::randomize() {
    for(Index l : m_site_selection) {
      m_current->set_occ(l, m_mtrand->randInt(m_max_allowed[l]));
    }
  }

}
