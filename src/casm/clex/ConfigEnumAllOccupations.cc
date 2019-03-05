#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/enumerator/Enumerator_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumAllOccupations_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumAllOccupations>();
  }
}

namespace CASM {

  const std::string ConfigEnumAllOccupations::enumerator_name = "ConfigEnumAllOccupations";

  std::string ConfigEnumAllOccupations::interface_help() {
    return
      "ConfigEnumAllOccupations: \n\n"

      "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
      "    Indicate supercells to enumerate all occupational configurations in. May \n"
      "    be a JSON array of supercell names, or a JSON object specifying          \n"
      "    supercells in terms of size and unit cell. By default, all existing      \n"
      "    supercells are used. See 'ScelEnum' description for details.         \n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  Examples:\n"
      "    To enumerate all occupations in supercells up to and including size 4:\n"
      "      casm enum --method ConfigEnumAllOccupations -i '{\"supercells\": {\"max\": 4}}' \n"
      "\n"
      "    To enumerate all occupations in all existing supercells:\n"
      "      casm enum --method ConfigEnumAllOccupations\n"
      "\n"
      "    To enumerate all occupations in all particular supercells:\n"
      "      casm enum --method ConfigEnumAllOccupations -i \n"
      "      '{ \n"
      "        \"supercells\": { \n"
      "          \"name\": [\n"
      "            \"SCEL1_1_1_1_0_0_0\",\n"
      "            \"SCEL2_1_2_1_0_0_0\",\n"
      "            \"SCEL4_1_4_1_0_0_0\"\n"
      "          ]\n"
      "        } \n"
      "      }' \n\n";
  }

  int ConfigEnumAllOccupations::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt,
    EnumeratorMap const *interface_map) {

    std::vector<ConfigEnumInput> in_configs = make_enumerator_input_configs(primclex, _kwargs, enum_opt, interface_map);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    return run(primclex, in_configs.begin(), in_configs.end(), filter_expr, CASM::dry_run(_kwargs, enum_opt));
  }

  namespace local_impl {
    std::vector<int> max_selected_occupation(ConfigEnumInput const &_config) {
      std::vector<int> max = _config.supercell().max_allowed_occupation();
      std::vector<int> maxselect;
      for(Index i : _config.sites())
        maxselect.push_back(max[i]);

      if(maxselect.empty())
        return max;
      else
        return maxselect;
    }
  }

  /// \brief Construct with a Supercell, using all permutations
  ConfigEnumAllOccupations::ConfigEnumAllOccupations(const ConfigEnumInput &_input_config) :
    m_selection(_input_config.sites().begin(), _input_config.sites().end()),
    m_counter(
      std::vector<int>(_input_config.sites().size(), 0),
      local_impl::max_selected_occupation(_input_config),
      std::vector<int>(_input_config.sites().size(), 1)),

    m_current(notstd::make_cloneable<Configuration>(_input_config.config())),
    m_subset_mode(false) {

    if(m_counter.size() != _input_config.sites().size())
      m_subset_mode = true;


    m_current->set_occupation(m_counter());
    reset_properties(*m_current);
    this->_initialize(&(*m_current));

    // Make sure that current() is a primitive canonical config
    if(!_check_current()) {
      increment();
    }

    // set step to 0
    if(valid()) {
      _set_step(0);
    }
    m_current->set_source(this->source(step()));
  }

  /// Implements _increment over all occupations
  void ConfigEnumAllOccupations::increment() {

    bool is_valid_config {false};

    while(!is_valid_config && ++m_counter) {
      for(Index l : m_selection) {
        m_current->set_occ(l, m_counter[l]);
      }
      is_valid_config = _check_current();
    }

    if(m_counter.valid()) {
      this->_increment_step();
    }
    else {
      this->_invalidate();
    }
    m_current->set_source(this->source(step()));
  }

  /// Returns true if current() is primitive and canonical
  bool ConfigEnumAllOccupations::_check_current() const {
    return current().is_primitive() && (!m_subset_mode && current().is_canonical());
  }

}
