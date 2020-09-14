#include "TestEnum_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/app/enum/EnumInterface_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_TestEnum_interface() {
    return new CASM::EnumInterface<CASM::TestEnum>();
  }
}

namespace CASM {

  const std::string TestEnum::enumerator_name = "TestEnum";

  std::string TestEnum::interface_help() {
    return
      "TestEnum: \n\n"

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
      "    Restricts enumeration to specified sites. Sites are specified in [b,i,j,k] convention,\n"
      "    where 'b' is sublattice index and [i,j,k] specifies linear combinations of primitive-\n"
      "    cell lattice vectors.\n"
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
      "          \"names\": [\n"
      "            \"SCEL1_1_1_1_0_0_0\",\n"
      "            \"SCEL2_1_2_1_0_0_0\",\n"
      "            \"SCEL4_1_4_1_0_0_0\"\n"
      "          ]\n"
      "        } \n"
      "      }' \n\n";
  }

  int TestEnum::run(
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
      std::vector<int> max = _config.configuration().supercell().max_allowed_occupation();

      std::vector<int> maxselect;
      for(Index i : _config.sites())
        maxselect.push_back(max[i]);

      //std::cout << "Revealing maxselect: " << maxselect << "\n";
      if(maxselect.empty())
        return max;
      else
        return maxselect;
    }
  }

  /// \brief Construct with a Supercell, using all permutations
  TestEnum::TestEnum(const ConfigEnumInput &_input_config) :
    m_selection(_input_config.sites().begin(), _input_config.sites().end()),
    m_counter(
      std::vector<int>(_input_config.sites().size(), 0),
      local_impl::max_selected_occupation(_input_config),
      std::vector<int>(_input_config.sites().size(), 1)),

    m_current(notstd::make_cloneable<Configuration>(_input_config.configuration())),
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
  void TestEnum::increment() {
    //std::cout << "Incrementing!\n";
    bool is_valid_config {false};

    while(!is_valid_config && ++m_counter) {
      for(Index l : m_selection) {
        m_current->set_occ(l, m_counter[l]);
      }
      //std::cout << "Set occupation to : " << m_current->occupation() << "\n";
      is_valid_config = _check_current();
    }
    if(m_counter.valid()) {
      //std::cout << "Escaped loop, valid state\n";

      this->_increment_step();
    }
    else {
      //std::cout << "Escaped loop, invalid state\n";
      this->_invalidate();
    }
    m_current->set_source(this->source(step()));
  }

  /// Returns true if current() is primitive and canonical
  bool TestEnum::_check_current() const {
    return current().is_primitive() && (m_subset_mode || current().is_canonical());
  }

}
