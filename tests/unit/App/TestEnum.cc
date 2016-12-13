#include "TestEnum.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/container/Enumerator_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_TestEnum_interface() {
    return new CASM::EnumInterface<CASM::TestEnum>();
  }
}

namespace CASM {

  const std::string TestEnum::enumerator_name = "TestEnum";

  const std::string TestEnum::interface_help =
    "TestEnum: \n\n"

    "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
    "    Indicate supercells to enumerate all occupational configurations in. May \n"
    "    be a JSON array of supercell names, or a JSON object specifying          \n"
    "    supercells in terms of size and unit cell. By default, all existing      \n"
    "    supercells are used. See 'ScelEnum' description for details.         \n\n"

    "  filter: string (optional, default=None)\n"
    "    A query command to use to filter which Configurations are kept.          \n"
    "\n"
    "  Examples:\n"
    "    To enumerate all occupations in supercells up to and including size 4:\n"
    "      casm enum --method TestEnum -i '{\"supercells\": {\"max\": 4}}' \n"
    "\n"
    "    To enumerate all occupations in all existing supercells:\n"
    "      casm enum --method TestEnum\n"
    "\n"
    "    To enumerate all occupations in all particular supercells:\n"
    "      casm enum --method TestEnum -i \n"
    "      '{ \n"
    "        \"supercells\": { \n"
    "          \"name\": [\n"
    "            \"SCEL1_1_1_1_0_0_0\",\n"
    "            \"SCEL2_1_2_1_0_0_0\",\n"
    "            \"SCEL4_1_4_1_0_0_0\"\n"
    "          ]\n"
    "        } \n"
    "      }' \n\n";

  int TestEnum::run(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    std::unique_ptr<ScelEnum> scel_enum = make_enumerator_scel_enum(primclex, _kwargs, enum_opt);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    auto lambda = [&](Supercell & scel) {
      return notstd::make_unique<TestEnum>(scel);
    };

    int returncode = insert_unique_canon_configs(
                       enumerator_name,
                       primclex,
                       scel_enum->begin(),
                       scel_enum->end(),
                       lambda,
                       filter_expr);

    return returncode;
  }


  /// \brief Construct with a Supercell, using all permutations
  TestEnum::TestEnum(Supercell &_scel) :
    m_counter(
      Array<int>(_scel.num_sites(), 0),
      _scel.max_allowed_occupation(),
      Array<int>(_scel.num_sites(), 1)) {

    m_current = notstd::make_cloneable<Configuration>(_scel, this->source(0), m_counter());
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
    _current().set_source(this->source(step()));
  }

  /// Implements _increment over all occupations
  void TestEnum::increment() {

    bool is_valid_config {false};

    while(!is_valid_config && ++m_counter) {

      _current().set_occupation(m_counter());
      is_valid_config = _check_current();
    }

    if(m_counter.valid()) {
      this->_increment_step();
    }
    else {
      this->_invalidate();
    }
    _current().set_source(this->source(step()));
  }

  /// Returns true if current() is primitive and canonical
  bool TestEnum::_check_current() const {
    return current().is_primitive() && current().is_canonical();
  }

}
