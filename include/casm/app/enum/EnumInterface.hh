#ifndef CASM_Enumerator
#define CASM_Enumerator

#include <string>

#include "casm/app/casm_functions.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/Configuration.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"

namespace CASM {

  namespace Completer {
    class EnumOptionBase;
    class EnumOption;
  }

  typedef InterfaceBase<Completer::EnumOption> EnumInterfaceBase;
  typedef InterfaceMap<Completer::EnumOption> EnumeratorMap;

  jsonParser &to_json(const Completer::EnumOption &enum_opt, jsonParser &json);

  // might be useful for other casm commands...
  template<typename OptionType>
  jsonParser make_json_input(const OptionType &opt);

  /// \brief Standardizes parsing casm enum input options to make ScelEnum JSON input
  jsonParser make_enumerator_scel_enum_input(
    jsonParser kwargs,
    const Completer::EnumOptionBase &enum_opt);

  /// \brief Standardizes parsing casm enum input options to make an ScelEnumProps
  ScelEnumProps make_enumerator_scel_enum_props(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOptionBase &enum_opt);

  /// \brief Standardizes parsing casm enum input options to make a SuperlatticeEnumerator
  std::unique_ptr<SuperlatticeEnumerator> make_enumerator_superlat_enum(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOptionBase &enum_opt);

  /// \brief Standardizes parsing casm enum input options to make an ScelEnum
  std::unique_ptr<ScelEnum> make_enumerator_scel_enum(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOptionBase &enum_opt);

  /// \brief Standardizes parsing casm enum filter expressions
  std::vector<std::string> make_enumerator_filter_expr(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  /// \brief Get dry-run value (default=false)
  bool dry_run(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  inline std::string dry_run_msg(bool dry_run) {
    return dry_run ? "(dry run) " : "";
  }

  struct InsertAlgorithmsOptions {

    std::string method;
    std::vector<std::string> filter_expr;
    bool dry_run = false;
    COORD_TYPE coord_type = COORD_TYPE::FRAC;
    bool primitive_only = true;
    int verbosity = 10;
  };



  /// \brief Standardizes insertion from enumerators that construct unique
  /// primitive canonical configurations
  template<typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    const PrimClex &primclex,
    ConfigEnumInput const &in_config,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool dry_run);

  /// \brief Standardizes insertion from enumerators that construct unique
  /// primitive canonical configurations
  template<typename InConfigIterator, typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    const PrimClex &primclex,
    InConfigIterator begin,
    InConfigIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool dry_run);

  /// \brief Standardizes insertion from enumerators that construct configurations
  template<typename ConfigEnumConstructor>
  int insert_configs(
    std::string method,
    const PrimClex &primclex,
    ConfigEnumInput const &in_config,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only,
    bool dry_run);

  /// \brief Standardizes insertion from enumerators that construct configurations
  template<typename InConfigIterator, typename ConfigEnumConstructor>
  int insert_configs(
    std::string method,
    const PrimClex &primclex,
    InConfigIterator begin,
    InConfigIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only,
    bool dry_run);

  /// \brief Standardizes insertion from enumerators that construct configurations
  template<typename LatticeIterator, typename ConfigEnumConstructor>
  int insert_configs_via_lattice_enum(
    std::string method,
    const PrimClex &primclex,
    LatticeIterator begin,
    LatticeIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only,
    bool dry_run);

  /// \brief Template class to be specialized for each enumerator that may be accessed via the API
  template<typename Derived>
  class EnumInterface : public EnumInterfaceBase {

  public:

    std::string help() const override {
      return Derived::interface_help();
    }

    std::string name() const override {
      return Derived::enumerator_name;
    }

    int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt, EnumeratorMap const *interface_map = nullptr) const override {
      return Derived::run(primclex, kwargs, enum_opt, interface_map);
    }

    std::unique_ptr<EnumInterfaceBase> clone() const {
      return std::unique_ptr<EnumInterfaceBase>(this->_clone());
    }

  private:

    EnumInterfaceBase *_clone() const override {
      return new EnumInterface<Derived>(*this);
    }


  };

  /// \brief Use to construct an InterfaceMap
  std::unique_ptr<InterfaceMap<Completer::EnumOption> > make_enumerator_map();

  std::unique_ptr<EnumeratorMap> make_standard_enumerator_map();

  std::vector<ConfigEnumInput> make_enumerator_input_configs(
    PrimClex const &primclex,
    jsonParser const &_kwargs,
    Completer::EnumOptionBase const &enum_opt,
    EnumeratorMap const *interface_map);

}

#endif
