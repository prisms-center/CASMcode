#ifndef CASM_ConfigEnumAllOccupations_impl
#define CASM_ConfigEnumAllOccupations_impl

#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/enumerator/Enumerator_impl.hh"

namespace CASM {

  template<typename InConfigIterator>
  int ConfigEnumAllOccupations::run(
    const PrimClex &primclex,
    InConfigIterator begin,
    InConfigIterator end,
    const std::vector<std::string> &filter_expr,
    bool dry_run) {

    auto lambda = [&](const ConfigEnumInput & in_config) {
      return notstd::make_unique<ConfigEnumAllOccupations>(in_config);
    };

    int returncode = insert_unique_canon_configs(
                       enumerator_name,
                       primclex,
                       begin,
                       end,
                       lambda,
                       filter_expr,
                       dry_run);

    return returncode;
  }
}

#endif
