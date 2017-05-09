#ifndef CASM_ConfigEnumAllOccupations_impl
#define CASM_ConfigEnumAllOccupations_impl

#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/container/Enumerator_impl.hh"

namespace CASM {

  template<typename ScelIterator>
  int ConfigEnumAllOccupations::run(
    const PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    const std::vector<std::string> &filter_expr) {

    auto lambda = [&](const Supercell & scel) {
      return notstd::make_unique<ConfigEnumAllOccupations>(scel);
    };

    int returncode = insert_unique_canon_configs(
                       enumerator_name,
                       primclex,
                       begin,
                       end,
                       lambda,
                       filter_expr);

    return returncode;
  }
}

#endif
