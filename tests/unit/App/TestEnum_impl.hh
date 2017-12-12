#ifndef CASM_TestEnum_impl
#define CASM_TestEnum_impl

#include "TestEnum.hh"
#include "casm/enumerator/Enumerator_impl.hh"

namespace CASM {

  template<typename ScelIterator>
  int TestEnum::run(
    const PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    const std::vector<std::string> &filter_expr,
    bool dry_run) {

    auto lambda = [&](const Supercell & scel) {
      return notstd::make_unique<TestEnum>(scel);
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
