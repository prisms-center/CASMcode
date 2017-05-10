#ifndef CASM_TestEnum_impl
#define CASM_TestEnum_impl

#include "TestEnum.hh"

namespace CASM {

  template<typename ScelIterator>
  int TestEnum::run(
    const PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    const std::vector<std::string> &filter_expr) {

    auto lambda = [&](const Supercell & scel) {
      return notstd::make_unique<TestEnum>(scel);
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
