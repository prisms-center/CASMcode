#ifndef CASM_TestEnum_impl
#define CASM_TestEnum_impl

#include "TestEnum.hh"
#include "casm/app/enum/EnumInterface_impl.hh"

namespace CASM {

  template<typename InConfigIterator>
  int TestEnum::run(
    const PrimClex &primclex,
    InConfigIterator begin,
    InConfigIterator end,
    const std::vector<std::string> &filter_expr,
    bool dry_run) {

    auto lambda = [&](const ConfigEnumInput & in_config) {
      return notstd::make_unique<TestEnum>(in_config);
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
