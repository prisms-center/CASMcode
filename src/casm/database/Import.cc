#include "casm/database/Import_impl.hh"

typedef std::back_insert_iterator<std::vector<CASM::fs::path> > vector_path_back_inserter;
typedef std::insert_iterator<std::set<CASM::fs::path> > set_path_inserter;

namespace CASM {
  namespace DB {

    template std::pair<vector_path_back_inserter, int>
    construct_pos_paths<vector_path_back_inserter>(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      vector_path_back_inserter result);

    template std::pair<set_path_inserter, int>
    construct_pos_paths<set_path_inserter>(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      set_path_inserter result);

  }
}
