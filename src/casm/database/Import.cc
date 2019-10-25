#include "casm/database/Import_impl.hh"

typedef std::back_insert_iterator<std::vector<CASM::fs::path> > vector_path_back_inserter;
typedef std::insert_iterator<std::set<CASM::fs::path> > set_path_inserter;

namespace CASM {
  namespace DB {

    jsonParser &to_json(ImportSettings const &_set, jsonParser &_json) {
      _json["data"] = _set.data;
      _json["additional_files"] = _set.additional_files;
      _json["overwrite"] = _set.overwrite;
      return _json;
    }

    jsonParser const &from_json(ImportSettings &_set, jsonParser const &_json) {
      _set.set_default();

      if(_json.contains("data"))
        _set.data = _json["data"].get<bool>();

      if(_json.contains("additional_files"))
        _set.additional_files = _json["additional_files"].get<bool>();

      if(_json.contains("overwrite"))
        _set.overwrite = _json["overwrite"].get<bool>();

      return _json;
    }

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
