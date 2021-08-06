#include "casm/database/Import_impl.hh"

typedef std::back_insert_iterator<std::vector<CASM::fs::path> >
    vector_path_back_inserter;
typedef std::insert_iterator<std::set<CASM::fs::path> > set_path_inserter;

namespace CASM {
namespace DB {

jsonParser &to_json(ImportSettings const &_set, jsonParser &_json) {
  //_json["import_configurations"] = _set.import_configurations;
  //_json["primitive_only"] = _set.primitive_only;
  _json["import_properties"] = _set.import_properties;
  _json["copy_structure_files"] = _set.copy_structure_files;
  _json["copy_additional_files"] = _set.copy_additional_files;
  _json["overwrite"] = _set.overwrite;
  //_json["output_as_json"] = _set.output_as_json;
  return _json;
}

jsonParser const &from_json(ImportSettings &_set, jsonParser const &_json) {
  _set = ImportSettings();

  // _json.get_if(_set.import_configurations, "import_configurations");
  //_json.get_if(_set.primitive_only, "primitive_only");
  _json.get_if(_set.import_properties, "import_properties");
  _json.get_if(_set.copy_structure_files, "copy_structure_files");
  _json.get_if(_set.copy_additional_files, "copy_additional_files");
  _json.get_if(_set.overwrite, "overwrite");
  //_json.get_if(_set.output_as_json, "output_as_json");

  return _json;
}

template std::pair<vector_path_back_inserter, int>
construct_pos_paths<vector_path_back_inserter>(
    const PrimClex &primclex, const Completer::ImportOption &import_opt,
    vector_path_back_inserter result);

template std::pair<set_path_inserter, int>
construct_pos_paths<set_path_inserter>(
    const PrimClex &primclex, const Completer::ImportOption &import_opt,
    set_path_inserter result);

}  // namespace DB
}  // namespace CASM
