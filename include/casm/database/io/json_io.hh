#ifndef CASM_database_json_io
#define CASM_database_json_io

#include <string>

#include "casm/global/enum.hh"

namespace CASM {

class jsonParser;
namespace DB {
template <typename T>
class Database;
template <typename T>
class Selection;
}  // namespace DB

namespace DB {
/// \brief Make a DB::Selection from JSON input
template <typename DataObject>
DB::Selection<DataObject> make_selection(DB::Database<DataObject> &db,
                                         const jsonParser &kwargs,
                                         std::string name_key,
                                         std::string sel_key,
                                         OnError on_error = OnError::THROW);
}  // namespace DB
}  // namespace CASM

#endif
