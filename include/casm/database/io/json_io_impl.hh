#ifndef CASM_database_json_io_impl
#define CASM_database_json_io_impl

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/Selection.hh"
#include "casm/database/io/json_io.hh"
#include "casm/misc/TypeInfo.hh"

namespace CASM {
  namespace DB {

    /// \brief Make a Selection from JSON input
    ///
    /// \param db The database to generate the selection from
    /// \param kwargs jsonParser with JSON input
    /// \param name_key Read object names from kwargs[name_key]. Expects array of string.
    /// \param sel_key Read selection name for kwargs[sel_key]. Expect string.
    /// \param on_error Indicates how to handle names that do not exist in database. OnError::WARN is
    /// not allowed and treated as OnError::THROW.
    ///
    /// Notes:
    /// - Reads selection from 'sel_key' first, using "NONE" if 'sel_key' does not
    ///   exist. Then reads names from 'name_key' array and adds them to the selection.
    /// - If 'name_key' and 'sel_key' do not exist, throw exception.
    ///
    template<typename DataObject>
    Selection<DataObject> make_selection(
      Database<DataObject> &db,
      const jsonParser &kwargs,
      std::string name_key,
      std::string sel_key,
      OnError on_error) {

      if(!kwargs.contains(name_key) && !kwargs.contains(sel_key)) {
        std::string msg = "Error in make_selection<" + type_name<DataObject>() + ">: One of " +
                          name_key + " or " + sel_key + " must be given.";
        throw std::runtime_error(msg);
      }

      std::vector<std::string> obj_names;
      std::string sel_name;
      kwargs.get_else(obj_names, name_key, std::vector<std::string>());
      kwargs.get_else(sel_name, sel_key, std::string("NONE"));

      DB::Selection<DataObject> sel(db, sel_name);
      for(const auto &name : obj_names) {

        // validate obj_names & handle errors
        if(!sel.db().count(name)) {
          if(on_error == OnError::THROW || on_error == OnError::WARN) {
            std::string msg = "Error in make_selection<" + type_name<DataObject>() + ">: \"" +
                              name + "\" is not in the database.";
            throw std::runtime_error(msg);
          }
          else if(on_error == OnError::CONTINUE) {
            continue;
          }
          else {
            std::string msg = "Error in make_selection<" + type_name<DataObject>() + ">: Unknown error";
            throw std::runtime_error(msg);
          }
        }
        sel.data()[name] = true;
      }
      return sel;
    }

  }
}

#endif
