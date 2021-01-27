#ifndef CASM_DBInterface_impl
#define CASM_DBInterface_impl

#include "casm/app/DBInterface.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/query.hh"
#include "casm/app/select.hh"
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/Selection.hh"

namespace CASM {
namespace DB {

template <typename DataObject>
InterfaceData<DataObject>::InterfaceData(const QueryCommand &cmd) {
  // set 'selected' column
  if (cmd.in_project() && cmd.vm().count("selection")) {
    m_sel.resize(1);
    m_sel[0].reset(new Selection<DataObject>(cmd.primclex().db<DataObject>(),
                                             cmd.opt().selection_path()));
    m_ss << cmd.opt().selection_path();
    cmd.primclex().settings().query_handler<DataObject>().set_selected(sel(0));
  }
  _make_dict(cmd);
}

template <typename DataObject>
InterfaceData<DataObject>::InterfaceData(const SelectCommand &cmd) {
  // set 'selected' column
  if (cmd.in_project() && cmd.vm().count("selections")) {
    m_sel.resize(cmd.opt().selection_paths().size());
    for (int i = 0; i < cmd.opt().selection_paths().size(); ++i) {
      if (i != 0) {
        m_ss << ", ";
      }
      m_ss << cmd.opt().selection_paths()[i];
      m_sel[i].reset(new Selection<DataObject>(cmd.primclex().db<DataObject>(),
                                               cmd.opt().selection_paths()[i]));
    }

    cmd.primclex().settings().query_handler<DataObject>().set_selected(sel(0));
  }
  _make_dict(cmd);
}

template <typename DataObject>
void InterfaceData<DataObject>::_make_dict(const APICommandBase &cmd) {
  // set 'selected' column
  if (cmd.in_project()) {
    m_dict = &cmd.primclex().settings().query_handler<DataObject>().dict();
  } else {
    m_standard_dict.reset(new DataFormatterDictionary<DataObject>());
    *m_standard_dict = make_dictionary<DataObject>();
    m_dict = m_standard_dict.get();
  }
}

}  // namespace DB
}  // namespace CASM

#endif
