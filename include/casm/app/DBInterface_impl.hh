#ifndef CASM_DBInterface_impl
#define CASM_DBInterface_impl

#include "casm/app/DBInterface.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/casm_io/DataFormatter.hh"

namespace CASM {
  namespace DB {

    template<typename DataObject>
    void InterfaceData<DataObject>::_make_dict(const APICommandBase &cmd) {

      // set 'selected' column
      if(cmd.in_project()) {
        m_dict = &cmd.primclex().settings().query_handler<DataObject>().dict();
      }
      else {
        m_standard_dict.reset(new DataFormatterDictionary<DataObject>());
        *m_standard_dict = make_dictionary<DataObject>();
        m_dict = m_standard_dict.get();
      }
    }

  }
}

#endif

