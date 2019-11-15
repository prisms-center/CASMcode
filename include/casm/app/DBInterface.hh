#ifndef CASM_DBInterface
#define CASM_DBInterface

#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"

namespace CASM {

  class APICommandBase;
  class QueryCommand;
  class SelectCommand;

  namespace DB {

    template<typename T> class Selection;

    template<typename CommandType, template<typename> class ImplType = CommandType::template Impl>
    struct ConstructImpl {

      typedef typename CommandType::ImplBase ImplBase;
      template<typename DataObject> using Impl = ImplType<DataObject>;

      ConstructImpl(
        std::unique_ptr<ImplBase> &_impl,
        const CommandType &_cmd) : impl(_impl), cmd(_cmd) {}

      template<typename T>
      void eval() {
        impl.reset(new Impl<T>(cmd));
      }

      std::unique_ptr<ImplBase> &impl;
      const CommandType &cmd;
    };


    template<typename DataObject>
    class InterfaceData {
    public:

      InterfaceData(const QueryCommand &cmd);
      InterfaceData(const SelectCommand &cmd);

      const DataFormatterDictionary<DataObject> &dict() const {
        return *m_dict;
      }

      std::string sel_str() const {
        return m_ss.str();
      }

      double sel_size() const {
        return m_sel.size();
      }

      Selection<DataObject> &sel(Index i = 0) {
        return *m_sel[i];
      }

      const Selection<DataObject> &sel(Index i = 0) const {
        return *m_sel[i];
      }

    private:


      void _make_dict(const APICommandBase &cmd);


      const DataFormatterDictionary<DataObject> *m_dict;

      std::unique_ptr<DataFormatterDictionary<DataObject>> m_standard_dict;

      std::vector<std::unique_ptr<Selection<DataObject>>> m_sel;

      std::stringstream m_ss;
    };
  }
}

#endif
