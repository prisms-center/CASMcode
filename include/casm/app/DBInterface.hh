#ifndef CASM_DBInterface
#define CASM_DBInterface

#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/DataFormatter.hh"

namespace CASM {

  class APICommandBase;

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

      template<typename CommandType>
      InterfaceData(const CommandType &cmd);

      /*
      InterfaceData(SelectCommand& cmd) {

        // set 'selected' column
        if(cmd.in_project() && cmd.vm.count("selections")) {
          m_sel.resize(cmd.opt().selection_paths().size());
          for(int i=0; i < cmd.opt().selection_paths().size(); ++i) {
            if(i!=0) {
              m_ss << ", ";
            }
            m_ss << cmd.opt().selection_paths()[i];
            m_sel[i].reset(
              new Selection<DataObject>(
                cmd.primclex().db<DataObject>(),
                cmd.opt().selection_paths()[i]));
          }

          cmd.primclex().settings().query_handler<DataObject>().set_selected(sel(0));
        }
      }

      InterfaceData(QueryCommand& cmd) {

        // set 'selected' column
        if(cmd.in_project() && cmd.vm().count("selection")) {
          m_sel.resize(1);
          m_sel[0].reset(
            new Selection<DataObject>(
              cmd.primclex().db<DataObject>(),
              cmd.opt().selection_path()));
          m_ss << cmd.opt().selection_path();
          cmd.primclex().settings().query_handler<DataObject>().set_selected(sel(0));
        }
      }

      InterfaceData(RmCommand& cmd) {
        // set 'selected' column
        if(cmd.in_project() && cmd.vm.count("selection")) {
            m_sel.resize(1);
            m_sel[0].reset(
              new Selection<DataObject>(
                cmd.primclex().db<DataObject>(),
                cmd.opt().selection_path()));
            m_ss << cmd.opt().selection_path();
          }
        }
      }
      */

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
