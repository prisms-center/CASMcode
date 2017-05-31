#include "casm/app/APICommand.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  APICommandBase::APICommandBase(const CommandArgs &_args) :
    Logging(_args),
    m_args(_args),
    m_in_project(!find_casmroot(_args.root).empty()) {}

  const CommandArgs &APICommandBase::args() const {
    return m_args;
  }

  fs::path APICommandBase::root() const {
    return m_args.root;
  }

  bool APICommandBase::in_project() const {
    return m_in_project;
  }

  PrimClex &APICommandBase::primclex() const {
    return primclex(log());
  }

  PrimClex &APICommandBase::primclex(Log &status_log) const {
    if(m_args.primclex) {
      return *m_args.primclex;
    }
    else if(!m_primclex) {
      m_primclex.reset(new PrimClex(root(), status_log));
    }
    return *m_primclex;
  }
}
