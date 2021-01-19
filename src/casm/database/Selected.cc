#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/Selected_impl.hh"

// explicit template instantiations
#define INST_Selected(r, data, type) template class Selected<type>;

namespace CASM {
namespace DB {

template <typename ObjType>
bool Selected<ObjType>::init(const ObjType &_tmplt) const {
  if (!m_selection) {
    if (m_selection_name.empty()) {
      m_selection_name = "MASTER";
    }
    m_selection = notstd::make_cloneable<Selection<ObjType> >(
        _tmplt.primclex().template db<ObjType>(), m_selection_name);
  } else if (m_selection_name.empty()) {
    m_selection_name = m_selection->name();
    if (m_selection_name.empty()) {
      m_selection_name = "unknown";
    }
  }
  return true;
}

template <typename ObjType>
std::unique_ptr<Selected<ObjType> > Selected<ObjType>::clone() const {
  return std::unique_ptr<Selected>(this->_clone());
}

template <typename ObjType>
Selected<ObjType> *Selected<ObjType>::_clone() const {
  return new Selected(*this);
}

template <typename ObjType>
std::string Selected<ObjType>::short_header(const ObjType &_obj) const {
  return this->name() + "(" + m_selection_name + ")";
}

template <typename ObjType>
bool Selected<ObjType>::evaluate(const ObjType &_obj) const {
  if (_obj.name().find("equiv") != std::string::npos) {
    return false;
  }
  return m_selection->is_selected(_obj.name());
}

template <typename ObjType>
bool Selected<ObjType>::parse_args(const std::string &args) {
  if ((m_selection && m_selection->data().size()) || m_selection_name.size()) {
    return false;
  }

  m_selection_name = args;
  return true;
}

BOOST_PP_SEQ_FOR_EACH(INST_Selected, _, CASM_DB_TYPES)
}  // namespace DB
}  // namespace CASM
