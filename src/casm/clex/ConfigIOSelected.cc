#include "casm/clex/ConfigIOSelected.hh"

namespace CASM {
  namespace ConfigIO {

    void Selected::init(const Configuration &_tmplt) const {
      if(m_selection.size() == 0) {
        if(m_selection_name.empty())
          m_selection_name = "MASTER";
        m_selection = ConstConfigSelection(_tmplt.get_primclex(), m_selection_name);
      }
      else if(m_selection_name.empty()) {
        m_selection_name = m_selection.name();
        if(m_selection_name.empty())
          m_selection_name = "unknown";
      }
    }

    std::unique_ptr<Selected> Selected::clone() const {
      return std::unique_ptr<Selected>(this->_clone());
    }

    Selected *Selected::_clone() const {
      return new Selected(*this);
    }

    std::string Selected::short_header(const Configuration &_config) const {
      return name() + "(" + m_selection_name + ")";
    }

    bool Selected::evaluate(const Configuration &_config) const {
      return m_selection.selected(_config.name());
    }

    bool Selected::parse_args(const std::string &args) {
      if(m_selection.size() || m_selection_name.size()) {
        return false;
      }

      m_selection_name = args;
      return true;
    }
  }
}
