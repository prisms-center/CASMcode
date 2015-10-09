#include "casm/clex/ConfigIOSelected.hh"

namespace CASM {
  namespace ConfigIO_impl {

    void SelectedConfigFormatter::init(const Configuration &_tmplt) const {
      if(m_selection.size() == 0) {
        if(m_selection_name.empty())
          m_selection_name = "MASTER";
        if(m_selection_name == "MASTER") {
          m_selection = ConstConfigSelection(_tmplt.get_primclex());
        }
        else if(fs::exists(m_selection_name)) {
          m_selection = ConstConfigSelection(_tmplt.get_primclex(), m_selection_name);
        }
        else {
          throw std::runtime_error("ERROR: Configuration selection " + m_selection_name + " does not refer to 'MASTER' list or a valid selection file.\n");
        }
      }
      else if(m_selection_name.empty()) {
        m_selection_name = m_selection.name();
        if(m_selection_name.empty())
          m_selection_name = "unknown";
      }
    }

    std::string SelectedConfigFormatter::short_header(const Configuration &_config) const {
      return name() + "(" + m_selection_name + ")";
    }

    void SelectedConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      _stream << m_selection.selected(_config.name());
    }

    void SelectedConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      _stream << std::string(name().size() + m_selection_name.size() + 1, ' ') <<  m_selection.selected(_config.name());
    }

    jsonParser &SelectedConfigFormatter::to_json(const Configuration &_config, jsonParser &json) const {
      return json = m_selection.selected(_config.name());
    }

    bool SelectedConfigFormatter::parse_args(const std::string &args) {
      if(m_selection.size() || m_selection_name.size()) {
        return false;
      }

      m_selection_name = args;
      return true;
    }
  }
}
