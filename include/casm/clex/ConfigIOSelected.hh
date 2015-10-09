#ifndef CONFIGIOSELECTED_HH
#define CONFIGIOSELECTED_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/clex/ConfigSelection.hh"

namespace CASM {
  namespace ConfigIO_impl {
    class SelectedConfigFormatter : public BaseDatumFormatter<Configuration> {
    public:
      SelectedConfigFormatter()
        : BaseDatumFormatter<Configuration>("selected_in", "Returns true if configuration is specified in given selection (default MASTER). Ex: 'selected_in(myselection.txt)'") {}

      SelectedConfigFormatter(const SelectedConfigFormatter &) = default;

      template<bool IsConst>
      SelectedConfigFormatter(ConfigSelection<IsConst> _selection):
        BaseDatumFormatter<Configuration>("selected_in", "Returns true if configuration is specified in given selection (default MASTER). Ex: 'selected_in(myselection.txt)'"),
        m_selection(_selection) {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new SelectedConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);

    private:
      mutable std::string m_selection_name;
      mutable ConstConfigSelection m_selection;
    };
  }
}
#endif
