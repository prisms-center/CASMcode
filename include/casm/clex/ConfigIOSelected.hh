#ifndef CONFIGIOSELECTED_HH
#define CONFIGIOSELECTED_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/clex/ConfigSelection.hh"
//#include "casm/app/casm_functions.hh"

namespace CASM {
  class Configuration;

  namespace ConfigIO {

    /// \brief Returns true if configuration is specified in given selection (default: MASTER)
    ///
    /// Ex: 'selected_in', 'selected_in(myselection.txt)'
    ///
    /// \ingroup ConfigIO
    ///
    class Selected : public BooleanAttribute<Configuration> {
    public:
      Selected()
        : BooleanAttribute<Configuration>("selected_in", "Returns true if configuration is specified in given selection (default MASTER). Ex: 'selected_in(myselection.txt)'") {}

      Selected(const Selected &) = default;

      template<bool IsConst>
      explicit Selected(ConfigSelection<IsConst> _selection):
        BooleanAttribute<Configuration>("selected_in", "Returns true if configuration is specified in given selection (default MASTER). Ex: 'selected_in(myselection.txt)'"),
        m_selection(_selection) {}


      // --- Required implementations -----------

      std::unique_ptr<Selected> clone() const;

      bool evaluate(const Configuration &_config) const override;


      // --- Specialized implementation -----------

      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;
      /*
            void inject(const Configuration &_config, DataStream &_stream, Index) const override;

            void print(const Configuration &_config, std::ostream &_stream, Index) const override;

            jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;
      */
      bool parse_args(const std::string &args) override;

    private:

      /// \brief Clone
      Selected *_clone() const override;

      mutable std::string m_selection_name;
      mutable ConstConfigSelection m_selection;
    };
  }
}
#endif
