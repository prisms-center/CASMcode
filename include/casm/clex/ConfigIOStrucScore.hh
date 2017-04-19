#ifndef CASM_ConfigIOStrucScore
#define CASM_ConfigIOStrucScore

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"

namespace CASM {

  class PrimClex;
  class Configuration;
  class ConfigMapper;

  namespace ConfigIO {

    /// \brief Evaluates the mapping of a configuration onto an arbitrary primitive structure
    ///
    /// \ingroup ConfigIO
    ///
    class StrucScore: public VectorXdAttribute<Configuration> {

    public:

      StrucScore();

      StrucScore(const StrucScore &RHS);

      // --- Required implementations -----------

      std::unique_ptr<StrucScore> clone() const {
        return std::unique_ptr<StrucScore>(this->_clone());
      }

      Eigen::VectorXd evaluate(const Configuration &_config) const override;


      // --- Specialized implementation -----------

      bool validate(const Configuration &_config) const override;

      std::string short_header(const Configuration &_config) const override;

      std::vector<std::string> col_header(const Configuration &_config) const override;

      bool parse_args(const std::string &args) override;

    protected:
      mutable std::unique_ptr<PrimClex> m_altprimclex;
      mutable std::unique_ptr<ConfigMapper> m_configmapper;
      fs::path m_prim_path;
      std::vector<std::string> m_prop_names;

    private:
      /// \brief Clone
      StrucScore *_clone() const override {
        return new StrucScore(*this);
      }
    };

  }
}
#endif

