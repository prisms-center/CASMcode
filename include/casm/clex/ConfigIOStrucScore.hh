#ifndef CASM_ConfigIOStrucScore
#define CASM_ConfigIOStrucScore

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {

  class Structure;
  class Configuration;
  class StrucMapper;

  namespace ConfigIO {

    /// \brief Evaluates the mapping of a configuration onto an arbitrary primitive structure
    ///
    /// \ingroup ConfigIO
    ///
    class StrucScore: public VectorXdAttribute<Configuration> {

    public:

      StrucScore();

      StrucScore(StrucScore const &RHS);

      // --- Required implementations -----------

      std::unique_ptr<StrucScore> clone() const {
        return std::unique_ptr<StrucScore>(this->_clone());
      }

      Eigen::VectorXd evaluate(Configuration const &_config) const override;


      // --- Specialized implementation -----------

      bool validate(Configuration const &_config) const override;

      std::string short_header(Configuration const &_config) const override;

      std::vector<std::string> col_header(Configuration const &_config) const override;

      bool parse_args(std::string const &args) override;

      bool init(Configuration const &tmplt) const override;

    protected:
      mutable std::unique_ptr<BasicStructure<Site> > m_altprim;
      mutable std::unique_ptr<StrucMapper> m_strucmapper;
      double m_lattice_weight;
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

