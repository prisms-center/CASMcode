#ifndef CONFIGIOSTRUCSCORE_HH
#define CONFIGIOSTRUCSCORE_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigMapping.hh"


namespace CASM {

  class Configuration;

  namespace ConfigIO {

    /// \brief Evaluates the mapping of a configuration onto an arbitrary primitive structure
    ///
    /// \ingroup ConfigIO
    ///
    class StrucScore: public VectorXdAttribute<Configuration> {
    public:
      StrucScore() :
        VectorXdAttribute<Configuration>("struc_score", "Evaluates the mapping of a configuration onto an arbitrary primitive structure, specified by its path. Allowed options are [ 'basis_score' (mean-square site displacement) | 'lattice_score' (lattice deformation metric having units Angstr.^2) | 'total_score' (w*lattice_score+(1.0-w)*basis_score) ].  The struc_score weighting parameter 'w' can be provided as an optional decimal parameter from 0.0 to 1.0 (default 0.5). Ex: struc_score(path/to/PRIM, basis_score, 0.4)"),
        m_configmapper(ConfigMapper(ConfigMapper::null_initializer)) {};

      StrucScore(const StrucScore &RHS) :
        VectorXdAttribute<Configuration>(RHS),
        m_altprimclex((RHS.m_altprimclex == nullptr) ? nullptr : new PrimClex(*RHS.m_altprimclex)),
        m_configmapper(RHS.m_configmapper),
        m_prim_path(RHS.m_prim_path),
        m_prop_names(RHS.m_prop_names) {
        m_configmapper.set_primclex(*m_altprimclex);
      }

      // --- Required implementations -----------

      std::unique_ptr<StrucScore> clone()const {
        return std::unique_ptr<StrucScore>(this->_clone());
      }

      Eigen::VectorXd evaluate(const Configuration &_config) const override;


      // --- Specialized implementation -----------

      bool validate(const Configuration &_config) const override;

      std::string short_header(const Configuration &_config) const override;

      std::vector<std::string> col_header(const Configuration &_config) const override;
      /*
            void inject(const Configuration &_config, DataStream &_stream, Index) const override;

            void print(const Configuration &_config, std::ostream &_stream, Index) const override;

            jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;
      */
      bool parse_args(const std::string &args) override;

    protected:
      mutable std::unique_ptr<PrimClex> m_altprimclex;
      mutable ConfigMapper m_configmapper;
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

