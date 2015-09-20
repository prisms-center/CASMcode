#ifndef CONFIGIOSTRUCSCORE_HH
#define CONFIGIOSTRUCSCORE_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/PrimClex.hh"


namespace CASM {

  class Configuration;

  namespace ConfigIO_impl {
    /*
     */

    class StrucScoreConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      StrucScoreConfigFormatter() :
        BaseDatumFormatter<Configuration>("struc_score", "Evaluates the mapping of a configuration onto an arbitrary primitive structure, specified by it's path. Allowed options are 'basis_score', which is the mean-square displacement and 'lattice_score' which is a lattice deformation metric having units Angstr.^2. Ex: struc_score(path/to/PRIM, basis_score)"),
        m_lattice_weight(0.5),
        m_altprimclex(Structure()) {};

      BaseDatumFormatter<Configuration> *clone()const {
        return new StrucScoreConfigFormatter(*this);
      }

      bool validate(const Configuration &_config) const override;

      std::string short_header(const Configuration &_config) const override;

      std::string long_header(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    protected:
      double m_lattice_weight;
      mutable PrimClex m_altprimclex;
      fs::path m_prim_path;
      std::vector<std::string> m_prop_names;

      std::vector<double> _evaluate(const Configuration &_config) const;
    };

  }
}
#endif

