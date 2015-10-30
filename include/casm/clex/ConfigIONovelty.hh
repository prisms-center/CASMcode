#ifndef CONFIGIONOVELTY_HH
#define CONFIGIONOVELTY_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/ConfigIO.hh"
namespace CASM {

  class Configuration;

  namespace ConfigIO_impl {
    /*
     */

    class NoveltyConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      NoveltyConfigFormatter() :
        BaseDatumFormatter<Configuration>("novelty", "Novelty of a configuration with respect to a set of configurations, measured using the Mahalanobis distance of its correlations. Accepts one argument, a configuration selection specifying the set used for measuring novelty. Ex: 'novelty(path/to/selection)'") {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new NoveltyConfigFormatter(*this);
      }

      /// Initialize the convex novelty and determine on-novelty configurations
      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    private:
      double _evaluate(const Configuration &_config) const;
      // specifies the dependent property to use (e.g., formation_energy or clex(formation_energy) )
      mutable std::string m_selection;
      mutable Eigen::MatrixXd m_gram_mat;
      mutable Eigen::VectorXd m_mean;
      mutable DataFormatter<Configuration> m_format;

    };
  }
}
#endif

