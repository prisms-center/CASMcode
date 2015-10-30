#ifndef CONFIGIONOVELTY_HH
#define CONFIGIONOVELTY_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/ConfigIO.hh"
namespace CASM {

  class Configuration;

  namespace ConfigIO_impl {
    /// \brief A DatumFormatter class to measure the 'novelty' of a configuration with respect to a population of configurations
    /// Larger numbers indicate a more novel configuration, and a very large number (>~100) indicates a configuration that
    /// is linearly independent from the population (in terms of its correlations)

    /// The novelty is based on the Mahalanobis distance.  Given a population of correlations, indexed by 'i'
    /// and basis functions, indexed by 'j', that depend on the DoFs of those configurations, we define a correlation matrix
    ///
    ///    corr_mat(i,j) = <function 'j' evaluated in configuration 'i'>
    ///
    /// and an average correlation, 'avg_corr', which is taken by averaging over the rows of corr_mat
    ///
    /// The covariance of the correlations is given as
    ///
    ///    covar = corr_mat.transpose()*corr_mat/Nconfig - avg_corr.transpose()*avg_corr
    ///
    /// The Mahalanobis distance, 'M', for a particular correlation vector, 'C' (which corresponds to a particular correlation) is
    ///
    ///    M = sqrt( (C-avg_corr) * inv(covar) * (C-avg_corr).transpose() )
    ///
    /// We call inv(covar) the 'Gram matrix', which is conventional terminology for scalar products.
    ///
    /// We use a slightly different definition for the novelty measure, 'N', which is
    ///
    ///    N = sqrt( (C-avg_corr) * inv(covar+epsilon*identity) * (C-avg_corr).transpose() / Ncorr )
    ///
    /// Which is regularized by adding a small matrix (epsilon*identity, where epsilon is ~1e-5/Ncorr) and divided by Ncorr, in order
    /// to get a number that does not depend strongly on the number of basis functions in the basis set

    class NoveltyConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      NoveltyConfigFormatter() :
        BaseDatumFormatter<Configuration>("novelty", "Novelty of a configuration with respect to a population of configurations, measured using the Mahalanobis distance of its correlations. Accepts one argument, a configuration selection specifying the population against which novelty is measured (default MASTER). Ex: novelty(path/to/selection)") {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new NoveltyConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    private:
      double _evaluate(const Configuration &_config) const;

      /// specifies which selection to use as the population
      mutable std::string m_selection;

      /// Gram matrix, which defind Mahalanobis scalar product
      mutable Eigen::MatrixXd m_gram_mat;

      /// The average correlation vector of the population
      mutable Eigen::VectorXd m_avg_corr;

      /// Formatter which is used to obtain correlations
      mutable DataFormatter<Configuration> m_format;

    };
  }
}
#endif

