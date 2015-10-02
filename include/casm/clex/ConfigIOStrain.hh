#ifndef CONFIGIOSTRAIN_HH
#define CONFIGIOSTRAIN_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/strain/StrainConverter.hh"


namespace CASM {

  class Configuration;

  namespace ConfigIO_impl {
    /*
     */

    class RelaxationStrainConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      RelaxationStrainConfigFormatter() :
        BaseDatumFormatter<Configuration>("relaxation_strain", "The strain of the configuration due to relaxation, measured relative to ideal lattice vectors. Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)]. Accepts strain convention as argument ('GL' [Green-Lagrange, Default], 'EA' [Euler-Almansi], 'B' [Biot], or 'H' [Hencky]). Accepts index as argument on interval [0,5]"),
        m_straincalc(true) {};

      BaseDatumFormatter<Configuration> *clone()const {
        return new RelaxationStrainConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      bool validate(const Configuration &_config) const override;

      std::string short_header(const Configuration &_config) const override;

      std::string long_header(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    protected:
      mutable StrainConverter m_straincalc;
      mutable std::string m_metric_name;

      Eigen::VectorXd _evaluate(const Configuration &_config) const;
    };

  }
}
#endif

