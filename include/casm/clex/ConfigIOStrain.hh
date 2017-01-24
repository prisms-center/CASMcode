#ifndef CONFIGIOSTRAIN_HH
#define CONFIGIOSTRAIN_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/strain/StrainConverter.hh"


namespace CASM {

  class Configuration;

  namespace ConfigIO {

    /// \brief The strain of the configuration due to relaxation, measured relative to ideal lattice vectors.
    ///
    /// Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].
    ///
    /// Accepts strain convention as argument:
    /// - 'GL' [Green-Lagrange, Default]
    /// - 'EA' [Euler-Almansi]
    /// - 'B' [Biot]
    /// - 'H' [Hencky]).
    ///
    /// Accepts index as argument on interval [0,5]
    ///
    /// Ex: 'relxation_strain', 'relaxation_strain(EA,2)'
    ///
    /// \ingroup ConfigIO
    ///
    class RelaxationStrain: public VectorXdAttribute<Configuration> {
    public:
      RelaxationStrain() :
        VectorXdAttribute<Configuration>("relaxation_strain", "The strain of the configuration due to relaxation, measured relative to ideal lattice vectors. Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)]. Accepts strain convention as argument ('GL' [Green-Lagrange, Default], 'EA' [Euler-Almansi], 'B' [Biot], or 'H' [Hencky]). Accepts index as argument on interval [0,5]"),
        m_straincalc(true) {};


      // --- Required implementations -----------

      std::unique_ptr<RelaxationStrain> clone()const {
        return std::unique_ptr<RelaxationStrain>(this->_clone());
      }

      Eigen::VectorXd evaluate(const Configuration &_config) const override;


      // --- Specialized implementation -----------

      void init(const Configuration &_tmplt) const override;

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
      mutable StrainConverter m_straincalc;
      mutable std::string m_metric_name;

    private:
      /// \brief Clone
      RelaxationStrain *_clone() const override {
        return new RelaxationStrain(*this);
      }

    };


    /// \brief The strain of the configuration due to relaxation, measured relative to ideal lattice vectors.
    ///
    /// Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].
    ///
    /// Accepts strain convention as argument:
    /// - 'GL' [Green-Lagrange, Default]
    /// - 'EA' [Euler-Almansi]
    /// - 'B' [Biot]
    /// - 'H' [Hencky]).
    ///
    /// Accepts index as argument on interval [0,5]
    ///
    /// Ex: 'relxation_strain', 'relaxation_strain(EA,2)'
    ///
    /// \ingroup ConfigIO
    ///
    class DoFStrain: public VectorXdAttribute<Configuration> {
    public:
      DoFStrain() :
        VectorXdAttribute<Configuration>("dof_strain", "The imposed strain of the configuration due to relaxation, measured relative to ideal lattice vectors. Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)]. Accepts strain convention as argument ('GL' [Green-Lagrange, Default], 'EA' [Euler-Almansi], 'B' [Biot], or 'H' [Hencky]). Accepts index as argument on interval [0,5]"),
        m_straincalc(true) {};


      // --- Required implementations -----------

      std::unique_ptr<DoFStrain> clone()const {
        return std::unique_ptr<DoFStrain>(this->_clone());
      }

      Eigen::VectorXd evaluate(const Configuration &_config) const override;


      // --- Specialized implementation -----------

      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;

      std::vector<std::string> col_header(const Configuration &_config) const override;
      /*
            void inject(const Configuration &_config, DataStream &_stream, Index) const override;

            void print(const Configuration &_config, std::ostream &_stream, Index) const override;

            jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;
      */
      bool parse_args(const std::string &args) override;

    protected:
      mutable StrainConverter m_straincalc;
      mutable std::string m_metric_name;

    private:
      /// \brief Clone
      DoFStrain *_clone() const override {
        return new DoFStrain(*this);
      }

    };

  }
}
#endif

