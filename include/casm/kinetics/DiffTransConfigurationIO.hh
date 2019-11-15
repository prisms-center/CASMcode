#ifndef CASM_DiffTransConfigurationIO
#define CASM_DiffTransConfigurationIO

#include "casm/clex/Supercell.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools.hh"
#include "casm/container/ContainerTraits.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"


namespace CASM {

  namespace Kinetics {

    namespace DiffTransConfigIO {

      namespace DiffTransConfigIO_impl {

        /// \brief Returns fraction of sites occupied by a species
        ///
        /// Fraction of sites occupied by a species, including vacancies. No argument
        /// prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.
        ///
        class MolDependent : public VectorXdAttribute<DiffTransConfiguration> {

        public:

          MolDependent(const std::string &_name, const std::string &_desc) :
            VectorXdAttribute<DiffTransConfiguration>(_name, _desc) {}


          // --- Specialized implementation -----------

          /// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
          bool parse_args(const std::string &args) override;

          /// \brief Adds index rules corresponding to the parsed args
          bool init(const DiffTransConfiguration &_tmplt) const override;

          /// \brief col_header returns: {'name(Au)', 'name(Pt)', ...}
          std::vector<std::string> col_header(const DiffTransConfiguration &_tmplt) const override;

        private:
          mutable std::vector<std::string> m_mol_names;

        };

      }

      /// \brief Calculate param composition of a DiffTransConfiguration
      ///
      class Comp : public VectorXdAttribute<DiffTransConfiguration> {

      public:

        static const std::string Name;

        static const std::string Desc;


        Comp() :
          Base1DDatumFormatter(Name, Desc) {}


        // --- Required implementations -----------

        /// \brief Returns the parametric composition
        Eigen::VectorXd evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<Comp> clone() const {
          return std::unique_ptr<Comp>(this->_clone());
        }


        // --- Specialized implementation -----------

        /// \brief Returns true if the PrimClex has composition axes
        bool validate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Expects arguments of the form 'comp' or 'comp(a)', 'comp(b)', etc.
        bool parse_args(const std::string &args) override;

        /// \brief col_header returns: {'comp(a)', 'comp(b)', ...'}
        std::vector<std::string> col_header(const DiffTransConfiguration &_tmplt) const override;

      private:

        /// \brief Clone using copy constructor
        Comp *_clone() const override {
          return new Comp(*this);
        }

      };

      /// \brief Calculate param composition of a DiffTransConfiguration
      ///
      class LocalComp : public StringAttribute<DiffTransConfiguration> {

      public:

        static const std::string Name;

        static const std::string Desc;


        LocalComp(): BaseValueFormatter(Name, Desc) {}


        // --- Required implementations -----------

        /// \brief Returns the local composition in terms off species count
        std::string evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<LocalComp> clone() const {
          return std::unique_ptr<LocalComp>(this->_clone());
        }


        // --- Specialized implementation -----------

        /// \brief Expects arguments of the form 'local_comp(radius)' etc.
        bool parse_args(const std::string &args) override;


      private:

        /// \brief Clone using copy constructor
        LocalComp *_clone() const override {
          return new LocalComp(*this);
        }

        double m_radius;

      };


      /// \brief Calculate number of each species per unit cell
      ///
      class CompN : public DiffTransConfigIO_impl::MolDependent {

      public:

        static const std::string Name;

        static const std::string Desc;


        CompN() :
          MolDependent(Name, Desc) {}


        // --- Required implementations -----------

        /// \brief Returns the parametric composition
        Eigen::VectorXd evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<CompN> clone() const {
          return std::unique_ptr<CompN>(this->_clone());
        }

      private:

        /// \brief Clone using copy constructor
        CompN *_clone() const override {
          return new CompN(*this);
        }

      };


      /// \brief Returns fraction of sites occupied by a species, including vacancies
      ///
      /// Fraction of sites occupied by a species, including vacancies. No argument
      /// prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.
      ///
      class SiteFrac : public DiffTransConfigIO_impl::MolDependent {

      public:

        static const std::string Name;

        static const std::string Desc;


        SiteFrac() : MolDependent(Name, Desc) {}


        // --- Required implementations -----------

        /// \brief Returns the site fraction
        Eigen::VectorXd evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<SiteFrac> clone() const {
          return std::unique_ptr<SiteFrac>(this->_clone());
        }

      private:

        /// \brief Clone using copy constructor
        SiteFrac *_clone() const override {
          return new SiteFrac(*this);
        }

      };


      /// \brief Returns fraction of all species that are a particular species, excluding vacancies
      ///
      /// Fraction of species that are a particular species, excluding vacancies.
      /// Without argument, all values are printed. Ex: atom_frac(Au), atom_frac(Pt), etc.
      ///
      class AtomFrac : public DiffTransConfigIO_impl::MolDependent {

      public:

        static const std::string Name;

        static const std::string Desc;


        AtomFrac() : MolDependent(Name, Desc) {}


        // --- Required implementations -----------

        /// \brief Returns the atom fraction
        Eigen::VectorXd evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<AtomFrac> clone() const {
          return std::unique_ptr<AtomFrac>(this->_clone());
        }

      private:
        /// \brief Clone using copy constructor
        AtomFrac *_clone() const override {
          return new AtomFrac(*this);
        }

      };

      /// \brief In the future, AtomFrac will actually be atoms only
      typedef AtomFrac SpeciesFrac;

      /// \brief Returns local correlation values
      ///
      /// Evaluated basis function values, normalized per primitive cell;
      ///
      class LocalCorr : public VectorXdAttribute<DiffTransConfiguration> {

      public:

        static const std::string Name;

        static const std::string Desc;


        LocalCorr() : VectorXdAttribute<DiffTransConfiguration>(Name, Desc), m_clex_name("") {}

        LocalCorr(const Clexulator &clexulator) :
          VectorXdAttribute<DiffTransConfiguration>(Name, Desc),
          m_clexulator(clexulator) {}


        // --- Required implementations -----------

        /// \brief Returns the atom fraction
        Eigen::VectorXd evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<LocalCorr> clone() const {
          return std::unique_ptr<LocalCorr>(this->_clone());
        }


        // --- Specialized implementation -----------

        /// \brief If not yet initialized, use the global clexulator from the PrimClex
        bool init(const DiffTransConfiguration &_tmplt) const override;

        /// \brief Expects 'corr', 'corr(clex_name)', 'corr(index_expression)', or
        /// 'corr(clex_name,index_expression)'
        bool parse_args(const std::string &args) override;


      private:

        /// \brief Clone using copy constructor
        LocalCorr *_clone() const override {
          return new LocalCorr(*this);
        }

        mutable Clexulator m_clexulator;
        mutable std::string m_clex_name;

      };




      /// \brief Returns predicted kra barrier
      ///
      /// Returns predicted kra energy (only kra energy for now)
      ///
      class LocalClex : public ScalarAttribute<DiffTransConfiguration> {

      public:

        static const std::string Name;

        static const std::string Desc;


        LocalClex();

        /// \brief Construct with Clexulator, ECI, and either 'formation_energy'
        LocalClex(const Clexulator &clexulator, const ECIContainer &eci);


        // --- Required implementations -----------

        /// \brief Returns the atom fraction
        double evaluate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief Clone using copy constructor
        std::unique_ptr<LocalClex> clone() const;


        // --- Specialized implementation -----------

        // validate: predicted property can't always be calculated if clexulator and
        // eci were obtained ok, checks to see if clexulator and eci match hop within dtconfig
        bool validate(const DiffTransConfiguration &dtconfig) const override;

        /// \brief If not yet initialized, use the global clexulator and eci from the PrimClex
        bool init(const DiffTransConfiguration &_tmplt) const override;

        /// \brief Expects 'local_clex', 'local_clex(kra)'
        bool parse_args(const std::string &args) override;

        /// \brief Short header returns: 'local_clex(kra)',  etc.
        std::string short_header(const DiffTransConfiguration &_tmplt) const override {
          return "local_clex(" + m_clex_name + ")";
        }

      private:

        /// \brief Clone using copy constructor
        LocalClex *_clone() const override;

        mutable std::string m_clex_name;
        mutable Clexulator m_clexulator;
        mutable ECIContainer m_eci;
      };






      template<typename ValueType>
      using GenericDiffTransConfigFormatter = GenericDatumFormatter<ValueType, DiffTransConfiguration>;

      DiffTransConfigIO::GenericDiffTransConfigFormatter<Index> multiplicity();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<Index> scel_size();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<Index> transformation_size();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<Index> suborbit_ind();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> dtconfigname();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> from_configname();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> to_configname();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> scelname();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> orbitname();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> bg_configname();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> calc_status();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> failure_type();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<std::string> subgroup_name();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<bool> is_calculated();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<bool> is_canonical();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<double> kra_barrier();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<double> min_perturb_radius();

      DiffTransConfigIO::GenericDiffTransConfigFormatter<double> max_perturb_radius();

    }

  }

  template<>
  StringAttributeDictionary<Kinetics::DiffTransConfiguration> make_string_dictionary<Kinetics::DiffTransConfiguration>();

  template<>
  BooleanAttributeDictionary<Kinetics::DiffTransConfiguration> make_boolean_dictionary<Kinetics::DiffTransConfiguration>();

  template<>
  IntegerAttributeDictionary<Kinetics::DiffTransConfiguration> make_integer_dictionary<Kinetics::DiffTransConfiguration>();

  template<>
  ScalarAttributeDictionary<Kinetics::DiffTransConfiguration> make_scalar_dictionary<Kinetics::DiffTransConfiguration>();

  template<>
  VectorXiAttributeDictionary<Kinetics::DiffTransConfiguration> make_vectorxi_dictionary<Kinetics::DiffTransConfiguration>();

  template<>
  VectorXdAttributeDictionary<Kinetics::DiffTransConfiguration> make_vectorxd_dictionary<Kinetics::DiffTransConfiguration>();

}

#endif
