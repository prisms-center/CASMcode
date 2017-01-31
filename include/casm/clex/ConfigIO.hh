#ifndef CONFIGIO_HH
#define CONFIGIO_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/container/ContainerTraits.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"

namespace CASM {

  /// \defgroup ConfigIO Configuration Queries
  ///
  /// \brief Data formatters that return Configuration properties
  ///
  /// \ingroup DataFormatter


  class Configuration;
  template<typename DataObject>
  class Norm;

  ///
  template<bool IsConst>
  class ConfigSelection;

  /**  \addtogroup ConfigIO
       @{
   */

  namespace ConfigIO_impl {

    /// \brief Returns fraction of sites occupied by a species
    ///
    /// Fraction of sites occupied by a species, including vacancies. No argument
    /// prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.
    ///
    class MolDependent : public VectorXdAttribute<Configuration> {

    public:

      MolDependent(const std::string &_name, const std::string &_desc) :
        VectorXdAttribute<Configuration>(_name, _desc) {}


      // --- Specialized implementation -----------

      /// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
      bool parse_args(const std::string &args) override;

      /// \brief Adds index rules corresponding to the parsed args
      void init(const Configuration &_tmplt) const override;

      /// \brief col_header returns: {'name(Au)', 'name(Pt)', ...}
      std::vector<std::string> col_header(const Configuration &_tmplt) const override;

    private:
      mutable std::vector<std::string> m_mol_names;

    };

  }

  /// Contains ConfigIO classes and functions
  namespace ConfigIO {

    class Selected;

    /// \brief Template alias for Configuration formatters of specified ValueType
    ///
    template<typename ValueType>
    using GenericConfigFormatter = GenericDatumFormatter<ValueType, Configuration>;


    template<typename ValueType>
    using ConstantValue = ConstantValueFormatter<ValueType, Configuration>;

    /// \brief Calculate param composition of a Configuration
    ///
    class Comp : public VectorXdAttribute<Configuration> {

    public:

      static const std::string Name;

      static const std::string Desc;


      Comp() :
        Base1DDatumFormatter(Name, Desc) {}


      // --- Required implementations -----------

      /// \brief Returns the parametric composition
      Eigen::VectorXd evaluate(const Configuration &config) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Comp> clone() const {
        return std::unique_ptr<Comp>(this->_clone());
      }


      // --- Specialized implementation -----------

      /// \brief Returns true if the PrimClex has composition axes
      bool validate(const Configuration &config) const override;

      /// \brief Expects arguments of the form 'comp' or 'comp(a)', 'comp(b)', etc.
      bool parse_args(const std::string &args) override;

      /// \brief col_header returns: {'comp(a)', 'comp(b)', ...'}
      std::vector<std::string> col_header(const Configuration &_tmplt) const override;

    private:

      /// \brief Clone using copy constructor
      Comp *_clone() const override {
        return new Comp(*this);
      }

    };


    /// \brief Calculate number of each species per unit cell
    ///
    class CompN : public ConfigIO_impl::MolDependent {

    public:

      static const std::string Name;

      static const std::string Desc;


      CompN() :
        MolDependent(Name, Desc) {}


      // --- Required implementations -----------

      /// \brief Returns the parametric composition
      Eigen::VectorXd evaluate(const Configuration &config) const override;

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
    class SiteFrac : public ConfigIO_impl::MolDependent {

    public:

      static const std::string Name;

      static const std::string Desc;


      SiteFrac() : MolDependent(Name, Desc) {}


      // --- Required implementations -----------

      /// \brief Returns the site fraction
      Eigen::VectorXd evaluate(const Configuration &config) const override;

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
    class AtomFrac : public ConfigIO_impl::MolDependent {

    public:

      static const std::string Name;

      static const std::string Desc;


      AtomFrac() : MolDependent(Name, Desc) {}


      // --- Required implementations -----------

      /// \brief Returns the atom fraction
      Eigen::VectorXd evaluate(const Configuration &config) const override;

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

    /// \brief Returns the site-specific magnetic moments
    ///
    /// Site-specific magnetic moments, should they exist
    ///

    class MagBase : public ConfigIO_impl::MolDependent {
    /* class MagBase : public VectorXdAttribute<Configuration> { */

    public:

      static const std::string Name;

      static const std::string Desc;


      MagBase() : MolDependent(Name, Desc) {}

      // --- Required implementations -----------

      /// \brief Returns the mag sites
      Eigen::VectorXd evaluate(const Configuration &config) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<MagBase> clone() const {
        return std::unique_ptr<MagBase>(this->_clone());
      }

      // --- Specialized implementation -----------

      /// \brief Returns true if the Configuration has relaxed_mag
      bool validate(const Configuration &config) const override;

    private:

      /// \brief Clone using copy constructor
      MagBase *_clone() const override {
        return new MagBase(*this);
      }

    };


    /// \brief Returns correlation values
    ///
    /// Evaluated basis function values, normalized per primitive cell;
    ///
    class Corr : public VectorXdAttribute<Configuration> {

    public:

      static const std::string Name;

      static const std::string Desc;


      Corr() : VectorXdAttribute<Configuration>(Name, Desc), m_clex_name("") {}

      Corr(const Clexulator &clexulator) :
        VectorXdAttribute<Configuration>(Name, Desc),
        m_clexulator(clexulator) {}


      // --- Required implementations -----------

      /// \brief Returns the atom fraction
      Eigen::VectorXd evaluate(const Configuration &config) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Corr> clone() const {
        return std::unique_ptr<Corr>(this->_clone());
      }


      // --- Specialized implementation -----------

      /// \brief If not yet initialized, use the global clexulator from the PrimClex
      void init(const Configuration &_tmplt) const override;

      /// \brief Expects 'corr', 'corr(clex_name)', 'corr(index_expression)', or
      /// 'corr(clex_name,index_expression)'
      bool parse_args(const std::string &args) override;


    private:

      /// \brief Clone using copy constructor
      Corr *_clone() const override {
        return new Corr(*this);
      }

      mutable Clexulator m_clexulator;
      mutable std::string m_clex_name;

    };

    /// \brief Returns predicted formation energy
    ///
    /// Returns predicted formation energy (only formation energy for now)
    ///
    class Clex : public ScalarAttribute<Configuration> {

    public:

      static const std::string Name;

      static const std::string Desc;


      Clex();

      /// \brief Construct with Clexulator, ECI, and either 'formation_energy' or 'formation_energy_per_species'
      Clex(const Clexulator &clexulator, const ECIContainer &eci, const Norm<Configuration> &norm);


      // --- Required implementations -----------

      /// \brief Returns the atom fraction
      double evaluate(const Configuration &config) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Clex> clone() const;


      // --- Specialized implementation -----------

      // validate: predicted property can always be calculated if clexulator and
      // eci were obtained ok

      /// \brief If not yet initialized, use the global clexulator and eci from the PrimClex
      void init(const Configuration &_tmplt) const override;

      /// \brief Expects 'clex', 'clex(formation_energy)', or 'clex(formation_energy_per_species)'
      bool parse_args(const std::string &args) override;

      /// \brief Short header returns: 'clex(formation_energy)', 'clex(formation_energy_per_species)', etc.
      std::string short_header(const Configuration &_tmplt) const override {
        return "clex(" + m_clex_name + ")";
      }

    private:

      /// \brief Returns the normalization
      double _norm(const Configuration &config) const;

      /// \brief Clone using copy constructor
      Clex *_clone() const override;

      mutable std::string m_clex_name;
      mutable Clexulator m_clexulator;
      mutable ECIContainer m_eci;
      mutable notstd::cloneable_ptr<Norm<Configuration> > m_norm;
    };

  }

  namespace ConfigIO {

    /// \brief Constructs DataFormmaterDictionary containing all Configuration DatumFormatters
    //void initialize_formatting_dictionary(DataFormatterDictionary<Configuration> &dict);


    ConfigIO::GenericConfigFormatter<std::string> configname();

    ConfigIO::GenericConfigFormatter<std::string> scelname();

    ConfigIO::GenericConfigFormatter<std::string> calc_status();

    ConfigIO::GenericConfigFormatter<std::string> failure_type();

    ConfigIO::GenericConfigFormatter<Index> scel_size();

    ConfigIO::GenericConfigFormatter<Index> multiplicity();

    ConfigIO::GenericConfigFormatter<std::string> subgroup_name();

    template<bool IsConst>
    ConfigIO::Selected selected_in(const ConfigSelection<IsConst> &_selection);

    ConfigIO::Selected selected_in();

    ConfigIO::GenericConfigFormatter<bool> is_calculated();

    ConfigIO::GenericConfigFormatter<bool> is_primitive();

    ConfigIO::GenericConfigFormatter<bool> is_canonical();

    ConfigIO::GenericConfigFormatter<double> relaxed_energy();

    ConfigIO::GenericConfigFormatter<double> relaxed_energy_per_species();

    ConfigIO::GenericConfigFormatter<double> reference_energy();

    ConfigIO::GenericConfigFormatter<double> reference_energy_per_species();

    ConfigIO::GenericConfigFormatter<double> formation_energy();

    ConfigIO::GenericConfigFormatter<double> formation_energy_per_species();

    Generic1DDatumFormatter<std::vector<double>, Configuration > relaxation_strain();

    ConfigIO::GenericConfigFormatter<double> rms_force();

    ConfigIO::GenericConfigFormatter<double> basis_deformation();

    ConfigIO::GenericConfigFormatter<double> lattice_deformation();

    ConfigIO::GenericConfigFormatter<double> volume_relaxation();

    ConfigIO::GenericConfigFormatter<double> relaxed_magmom();

    ConfigIO::GenericConfigFormatter<double> relaxed_magmom_per_species();

  }

  template<>
  StringAttributeDictionary<Configuration> make_string_dictionary<Configuration>();

  template<>
  BooleanAttributeDictionary<Configuration> make_boolean_dictionary<Configuration>();

  template<>
  IntegerAttributeDictionary<Configuration> make_integer_dictionary<Configuration>();

  template<>
  ScalarAttributeDictionary<Configuration> make_scalar_dictionary<Configuration>();

  template<>
  VectorXdAttributeDictionary<Configuration> make_vectorxd_dictionary<Configuration>();

  /** @} */

}

#endif

