#ifndef CONFIGIO_HH
#define CONFIGIO_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/container/ContainerTraits.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"

namespace CASM {
  
  /// \defgroup ConfigIO
  ///    
  /// \brief DatumFormatters that act on Configuration
  ///    
  /// \ingroup Configuration
  
  
  class Configuration;
  
  /*
  /// \brief Provides DataFormatterDictionary and constructs DataFormatter for Configuration
  ///
  /// \ingroup ConfigIO
  ///
  class ConfigIOParser : public DataFormatterParser<Configuration> {
    static int hack;
  };
  */
  
  typedef DataFormatterParser<Configuration> ConfigIOParser;
  

  /// 
  template<bool IsConst>
  class ConfigSelection;
  
  namespace ConfigIO_impl {
    
    /// \brief Returns fraction of sites occupied by a species
    ///
    /// Fraction of sites occupied by a species, including vacancies. No argument 
    /// prints all available values. Ex: site_frac(Au), site_frac(Pt), etc. 
    /// 
    /// \ingroup ConfigIO
    ///
    class MolDependent : public VectorXdAttribute<Configuration> {
    
    public:
      
      MolDependent(const std::string& _name, const std::string& _desc) : 
        VectorXdAttribute<Configuration>(_name, _desc) {}
      
      
      // --- Specialized implementation -----------
      
      /// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
      bool parse_args(const std::string &args) override;
      
      /// \brief Adds index rules corresponding to the parsed args
      void init(const Configuration &_tmplt) const override;
      
      /// \brief Long header returns: 'name(Au)   name(Pt)   ...'
      std::string long_header(const Configuration &_tmplt) const override;
      
    private:
      mutable std::vector<std::string> m_mol_names;
      
    };
    
/*    
    /// \brief Contains some other DatumFormatter from a specified DataFormatterDictionary
    ///
    /// May be one of:
    /// - calc, ref, form, corr, clex
    ///
    /// Examples: 
    /// - calc(relaxed_energy)
    /// - ref(relaxed_energy)
    /// - form(relaxed_energy)
    /// - corr(relaxed_energy)
    /// - clex(relaxed_energy)
    /// 
    /// \ingroup ConfigIO
    ///
    template<typename DataObject> 
    class ScalarPropertyFormatter : public BaseDatumFormatter<Configuration> {
    
    public:
    
      
      /// \brief Construct with DataFormatterDictionary of allowed Formatter
      VariableFormatter(std::string _name,
                        std::string _desc,
                        const Project& proj) : 
        BaseDatumFormatter(Name, Desc),
        m_dict(make_property_dictionary<DataObject>(proj)) {}
      
      
      /// \brief Clone using copy constructor
      virtual BaseDatumFormatter<Configuration> *clone() const {
        return new VariableFormatter(*this);
      }
      
      /// \brief The name given to this VariableFormatter, ex: 'clex' or 'ref'
      std::string var_name() const {
        return m_var_name;
      }
      
      virtual bool parse_args(const std::string &args) override {
        
        // if this already contains a formatter, try to parse args into that
        if(m_formatter.unique().get() != nullptr) {
          return m_formatter->parse_args(args);
        }
        
        std::vector<std::string> tag_names;
        std::vector<std::string> sub_exprs;
        split_formatter_expression(args, tag_names, sub_exprs);
        
        // we expect exactly 1 DatumFormatter to be specified
        if(tag_names.size() != 1) {
          throw std::runtime_error(
            std::string("Error in VariableFormatter '") + var_name() + 
            "': Expected one argument, received: '" + args + "'");
        }
        
        m_formatter = m_dict.lookup(tag_names[0]);
        
        m_formatter->parse_args(sub_exprs[0]);
        
        return true;
      }
      
      /// \brief Call init on the stored formatter
      void init(const Configuration &_config) const override {
        m_formatter->init(_config);
      }
      
      /// \brief The short header includes the variable name and stored formatter
      ///
      /// Ex: 'clex(formation_energy)' or 'ref(relaxed_energy)'
      ///
      std::string short_header(const Configuration &_config) const override {
        return var_name() + "(" + m_formatter->name() + ")";
      }
      
      /// \brief Enclose each token in the stored formatter's long_header with the variable name
      ///
      /// Ex: 'clex(relaxation_strain(0))    clex(relaxation_strain(1))   ....'
      ///
      std::string long_header(const Configuration &_config) const override {
        std::istringstream ss(m_formatter->long_header(_config));
        std::vector<std::string> tokenized {
          std::istream_iterator<std::string>(ss), 
          std::istream_iterator<std::string>()};
        
        std::string result = "";
        for(int i=0; i<tokenized.size(); ++i) {
          result += var_name() + "(" + tokenized[i] + ")";
          if(i != tokenized.size()-1) {
            result += "   ";
          }
        }
        
        return result;
      }

      void inject(const Configuration &_config, DataStream &_stream, Index pass_index) const override {
        m_formatter->inject(_config, _stream, pass_index);
      }

      void print(const Configuration &_config, std::ostream &_stream, Index pass_index) const override {
        m_formatter->print(_config, _stream, pass_index);
      }

      jsonParser &to_json(const Configuration &_config, jsonParser &json) const override {
        return m_formatter->to_json(_config, json);
      }

      
    private:
      
      mutable std::string m_var_name;
      
      notstd::cloneable_ptr<BaseDatumFormatter<Configuration> > m_formatter;
      
      DataFormatterDictionary<Configuration> m_dict;
      
    };
*/
    
  }
  
  /// Contains ConfigIO classes and functions
  namespace ConfigIO {

    class Selected;

    /// \brief Template alias for Configuration formatters of specified ValueType 
    ///
    /// \ingroup ConfigIO
    ///
    template<typename ValueType>
    using GenericConfigFormatter = GenericDatumFormatter<ValueType, Configuration>;
    
    
    template<typename ValueType>
    using ConstantValue = ConstantValueFormatter<ValueType, Configuration>;
    
    /// \brief Calculate param composition of a Configuration
    ///
    /// \ingroup ConfigIO
    class Comp : public VectorXdAttribute<Configuration> {
    
    public:
      
      static const std::string Name;
      
      static const std::string Desc;
      
      
      Comp() : 
        Base1DDatumFormatter(Name, Desc) {}
      
      
      // --- Required implementations -----------
      
      /// \brief Returns the parametric composition
      Eigen::VectorXd evaluate(const Configuration& config) const override;
      
      /// \brief Clone using copy constructor
      std::unique_ptr<Comp> clone() const {
        return std::unique_ptr<Comp>(this->_clone());
      }
      
      
      // --- Specialized implementation -----------
      
      /// \brief Returns true if the PrimClex has composition axes
      bool validate(const Configuration& config) const override;
      
      /// \brief Expects arguments of the form 'comp' or 'comp(a)', 'comp(b)', etc.
      bool parse_args(const std::string &args) override;
      
      /// \brief Long header returns: 'comp(a)   comp(b)   ...'
      std::string long_header(const Configuration &_tmplt) const override;
      
    private:
      
      /// \brief Clone using copy constructor
      Comp* _clone() const override {
        return new Comp(*this);
      }
      
    };
    
    
    /// \brief Calculate number of each species per unit cell
    ///
    /// \ingroup ConfigIO
    class CompN : public ConfigIO_impl::MolDependent {
    
    public:
      
      static const std::string Name;
      
      static const std::string Desc;
      
      
      CompN() : 
        MolDependent(Name, Desc) {}
      
      
      // --- Required implementations -----------
      
      /// \brief Returns the parametric composition
      Eigen::VectorXd evaluate(const Configuration& config) const override;
      
      /// \brief Clone using copy constructor
      std::unique_ptr<CompN> clone() const {
        return std::unique_ptr<CompN>(this->_clone());
      }
    
    private:
      
      /// \brief Clone using copy constructor
      CompN* _clone() const override {
        return new CompN(*this);
      }
    
    };
    
    
    /// \brief Returns fraction of sites occupied by a species, including vacancies
    ///
    /// Fraction of sites occupied by a species, including vacancies. No argument 
    /// prints all available values. Ex: site_frac(Au), site_frac(Pt), etc. 
    /// 
    /// \ingroup ConfigIO
    ///
    class SiteFrac : public ConfigIO_impl::MolDependent {
    
    public:
      
      static const std::string Name;
      
      static const std::string Desc;
      
      
      SiteFrac() : MolDependent(Name, Desc) {}
      
      
      // --- Required implementations -----------
      
      /// \brief Returns the site fraction
      Eigen::VectorXd evaluate(const Configuration& config) const override;
      
      /// \brief Clone using copy constructor
      std::unique_ptr<SiteFrac> clone() const {
        return std::unique_ptr<SiteFrac>(this->_clone());
      }
    
    private:
      
      /// \brief Clone using copy constructor
      SiteFrac* _clone() const override {
        return new SiteFrac(*this);
      }
      
    };
    
   
    /// \brief Returns fraction of all species that are a particular species, excluding vacancies
    ///
    /// Fraction of species that are a particular species, excluding vacancies. 
    /// Without argument, all values are printed. Ex: atom_frac(Au), atom_frac(Pt), etc. 
    /// 
    /// \ingroup ConfigIO
    ///
    class AtomFrac : public ConfigIO_impl::MolDependent {
    
    public:
      
      static const std::string Name;
      
      static const std::string Desc;
      
      
      AtomFrac() : MolDependent(Name, Desc) {}
      
      
      // --- Required implementations -----------
      
      /// \brief Returns the atom fraction
      Eigen::VectorXd evaluate(const Configuration& config) const override;
      
      /// \brief Clone using copy constructor
      std::unique_ptr<AtomFrac> clone() const {
        return std::unique_ptr<AtomFrac>(this->_clone());
      }
    
    private:
      /// \brief Clone using copy constructor
      AtomFrac* _clone() const override {
        return new AtomFrac(*this);
      }
      
    };
    
    /// \brief In the future, AtomFrac will actually be atoms only
    typedef AtomFrac SpeciesFrac;
    
    
    /// \brief Returns average correlation values, normalized per primitive cell
    ///
    /// Average correlation values, normalized per primitive cell; accepts range
    /// as argument. Ex: corr, corr(ind1:ind2)"
    /// 
    /// \ingroup ConfigIO
    ///
    class Corr : public VectorXdAttribute<Configuration> {
    
    public:
      
      static const std::string Name;
      
      static const std::string Desc;
      
      
      Corr() : VectorXdAttribute<Configuration>(Name, Desc) {}
      
      Corr(const Clexulator &clexulator) : 
        VectorXdAttribute<Configuration>(Name, Desc), 
        m_clexulator(clexulator) {}
      
      
      // --- Required implementations -----------
      
      /// \brief Returns the atom fraction
      Eigen::VectorXd evaluate(const Configuration& config) const override;
      
      /// \brief Clone using copy constructor
      std::unique_ptr<Corr> clone() const {
        return std::unique_ptr<Corr>(this->_clone());
      }
      
      
      // --- Specialized implementation -----------
      
      /// \brief If not yet initialized, use the global clexulator from the PrimClex
      void init(const Configuration &_tmplt) const override;
    
    private:
      
      /// \brief Clone using copy constructor
      Corr* _clone() const override {
        return new Corr(*this);
      }

      mutable Clexulator m_clexulator;
      
    };

  }
  
  namespace ConfigIO {
    
    /// \brief Constructs DataFormmaterDictionary containing all Configuration DatumFormatters
    //void initialize_formatting_dictionary(DataFormatterDictionary<Configuration> &dict);
    

    ConfigIO::GenericConfigFormatter<std::string> configname();

    ConfigIO::GenericConfigFormatter<std::string> scelname();

    ConfigIO::GenericConfigFormatter<Index> scel_size();

    ConfigIO::GenericConfigFormatter<bool> selected();

    template<bool IsConst>
    ConfigIO::Selected selected_in(const ConfigSelection<IsConst> &_selection);

    ConfigIO::Selected selected_in();

    ConfigIO::GenericConfigFormatter<bool> is_calculated();

    ConfigIO::GenericConfigFormatter<double> formation_energy();
    
    ConfigIO::GenericConfigFormatter<double> formation_energy_per_species();

    Generic1DDatumFormatter<std::vector<double>, Configuration > relaxation_strain();

    ConfigIO::GenericConfigFormatter<double> rms_force();

    ConfigIO::GenericConfigFormatter<double> basis_deformation();

    ConfigIO::GenericConfigFormatter<double> lattice_deformation();

    ConfigIO::GenericConfigFormatter<double> volume_relaxation();


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

}

#endif

