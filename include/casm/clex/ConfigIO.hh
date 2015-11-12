#ifndef CONFIGIO_HH
#define CONFIGIO_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"

namespace CASM {

  class Configuration;
  class ConfigIOParser : public DataFormatterParser<Configuration> {
    static int hack;
  };

  template<bool IsConst>
  class ConfigSelection;

  namespace ConfigIO_impl {

    class SelectedConfigFormatter;

    template<typename ValueType>
    using GenericConfigFormatter = GenericDatumFormatter<ValueType, Configuration>;

    void init_parser();

    /*
     */

    class CompConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      CompConfigFormatter() : BaseDatumFormatter<Configuration>("comp", "Parametric composition parameters, individual label as argument. Without argument, all values are printed. Ex: comp(a), comp(b), etc.") {}
      BaseDatumFormatter<Configuration> *clone() const {
        return new CompConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string long_header(const Configuration &_config) const override;

      bool multi_parsable()const {
        return true;
      }

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    };

    /*
     */

    class SiteFracConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      SiteFracConfigFormatter() : BaseDatumFormatter<Configuration>("site_frac", "Fraction of sites occupied by a species, including vacancies. No argument prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.") {}
      BaseDatumFormatter<Configuration> *clone() const {
        return new SiteFracConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string long_header(const Configuration &_config) const override;

      bool multi_parsable()const {
        return true;
      }

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    private:
      mutable std::vector<std::string> m_mol_names;
    };

    /*
     */

    class AtomFracConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      AtomFracConfigFormatter() : BaseDatumFormatter<Configuration>("atom_frac", "Fraction of atoms that are a particular species, excluding vacancies. Without argument, all values are printed. Ex: atom_frac(Au), atom_frac(Pt), etc.") {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new AtomFracConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string long_header(const Configuration &_config) const override;

      bool multi_parsable()const {
        return true;
      }

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    private:
      mutable std::vector<std::string> m_mol_names;
    };

    /*
     */

    class CorrConfigFormatter: public BaseDatumFormatter<Configuration> {

    public:
      CorrConfigFormatter() :
        BaseDatumFormatter<Configuration>("corr",
                                          "Average correlation values, normalized per primitive cell; accepts range as argument, for example corr(ind1:ind2)") {}

      CorrConfigFormatter(const Clexulator &clexulator) :
        BaseDatumFormatter<Configuration>("corr",
                                          "Average correlation values, normalized per primitive cell; accepts range as argument, for example corr(ind1:ind2)"),
        m_clexulator(clexulator) {}

      BaseDatumFormatter<Configuration> *clone() const {
        return new CorrConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string long_header(const Configuration &_config) const override;

      bool multi_parsable() const {
        return true;
      }

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args) {
        _parse_index_expression(args);
        return true;
      }

    private:

      mutable Clexulator m_clexulator;

    };

    /*
     */

    class ClexConfigFormatter: public BaseDatumFormatter<Configuration> {
    public:
      ClexConfigFormatter() : BaseDatumFormatter<Configuration>("clex",
                                                                  "Evaluated cluster expansion, using the current ECI. Example: clex(formation_energy)") {}
      BaseDatumFormatter<Configuration> *clone() const {
        return new ClexConfigFormatter(*this);
      }

      void init(const Configuration &_tmplt) const override;

      std::string short_header(const Configuration &_config) const override;

      void inject(const Configuration &_config, DataStream &_stream, Index) const override;

      void print(const Configuration &_config, std::ostream &_stream, Index) const override;

      jsonParser &to_json(const Configuration &_config, jsonParser &json)const override;

      bool parse_args(const std::string &args);
    private:
      mutable std::string m_clex_name;
      mutable Clexulator m_clexulator;
      mutable ECIContainer m_eci;

    };

  }
  namespace ConfigIO {
    void initialize_formatting_dictionary(DataFormatterDictionary<Configuration> &dict);


    inline
    ConfigIO_impl::AtomFracConfigFormatter atom_frac(const std::string &args = "") {
      ConfigIO_impl::AtomFracConfigFormatter tfmt;
      tfmt.parse_args(args);
      return tfmt;
    }

    inline
    ConfigIO_impl::CompConfigFormatter param_composition(const std::string &args = "") {
      ConfigIO_impl::CompConfigFormatter tfmt;
      tfmt.parse_args(args);
      return tfmt;
    }

    inline
    ConfigIO_impl::CorrConfigFormatter corr(const Clexulator &clex) {
      return ConfigIO_impl::CorrConfigFormatter(clex);
    }

    template< typename ValueType >
    ConstantValueFormatter<ValueType, Configuration> constant_value(const std::string &header, ValueType value) {
      return ConstantValueFormatter<ValueType, Configuration>(header, value);
    }

    ConfigIO_impl::GenericConfigFormatter<std::string> configname();

    ConfigIO_impl::GenericConfigFormatter<std::string> scelname();

    ConfigIO_impl::GenericConfigFormatter<Index> scel_size();

    ConfigIO_impl::GenericConfigFormatter<bool> selected();

    template<bool IsConst>
    ConfigIO_impl::SelectedConfigFormatter selected_in(const ConfigSelection<IsConst> &_selection);

    ConfigIO_impl::SelectedConfigFormatter selected_in();

    ConfigIO_impl::GenericConfigFormatter<bool> is_calculated();

    ConfigIO_impl::GenericConfigFormatter<double> formation_energy();

    Generic1DDatumFormatter<std::vector<double>, Configuration >relaxation_strain();

    ConfigIO_impl::GenericConfigFormatter<double> rms_force();

    ConfigIO_impl::GenericConfigFormatter<double> basis_deformation();

    ConfigIO_impl::GenericConfigFormatter<double> lattice_deformation();

    ConfigIO_impl::GenericConfigFormatter<double> dist_from_hull();

    ConfigIO_impl::GenericConfigFormatter<double> volume_relaxation();


  }


}
#endif
