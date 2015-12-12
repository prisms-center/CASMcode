#include <functional>
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigIOSelected.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {
  int ConfigIOParser::hack = ConfigIOParser::init(std::function<void(DataFormatterDictionary<Configuration>&) >(ConfigIO::initialize_formatting_dictionary));

  namespace ConfigIO_impl {


    //"DFT formation energy, normalized per primitive cell and measured relative to current reference states"
    double get_config_formation_energy(const Configuration &_config) {
      return _config.delta_properties().contains("relaxed_energy") ? _config.delta_properties()["relaxed_energy"].get<double>() : NAN;
    }

    //"DFT formation energy, normalized per atom and measured relative to current reference states"
    double get_config_formation_energy_per_species(const Configuration &_config) {
      return _config.delta_properties().contains("relaxed_energy") ? formation_energy_per_species(_config) : NAN;
    }

    bool has_config_formation_energy(const Configuration &_config) {
      return _config.delta_properties().contains("relaxed_energy");
    }

    //"Green-Lagrange strain of dft-relaxed configuration, relative to the ideal crystal.
    //Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval [0,5]"
    //std::vector<double> get_config_relaxation_strain(const Configuration &_config) {
    //return _config.calc_properties().contains("relaxation_strain") ? _config.calc_properties()["relaxation_strain"].get<std::vector<double> >() : std::vector<double>(6, NAN);
    //}

    //bool has_config_relaxation_strain(const Configuration &_config) {
    //return _config.calc_properties().contains("relaxation_strain");
    //}

    //"Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)"
    double get_config_rms_force(const Configuration &_config) {
      return _config.calc_properties().contains("rms_force") ? _config.calc_properties()["rms_force"].get<double>() : NAN;
    }

    bool has_config_rms_force(const Configuration &_config) {
      return _config.calc_properties().contains("rms_force");
    }

    //"Cost function that describes the degree to which basis sites have relaxed."
    double get_config_basis_deformation(const Configuration &_config) {
      return _config.calc_properties().contains("basis_deformation") ? _config.calc_properties()["basis_deformation"].get<double>() : NAN;
    }

    bool has_config_basis_deformation(const Configuration &_config) {
      return _config.calc_properties().contains("basis_deformation");
    }

    //"Cost function that describes the degree to which lattice has relaxed."
    double get_config_lattice_deformation(const Configuration &_config) {
      return _config.calc_properties().contains("lattice_deformation") ? _config.calc_properties()["lattice_deformation"].get<double>() : NAN;
    }

    bool has_config_lattice_deformation(const Configuration &_config) {
      return _config.calc_properties().contains("lattice_deformation");
    }

    //"Change in volume due to relaxation, expressed as the ratio V/V_0."
    double get_config_volume_relaxation(const Configuration &_config) {
      return _config.calc_properties().contains("volume_relaxation") ? _config.calc_properties()["volume_relaxation"].get<double>() : NAN;
    }

    bool has_config_volume_relaxation(const Configuration &_config) {
      return _config.calc_properties().contains("volume_relaxation");
    }

    //"Calculation status."
    std::string get_config_status(const Configuration &_config) {
      return _config.status();
    }

    bool has_config_status(const Configuration &_config) {
      return !_config.status().empty();
    }

    //"Reason for calculation failure."
    std::string get_config_failure_type(const Configuration &_config) {
      return _config.failure_type();
    }

    bool has_config_failure_type(const Configuration &_config) {
      return !_config.failure_type().empty();
    }

    //"Distance from DFT hull, extracted from config_list database."
    double get_config_dist_from_hull(const Configuration &_config) {
      return _config.generated_properties().contains("dist_from_hull") ? _config.generated_properties()["dist_from_hull"].get<double>() : NAN;
    }

    bool has_config_dist_from_hull(const Configuration &_config) {
      return _config.generated_properties().contains("dist_from_hull");
    }


    double get_config_deformation_change(const Configuration &_config) {
      return std::abs(_config.deformation().determinant()) - 1.0;
    }


    //****************************************************************************************

    std::string CorrConfigFormatter::long_header(const Configuration &_tmplt) const {

      Correlation corr = correlations(_tmplt, m_clexulator);

      std::stringstream ss, word_ss;
      if(_index_rules().size() == 0) {
        for(Index lin_ind = 0; lin_ind < corr.size(); lin_ind++) {
          word_ss.str(std::string());
          word_ss.clear();
          word_ss << "corr(" << lin_ind << ")";
          ss << ' ' << std::setw(16) << word_ss.str();
        }
      }
      else {
        for(Index i = 0; i < _index_rules().size(); i++) {
          word_ss.str(std::string());
          word_ss.clear();
          word_ss << "corr(" << _index_rules()[i][0] << ')';
          ss << ' ' << std::setw(16) << word_ss.str();
        }
      }

      return ss.str();
    }

    //****************************************************************************************

    void CorrConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {

      Correlation corr = correlations(_config, m_clexulator);

      //Cases
      if(_index_rules().size() == 0) {
        for(Index nc = 0; nc < corr.size(); nc++) {
          _stream << corr[nc];
        }
      }
      else if(_index_rules()[0].size() == 1) {
        IndexContainer::const_iterator it(_index_rules().cbegin()), it_end(_index_rules().cend());
        for(; it != it_end; ++it) {
          if((*it)[0] < corr.size())
            _stream <<  corr[(*it)[0]];
          else
            _stream <<  double(NAN) << DataStream::failbit;
        }
      }

    }

    //****************************************************************************************

    void CorrConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {

      Correlation corr = correlations(_config, m_clexulator);

      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      //Cases
      if(_index_rules().size() == 0) {
        for(Index nc = 0; nc < corr.size(); nc++) {
          _stream << ' ' << std::setw(16) << corr[nc];
        }
      }
      else if(_index_rules()[0].size() == 1) {
        IndexContainer::const_iterator it(_index_rules().cbegin()), it_end(_index_rules().cend());
        for(; it != it_end; ++it) {
          if((*it)[0] < corr.size())
            _stream <<  ' ' << std::setw(16) << corr[(*it)[0]];
          else
            _stream <<  std::setw(17) << "unknown";
        }
      }

    }

    //****************************************************************************************

    jsonParser &CorrConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      json = correlations(_config, m_clexulator);
      return json;
    }

    //****************************************************************************************
    void CorrConfigFormatter::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        m_clexulator = _tmplt.get_primclex().global_clexulator();
      }
    };

    //****************************************************************************************
    bool ClexConfigFormatter::parse_args(const std::string &args) {
      if(m_clex_name.size())
        return false;
      m_clex_name = args;

      if(m_clex_name == "formation_energy") {
        m_inject = [](const Configuration & _config, DataStream & _stream, Index) {
          _stream << clex_formation_energy(_config);
        };

        m_print = [](const Configuration & _config, std::ostream & _stream, Index) {
          _stream << clex_formation_energy(_config);
        };

        m_to_json = [](const Configuration & _config, jsonParser & json)->jsonParser& {
          json = clex_formation_energy(_config);
          return json;
        };
      }
      else if(m_clex_name == "formation_energy_per_atom") {
        m_inject = [](const Configuration & _config, DataStream & _stream, Index) {
          _stream << clex_formation_energy_per_species(_config);
        };

        m_print = [](const Configuration & _config, std::ostream & _stream, Index) {
          _stream << clex_formation_energy_per_species(_config);
        };

        m_to_json = [](const Configuration & _config, jsonParser & json)->jsonParser& {
          json = clex_formation_energy_per_species(_config);
          return json;
        };
      }
      else {
        throw std::runtime_error(
          std::string("ERROR parsing '") + "clex(" + m_clex_name + ")' unknown argument." +
          " Options are 'formation_energy' or 'formation_energy_per_atom'.");
      }

      return true;
    }

    //****************************************************************************************

    void ClexConfigFormatter::init(const Configuration &_tmplt) const {

    };

    //****************************************************************************************

    std::string ClexConfigFormatter::short_header(const Configuration &_tmplt) const {
      return "clex(" + m_clex_name + ")";
    }

    //****************************************************************************************

    void ClexConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index i) const {
      m_inject(_config, _stream, i);
    }

    //****************************************************************************************

    void ClexConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index i) const {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      m_print(_config, _stream, i);

    }

    //****************************************************************************************

    jsonParser &ClexConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      return m_to_json(_config, json);
    }

    //****************************************************************************************
    bool SiteFracConfigFormatter::parse_args(const std::string &args) {
      if(args.size() > 0)
        m_mol_names.push_back(args);
      return true;
    }
    //****************************************************************************************
    void SiteFracConfigFormatter::init(const Configuration &_tmplt) const {
      Array<Molecule> struc_molecule = _tmplt.get_primclex().get_prim().get_struc_molecule();

      if(m_mol_names.size() == 0) {
        for(Index i = 0; i < struc_molecule.size(); i++) {
          _add_rule(std::vector<Index>({i}));
          m_mol_names.push_back(struc_molecule[i].name);
        }
      }
      else {
        for(Index n = 0; n < m_mol_names.size(); n++) {
          Index i = 0;
          for(i = 0; i < struc_molecule.size(); i++) {
            if(struc_molecule[i].name == m_mol_names[n]) {
              _add_rule(std::vector<Index>({i}));
              break;
            }
          }
          if(i == struc_molecule.size())
            throw std::runtime_error(std::string("Format tag: 'site_frac(") + m_mol_names[n] + ")' does not correspond to a viable composition.\n");
        }
      }
    }

    //****************************************************************************************

    std::string SiteFracConfigFormatter::long_header(const Configuration &_tmplt) const {
      std::string t_header;
      for(Index c = 0; c < m_mol_names.size(); c++) {
        t_header += name() + "(" + m_mol_names[c] + ")";
        if(c != m_mol_names.size() - 1) {
          t_header += "   ";
        }
      }
      return t_header;
    }

    //****************************************************************************************

    void SiteFracConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      auto comp = _config.get_true_composition();

      for(Index c = 0; c < _index_rules().size(); c++) {
        _stream << comp[_index_rules()[c][0]];
      }
      return;
    }

    //****************************************************************************************

    void SiteFracConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      auto comp = _config.get_true_composition();

      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);


      for(Index c = 0; c < _index_rules().size(); c++) {
        _stream << comp[_index_rules()[c][0]];
        if(c != _index_rules().size() - 1) {
          _stream << "      ";
        }
      }
      return;
    }

    //****************************************************************************************

    jsonParser &SiteFracConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      json = _config.get_true_composition();
      return json;
    }

    //****************************************************************************************

    bool AtomFracConfigFormatter::parse_args(const std::string &args) {
      if(args.size() > 0)
        m_mol_names.push_back(args);
      return true;
    }

    //****************************************************************************************

    void AtomFracConfigFormatter::init(const Configuration &_tmplt) const {
      Array<Molecule> struc_molecule = _tmplt.get_primclex().get_prim().get_struc_molecule();
      if(m_mol_names.size() == 0) {
        for(Index i = 0; i < struc_molecule.size(); i++) {
          _add_rule(std::vector<Index>({i}));
          m_mol_names.push_back(struc_molecule[i].name);
        }
      }
      else {
        for(Index n = 0; n < m_mol_names.size(); n++) {
          Index i = 0;
          for(i = 0; i < struc_molecule.size(); i++) {
            if(struc_molecule[i].name == m_mol_names[n]) {
              _add_rule(std::vector<Index>({i}));
              break;
            }
          }
          if(i == struc_molecule.size())
            throw std::runtime_error(std::string("Format tag: 'atom_frac(") + m_mol_names[n] + ")' does not correspond to a viable composition.\n");
        }
      }
    }

    //****************************************************************************************

    std::string AtomFracConfigFormatter::long_header(const Configuration &_tmplt) const {
      std::string t_header;
      for(Index c = 0; c < m_mol_names.size(); c++) {
        t_header += name() + "(" + m_mol_names[c] + ")";
        if(c != m_mol_names.size() - 1) {
          t_header += "   ";
        }
      }
      return t_header;
    }

    //****************************************************************************************

    void AtomFracConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      auto comp = _config.get_composition();

      for(Index c = 0; c < _index_rules().size(); c++) {
        _stream << comp[_index_rules()[c][0]];
      }
      return;
    }

    //****************************************************************************************

    void AtomFracConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      auto comp = _config.get_composition();

      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);


      for(Index c = 0; c < _index_rules().size(); c++) {
        _stream << comp[_index_rules()[c][0]];
        if(c != _index_rules().size() - 1) {
          _stream << "      ";
        }
      }
      return;
    }

    //****************************************************************************************

    jsonParser &AtomFracConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      json = _config.get_composition();
      return json;
    }
    //****************************************************************************************

    bool CompConfigFormatter::parse_args(const std::string &args) {
      if(args.size() == 1) {
        _add_rule(std::vector<Index>({(Index)(args[0] - 'a')}));
      }
      else if(args.size() > 1) {
        throw std::runtime_error(std::string("Format tag: 'comp(") + args + ") is invalid.\n");
        return false;
      }
      return true;
    }


    //****************************************************************************************

    void CompConfigFormatter::init(const Configuration &_tmplt) const {

      Index Nind = _tmplt.get_param_composition().size();
      if(_index_rules().size() == 0) {
        for(Index i = 0; i < Nind; i++)
          _add_rule(std::vector<Index>({i}));
      }
      else {
        for(Index i = 0; i < _index_rules().size(); i++) {
          if(_index_rules()[i][0] < 0 || _index_rules()[i][0] >= Nind)
            throw std::runtime_error(std::string("Format tag: 'comp(") + std::to_string((char)('a' + _index_rules()[i][0])) + ") is out of bounds.\n");

        }
      }
    }

    //****************************************************************************************

    std::string CompConfigFormatter::long_header(const Configuration &_tmplt) const {
      std::string t_header;
      for(Index c = 0; c < _index_rules().size(); c++) {
        t_header += name() + "(";
        t_header.push_back((char)('a' + _index_rules()[c][0]));
        t_header.push_back(')');
        if(c != _index_rules().size() - 1) {
          t_header += "        ";
        }
      }
      return t_header;
    }

    //****************************************************************************************

    void CompConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      Eigen::VectorXd p_comp = _config.get_param_composition();

      for(Index c = 0; c < _index_rules().size(); c++) {
        _stream << p_comp(_index_rules()[c][0]);
      }
      return;
    }

    //****************************************************************************************

    void CompConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      Eigen::VectorXd p_comp = _config.get_param_composition();

      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      for(Index c = 0; c < _index_rules().size(); c++) {
        _stream << p_comp(_index_rules()[c][0]);
        if(c != _index_rules().size() - 1) {
          _stream << "     ";
        }
      }
      return;
    }

    //****************************************************************************************

    jsonParser &CompConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      json = _config.get_param_composition();
      return json;
    }

    //****************************************************************************************
    /*End ConfigIO_impl*/
  }

  namespace ConfigIO {
    template<>
    ConfigIO_impl::SelectedConfigFormatter selected_in(const ConfigSelection<true> &_selection) {
      return ConfigIO_impl::SelectedConfigFormatter(_selection);
    }

    template<>
    ConfigIO_impl::SelectedConfigFormatter selected_in(const ConfigSelection<false> &_selection) {
      return ConfigIO_impl::SelectedConfigFormatter(_selection);
    }

    ConfigIO_impl::SelectedConfigFormatter selected_in() {
      return ConfigIO_impl::SelectedConfigFormatter();
    }


    ConfigIO_impl::GenericConfigFormatter<std::string> configname() {
      return ConfigIO_impl::GenericConfigFormatter<std::string>("configname",
                                                                "Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
      [](const Configuration & config)->std::string{
        return config.name();
      });
    }

    ConfigIO_impl::GenericConfigFormatter<std::string> scelname() {
      return ConfigIO_impl::GenericConfigFormatter<std::string>("scelname",
                                                                "Supercell name, in the form 'SCEL#_#_#_#_#_#_#'",
      [](const Configuration & config)->std::string{
        return config.get_supercell().get_name();
      });
    }

    ConfigIO_impl::GenericConfigFormatter<Index> scel_size() {
      return ConfigIO_impl::GenericConfigFormatter<Index>("scel_size",
                                                          "Supercell volume, given as the integer number of unit cells",
      [](const Configuration & config)->Index{
        return config.get_supercell().volume();
      });
    }


    ConfigIO_impl::GenericConfigFormatter<bool> selected() {
      return ConfigIO_impl::GenericConfigFormatter<bool>("selected",
                                                         "Specifies whether configuration is selected (1/true) or not (0/false)",
      [](const Configuration & config)->bool{
        return config.selected();
      });
    }

    ConfigIO_impl::GenericConfigFormatter<double> formation_energy() {
      return ConfigIO_impl::GenericConfigFormatter<double>("formation_energy",
                                                           "DFT formation energy, normalized per primitive cell and measured relative to current reference states",
                                                           ConfigIO_impl::get_config_formation_energy,
                                                           ConfigIO_impl::has_config_formation_energy);
    }

    ConfigIO_impl::GenericConfigFormatter<double> formation_energy_per_species() {
      return ConfigIO_impl::GenericConfigFormatter<double>("formation_energy_per_atom",
                                                           "DFT formation energy, normalized per atom and measured relative to current reference states",
                                                           ConfigIO_impl::get_config_formation_energy_per_species,
                                                           ConfigIO_impl::has_config_formation_energy);
    }

    /*Generic1DDatumFormatter<std::vector<double>, Configuration >relaxation_strain() {
      return Generic1DDatumFormatter<std::vector<double>, Configuration >("relaxation_strain",
                                                                          "Green-Lagrange strain of dft-relaxed configuration, relative to the ideal crystal.  Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval [0,5]",
                                                                          ConfigIO_impl::get_config_relaxation_strain,
                                                                          ConfigIO_impl::has_config_relaxation_strain,
      [](const std::vector<double> &cont)->Index{
        return 6;
      });
      }*/

    ConfigIO_impl::GenericConfigFormatter<bool> is_calculated() {
      return ConfigIO_impl::GenericConfigFormatter<bool>("is_calculated",
                                                         "True (1) if all current properties have been been calculated for the configuration",
                                                         CASM::is_calculated);
    }

    ConfigIO_impl::GenericConfigFormatter<double> rms_force() {
      return ConfigIO_impl::GenericConfigFormatter<double>("rms_force",
                                                           "Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)",
                                                           ConfigIO_impl::get_config_rms_force,
                                                           ConfigIO_impl::has_config_rms_force);
    }

    ConfigIO_impl::GenericConfigFormatter<double> basis_deformation() {
      return ConfigIO_impl::GenericConfigFormatter<double>("basis_deformation",
                                                           "Cost function that describes the degree to which basis sites have relaxed",
                                                           ConfigIO_impl::get_config_basis_deformation,
                                                           ConfigIO_impl::has_config_basis_deformation);
    }

    ConfigIO_impl::GenericConfigFormatter<double> lattice_deformation() {
      return ConfigIO_impl::GenericConfigFormatter<double>("lattice_deformation",
                                                           "Cost function that describes the degree to which lattice has relaxed.",
                                                           ConfigIO_impl::get_config_lattice_deformation,
                                                           ConfigIO_impl::has_config_lattice_deformation);
    }

    ConfigIO_impl::GenericConfigFormatter<double> volume_relaxation() {
      return ConfigIO_impl::GenericConfigFormatter<double>("volume_relaxation",
                                                           "Change in volume due to relaxation, expressed as the ratio V/V_0.",
                                                           ConfigIO_impl::get_config_volume_relaxation,
                                                           ConfigIO_impl::has_config_volume_relaxation);
    }

    ConfigIO_impl::GenericConfigFormatter<std::string> status() {
      return ConfigIO_impl::GenericConfigFormatter<std::string>("status",
                                                           "Status of calculation.",
                                                           ConfigIO_impl::get_config_status,
                                                           ConfigIO_impl::has_config_status);
    }

    ConfigIO_impl::GenericConfigFormatter<std::string> failure_type() {
      return ConfigIO_impl::GenericConfigFormatter<std::string>("failure_type",
                                                           "Reason for calculation failure.",
                                                           ConfigIO_impl::get_config_failure_type,
                                                           ConfigIO_impl::has_config_failure_type);
    }

    void initialize_formatting_dictionary(DataFormatterDictionary<Configuration> &dict) {
      dict // <-- add to dict

      //Self-contained formatters
      .add_formatter(ConfigIO_impl::CorrConfigFormatter())
      .add_formatter(ConfigIO_impl::ClexConfigFormatter())
      .add_formatter(ConfigIO_impl::CompConfigFormatter())
      .add_formatter(ConfigIO_impl::AtomFracConfigFormatter())
      .add_formatter(ConfigIO_impl::SiteFracConfigFormatter())
      .add_formatter(ConfigIO_impl::StrucScoreConfigFormatter())
      .add_formatter(ConfigIO_impl::OnHullConfigFormatter())
      .add_formatter(ConfigIO_impl::HullDistConfigFormatter())
      .add_formatter(ConfigIO_impl::OnClexHullConfigFormatter())
      .add_formatter(ConfigIO_impl::ClexHullDistConfigFormatter())
      .add_formatter(ConfigIO_impl::RelaxationStrainConfigFormatter())
      .add_formatter(ConfigIO_impl::NoveltyConfigFormatter())

      //Generic formatters
      .add_formatter(ConfigIO::selected())
      .add_formatter(ConfigIO::configname())
      .add_formatter(ConfigIO::scelname())
      .add_formatter(ConfigIO::scel_size())
      .add_formatter(ConfigIO::formation_energy())
      .add_formatter(ConfigIO::formation_energy_per_species())
      .add_formatter(ConfigIO::is_calculated())
      .add_formatter(ConfigIO::rms_force())
      .add_formatter(ConfigIO::basis_deformation())
      .add_formatter(ConfigIO::lattice_deformation())
      .add_formatter(ConfigIO::volume_relaxation())
      .add_formatter(ConfigIO::status())
      .add_formatter(ConfigIO::failure_type())
      .add_formatter(ConfigIO::selected_in())

      //Formatter operators
      .add_formatter(format_operator_add<Configuration>())
      .add_formatter(format_operator_sub<Configuration>())
      .add_formatter(format_operator_mult<Configuration>())
      .add_formatter(format_operator_div<Configuration>())
      .add_formatter(format_operator_exp<Configuration>())
      .add_formatter(format_operator_sq<Configuration>())
      .add_formatter(format_operator_sqrt<Configuration>())
      .add_formatter(format_operator_neg<Configuration>())
      .add_formatter(format_operator_and<Configuration>())
      .add_formatter(format_operator_or<Configuration>())
      .add_formatter(format_operator_xor<Configuration>())
      .add_formatter(format_operator_not<Configuration>())
      .add_formatter(format_operator_min<Configuration>())
      .add_formatter(format_operator_max<Configuration>())
      .add_formatter(format_operator_imin<Configuration>())
      .add_formatter(format_operator_imax<Configuration>())
      .add_formatter(format_operator_eq<Configuration>())
      .add_formatter(format_operator_lt<Configuration>())
      .add_formatter(format_operator_le<Configuration>())
      .add_formatter(format_operator_gt<Configuration>())
      .add_formatter(format_operator_ge<Configuration>())
      .add_formatter(format_operator_regex_match<Configuration>())
      .add_formatter(format_operator_regex_search<Configuration>());

    }

    /*End ConfigIO*/
  }

  //****************************************************************************************
  //****************************************************************************************
  /*
    DataFormatterDictionary<Configuration> &ConfigIOParser::dictionary() {
    static DataFormatterDictionary<Configuration>
    m_dict = std::function<void(DataFormatterDictionary<Configuration>&) >(ConfigIO_impl::initialize_config_formatter);
    return m_dict;
    }
  */
  //  ConfigIOParser::init(std::function<void(DataFormatterDictionary<Configuration>&) >(ConfigIO_impl::initialize_config_formatter));
}

