#include <functional>
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"

namespace CASM {
  int ConfigIOParser::hack = ConfigIOParser::init(std::function<void(DataFormatterDictionary<Configuration>&) >(ConfigIO::initialize_formatting_dictionary));

  namespace ConfigIO_impl {


    //"DFT formation energy, normalized per primitive cell and measured relative to current reference states"
    double get_config_formation_energy(const Configuration &_config) {
      return _config.delta_properties().contains("relaxed_energy") ? _config.delta_properties()["relaxed_energy"].get<double>() : NAN;
    }

    bool has_config_formation_energy(const Configuration &_config) {
      return _config.delta_properties().contains("relaxed_energy");
    }

    //"Green-Lagrange strain of dft-relaxed configuration, relative to the ideal crystal.
    //Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval [0,5]"
    std::vector<double> get_config_relaxation_strain(const Configuration &_config) {
      return _config.calc_properties().contains("relaxation_strain") ? _config.calc_properties()["relaxation_strain"].get<std::vector<double> >() : std::vector<double>(6, NAN);
    }

    bool has_config_relaxation_strain(const Configuration &_config) {
      return _config.calc_properties().contains("relaxation_strain");
    }

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

    //"Change in volume due to relaxation, expressed as the ratio V/V_0."
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
      return true;
    }

    //****************************************************************************************

    void ClexConfigFormatter::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        m_clexulator = _tmplt.get_primclex().global_clexulator();
      }

      m_eci = _tmplt.get_primclex().global_eci(m_clex_name);
    };

    //****************************************************************************************

    std::string ClexConfigFormatter::short_header(const Configuration &_tmplt) const {
      return "clex(" + m_clex_name + ")";
    }

    //****************************************************************************************

    void ClexConfigFormatter::inject(const Configuration &_config, DataStream &_stream, Index) const {
      _stream << m_eci *correlations(_config, m_clexulator);
    }

    //****************************************************************************************

    void ClexConfigFormatter::print(const Configuration &_config, std::ostream &_stream, Index) const {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      _stream << m_eci *correlations(_config, m_clexulator);

    }

    //****************************************************************************************

    jsonParser &ClexConfigFormatter::to_json(const Configuration &_config, jsonParser &json)const {
      json = m_eci * correlations(_config, m_clexulator);
      return json;
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

    Generic1DDatumFormatter<std::vector<double>, Configuration >relaxation_strain() {
      return Generic1DDatumFormatter<std::vector<double>, Configuration >("relaxation_strain",
                                                                          "Green-Lagrange strain of dft-relaxed configuration, relative to the ideal crystal.  Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval [0,5]",
                                                                          ConfigIO_impl::get_config_relaxation_strain,
                                                                          ConfigIO_impl::has_config_relaxation_strain,
      [](const std::vector<double> &cont)->Index{
        return 6;
      });
    }

    ConfigIO_impl::GenericConfigFormatter<bool> is_calculated() {
      return ConfigIO_impl::GenericConfigFormatter<bool>("is_calculated",
                                                         "True (1) if all current properties have been been calculated for the configuration",
      [](const Configuration & config)->bool{
        return std::all_of(config.get_primclex().get_curr_property().begin(),
        config.get_primclex().get_curr_property().end(),
        [&](const std::string & key) {
          return config.calc_properties().contains(key);
        });
      });
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

    ConfigIO_impl::GenericConfigFormatter<double> dist_from_hull() {
      return ConfigIO_impl::GenericConfigFormatter<double>("dist_from_hull",
                                                           "Distance from DFT hull, extracted from config_list database",
                                                           ConfigIO_impl::get_config_dist_from_hull,
                                                           ConfigIO_impl::has_config_dist_from_hull);
    }

    ConfigIO_impl::GenericConfigFormatter<double> volume_relaxation() {
      return ConfigIO_impl::GenericConfigFormatter<double>("volume_relaxation",
                                                           "Change in volume due to relaxation, expressed as the ratio V/V_0.",
                                                           ConfigIO_impl::get_config_volume_relaxation,
                                                           ConfigIO_impl::has_config_volume_relaxation);
    }

    void initialize_formatting_dictionary(DataFormatterDictionary<Configuration> &dict) {
      dict // <-- add to dict
      //Self-contained formatters
      .add_formatter(ConfigIO_impl::CorrConfigFormatter())
      .add_formatter(ConfigIO_impl::ClexConfigFormatter())
      .add_formatter(ConfigIO_impl::CompConfigFormatter())
      .add_formatter(ConfigIO_impl::AtomFracConfigFormatter())
      .add_formatter(ConfigIO_impl::SiteFracConfigFormatter())

      //Generic formatters
      .add_formatter(ConfigIO::selected())
      .add_formatter(ConfigIO::configname())
      .add_formatter(ConfigIO::scelname())
      .add_formatter(ConfigIO::scel_size())
      .add_formatter(ConfigIO::formation_energy())
      .add_formatter(ConfigIO::is_calculated())
      .add_formatter(ConfigIO::relaxation_strain())
      .add_formatter(ConfigIO::rms_force())
      .add_formatter(ConfigIO::basis_deformation())
      .add_formatter(ConfigIO::lattice_deformation())
      .add_formatter(ConfigIO::volume_relaxation())

      .add_formatter(ConfigIO_impl::StrucScoreConfigFormatter())
      //hull formatters with specialized naming to disambiguate
      .add_formatter(ConfigIO_impl::OnHullConfigFormatter("on_hull",
                                                          std::string("Whether configuration is the formation_energy convex hull (i.e., is a groundstate).")
                                                          + " Accepts argument $selection (one of: <filename>, 'all', MASTER <--default)"
                                                          /* and $composition_type (one of: comp, atom_frac, site_frac).  "*/
                                                          + "  Ex: on_hull(MASTER" +/*,atom_frac*/ ").", "formation_energy"))
      .add_formatter(ConfigIO_impl::OnHullConfigFormatter("on_clex_hull",
                                                          std::string("Whether configuration is on the *cluster-expanded* formation_energy convex hull")
                                                          + " (i.e., is a *predicted* groundstate)."
                                                          + " Accepts argument $selection (one of: <filename>, 'all', MASTER <--default)"
                                                          /*and $composition_type (one of: comp, atom_frac, site_frac).  "*/
                                                          + " Ex: on_hull(MASTER).", "clex(formation_energy)"))
      .add_formatter(ConfigIO_impl::HullDistConfigFormatter("hull_dist",
                                                            std::string("Distance, in eV, of a configuration's formation_energy from the convex hull.")
                                                            + " Accepts argument $selection (one of: <filename>, 'all', MASTER <--default)"
                                                            /*and $composition_type (one of: comp, atom_frac, site_frac).  "*/
                                                            + " Ex: on_hull(MASTER).", "formation_energy"))
      .add_formatter(ConfigIO_impl::HullDistConfigFormatter("clex_hull_dist",
                                                            std::string("Distance, in eV, of a configuration's *cluster-expanded* formation_energy")
                                                            + " from the convex hull of *cluster-expanded* formation energies."
                                                            + " Accepts argument $selection (one of: <filename>, 'all', MASTER <--default)"
                                                            /*and $composition_type (one of: comp, atom_frac, site_frac).  "*/
                                                            + " Ex: on_hull(MASTER).", "clex(formation_energy)"))
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
      .add_formatter(format_operator_not<Configuration>())
      .add_formatter(format_operator_min<Configuration>())
      .add_formatter(format_operator_max<Configuration>())
      .add_formatter(format_operator_imin<Configuration>())
      .add_formatter(format_operator_imax<Configuration>())
      .add_formatter(format_operator_eq<Configuration>())
      .add_formatter(format_operator_lt<Configuration>())
      .add_formatter(format_operator_le<Configuration>())
      .add_formatter(format_operator_gt<Configuration>())
      .add_formatter(format_operator_ge<Configuration>());

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

