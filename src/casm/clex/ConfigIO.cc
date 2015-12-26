#include "casm/clex/ConfigIO.hh"

#include <functional>
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigIOSelected.hh"

namespace CASM {
  //int ConfigIOParser::hack = ConfigIOParser::init(std::function<void(DataFormatterDictionary<Configuration>&) >(ConfigIO::initialize_formatting_dictionary));
  
  namespace ConfigIO_impl {
    
    /// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
    bool MolDependent::parse_args(const std::string &args) {
      if(args.size() > 0)
        m_mol_names.push_back(args);
      return true;
    }
    
    /// \brief Adds index rules corresponding to the parsed args
    void MolDependent::init(const Configuration &_tmplt) const {
      auto struc_molecule = _tmplt.get_primclex().get_prim().get_struc_molecule();

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
            throw std::runtime_error(std::string("Format tag: '") + name() + "(" + 
              m_mol_names[n] + ")' does not correspond to a viable composition.\n");
        }
      }
    }
    
    /// \brief Long header returns: 'name(Au)   name(Pt)   ...'
    std::string MolDependent::long_header(const Configuration &_tmplt) const {
      std::string t_header;
      for(Index c = 0; c < m_mol_names.size(); c++) {
        t_header += name() + "(" + m_mol_names[c] + ")";
        if(c != m_mol_names.size() - 1) {
          t_header += "   ";
        }
      }
      return t_header;
    }
  }
  
  namespace ConfigIO {


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

    double get_config_deformation_change(const Configuration &_config) {
      return std::abs(_config.deformation().determinant()) - 1.0;
    }
   
    
    // --- Comp implementations -----------
    
    const std::string Comp::Name = "comp";
    
    const std::string Comp::Desc = 
      "Parametric composition parameters, individual label as argument. "
      "Without argument, all values are printed. Ex: comp(a), comp(b), etc.";
    
    /// \brief Returns the parametric composition
    Eigen::VectorXd Comp::evaluate(const Configuration& config) const {
      return comp(config);
    }
    
    /// \brief Returns true if the PrimClex has composition axes
    bool Comp::validate(const Configuration& config) const {
      return config.get_primclex().has_composition_axes();
    }
    
    /// \brief Expects arguments of the form 'comp(a)', 'comp(b)', etc.
    bool Comp::parse_args(const std::string &args) {
      if(args.size() == 1) {
        _add_rule(std::vector<Index>({(Index)(args[0] - 'a')}));
      }
      else if(args.size() > 1) {
        throw std::runtime_error(std::string("Format tag: 'comp(") + args + ") is invalid.\n");
        return false;
      }
      return true;
    }
    
    /// \brief Long header returns: 'comp(a)   comp(b)   ...'
    std::string Comp::long_header(const Configuration &_tmplt) const {
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
    
    
    // --- CompN implementations -----------
    
    const std::string CompN::Name = "comp_n";
    
    const std::string CompN::Desc = 
      "Number of each species per unit cell, including vacancies. "
      "No argument prints all available values. Ex: comp_n, comp_n(Au), comp_n(Pt), etc.";
    
    /// \brief Returns the number of each species per unit cell
    Eigen::VectorXd CompN::evaluate(const Configuration& config) const {
      return comp_n(config);
    }
    
    
    // --- SiteFrac implementations -----------
    
    const std::string SiteFrac::Name = "site_frac";
    
    const std::string SiteFrac::Desc = 
      "Fraction of sites occupied by a species, including vacancies. "
      "No argument prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.";
    
    /// \brief Returns the site fraction
    Eigen::VectorXd SiteFrac::evaluate(const Configuration& config) const {
      return site_frac(config);
    }
    
    
    // --- AtomFrac implementations -----------
    
    const std::string AtomFrac::Name = "atom_frac";
    
    const std::string AtomFrac::Desc =
      "Fraction of atoms that are a particular species, excluding vacancies.  "
      "Without argument, all values are printed. Ex: atom_frac(Au), atom_frac(Pt), etc.";
    
    /// \brief Returns the site fraction
    Eigen::VectorXd AtomFrac::evaluate(const Configuration& config) const {
      return species_frac(config);
    }
    
    
    // --- Corr implementations -----------
    
    const std::string Corr::Name ="corr";
    
    const std::string Corr::Desc =
      "Average correlation values, normalized per primitive cell; "
      "accepts range as argument, for example corr(ind1:ind2)";
    
    /// \brief Returns the atom fraction
    Eigen::VectorXd Corr::evaluate(const Configuration& config) const {
      return correlations(config, m_clexulator);
    }
      
    /// \brief If not yet initialized, use the global clexulator from the PrimClex
    void Corr::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        m_clexulator = _tmplt.get_primclex().global_clexulator();
      }
    }
        
    
    /*End ConfigIO*/
  }

  namespace ConfigIO {
    
    template<>
    ConfigIO::Selected selected_in(const ConfigSelection<true> &_selection) {
      return ConfigIO::Selected(_selection);
    }

    template<>
    ConfigIO::Selected selected_in(const ConfigSelection<false> &_selection) {
      return ConfigIO::Selected(_selection);
    }

    ConfigIO::Selected selected_in() {
      return ConfigIO::Selected();
    }


    ConfigIO::GenericConfigFormatter<std::string> configname() {
      return ConfigIO::GenericConfigFormatter<std::string>("configname",
                                                                "Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
      [](const Configuration & config)->std::string{
        return config.name();
      });
    }

    ConfigIO::GenericConfigFormatter<std::string> scelname() {
      return ConfigIO::GenericConfigFormatter<std::string>("scelname",
                                                                "Supercell name, in the form 'SCEL#_#_#_#_#_#_#'",
      [](const Configuration & config)->std::string{
        return config.get_supercell().get_name();
      });
    }

    ConfigIO::GenericConfigFormatter<Index> scel_size() {
      return ConfigIO::GenericConfigFormatter<Index>("scel_size",
                                                          "Supercell volume, given as the integer number of unit cells",
      [](const Configuration & config)->Index{
        return config.get_supercell().volume();
      });
    }


    ConfigIO::GenericConfigFormatter<bool> selected() {
      return ConfigIO::GenericConfigFormatter<bool>("selected",
                                                         "Specifies whether configuration is selected (1/true) or not (0/false)",
      [](const Configuration & config)->bool{
        return config.selected();
      });
    }

    ConfigIO::GenericConfigFormatter<double> formation_energy() {
      return ConfigIO::GenericConfigFormatter<double>(
        "formation_energy",
        "DFT formation energy, normalized per primitive cell and measured "
        "relative to current reference states",
        ConfigIO::get_config_formation_energy,
        ConfigIO::has_config_formation_energy);
    }
    
    ConfigIO::GenericConfigFormatter<double> formation_energy_per_species() {
      return ConfigIO::GenericConfigFormatter<double>(
        "formation_energy_per_atom",
        "DFT formation energy, normalized per atom and measured relative to "
        "current reference states",
        ConfigIO::get_config_formation_energy_per_species,
        ConfigIO::has_config_formation_energy);
    }
    
    /*Generic1DDatumFormatter<std::vector<double>, Configuration >relaxation_strain() {
      return Generic1DDatumFormatter<std::vector<double>, Configuration >("relaxation_strain",
                                                                          "Green-Lagrange strain of dft-relaxed configuration, relative to the ideal crystal.  Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval [0,5]",
                                                                          ConfigIO::get_config_relaxation_strain,
                                                                          ConfigIO::has_config_relaxation_strain,
      [](const std::vector<double> &cont)->Index{
        return 6;
      });
      }*/

    ConfigIO::GenericConfigFormatter<bool> is_calculated() {
      return ConfigIO::GenericConfigFormatter<bool>("is_calculated",
                                                         "True (1) if all current properties have been been calculated for the configuration",
                                                         CASM::is_calculated);
    }

    ConfigIO::GenericConfigFormatter<double> rms_force() {
      return ConfigIO::GenericConfigFormatter<double>("rms_force",
                                                           "Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)",
                                                           ConfigIO::get_config_rms_force,
                                                           ConfigIO::has_config_rms_force);
    }

    ConfigIO::GenericConfigFormatter<double> basis_deformation() {
      return ConfigIO::GenericConfigFormatter<double>("basis_deformation",
                                                           "Cost function that describes the degree to which basis sites have relaxed",
                                                           ConfigIO::get_config_basis_deformation,
                                                           ConfigIO::has_config_basis_deformation);
    }

    ConfigIO::GenericConfigFormatter<double> lattice_deformation() {
      return ConfigIO::GenericConfigFormatter<double>("lattice_deformation",
                                                           "Cost function that describes the degree to which lattice has relaxed.",
                                                           ConfigIO::get_config_lattice_deformation,
                                                           ConfigIO::has_config_lattice_deformation);
    }

    ConfigIO::GenericConfigFormatter<double> volume_relaxation() {
      return ConfigIO::GenericConfigFormatter<double>("volume_relaxation",
                                                           "Change in volume due to relaxation, expressed as the ratio V/V_0.",
                                                           ConfigIO::get_config_volume_relaxation,
                                                           ConfigIO::has_config_volume_relaxation);
    }

    /*End ConfigIO*/
  }
  
  template<>
  StringAttributeDictionary<Configuration> make_string_dictionary<Configuration>() {
    
    using namespace ConfigIO;
    StringAttributeDictionary<Configuration> dict;
    
   dict.insert(
      configname(),
      scelname()
    );
    
    return dict;
  }
  
  template<>
  BooleanAttributeDictionary<Configuration> make_boolean_dictionary<Configuration>() {
    
    using namespace ConfigIO;
    BooleanAttributeDictionary<Configuration> dict;
    
    dict.insert(
      OnHull(),
      OnClexHull(),
      selected(),
      is_calculated(),
      selected_in()
    );
    
    return dict;
  }

  template<>
  IntegerAttributeDictionary<Configuration> make_integer_dictionary<Configuration>() {
    
    using namespace ConfigIO;
    IntegerAttributeDictionary<Configuration> dict;
    
    dict.insert(
      scel_size()
    );
    
    return dict;
  }
  
  template<>
  ScalarAttributeDictionary<Configuration> make_scalar_dictionary<Configuration>() {
    
    using namespace ConfigIO;
    ScalarAttributeDictionary<Configuration> dict;
    
    dict.insert(
      //Clex(),
      HullDist(),
      ClexHullDist(),
      Novelty(),
      formation_energy(),
      formation_energy_per_species(),
      rms_force(),
      basis_deformation(),
      lattice_deformation(),
      volume_relaxation()
    );
    
    return dict;
  }
  
  template<>
  VectorXdAttributeDictionary<Configuration> make_vectorxd_dictionary<Configuration>() {
    
    using namespace ConfigIO;
    VectorXdAttributeDictionary<Configuration> dict;
    
    dict.insert(
      Corr(), 
      CompN(),
      Comp(),
      AtomFrac(),
      SiteFrac(),
      StrucScore()
    );
    
    return dict;
  }
  
}

