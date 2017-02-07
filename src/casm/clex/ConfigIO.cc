#include "casm/clex/ConfigIO.hh"

#include <functional>
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Norm.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigIOSelected.hh"

namespace CASM {

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

    /// \brief col_header returns: {'name(Au)', 'name(Pt)', ...}
    std::vector<std::string> MolDependent::col_header(const Configuration &_tmplt) const {
      std::vector<std::string> col;
      for(Index c = 0; c < m_mol_names.size(); c++) {
        col.push_back(name() + "(" + m_mol_names[c] + ")");
      }
      return col;
    }
  }

  namespace ConfigIO {

    // --- Comp implementations -----------

    const std::string Comp::Name = "comp";

    const std::string Comp::Desc =
      "Parametric composition parameters, individual label as argument. "
      "Without argument, all values are printed. Ex: comp(a), comp(b), etc.";

    /// \brief Returns the parametric composition
    Eigen::VectorXd Comp::evaluate(const Configuration &config) const {
      return comp(config);
    }

    /// \brief Returns true if the PrimClex has composition axes
    bool Comp::validate(const Configuration &config) const {
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

    /// \brief col_header returns: {'comp(a)', 'comp(b)', ...}
    std::vector<std::string> Comp::col_header(const Configuration &_tmplt) const {
      std::vector<std::string> col;
      for(Index c = 0; c < _index_rules().size(); c++) {
        col.push_back(name() + "(" + (char)('a' + _index_rules()[c][0]) + ")");
      }
      return col;
    }


    // --- CompN implementations -----------

    const std::string CompN::Name = "comp_n";

    const std::string CompN::Desc =
      "Number of each species per unit cell, including vacancies. "
      "No argument prints all available values. Ex: comp_n, comp_n(Au), comp_n(Pt), etc.";

    /// \brief Returns the number of each species per unit cell
    Eigen::VectorXd CompN::evaluate(const Configuration &config) const {
      return comp_n(config);
    }


    // --- SiteFrac implementations -----------

    const std::string SiteFrac::Name = "site_frac";

    const std::string SiteFrac::Desc =
      "Fraction of sites occupied by a species, including vacancies. "
      "No argument prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.";

    /// \brief Returns the site fraction
    Eigen::VectorXd SiteFrac::evaluate(const Configuration &config) const {
      return site_frac(config);
    }


    // --- AtomFrac implementations -----------

    const std::string AtomFrac::Name = "atom_frac";

    const std::string AtomFrac::Desc =
      "Fraction of atoms that are a particular species, excluding vacancies.  "
      "Without argument, all values are printed. Ex: atom_frac(Au), atom_frac(Pt), etc.";

    /// \brief Returns the site fraction
    Eigen::VectorXd AtomFrac::evaluate(const Configuration &config) const {
      return species_frac(config);
    }

    // --- MagBase implementations ---

    const std::string MagBase::Name = "relaxed_mag";

    const std::string MagBase::Desc =
      "Relaxed magnetic moment on each basis site. " ;

    /// \brief Returns true if the Configuration has relaxed_mag
    bool MagBase::validate(const Configuration &config) const {
      return config.calc_properties().contains("relaxed_mag");
    }

    /// \brief Returns the atom fraction
    Eigen::VectorXd MagBase::evaluate(const Configuration &config) const {
      return relaxed_mag(config);
    }

    // --- Corr implementations -----------

    const std::string Corr::Name = "corr";

    const std::string Corr::Desc =
      "Correlation values (evaluated basis functions, normalized per primitive cell). "
      "If no arguements, prints all correlations, using the basis set for the default "
      "cluster expansion as listed by 'casm settings -l'. "
      "If one argument, accepts either: "
      "1) a cluster expansion name, for example 'corr(formation_energy)', and "
      "evaluates all basis functions, or "
      "2) an integer index or range of indices of basis functions to evaluate, "
      "for example 'corr(6)', or 'corr(0:6)'. "
      "If two arguments, accepts cluster expansion name and an integer index or "
      "range of basis functions to evaluate, for example 'corr(formation_energy,6)' "
      "or 'corr(formation_energy,0:6)'.";

    /// \brief Returns the atom fraction
    Eigen::VectorXd Corr::evaluate(const Configuration &config) const {
      return correlations(config, m_clexulator);
    }

    /// \brief If not yet initialized, use the default clexulator from the PrimClex
    void Corr::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        const PrimClex &primclex = _tmplt.get_primclex();
        ClexDescription desc = m_clex_name.empty() ?
                               primclex.settings().default_clex() : primclex.settings().clex(m_clex_name);
        m_clexulator = primclex.clexulator(desc);
      }

      VectorXdAttribute<Configuration>::init(_tmplt);

    }

    /// \brief Expects 'corr', 'corr(clex_name)', 'corr(index_expression)', or
    /// 'corr(clex_name,index_expression)'
    bool Corr::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

      if(!splt_vec.size()) {
        return true;
      }
      else if(splt_vec.size() == 1) {
        if((splt_vec[0].find_first_not_of("0123456789") == std::string::npos) ||
           (splt_vec[0].find(':') != std::string::npos)) {
          _parse_index_expression(splt_vec[0]);
        }
        else {
          m_clex_name = splt_vec[0];
        }
      }
      else if(splt_vec.size() == 2) {
        m_clex_name = splt_vec[0];
        _parse_index_expression(splt_vec[1]);
      }
      else {
        std::stringstream ss;
        ss << "Too many arguments for 'clex'.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }
      return true;
    }

    // --- Clex implementations -----------

    const std::string Clex::Name = "clex";

    const std::string Clex::Desc =
      "Predicted property value."
      " Accepts arguments ($clex_name,$norm)."
      " ($clex_name is a cluster expansion name as listed by 'casm settings -l', default=the default clex)"
      " ($norm is the normalization, either 'per_species', or 'per_unitcell' <--default)";

    Clex::Clex() :
      ScalarAttribute<Configuration>(Name, Desc) {
      parse_args("");
    }

    Clex::Clex(const Clexulator &clexulator, const ECIContainer &eci, const Norm<Configuration> &norm) :
      ScalarAttribute<Configuration>(Name, Desc),
      m_clexulator(clexulator),
      m_eci(eci),
      m_norm(norm) {
    }

    /// \brief Returns the atom fraction
    double Clex::evaluate(const Configuration &config) const {
      return m_eci * correlations(config, m_clexulator) / _norm(config);
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<Clex> Clex::clone() const {
      return std::unique_ptr<Clex>(this->_clone());
    }

    /// \brief If not yet initialized, use the default cluster expansion from the PrimClex
    void Clex::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        const PrimClex &primclex = _tmplt.get_primclex();
        ClexDescription desc = m_clex_name.empty() ?
                               primclex.settings().default_clex() : primclex.settings().clex(m_clex_name);
        m_clexulator = primclex.clexulator(desc);
        m_eci = primclex.eci(desc);
        if(m_eci.index().back() >= m_clexulator.corr_size()) {
          Log &err_log = default_err_log();
          err_log.error<Log::standard>("bset and eci mismatch");
          err_log << "basis set size: " << m_clexulator.corr_size() << std::endl;
          err_log << "max eci index: " << m_eci.index().back() << std::endl;
          throw std::runtime_error("Error: bset and eci mismatch");
        }
      }
    }

    /// \brief Expects 'clex', 'clex(formation_energy)', or 'clex(formation_energy,per_species)'
    bool Clex::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

      m_clex_name = "";
      if(splt_vec.size()) {
        m_clex_name = splt_vec[0];
      }

      m_norm = notstd::make_cloneable<Norm<Configuration> >();
      if(splt_vec.size() == 2) {
        if(splt_vec[1] == "per_unitcell") {
          m_norm = notstd::make_cloneable<Norm<Configuration> >();
        }
        else if(splt_vec[1] == "per_species") {
          m_norm = notstd::make_cloneable<NormPerSpecies>();
        }
        else {
          std::stringstream ss;
          ss << "Error parsing second argument for 'clex'.  Received: " << args << "\n";
          throw std::runtime_error(ss.str());
        }
      }

      if(splt_vec.size() > 2) {
        std::stringstream ss;
        ss << "Too many arguments for 'clex'.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }

      return true;
    }

    /// \brief Returns the normalization
    double Clex::_norm(const Configuration &config) const {
      return (*m_norm)(config);
    }

    /// \brief Clone using copy constructor
    Clex *Clex::_clone() const {
      return new Clex(*this);
    }


    /*End ConfigIO*/
  }

  namespace ConfigIO {

    template<>
    Selected selected_in(const ConfigSelection<true> &_selection) {
      return Selected(_selection);
    }

    template<>
    Selected selected_in(const ConfigSelection<false> &_selection) {
      return Selected(_selection);
    }

    Selected selected_in() {
      return Selected();
    }


    GenericConfigFormatter<std::string> configname() {
      return GenericConfigFormatter<std::string>("configname",
                                                 "Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
      [](const Configuration & config)->std::string {
        return config.name();
      });
    }

    GenericConfigFormatter<std::string> scelname() {
      return GenericConfigFormatter<std::string>("scelname",
                                                 "Supercell name, in the form 'SCEL#_#_#_#_#_#_#'",
      [](const Configuration & config)->std::string {
        return config.get_supercell().get_name();
      });
    }



    GenericConfigFormatter<std::string> calc_status() {
      return GenericConfigFormatter<std::string>("calc_status",
                                                 "Status of calculation.",
                                                 CASM::calc_status,
                                                 CASM::has_calc_status);
    }

    GenericConfigFormatter<std::string> failure_type() {
      return GenericConfigFormatter<std::string>("failure_type",
                                                 "Reason for calculation failure.",
                                                 CASM::failure_type,
                                                 CASM::has_failure_type);
    }


    GenericConfigFormatter<Index> scel_size() {
      return GenericConfigFormatter<Index>("scel_size",
                                           "Supercell volume, given as the integer number of unit cells",
      [](const Configuration & config)->Index {
        return config.get_supercell().volume();
      });
    }

    GenericConfigFormatter<Index> multiplicity() {
      return GenericConfigFormatter<Index>("multiplicity",
                                           "Symmetric multiplicity of the configuration, excluding translational equivalents.",
      [](const Configuration & config)->Index {
        return config.get_prim().factor_group().size() / config.factor_group().size();
      });
    }

    GenericConfigFormatter<std::string> pointgroup_name() {
      return GenericConfigFormatter<std::string>("pointgroup_name",
                                                 "Name of the configuration's point group.",
      [](const Configuration & config)->std::string{
        return config.point_group().get_name();
      });
    }

    /*
    GenericConfigFormatter<bool> selected() {
      return GenericConfigFormatter<bool>("selected",
                                          "Specifies whether configuration is selected (1/true) or not (0/false)",
      [](const Configuration & config)->bool{
        return config.selected();
      });
    }
    */

    GenericConfigFormatter<double> relaxed_energy() {
      return GenericConfigFormatter<double>(
               "relaxed_energy",
               "DFT relaxed energy, normalized per primitive cell",
               CASM::relaxed_energy,
               has_relaxed_energy);
    }

    GenericConfigFormatter<double> relaxed_energy_per_species() {
      return GenericConfigFormatter<double>(
               "relaxed_energy_per_atom",
               "DFT relaxed energy, normalized per atom",
               CASM::relaxed_energy_per_species,
               has_relaxed_energy);
    }

    GenericConfigFormatter<double> reference_energy() {
      return GenericConfigFormatter<double>(
               "reference_energy",
               "reference energy, normalized per primitive cell, as determined by current reference states",
               CASM::reference_energy,
               has_reference_energy);
    }

    GenericConfigFormatter<double> reference_energy_per_species() {
      return GenericConfigFormatter<double>(
               "reference_energy_per_atom",
               "reference energy, normalized per atom, as determined by current reference states",
               CASM::reference_energy_per_species,
               has_reference_energy);
    }

    GenericConfigFormatter<double> formation_energy() {
      return GenericConfigFormatter<double>(
               "formation_energy",
               "DFT formation energy, normalized per primitive cell and measured "
               "relative to current reference states",
               CASM::formation_energy,
               has_formation_energy);
    }

    GenericConfigFormatter<double> formation_energy_per_species() {
      return GenericConfigFormatter<double>(
               "formation_energy_per_atom",
               "DFT formation energy, normalized per atom and measured relative to "
               "current reference states",
               CASM::formation_energy_per_species,
               has_formation_energy);
    }

    /*Generic1DDatumFormatter<std::vector<double>, Configuration >relaxation_strain() {
      return Generic1DDatumFormatter<std::vector<double>, Configuration >("relaxation_strain",
                                                                          "Green-Lagrange strain of dft-relaxed configuration, relative to the ideal crystal.  Ordered as [E(0,0), E(1,1), E(2,2), E(1,2), E(0,2), E(0,1)].  Accepts index as argument on interval [0,5]",
                                                                          CASM::relaxation_strain,
                                                                          has_relaxation_strain,
      [](const std::vector<double> &cont)->Index{
        return 6;
      });
      }*/

    GenericConfigFormatter<bool> is_calculated() {
      return GenericConfigFormatter<bool>("is_calculated",
                                          "True (1) if all current properties have been been calculated for the configuration",
                                          CASM::is_calculated);
    }

    GenericConfigFormatter<bool> is_primitive() {
      return GenericConfigFormatter<bool>("is_primitive",
                                          "True (1) if the configuration cannot be described within a smaller supercell",
                                          CASM::is_primitive);
    }

    GenericConfigFormatter<bool> is_canonical() {
      return GenericConfigFormatter<bool>("is_canonical",
                                          "True (1) if the configuration cannot be transfromed by symmetry to a configuration with higher lexicographic order",
                                          CASM::is_canonical);
    }

    GenericConfigFormatter<double> rms_force() {
      return GenericConfigFormatter<double>("rms_force",
                                            "Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)",
                                            CASM::rms_force,
                                            has_rms_force);
    }

    GenericConfigFormatter<double> basis_deformation() {
      return GenericConfigFormatter<double>("basis_deformation",
                                            "Cost function that describes the degree to which basis sites have relaxed",
                                            CASM::basis_deformation,
                                            has_basis_deformation);
    }

    GenericConfigFormatter<double> lattice_deformation() {
      return GenericConfigFormatter<double>("lattice_deformation",
                                            "Cost function that describes the degree to which lattice has relaxed.",
                                            CASM::lattice_deformation,
                                            has_lattice_deformation);
    }

    GenericConfigFormatter<double> volume_relaxation() {
      return GenericConfigFormatter<double>("volume_relaxation",
                                            "Change in volume due to relaxation, expressed as the ratio V/V_0.",
                                            CASM::volume_relaxation,
                                            has_volume_relaxation);
    }

    GenericConfigFormatter<double> relaxed_magmom() {
      return GenericConfigFormatter<double>("relaxed_magmom",
                                            "Relaxed magnetic moment, normalized per primative cell.",
                                            CASM::relaxed_magmom,
                                            has_relaxed_magmom);
    }

    GenericConfigFormatter<double> relaxed_magmom_per_species() {
      return GenericConfigFormatter<double>("relaxed_magmom_per_atom",
                                            "Relaxed magnetic moment, normalized per atom.",
                                            CASM::relaxed_magmom_per_species,
                                            has_relaxed_magmom);
    }



    /*End ConfigIO*/
  }

  template<>
  StringAttributeDictionary<Configuration> make_string_dictionary<Configuration>() {

    using namespace ConfigIO;
    StringAttributeDictionary<Configuration> dict;

    dict.insert(
      configname(),
      scelname(),
      calc_status(),
      failure_type(),
      pointgroup_name()
    );

    return dict;
  }

  template<>
  BooleanAttributeDictionary<Configuration> make_boolean_dictionary<Configuration>() {

    using namespace ConfigIO;
    BooleanAttributeDictionary<Configuration> dict;

    dict.insert(
      is_calculated(),
      is_canonical(),
      is_primitive(),
      //selected(),
      selected_in(),
      OnClexHull(),
      OnHull()
    );

    return dict;
  }

  template<>
  IntegerAttributeDictionary<Configuration> make_integer_dictionary<Configuration>() {

    using namespace ConfigIO;
    IntegerAttributeDictionary<Configuration> dict;

    dict.insert(
      scel_size(),
      multiplicity()
    );

    return dict;
  }

  template<>
  ScalarAttributeDictionary<Configuration> make_scalar_dictionary<Configuration>() {

    using namespace ConfigIO;
    ScalarAttributeDictionary<Configuration> dict;

    dict.insert(
      Clex(),
      HullDist(),
      ClexHullDist(),
      Novelty(),
      relaxed_energy(),
      relaxed_energy_per_species(),
      reference_energy(),
      reference_energy_per_species(),
      formation_energy(),
      formation_energy_per_species(),
      rms_force(),
      basis_deformation(),
      lattice_deformation(),
      volume_relaxation(),
      relaxed_magmom(),
      relaxed_magmom_per_species()
    );

    return dict;
  }

  template<>
  VectorXdAttributeDictionary<Configuration> make_vectorxd_dictionary<Configuration>() {

    using namespace ConfigIO;
    VectorXdAttributeDictionary<Configuration> dict;

    dict.insert(
      AtomFrac(),
      Comp(),
      CompN(),
      Corr(),
      RelaxationStrain(),
      DoFStrain(),
      SiteFrac(),
      StrucScore(),
      MagBase()
    );

    return dict;
  }

}

