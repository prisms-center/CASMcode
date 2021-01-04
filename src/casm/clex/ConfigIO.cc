#include "casm/clex/ConfigIO.hh"

#include <functional>
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Norm.hh"
#include "casm/clex/Calculable.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/MappedProperties.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/Selected.hh"
#include "casm/database/ConfigDatabase.hh"

namespace CASM {

  template class BaseDatumFormatter<Configuration>;
  template class DataFormatterOperator<bool, std::string, Configuration>;
  template class DataFormatterOperator<bool, bool, Configuration>;
  template class DataFormatterOperator<bool, double, Configuration>;
  template class DataFormatterOperator<double, double, Configuration>;
  template class DataFormatterOperator<Index, double, Configuration>;
  template class DataFormatter<Configuration>;
  template bool DataFormatter<Configuration>::evaluate_as_scalar<bool>(Configuration const &) const;
  template double DataFormatter<Configuration>::evaluate_as_scalar<double>(Configuration const &) const;
  template class DataFormatterDictionary<Configuration>;

  namespace ConfigIO_impl {

    /// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
    bool MolDependent::parse_args(const std::string &args) {
      if(args.size() > 0)
        m_mol_names.push_back(args);
      return true;
    }

    /// \brief Adds index rules corresponding to the parsed args
    bool MolDependent::init(const Configuration &_tmplt) const {
      auto struc_mol = xtal::struc_molecule_name(_tmplt.primclex().prim());

      if(m_mol_names.size() == 0) {
        for(Index i = 0; i < struc_mol.size(); i++) {
          _add_rule(std::vector<Index>({i}));
          m_mol_names.push_back(struc_mol[i]);
        }
      }
      else {
        for(Index n = 0; n < m_mol_names.size(); n++) {
          Index i = 0;
          for(i = 0; i < struc_mol.size(); i++) {
            if(struc_mol[i] == m_mol_names[n]) {
              _add_rule(std::vector<Index>({i}));
              break;
            }
          }
          if(i == struc_mol.size())
            throw std::runtime_error(std::string("Format tag: '") + name() + "(" +
                                     m_mol_names[n] + ")' does not correspond to a viable composition.\n");
        }
      }
      return true;
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
      return config.primclex().has_composition_axes();
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

    // --- Corr implementations -----------

    const std::string Corr::Name = "corr";

    const std::string Corr::Desc =
      "Correlation values (evaluated basis functions, normalized per primitive cell). "
      "If no arguments, prints all correlations, using the basis set for the default "
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
    bool Corr::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        const PrimClex &primclex = _tmplt.primclex();
        ClexDescription desc = m_clex_name.empty() ?
                               primclex.settings().default_clex() : primclex.settings().clex(m_clex_name);
        m_clexulator = primclex.clexulator(desc.bset);
      }

      VectorXdAttribute<Configuration>::init(_tmplt);
      return true;
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
        ss << "Too many arguments for 'corr'.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }
      return true;
    }

    // --- GradCorr implementations -----------

    const std::string GradCorr::Name = "gradcorr";

    const std::string GradCorr::Desc =
      "Gradiant of correlation values (evaluated basis functions), with respect to a specified "
      "degree of freedom (DoF). For each configuration, output is a (D*N x C) matrix, where 'D' "
      "is DoF dimension, 'N' is either 1 (for global DoF) or number of sites in the configuration "
      "(for site DoF), and 'C' is number of basis functions. Gradient components are ordered such "
      "that components corresponding to particular site are listed in consecutive rows. Requires "
      "at least one argument, specifying the DoF with repect to which gradient is taken [e.g., 'gradcorr(disp)']. "
      "Basis functions are the basis set for the default cluster expansion, as listed by 'casm settings -l', "
      "unless otherwise specified. Accepts up to three additional arguments:\n"
      "1) a cluster expansion name, e.g. 'gradcorr(disp,formation_energy)' and/or\n"
      "2) a pair of indices, or index ranges, e.g. 'gradcorr(disp,5,2)', 'gradcorr(disp,:,3:5)', 'gradcorr(disp,0:4,3)'";

    /// \brief Returns the atom fraction
    Eigen::MatrixXd GradCorr::evaluate(const Configuration &config) const {
      return gradcorrelations(config, m_clexulator, m_key);
    }

    /// \brief If not yet initialized, use the default clexulator from the PrimClex
    bool GradCorr::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        const PrimClex &primclex = _tmplt.primclex();
        ClexDescription desc = m_clex_name.empty() ?
                               primclex.settings().default_clex() : primclex.settings().clex(m_clex_name);
        m_clexulator = primclex.clexulator(desc.bset);
      }

      MatrixXdAttribute<Configuration>::init(_tmplt);
      return true;
    }

    /// \brief Expects 'corr', 'corr(clex_name)', 'corr(index_expression)', or
    /// 'corr(clex_name,index_expression)'
    bool GradCorr::parse_args(const std::string &args) {
      //std::cout << "parsing args: " << args << "\n";
      std::vector<std::string> split_vec;
      boost::split(split_vec, args, boost::is_any_of(","), boost::token_compress_on);
      //std::cout << "after split: " << split_vec << "\n";
      if(!split_vec.size()) {
        throw std::runtime_error("'gradcorr' query requires at least one argument, corresponding to the independent variable wrt which gradient is to be computed.");
        return false;
      }
      else if(split_vec.size() > 4) {
        std::stringstream ss;
        ss << "Too many arguments for 'gradcorr'.  Received: " << args << "\n";
        throw std::runtime_error(ss.str());
      }

      boost::erase_all(split_vec[0], "'");
      //std::cout << "Now split_vec[0] is " << split_vec[0] << "\n";
      m_key = split_vec[0];

      for(Index i = 1; i < split_vec.size(); ++i) {
        if((split_vec[i].find_first_not_of("0123456789") == std::string::npos) ||
           (split_vec[i].find(':') != std::string::npos)) {
          _parse_index_expression(split_vec[i] + "," + split_vec[i + 1]);
          ++i;
        }
        else {
          m_clex_name = split_vec[i];
        }
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
      m_norm(norm.clone()) {
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
    bool Clex::init(const Configuration &_tmplt) const {
      if(!m_clexulator.initialized()) {
        const PrimClex &primclex = _tmplt.primclex();
        ClexDescription desc = m_clex_name.empty() ?
                               primclex.settings().default_clex() : primclex.settings().clex(m_clex_name);
        m_clexulator = primclex.clexulator(desc.bset);
        m_eci = primclex.eci(desc);
        if(m_eci.index().back() >= m_clexulator.corr_size()) {
          Log &err_log = CASM::err_log();
          err_log.error<Log::standard>("bset and eci mismatch");
          err_log << "basis set size: " << m_clexulator.corr_size() << std::endl;
          err_log << "max eci index: " << m_eci.index().back() << std::endl;
          throw std::runtime_error("Error: bset and eci mismatch");
        }
      }
      return true;
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
        return config.supercell().name();
      });
    }



    GenericConfigFormatter<std::string> calc_status() {
      return GenericConfigFormatter<std::string>("calc_status",
                                                 "Status of calculation.",
                                                 [](const Configuration & config)->std::string{return CASM::calc_status<Configuration>(config);},
                                                 [](const Configuration & config)->bool{return CASM::has_calc_status<Configuration>(config);});
    }

    GenericConfigFormatter<std::string> failure_type() {
      return GenericConfigFormatter<std::string>("failure_type",
                                                 "Reason for calculation failure.",
                                                 [](const Configuration & config)->std::string{return CASM::failure_type<Configuration>(config);},
                                                 [](const Configuration & config)->bool{return CASM::has_failure_type<Configuration>(config);});
    }


    GenericConfigFormatter<Index> scel_size() {
      return GenericConfigFormatter<Index>("scel_size",
                                           "Supercell volume, given as the integer number of primitive cells",
      [](const Configuration & config)->Index {
        return config.supercell().volume();
      });
    }

    GenericConfigFormatter<Index> multiplicity() {
      return GenericConfigFormatter<Index>("multiplicity",
                                           "Symmetric multiplicity of the configuration, excluding translational equivalents.",
      [](const Configuration & config)->Index {
        return config.multiplicity();
      });
    }

    GenericConfigFormatter<std::string> point_group_name() {
      return GenericConfigFormatter<std::string>("point_group_name",
                                                 "Name of the configuration's point group.",
      [](const Configuration & config)->std::string {
        return config.point_group_name();
      });
    }

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
                                          [](const Configuration & config)->bool{return CASM::is_calculated(config);});
    }

    GenericConfigFormatter<bool> is_primitive() {
      return GenericConfigFormatter<bool>("is_primitive",
                                          "True (1) if the configuration cannot be described within a smaller supercell",
                                          CASM::is_primitive);
    }

    GenericConfigFormatter<bool> is_canonical() {
      return GenericConfigFormatter<bool>("is_canonical",
                                          "True (1) if the configuration cannot be transformed by symmetry to a configuration with higher lexicographic order",
                                          CASM::is_canonical);
    }

    GenericConfigFormatter<double> rms_force() {
      return GenericConfigFormatter<double>("rms_force",
                                            "Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)",
                                            CASM::rms_force,
                                            has_rms_force);
    }

    GenericConfigFormatter<double> atomic_deformation() {
      return GenericConfigFormatter<double>("atomic_deformation",
                                            "Cost function that describes the degree to which basis sites have relaxed",
                                            CASM::atomic_deformation,
                                            has_atomic_deformation);
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

    ConfigIO::GenericConfigFormatter<jsonParser> config() {
      return GenericConfigFormatter<jsonParser>(
               "config",
               "Structure resulting from application of DoF, formatted as JSON",
      [](Configuration const & configuration) {
        jsonParser json = jsonParser::object();
        to_json(make_simple_structure(configuration), json);
        return json;
      });
    }

    ConfigIO::GenericConfigFormatter<jsonParser> dof() {
      return GenericConfigFormatter<jsonParser>(
               "dof",
               "All degrees of freedom (DoF), formatted as JSON",
      [](Configuration const & configuration) {
        return jsonParser {configuration.configdof()};
      });
    }

    ConfigIO::GenericConfigFormatter<jsonParser> properties() {
      return GenericConfigFormatter<jsonParser>(
               "properties",
               "All mapped properties, by calctype, formatted as JSON",
      [](Configuration const & configuration) {
        return jsonParser {configuration.calc_properties_map()};
      });
    }

    ConfigIO::GenericConfigFormatter<std::string> poscar() {
      return GenericConfigFormatter<std::string>(
               "poscar",
               "Structure resulting from application of DoF, formatted as VASP POSCAR",
      [](Configuration const & configuration) {
        return pos_string(configuration);
      });

    }

    /*End ConfigIO*/
  }

  template<>
  StringAttributeDictionary<Configuration> make_string_dictionary<Configuration>() {

    using namespace ConfigIO;
    StringAttributeDictionary<Configuration> dict;

    dict.insert(
      name<Configuration>(),
      configname(),
      alias<Configuration>(),
      alias_or_name<Configuration>(),
      scelname(),
      calc_status(),
      failure_type(),
      point_group_name(),
      poscar()
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
      DB::Selected<Configuration>(),
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
      atomic_deformation(),
      lattice_deformation(),
      volume_relaxation(),
      relaxed_magmom(),
      relaxed_magmom_per_species()
    );

    return dict;
  }

  template<>
  VectorXiAttributeDictionary<Configuration> make_vectorxi_dictionary<Configuration>() {
    using namespace ConfigIO;
    VectorXiAttributeDictionary<Configuration> dict;
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
      StrucScore()
    );

    return dict;
  }

  template<>
  MatrixXdAttributeDictionary<Configuration> make_matrixxd_dictionary<Configuration>() {

    using namespace ConfigIO;
    MatrixXdAttributeDictionary<Configuration> dict;

    dict.insert(
      GradCorr()
    );

    return dict;
  }

  template<>
  DataFormatterDictionary<Configuration, BaseValueFormatter<jsonParser, Configuration>> make_json_dictionary<Configuration>() {

    using namespace ConfigIO;
    DataFormatterDictionary<Configuration, BaseValueFormatter<jsonParser, Configuration>> dict;

    dict.insert(
      config(),
      dof(),
      properties()
    );

    return dict;
  }

}
