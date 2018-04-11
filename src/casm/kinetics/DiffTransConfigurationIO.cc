#include "casm/kinetics/DiffTransConfigurationIO.hh"

#include "casm/kinetics/DiffTransConfiguration_impl.hh"
#include "casm/casm_io/DataFormatterTools_impl.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/Selected.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include  "casm/clusterography/ClusterOrbits.hh"
namespace CASM {

  template class BaseDatumFormatter<Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<bool, std::string, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<bool, bool, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<bool, double, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<double, double, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<Index, double, Kinetics::DiffTransConfiguration>;
  template class DataFormatter<Kinetics::DiffTransConfiguration>;
  template bool DataFormatter<Kinetics::DiffTransConfiguration>::evaluate_as_scalar<bool>(Kinetics::DiffTransConfiguration const &) const;
  template double DataFormatter<Kinetics::DiffTransConfiguration>::evaluate_as_scalar<double>(Kinetics::DiffTransConfiguration const &) const;
  template class DataFormatterDictionary<Kinetics::DiffTransConfiguration>;

  namespace Kinetics {
    namespace DiffTransConfigIO {

      namespace DiffTransConfigIO_impl {

        /// \brief Expects arguments of the form 'name' or 'name(Au)', 'name(Pt)', etc.
        bool MolDependent::parse_args(const std::string &args) {
          if(args.size() > 0)
            m_mol_names.push_back(args);
          return true;
        }

        /// \brief Adds index rules corresponding to the parsed args
        void MolDependent::init(const DiffTransConfiguration &_tmplt) const {
          auto struc_molecule = _tmplt.primclex().prim().struc_molecule();

          if(m_mol_names.size() == 0) {
            for(Index i = 0; i < struc_molecule.size(); i++) {
              _add_rule(std::vector<Index>({i}));
              m_mol_names.push_back(struc_molecule[i].name());
            }
          }
          else {
            for(Index n = 0; n < m_mol_names.size(); n++) {
              Index i = 0;
              for(i = 0; i < struc_molecule.size(); i++) {
                if(struc_molecule[i].name() == m_mol_names[n]) {
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
        std::vector<std::string> MolDependent::col_header(const DiffTransConfiguration &_tmplt) const {
          std::vector<std::string> col;
          for(Index c = 0; c < m_mol_names.size(); c++) {
            col.push_back(name() + "(" + m_mol_names[c] + ")");
          }
          return col;
        }
      }



      // --- Comp implementations -----------

      const std::string Comp::Name = "comp";

      const std::string Comp::Desc =
        "Parametric composition parameters, individual label as argument. "
        "Without argument, all values are printed. Ex: comp(a), comp(b), etc.";

      /// \brief Returns the parametric composition
      Eigen::VectorXd Comp::evaluate(const DiffTransConfiguration &dtconfig) const {
        return comp(dtconfig.from_config());
      }

      /// \brief Returns true if the PrimClex has composition axes
      bool Comp::validate(const DiffTransConfiguration &dtconfig) const {
        return dtconfig.primclex().has_composition_axes();
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
      std::vector<std::string> Comp::col_header(const DiffTransConfiguration &_tmplt) const {
        std::vector<std::string> col;
        for(Index c = 0; c < _index_rules().size(); c++) {
          col.push_back(name() + "(" + (char)('a' + _index_rules()[c][0]) + ")");
        }
        return col;
      }


      // --- LocalComp implementations -----------

      const std::string LocalComp::Name = "local_comp";

      const std::string LocalComp::Desc =
        "Number of each species within radius, including vacancies. ";

      /// \brief Returns the number of each species per unit cell
      std::string LocalComp::evaluate(const DiffTransConfiguration &dtconfig) const {
        std::vector<UnitCellCoord> candidate_sites;
        neighborhood<Site>(dtconfig.diff_trans(), m_radius, all_sites_filter, std::back_inserter(candidate_sites), dtconfig.primclex().crystallography_tol());
        std::map<std::string, int> species_map;
        for(auto &site : candidate_sites) {
          Index l = dtconfig.from_config().linear_index(site);
          Index occ_idx = dtconfig.from_config().occ(l);
          std::string occ_name = dtconfig.from_config().prim().basis[ dtconfig.from_config().sublat(l) ].site_occupant()[occ_idx].name();
          if(species_map.find(occ_name) == species_map.end()) {
            species_map.emplace(occ_name, 1);
          }
          else {
            species_map[occ_name]++;
          }
        }
        std::string ret_string;
        for(auto it = species_map.begin(); it != species_map.end(); ++it) {
          if(ret_string != "") {
            ret_string = ret_string + "_";
          }
          ret_string = ret_string + it->first + "_" + std::to_string(it->second) ;
        }
        return ret_string;
      }

      bool LocalComp::parse_args(const std::string &args) {
        try {
          double rad = std::stod(args);
          if(rad >= 0) {
            m_radius = rad;
          }
          else {
            throw std::runtime_error(std::string("Format tag: 'local_comp(") + args + ") is invalid.\n");
            return false;
          }
        }
        catch(const std::invalid_argument &e) {
          throw std::runtime_error("Argument could not be converted to double. Received: " + args);
        }
        return true;
      }

      // --- CompN implementations -----------

      const std::string CompN::Name = "comp_n";

      const std::string CompN::Desc =
        "Number of each species per unit cell, including vacancies. "
        "No argument prints all available values. Ex: comp_n, comp_n(Au), comp_n(Pt), etc.";

      /// \brief Returns the number of each species per unit cell
      Eigen::VectorXd CompN::evaluate(const DiffTransConfiguration &dtconfig) const {
        return comp_n(dtconfig.from_config());
      }


      // --- SiteFrac implementations -----------

      const std::string SiteFrac::Name = "site_frac";

      const std::string SiteFrac::Desc =
        "Fraction of sites occupied by a species, including vacancies. "
        "No argument prints all available values. Ex: site_frac(Au), site_frac(Pt), etc.";

      /// \brief Returns the site fraction
      Eigen::VectorXd SiteFrac::evaluate(const DiffTransConfiguration &dtconfig) const {
        return site_frac(dtconfig.from_config());
      }


      // --- AtomFrac implementations -----------

      const std::string AtomFrac::Name = "atom_frac";

      const std::string AtomFrac::Desc =
        "Fraction of atoms that are a particular species, excluding vacancies.  "
        "Without argument, all values are printed. Ex: atom_frac(Au), atom_frac(Pt), etc.";

      /// \brief Returns the site fraction
      Eigen::VectorXd AtomFrac::evaluate(const DiffTransConfiguration &dtconfig) const {
        return species_frac(dtconfig.from_config());
      }




      // --- LocalCorr implementations -----------

      const std::string LocalCorr::Name = "local_corr";

      const std::string LocalCorr::Desc =
        "Local Correlation values (evaluated basis functions). "
        "If no arguments, prints all local correlations, using the basis set for the default "
        "cluster expansion for this diff_trans_config as listed by 'casm settings -l'. "
        "If one argument, accepts either: "
        "1) a cluster expansion name, for example 'local_corr(kra_barrier)', and "
        "evaluates all basis functions, or "
        "2) an integer index or range of indices of basis functions to evaluate, "
        "for example 'local_corr(6)', or 'local_corr(0:6)'. "
        "If two arguments, accepts cluster expansion name and an integer index or "
        "range of basis functions to evaluate, for example 'local_corr(kra_barrier,6)' "
        "or 'local_corr(kra_barrier,0:6)'.";

      /// \brief Returns the correlations
      Eigen::VectorXd LocalCorr::evaluate(const DiffTransConfiguration &dtconfig) const {
        return correlations(dtconfig, m_clexulator);
      }

      /// \brief If not yet initialized, use the default clexulator from the PrimClex
      void LocalCorr::init(const DiffTransConfiguration &_tmplt) const {
        if(!m_clexulator.initialized()) {
          const PrimClex &primclex = _tmplt.primclex();
          //Need to get default clex for a given hop based on the orbit name of the diff_trans_config
          //unless each hop has its own basis set folder similar to default
          //The fact that each hop has its own basis set folder might cause problems for selecting defaults
          //for selections in which diff_trans_configs contain different orbit names.
          //The evaluation should return a characteristic value (negative value) for mismatch hops to corr request
          ClexDescription desc = m_clex_name.empty() ?
                                 primclex.settings().default_clex() : primclex.settings().clex(m_clex_name);
          m_clexulator = primclex.clexulator(desc);
        }

        VectorXdAttribute<DiffTransConfiguration>::init(_tmplt);

      }

      /// \brief Expects 'local_corr', 'local_corr(clex_name)', 'local_corr(index_expression)', or
      /// 'local_corr(clex_name,index_expression)'
      bool LocalCorr::parse_args(const std::string &args) {
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
          ss << "Too many arguments for 'local_corr'.  Received: " << args << "\n";
          throw std::runtime_error(ss.str());
        }
        return true;
      }

      // --- LocalClex implementations -----------

      const std::string LocalClex::Name = "local_clex";

      const std::string LocalClex::Desc =
        "Predicted local property value."
        " Accepts arguments ($clex_name)."
        " ($clex_name is a cluster expansion name as listed by 'casm settings -l', default=the default clex)";

      LocalClex::LocalClex() :
        ScalarAttribute<DiffTransConfiguration>(Name, Desc) {
        parse_args("");
      }

      LocalClex::LocalClex(const Clexulator &clexulator, const ECIContainer &eci) :
        ScalarAttribute<DiffTransConfiguration>(Name, Desc),
        m_clexulator(clexulator),
        m_eci(eci) {
      }

      /// \brief Returns the atom fraction
      double LocalClex::evaluate(const DiffTransConfiguration &dtconfig) const {
        return m_eci * correlations(dtconfig, m_clexulator);
      }

      /// \brief Clone using copy constructor
      std::unique_ptr<LocalClex> LocalClex::clone() const {
        return std::unique_ptr<LocalClex>(this->_clone());
      }

      /// \brief Checks to see if clexulator and eci match orbit name of dtconfig
      bool LocalClex::validate(const DiffTransConfiguration &dtconfig) const {
        //actually check here
        return false;
      }

      /// \brief If not yet initialized, use the default cluster expansion from the PrimClex
      void LocalClex::init(const DiffTransConfiguration &_tmplt) const {
        if(!m_clexulator.initialized()) {
          const PrimClex &primclex = _tmplt.primclex();
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

      /// \brief Expects 'clex', 'clex(formation_energy)'
      bool LocalClex::parse_args(const std::string &args) {
        std::vector<std::string> splt_vec;
        boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

        m_clex_name = "";
        if(splt_vec.size()) {
          m_clex_name = splt_vec[0];
        }

        if(splt_vec.size() > 1) {
          std::stringstream ss;
          ss << "Too many arguments for 'local_clex'.  Received: " << args << "\n";
          throw std::runtime_error(ss.str());
        }

        return true;
      }

      /// \brief Clone using copy constructor
      LocalClex *LocalClex::_clone() const {
        return new LocalClex(*this);
      }


      GenericDiffTransConfigFormatter<std::string> from_configname() {
        return GenericDiffTransConfigFormatter<std::string>("from_configname",
                                                            "canonical from Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
        [](const DiffTransConfiguration & dtconfig)->std::string {
          return dtconfig.from_configname();
        });
      }

      GenericDiffTransConfigFormatter<std::string> to_configname() {
        return GenericDiffTransConfigFormatter<std::string>("to_configname",
                                                            "canonical to Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
        [](const DiffTransConfiguration & dtconfig)->std::string {
          return dtconfig.to_configname();
        });
      }

      GenericDiffTransConfigFormatter<std::string> scelname() {
        return GenericDiffTransConfigFormatter<std::string>("scelname",
                                                            "canonical Supercell name, in the form 'SCEL#_#_#_#_#_#_#'",
        [](const DiffTransConfiguration & dtconfig)->std::string {
          return dtconfig.from_config().supercell().name();
        });
      }

      GenericDiffTransConfigFormatter<std::string> orbitname() {
        return GenericDiffTransConfigFormatter<std::string>("orbitname",
                                                            "canonical orbit name, in the form 'diff_trans/#'",
        [](const DiffTransConfiguration & dtconfig)->std::string {
          return dtconfig.orbit_name();
        });
      }

      GenericDiffTransConfigFormatter<std::string> suborbit_ind() {
        return GenericDiffTransConfigFormatter<std::string>("suborbit_ind",
                                                            "suborbit index of $orbitname within $bg_configname tiled into $scelname, distinguishes unique types of hops given a decorated background",
        [](const DiffTransConfiguration & dtconfig)->std::string {
          return std::to_string(dtconfig.suborbit_ind());
        });
      }

      GenericDiffTransConfigFormatter<std::string> bg_configname() {
        return GenericDiffTransConfigFormatter<std::string>("bg_configname",
                                                            "canonical Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#', that represents the background configuration this was generated from",
        [](const DiffTransConfiguration & dtconfig)->std::string {
          return dtconfig.bg_configname();
        });
      }
      GenericDiffTransConfigFormatter<std::string> calc_status() {
        return GenericDiffTransConfigFormatter<std::string>("calc_status",
                                                            "Status of calculation.",
                                                            [](const DiffTransConfiguration & config)->std::string{return CASM::calc_status<DiffTransConfiguration>(config);},
                                                            [](const DiffTransConfiguration & config)->bool{return CASM::has_calc_status<DiffTransConfiguration>(config);});
      }

      GenericDiffTransConfigFormatter<std::string> failure_type() {
        return GenericDiffTransConfigFormatter<std::string>("failure_type",
                                                            "Reason for calculation failure.",
                                                            [](const DiffTransConfiguration & config)->std::string{return CASM::failure_type<DiffTransConfiguration>(config);},
                                                            [](const DiffTransConfiguration & config)->bool{return CASM::has_failure_type<DiffTransConfiguration>(config);});
      }


      GenericDiffTransConfigFormatter<double> min_perturb_radius() {
        return GenericDiffTransConfigFormatter<double>("min_perturb_radius",
                                                       "gives the distance from the diffusion hop of the closest perturbed site from the hop. A value of 0 indicates no perturbation",
                                                       CASM::Kinetics::min_perturb_rad);
      }
      GenericDiffTransConfigFormatter<double> max_perturb_radius() {
        return GenericDiffTransConfigFormatter<double>("max_perturb_radius",
                                                       "gives the distance from the diffusion hop of the furthest perturbed site from the hop. A value of 0 indicates no perturbation",
                                                       CASM::Kinetics::max_perturb_rad);
      }


      GenericDiffTransConfigFormatter<double> kra_barrier() {
        return GenericDiffTransConfigFormatter<double>("kra",
                                                       "kinetically resolved activation barrier, this is the distance from the average energy of the two endpoints to the highest point on the diffusion energy landscape",
                                                       CASM::Kinetics::kra,
                                                       CASM::Kinetics::has_kra);
      }
      GenericDiffTransConfigFormatter<bool> is_calculated() {
        return GenericDiffTransConfigFormatter<bool>("is_calculated",
                                                     "True (1) if all current properties have been been calculated for the diff_trans_config",
                                                     [](const DiffTransConfiguration & config)->bool{return CASM::is_calculated(config);});
      }
      GenericDiffTransConfigFormatter<bool> is_canonical() {
        return GenericDiffTransConfigFormatter<bool>("is_canonical",
                                                     "True (1) if the diff_trans_config cannot be transformed by symmetry to a diff_trans_config with higher lexicographic order",
                                                     [](const DiffTransConfiguration & config)->bool{return config.is_canonical();});
      }
    }
  }

  template<>
  StringAttributeDictionary<Kinetics::DiffTransConfiguration> make_string_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    StringAttributeDictionary<Kinetics::DiffTransConfiguration> dict;

    dict.insert(
      name<Kinetics::DiffTransConfiguration>(),
      from_configname(),
      to_configname(),
      bg_configname(),
      scelname(),
      orbitname(),
      suborbit_ind(),
      calc_status(),
      failure_type(),
      LocalComp()
    );
    return dict;
  }

  template<>
  BooleanAttributeDictionary<Kinetics::DiffTransConfiguration> make_boolean_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    BooleanAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    dict.insert(
      DB::Selected<Kinetics::DiffTransConfiguration>(),
      is_calculated(),
      is_canonical()
    );
    return dict;
  }

  template<>
  IntegerAttributeDictionary<Kinetics::DiffTransConfiguration> make_integer_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    IntegerAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  ScalarAttributeDictionary<Kinetics::DiffTransConfiguration> make_scalar_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    ScalarAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    dict.insert(
      LocalClex(),
      kra_barrier(),
      max_perturb_radius(),
      min_perturb_radius()
    );
    return dict;
  }

  template<>
  VectorXiAttributeDictionary<Kinetics::DiffTransConfiguration> make_vectorxi_dictionary<Kinetics::DiffTransConfiguration>() {
    using namespace Kinetics::DiffTransConfigIO;
    VectorXiAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  VectorXdAttributeDictionary<Kinetics::DiffTransConfiguration> make_vectorxd_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    VectorXdAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    dict.insert(
      LocalCorr(),
      AtomFrac(),
      Comp(),
      CompN(),
      SiteFrac()
    );
    return dict;
  }

}
