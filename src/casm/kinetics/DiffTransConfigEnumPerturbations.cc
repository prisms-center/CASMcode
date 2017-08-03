#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/symmetry/SubOrbits_impl.hh"
#include "casm/database/Selection.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"
#include "casm/app/QueryHandler.hh"


extern "C" {
  CASM::EnumInterfaceBase *make_DiffTransConfigEnumPerturbations_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffTransConfigEnumPerturbations>();
  }
}
namespace {
  class tmpBuff : public std::streambuf {
  public:
    int overflow(int c) {
      return c;
    }
  };
}

namespace CASM {

  namespace Kinetics {

    Perturbation::Perturbation(std::set<OccupationTransformation> &from_set) {
      for(const OccupationTransformation &item : from_set) {
        this->insert(item);
      }
    }

    Perturbation &Perturbation::apply_sym(const SymOp &op) {
      std::set<OccupationTransformation> tmp;
      for(const OccupationTransformation &occ_trans : *this) {
        tmp.insert(copy_apply(op, occ_trans));
      }
      this->clear();
      for(const OccupationTransformation &item : tmp) {
        this->insert(item);
      }
      return *this;
    }

    DiffTransConfigEnumPerturbations::DiffTransConfigEnumPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
      const jsonParser &local_bspecs // or iterators over IntegralClusters
    ) :
      m_background_config(background_config), m_diff_trans_orbit(diff_trans_orbit),
      m_local_bspecs(local_bspecs), m_base_config(background_config, diff_trans_orbit.prototype()) {
      //Enumerate unique diffusion transformation set for this background config
      _init_unique_difftrans();

      if(!m_unique_difftrans.size()) {
        _invalidate();
      }
      //Pick the first base DiffTransConfiguration
      _init_base_config();

      //Initialize perturbation set for m_base_dtc
      _init_perturbations();

      if(!m_perturbations.size()) {
        increment();
      }
      //Set current DiffTransConfiguration
      _set_current();
      while(!m_current->is_valid_neb() && m_unique_difftrans_it != m_unique_difftrans.end()) {
        increment();
      }
      if(m_current->is_valid_neb()) {
        this-> _initialize(&(*m_current));
        _set_step(0);
      }
    }

    const std::string DiffTransConfigEnumPerturbations::enumerator_name = "DiffTransConfigEnumPerturbations";
    const std::string DiffTransConfigEnumPerturbations::interface_help =
      "DiffTransConfigEnumPerturbations: \n\n"

      "  orbits: JSON array of strings \n"
      "    Indicate which diffusion transformation orbits are of interest. The \n"
      "    JSON array \"orbits\" should be the names of the orbits of interest.\n"
      "              \n"
      "  background_configs: JSON array of strings \n "
      "    Indicate which configurations will be the background structures for the transformations to occur in.\n"
      "    The JSON array of strings \"background_configs\" should be names of Configurations\n"
      "    that exist in your CASM project.\n\n"
      ""
      "  local_bspecs: JSON object (optional,default= {}) \n "
      "    Indicate the local bspecs hat indicate the clusters around \n"
      "    the transformation that should be perturbed. The string \"local_bspecs\" should be a local bspecs style\n"
      "    initialization used in casm bset enumeration.\n "
      "    This option takes precedence over the following option.\n\n"
      ""
      ""
      "  local_bspecs_filepath: string (optional,default=\"\") \n "
      "    Indicate the local bspecs file that indicates the clusters around \n"
      "    the transformation that should be perturbed. The string \"local_bspecs_filepath\" should be the file path\n"
      "    to the JSON file containing the local bspecs.\n "
      ""
      "  Example:\n"
      "  {\n"
      "   \"orbits\":[\"diff_trans/0\",\"diff_trans/1\"],\n"
      "   \"background_configs\":[\"SCEL8_2_2_2_0_0_0/2\"],\n"
      "    \"local_bspecs\":{\n"
      "       \"basis_functions\" : {\n"
      "        \"site_basis_functions\" : \"occupation\"\n"
      "      },\n"
      "      \"orbit_branch_specs\" : { \n"
      "       \"1\" : {\"cutoff_radius\" : 6.0},\n"
      "       \"2\" : {\"max_length\" : 6.01,\"cutoff_radius\" : 6.0},\n"
      "       \"3\" : {\"max_length\" : 4.01,\"cutoff_radius\" : 5.0}\n"
      "      }\n"
      "    }\n"
      "  }\n\n"
      ;

    void DiffTransConfigEnumPerturbations::increment() {
      bool is_valid_config {false};
      while(!is_valid_config) {
        if(m_perturb_it == m_perturbations.end()) {
          _increment_base_config();
          if(m_unique_difftrans_it == m_unique_difftrans.end()) {
            return;
          }
          _init_perturbations();
          m_perturb_it = m_perturbations.begin();
        }

        _set_current();
        m_perturb_it++;
        is_valid_config = _check_current();
      }
    };

    int DiffTransConfigEnumPerturbations::run(const PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt) {

      jsonParser orbitnames;
      if(!_kwargs.get_if(orbitnames, "orbits")) {
        default_err_log() << "DiffTransConfigEnumPerturbations currently has no default and requires a correct JSON with a orbits tag within it" << std::endl;
        default_err_log() << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser confignames;
      if(!_kwargs.get_if(confignames, "background_configs")) {
        default_err_log() << "DiffTransConfigEnumPerturbations currently has no default and requires a correct JSON with a background_configs tag within it" << std::endl;
        default_err_log() << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser kwargs;
      jsonParser local_bspecs;
      if(_kwargs.get_if(kwargs, "local_bspecs")) {
        local_bspecs = kwargs;
      }
      else if(_kwargs.get_if(kwargs, "local_bspecs_filepath")) {
        //look for filepath given and load it
        fs::path local_bspecs_path = _kwargs.get<std::string>();
        jsonParser tmp {local_bspecs_path};
        local_bspecs = tmp;
      }
      else {
        //look for default local bspecs file.
        fs::path local_bspecs_path = "default/path/to/local_bspecs.json";
        jsonParser tmp {local_bspecs_path};
        local_bspecs = tmp;
      }


      Log &log = primclex.log();
      auto &db_diff_trans_configs = primclex.db<DiffTransConfiguration>();

      Index Ninit = db_diff_trans_configs.size();
      log << "# DiffTransConfiguration in this project: " << Ninit << "\n" << std::endl;

      log.begin(enumerator_name);
      for(const auto &configname_parser : confignames) {
        std::string configname = configname_parser.get<std::string>();
        log << "Searching in " << configname << "..." << std::endl;

        Configuration bg_config = *primclex.db<Configuration>().find(configname);
        for(const auto &orbitname_parser : orbitnames) {
          std::string orbitname = orbitname_parser.get<std::string>();
          log << "\tUsing " << orbitname << "..." << std::flush;
          Index Ninit_spec = db_diff_trans_configs.size() ;

          PrimPeriodicDiffTransOrbit dtorbit = *primclex.db<PrimPeriodicDiffTransOrbit>().find(orbitname);
          /// check if configuration is big enough for local bspecs here
          /// give warning if not

          tmpBuff streambuff;
          std::ostream dead(&streambuff);
          std::vector<LocalIntegralClusterOrbit> local_orbits;
          make_local_orbits(
            dtorbit.prototype(),
            local_bspecs,
            alloy_sites_filter,
            primclex.crystallography_tol(),
            std::back_inserter(local_orbits),
            dead);
          if(has_local_neighborhood_overlap(local_orbits, bg_config.supercell())) {
            log << "WARNING!!! CHOICE OF BACKGROUND CONFIGURATION " << configname <<
                "\nRESULTS IN AN OVERLAP IN THE LOCAL CLUSTERS OF " << orbitname <<
                "\nWITH ITS PERIODIC IMAGES. CONSIDER CHOOSING \n" <<
                "A LARGER BACKGROUND CONFIGURATION." << std::endl;
          }
          DiffTransConfigEnumPerturbations enumerator(bg_config, dtorbit, local_bspecs);
          auto enum_it = enumerator.begin();
          std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);
          if(!filter_expr.empty()) {
            try {
              DataFormatter<DiffTransConfiguration> filter = primclex.settings().query_handler<DiffTransConfiguration>().dict().parse(filter_expr);
              while(enum_it != enumerator.end()) {
                ValueDataStream<bool> _stream;
                _stream << filter(*enum_it);
                if(_stream.value()) {
                  db_diff_trans_configs.insert(enum_it->canonical_form());
                }
                ++enum_it;
              }
            }
            catch(std::exception &e) {
              primclex.err_log() << "Cannot filter difftransconfigs using the expression provided: \n" << e.what() << "\nExiting...\n";
              return ERR_INVALID_ARG;
            }
          }
          else {
            primclex.db<DiffTransConfiguration>().insert(enumerator.begin(), enumerator.end());
          }

          Index Nfinal_spec = db_diff_trans_configs.size();
          log << "Found " << Nfinal_spec - Ninit_spec << " new difftransconfigs" << std::endl;
        }
      }

      log << "  DONE." << std::endl << std::endl;

      Index Nfinal = db_diff_trans_configs.size();

      log << "# new DiffTransConfiguration: " << Nfinal - Ninit << "\n";
      log << "# DiffTransConfiguration in this project: " << Nfinal << "\n" << std::endl;

      log << "Writing DiffTransConfiguration database..." << std::endl;
      db_diff_trans_configs.commit();
      log << "  DONE" << std::endl;
      return 0;
    }

    ///------------------------------------Internal functions-------------------------------///

    ///sets the first m_base_config as the orbit prototype placed in m_background_config
    void DiffTransConfigEnumPerturbations::_init_base_config() {
      std::set<Index> unique_indeces;
      ScelPeriodicDiffTransSymCompare symcompare(m_background_config.supercell().prim_grid(),
                                                 m_background_config.supercell().crystallography_tol());
      DiffusionTransformation prepped = symcompare.prepare(*(m_unique_difftrans.begin()));
      Configuration bg_config = make_attachable(prepped, m_background_config);
      DiffTransConfiguration tmp(bg_config, prepped);
      m_base_config = tmp;//.canonical_form();
      m_current = notstd::make_cloneable<DiffTransConfiguration>(m_base_config);
      return;
    }

    /// Uses make_suborbit_generators to initialize set of the unique diffusion transformations in
    /// this configuration
    void DiffTransConfigEnumPerturbations::_init_unique_difftrans() {
      std::vector<PrimPeriodicDiffTransOrbit> orbit_vec;
      orbit_vec.push_back(m_diff_trans_orbit);
      std::vector<DiffusionTransformation> subprototypes;
      make_suborbit_generators(orbit_vec.begin(), orbit_vec.end(), m_background_config, std::back_inserter(subprototypes));
      m_unique_difftrans.insert(subprototypes.begin(), subprototypes.end());
      m_unique_difftrans_it = m_unique_difftrans.begin();
      return;
    }

    /// Uses m_local_bspecs to initialize set of the possible perturbations in
    /// this configuration
    void DiffTransConfigEnumPerturbations::_init_perturbations() {
      m_perturbations.clear();
      const SymGroup &scel_grp = m_base_config.from_config().supercell().factor_group();
      Kinetics::ScelPeriodicDiffTransSymCompare dt_sym_compare(m_base_config.from_config().supercell().prim_grid(),
                                                               m_base_config.from_config().primclex().crystallography_tol());
      SymGroup generating_group = make_invariant_subgroup(m_base_config.diff_trans(), scel_grp, dt_sym_compare);
      std::vector<LocalIntegralClusterOrbit> local_orbits;
      tmpBuff streambuff;
      std::ostream dead(&streambuff);
      make_local_orbits(
        m_base_config.diff_trans(),
        m_local_bspecs,
        alloy_sites_filter,
        m_base_config.from_config().primclex().crystallography_tol(),
        std::back_inserter(local_orbits),
        dead,
        generating_group);

      for(const auto &orbit : local_orbits) {
        Eigen::VectorXi max_count(orbit.prototype().size());
        for(int i = 0; i < max_count.size(); ++i) {
          max_count(i) = m_base_config.from_config().supercell().max_allowed_occupation()
                         [m_base_config.from_config().supercell().linear_index(orbit.prototype()[i])];
        }
        EigenCounter<Eigen::VectorXi> occ_counter(Eigen::VectorXi::Zero(orbit.prototype().size()),
                                                  max_count, Eigen::VectorXi::Constant(orbit.prototype().size(), 1));
        do {
          std::set<OccupationTransformation> my_occ_transforms;
          for(int i = 0; i < orbit.prototype().size(); ++i) {
            my_occ_transforms.emplace(orbit.prototype()[i],
                                      m_base_config.from_config().occ(m_base_config.from_config().supercell().linear_index(orbit.prototype()[i])),
                                      occ_counter()[i]);
          }
          Perturbation proto_perturb(my_occ_transforms);
          Perturbation max_perturbation = proto_perturb;
          for(auto &op : generating_group) {
            //m_perturbations.insert(copy_apply(op, proto_perturb));
            if(copy_apply(op, proto_perturb) > max_perturbation) {
              max_perturbation = copy_apply(op, proto_perturb);
            }
          }
          m_perturbations.insert(max_perturbation);
        }
        while(++occ_counter);

      }
      m_perturb_it = m_perturbations.begin();
      return;
    }

    /// Applies current perturbation to m_base_config and stores result in m_current
    void DiffTransConfigEnumPerturbations::_set_current() {
      /// apply_perturbation(perturb,m_base_dtc);
      Configuration tmp {m_base_config.from_config()};
      //Need to make sure perturbations doesn't attempt to alter any sites of the hop
      std::set<Index> unique_indeces;
      for(auto &traj : m_base_config.diff_trans().specie_traj()) {
        unique_indeces.insert(tmp.supercell().linear_index(traj.from.uccoord));
      }
      for(const auto &occ_trans : *m_perturb_it) {
        if(unique_indeces.find(tmp.supercell().linear_index(occ_trans.uccoord)) != unique_indeces.end()) {
          ++m_perturb_it;
          increment();
          return;
        }
        occ_trans.apply_to(tmp);
      }
      DiffTransConfiguration ret_dtc(tmp, m_base_config.diff_trans());
      ret_dtc.set_orbit_name(m_diff_trans_orbit.name());
      *m_current = ret_dtc.canonical_form();
      return;
    }

    /// Checks if current is acceptable to insert
    bool DiffTransConfigEnumPerturbations::_check_current() {
      return current().is_canonical();
    }

    /// Moves to next unique diffusion transformation and places in m_background_config
    void DiffTransConfigEnumPerturbations::_increment_base_config() {
      ++m_unique_difftrans_it;
      if(m_unique_difftrans_it != m_unique_difftrans.end()) {
        ScelPeriodicDiffTransSymCompare symcompare(m_background_config.supercell().prim_grid(),
                                                   m_background_config.supercell().crystallography_tol());
        DiffusionTransformation prepped = symcompare.prepare(*m_unique_difftrans_it);
        Configuration bg_config = make_attachable(prepped, m_background_config);
        DiffTransConfiguration tmp(bg_config, prepped);
        if(!tmp.is_valid_neb()) {
          _increment_base_config();
          return;
        }
        m_base_config = tmp;//.canonical_form();
      }
      else {
        _invalidate();
      }
      return;
    }

    bool has_local_neighborhood_overlap(std::vector<LocalIntegralClusterOrbit> &local_orbits, const Supercell &scel) {
      std::set<int> present;
      std::set<UnitCellCoord> coords;
      for(auto &orbit : local_orbits) {
        for(auto &cluster : orbit) {
          for(int i = 0; i < cluster.size(); ++i) {
            coords.insert(cluster[i]);
          }
        }
      }
      for(auto &coord : coords) {
        if(!present.insert(scel.linear_index(coord)).second) {
          return true;
        }
      }
      //If no set insertion collision then no problems
      return false;
    }

    std::vector<Supercell> viable_supercells(std::vector<LocalIntegralClusterOrbit> &local_orbits, std::vector<Supercell> scel_options) {
      std::vector<Supercell> results;
      for(auto &scel : scel_options) {
        if(!has_local_neighborhood_overlap(local_orbits, scel)) {
          results.push_back(scel);
        }
      }
      return results;
    }

  }
}
