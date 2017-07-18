#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_DiffusionTransformationEnum_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffusionTransformationEnum>();
  }
}

namespace CASM {

  namespace Kinetics {

    namespace {

      std::ostream &operator<<(std::ostream &sout, const IntegralCluster &clust) {
        SitesPrinter printer;
        printer.print(clust, std::cout);
        sout << std::endl;
        return sout;
      }
    }

    /// \brief Construct with an IntegralCluster
    DiffusionTransformationEnum::DiffusionTransformationEnum(const IntegralCluster &clust) :
      m_cluster(clust),
      m_current(new DiffusionTransformation(clust.prim())) {

      // initialize to/from counter
      _init_occ_counter();

      // initialize the specie trajectory for the first valid set of to/from occupation values
      m_from_loc = _init_from_loc(m_occ_counter());
      m_to_loc = _init_to_loc(m_occ_counter());

      // set initial DiffTrans
      _set_current();

      if(!m_current->is_valid()) {
        increment();
      }

      if(!m_occ_counter.valid()) {
        _invalidate();
      }
      else {
        this->_initialize(&(*m_current));
        _set_step(0);
      }
    }

    const std::string DiffusionTransformationEnum::enumerator_name = "DiffusionTransformationEnum";
    const std::string DiffusionTransformationEnum::interface_help =
      "DiffusionTransformationEnum: \n\n"

      "  bspecs: JSON object \n"
      "    Indicate clusters to enumerate all occupational diffusion transformations. The \n"
      "    JSON item \"cspecs\" should be a cspecs style initialization of cluster number and sizes.\n"
      "    See below.          \n\n"
      ""
      "  require: JSON array of strings (optional,default=[]) \n "
      "    Indicate required species to enforce that a given species must be a part of the diffusion \n"
      "    transformation. The JSON array \"require\" should be an array of species names.\n"
      "    i.e. \"require\": [\"Va\",\"O\"] \n\n"
      "  exclude: JSON array of strings (optional,default=[]) \n "
      "    Indicate excluded species to enforce that a given species must not be a part of the diffusion \n"
      "    transformation. The JSON array \"exclude\" should be an array of species names.\n"
      "    i.e. \"exclude\": [\"Al\",\"Ti\"] \n\n"
      "  Example:\n"
      "  {\n"
      "   \"require\":[\"Va\"],\n"
      "   \"exclude\":[],\n"
      "    \"cspecs\":{\n"
      "      \"orbit_branch_specs\" : { \n"
      "       \"2\" : {\"max_length\" : 5.01},\n"
      "       \"3\" : {\"max_length\" : 5.01}\n"
      "      }\n"
      "    }\n"
      "  }\n\n"
      ;

    /// Implements increment
    void DiffusionTransformationEnum::increment() {

      // get the next valid DiffTrans
      // to do so, get the next valid specie trajectory
      do {

        // by taking a permutation of possible 'to' specie position
        bool valid_perm = std::next_permutation(m_to_loc.begin(), m_to_loc.end());

        // if no more possible specie trajectory,
        if(!valid_perm) {

          // get next valid from/to occupation values
          do {
            m_occ_counter++;
            _update_current_occ_transform();
          }
          while(m_occ_counter.valid() && !m_current->is_valid_occ_transform());

          // if no more possible from/to occupation values, return
          if(!m_occ_counter.valid()) {
            _invalidate();
            return;
          }
          else {
            m_from_loc = _init_from_loc(m_occ_counter());
            m_to_loc = _init_to_loc(m_occ_counter());
            _update_current_occ_transform();
            _set_current_loc();
          }
        }
        _update_current_to_loc();
      }
      while(!m_current->is_valid_specie_traj());

      _increment_step();
    }

    /// Implements run
    int DiffusionTransformationEnum::run(const PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt) {

      jsonParser kwargs;
      if(!_kwargs.contains("cspecs")) {
        primclex.err_log() << "DiffusionTransformationEnum currently has no default and requires a correct JSON with a bspecs tag within it" << std::endl;
        throw std::runtime_error("Error in DiffusionTransformationEnum: cspecs not found");
      }

      std::vector<std::string> require;
      std::vector<std::string> exclude;
      if(_kwargs.contains("require")) {
        for(auto it = _kwargs["require"].begin(); it != _kwargs["require"].end(); ++it) {
          require.push_back(from_json<std::string>(*it));
        }
      }
      if(_kwargs.contains("exclude")) {
        for(auto it = _kwargs["exclude"].begin(); it != _kwargs["exclude"].end(); ++it) {
          exclude.push_back(from_json<std::string>(*it));
        }
      }
      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
      std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

      auto end = make_prim_periodic_orbits(
                   primclex.prim(), _kwargs["cspecs"], alloy_sites_filter, primclex.crystallography_tol(), std::back_inserter(orbits), primclex.log());

      Log &log = primclex.log();
      auto &db_orbits = primclex.db<PrimPeriodicDiffTransOrbit>();

      Index Ninit = db_orbits.size();
      log << "# diffusion transformations in this project: " << Ninit << "\n" << std::endl;

      log.begin(enumerator_name);


      std::vector< PrimPeriodicDiffTransOrbit > diff_trans_orbits;
      auto end2 = make_prim_periodic_diff_trans_orbits(
                    orbits.begin(), orbits.end(), primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));

      for(auto &diff_trans_orbit : diff_trans_orbits) {
        auto speciemap = diff_trans_orbit.prototype().specie_count();
        auto it = speciemap.begin();
        for(; it != speciemap.end(); ++it) {
          if(std::find(require.begin(), require.end(), it->first.name()) != require.end() && it->second == 0) {
            break;
          }
          if(std::find(exclude.begin(), exclude.end(), it->first.name()) != exclude.end() && it->second != 0) {
            break;
          }
        }
        if(it == speciemap.end()) {
          //insert current into database
          db_orbits.insert(diff_trans_orbit);
        }
      }

      log << "  DONE." << std::endl << std::endl;

      Index Nfinal = db_orbits.size();

      log << "# new diffusion transformations: " << Nfinal - Ninit << "\n";
      log << "# diffusion transformations in this project: " << Nfinal << "\n" << std::endl;

      log << "Writing diffusion transformation database..." << std::endl;
      db_orbits.commit();
      log << "  DONE" << std::endl;

      return 0;
    }

    // -- Unique -------------------

    const Structure &DiffusionTransformationEnum::prim() const {
      return m_cluster.prim();
    }

    const IntegralCluster &DiffusionTransformationEnum::cluster() const {
      return m_cluster;
    }

    /// \brief The occ_counter contains the from/to occupation values for each site
    ///
    /// - Layout is: [from values | to values ]
    ///
    void DiffusionTransformationEnum::_init_occ_counter() {
      Index N = cluster().size();
      std::vector<Index> init_occ(N * 2, 0);
      std::vector<Index> final_occ(N * 2);
      for(int i = 0; i < N; i++) {
        final_occ[i] = cluster()[i].site().site_occupant().size() - 1;
        final_occ[i + N] = final_occ[i];
      }
      std::vector<Index> incr(N * 2, 1);

      m_occ_counter = Counter<std::vector<Index> >(init_occ, final_occ, incr);
    }

    /// \brief Returns container of 'from' specie locations
    std::vector<SpecieLocation> DiffusionTransformationEnum::_init_from_loc(const std::vector<Index> &occ_values) {
      return _init_loc(occ_values, 0);
    }

    /// \brief Returns container of 'to' specie locations
    std::vector<SpecieLocation> DiffusionTransformationEnum::_init_to_loc(const std::vector<Index> &occ_values) {
      return _init_loc(occ_values, cluster().size());
    }

    /// \brief Returns container of 'from' or 'to' specie locations
    ///
    /// - offset == 0 for 'from', N for 'to' specie locations
    ///
    std::vector<SpecieLocation> DiffusionTransformationEnum::_init_loc(const std::vector<Index> &occ_values, Index offset) {

      Index N = cluster().size();
      std::vector<SpecieLocation> loc;
      // for each 'from' occupant
      for(Index i = 0; i < N; ++i) {
        Index occ = occ_values[i + offset];
        UnitCellCoord uccoord = cluster()[i];
        Index mol_size = uccoord.site().site_occupant()[occ].size();
        // for each specie
        for(Index j = 0; j < mol_size; ++j) {
          loc.emplace_back(uccoord, occ, j);
        }
      }
      return loc;
    }

    /// \brief Uses m_cluster, m_occ_counter, m_from_loc, and m_to_loc to set m_current
    void DiffusionTransformationEnum::_set_current() {
      m_current->occ_transform().clear();
      m_current->occ_transform().reserve(cluster().size());
      for(const auto &uccoord : cluster()) {
        m_current->occ_transform().emplace_back(uccoord, 0, 0);
      }
      _update_current_occ_transform();
      _set_current_loc();
      _update_current_to_loc();
    }

    void DiffusionTransformationEnum::_update_current_occ_transform() {
      auto N = cluster().size();
      Index i = 0;
      for(auto &t : m_current->occ_transform()) {
        t.from_value = m_occ_counter()[i];
        t.to_value = m_occ_counter()[i + N];
        ++i;
      }
    }

    void DiffusionTransformationEnum::_set_current_loc() {
      m_current->specie_traj().clear();
      m_current->specie_traj().reserve(m_from_loc.size());
      for(const auto &t : m_from_loc) {
        m_current->specie_traj().emplace_back(t, t);
      }
    }

    void DiffusionTransformationEnum::_update_current_to_loc() {
      auto it = m_to_loc.begin();
      for(auto &t : m_current->specie_traj()) {
        t.to = *it++;
      }
    }

  }
}
