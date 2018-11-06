#include "casm/clex/ConfigEnumNormalCoords.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/clex/ScelEnum.hh"
//#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
//#include "casm/misc/CASM_math.hh"
//#include "casm/misc/CASM_Eigen_math.hh"
//#include "casm/misc/algorithm.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumNormalCoords_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumNormalCoords>();
  }
}

namespace CASM {

  struct MakeConfigInvariantSubgroup {

    MakeConfigInvariantSubgroup() {}

    template<typename PermuteOutputIterator>
    PermuteOutputIterator operator()(const Configuration &config, PermuteIterator begin, PermuteIterator end, PermuteOutputIterator result) {
      ConfigIsEquivalent f(config, config.crystallography_tol());
      return std::copy_if(begin, end, result, f);
    }

  };

  const std::string ConfigEnumNormalCoords::enumerator_name = "ConfigEnumNormalCoords";

  std::string ConfigEnumNormalCoords::interface_help() {
    return
      "ConfigEnumNormalCoords: \n\n"

      "  confignames: Array of strings (optional) \n"
      "    Names of configurations to be used as reference states. Normal coordinates are enum-\n"
      "    erated after zeroing the DoF values at selected sites of the specified configurations\n"
      "    and calculating the resulting symmetry of the selected sites.\n"
      "    Ex: \"confignames\" : [\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]\n\n"

      "  scelnames: Array of strings (optional) \n"
      "    Names of supercells used as reference states. Normal coordinates are enumerated starting\n"
      "    from the fully zeroed configuration of the specified supercells.\n"
      "    Ex: \"scelnames\" : [\"SCEL1_1_1_1_0_0_0\",\"SCEL2_2_1_1_0_0_0\"]\n\n"

      "  dof: string (required) \n"
      "    Name of site degree of freecom  for which normal coordinates are to be generated.\n"
      "    Must be one of the degrees of freedom under consideration in the current project,\n"
      "    as determined by prim.json\n\n"

      "  sublats: array of integers (optional, default none) \n"
      "    Restricts normal coordinate determination to specified sublattices. Each sublat-\n"
      "    tice index specifies the correspondings basis site in prim.json, indexed from 0.\n"
      "    Ex: \"sublats\" : [0,2]\n\n"

      "  sites: array of 4-entry integer arrays (optional, default none) \n"
      "    Restricts normal coordinate determination to specified sites. Sites are specified\n"
      "    in [b,i,j,k] convention, where 'b' is sublattice index and [i,j,k] speifies line-\n"
      "    ar combinations of primitive-cell lattice vectors.\n"
      "    Ex: \"sites\" : [[0,0,0,0],\n"
      "                     [2,0,0,0]]\n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  Examples:\n"
      "    To enumerate all strain perturbations of a particular configuration:\n"
      "      casm enum --method ConfigEnumNormalCoords -i \n"
      "      '{ \n"
      "        \"config\": \"SCEL4_1_4_1_0_0_0/3\",\n"
      "        \"analysis\": true,\n"
      "        } \n"
      "      }' \n\n";
  }

  int ConfigEnumNormalCoords::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    std::vector<Index> sublats;
    std::vector<UnitCellCoord> sites;
    DoFKey dof;
    std::vector<std::string> confignames;
    std::vector<std::string> filter_expr;
    try {
      if(!_kwargs.contains("dof")) {
        throw std::runtime_error("Field \"dof\" is required.\n");
      }
      from_json(dof, _kwargs["dof"]);
      DoF::traits(dof);
      if(_kwargs.contains("confignames")) {
        if(!_kwargs["confignames"].is_array())
          throw std::runtime_error("Field \"confignames\" must specify array of string values.");
        confignames = _kwargs["confignames"].get<std::vector<std::string> >();
      }

      if(_kwargs.contains("sublats"))
        from_json(sublats, _kwargs["sublats"]);

      if(_kwargs.contains("sites"))
        from_json(sites, _kwargs["sites"], primclex.prim());

      filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    }
    catch(std::exception &e) {
      throw std::runtime_error(std::string("Error parsing JSON arguments for ConfigEnumNormalCoords:") + e.what());
    }

    std::unique_ptr<ScelEnum> scel_enum = make_enumerator_scel_enum(primclex, _kwargs, enum_opt);


    for(Supercell const &scel : *scel_enum) {
      Index result = run(primclex,
                         Configuration::zeros(scel),
                         dof,
                         sublats,
                         sites,
                         filter_expr,
                         CASM::dry_run(_kwargs, enum_opt));
      if(result)
        return result;
    }

    for(std::string const &configname : enum_opt.config_strs())
      confignames.push_back(configname);

    std::cout << "confignames is: " << confignames << "\n";
    for(std::string const &configname : confignames) {
      auto it = primclex.const_db<Configuration>().find(configname);
      if(it == primclex.const_db<Configuration>().end())
        throw std::runtime_error("Specified configuration " + configname + " does not exist in database.\n");

      Index result = run(primclex,
                         *it,
                         dof,
                         sublats,
                         sites,
                         filter_expr,
                         CASM::dry_run(_kwargs, enum_opt));
      if(result)
        return result;
    }
    return 0;
  }

  int ConfigEnumNormalCoords::run(PrimClex const &_primclex,
                                  Configuration config,
                                  DoFKey const &_dof,
                                  std::vector<Index> const &_sublats,
                                  std::vector<UnitCellCoord> const &_sites,
                                  std::vector<std::string> const &_filter_expr,
                                  bool dry_run) {

    std::set<Index> pert_sites;
    for(Index b : _sublats) {
      Index V = config.supercell().volume();
      for(Index i = b * V; i < (b + 1)*V; ++i) {
        pert_sites.insert(i);
      }
    }

    for(UnitCellCoord const &ucc : _sites)
      pert_sites.insert(config.supercell().prim_grid().find(ucc));

    if(pert_sites.size() == 0) {
      config.configdof().local_dof(_dof).values().setZero();
    }
    else {
      for(Index s : pert_sites) {
        config.configdof().local_dof(_dof).site_value(s).setZero();
      }
    }

    std::vector<PermuteIterator> pg;
    for(PermuteIterator const &perm_it : config.factor_group()) {
      bool add_it = true;

      for(Index s : pert_sites) {
        if(pert_sites.count(perm_it.permute_ind(s)) == 0) {
          add_it = false;
          break;
        }
      }
      if(add_it)
        pg.push_back(perm_it);
    }


    // if pert_sites is emptywe will fill it now to contain all the sites.
    if(pert_sites.size() == 0) {
      for(Index s = 0; s < config.size(); s++) {
        pert_sites.insert(s);
      }
    }
    Eigen::MatrixXd norm_coords = collective_dof_normal_coords(pert_sites.begin(),
                                                               pert_sites.end(),
                                                               config.supercell().sym_info(),
                                                               _dof,
                                                               pg).transpose();

    std::cout << "norm_coords are \n" << norm_coords << "\n";

    auto constructor = [&](const Supercell & scel) {
      return notstd::make_unique<ConfigEnumNormalCoords>(config,
                                                         _dof,
                                                         pert_sites,
                                                         norm_coords);
    };

    int returncode = insert_configs(enumerator_name,
                                    _primclex,
                                    config.supercell(),
                                    constructor,
                                    _filter_expr,
                                    false,
                                    dry_run);

    return returncode;

  }

  ConfigEnumNormalCoords::ConfigEnumNormalCoords(const Configuration &_init,
                                                 DoFKey const &_dof,
                                                 std::set<Index> const &_sites,
                                                 Eigen::Ref<const Eigen::MatrixXd> const &_coords):
    m_current(_init),
    m_dof_key(_dof),
    m_sites(_sites),
    m_coords(_coords) {


    reset_properties(m_current);
    this->_initialize(&m_current);

    if(m_coords.cols() == 0) {
      this->_invalidate();
    }
    else {
      _set_dof();
    }
    m_current.set_source(this->source(step()));
  }

  void ConfigEnumNormalCoords::_set_dof() {
    Eigen::MatrixXd vals = m_current.configdof().local_dof(m_dof_key).values();
    Index l = 0;
    for(Index s : m_sites) {
      for(Index i = 0; i < vals.rows(); ++i, ++l) {
        vals(i, s) = m_coords(l, step());
      }
    }
    m_current.configdof().set_local_dof(m_dof_key, vals);
  }

  // Implements _increment
  void ConfigEnumNormalCoords::increment() {
    //bool is_valid_config(false);
    //std::cout << "Incrementing...\n";
    if(step() == m_coords.cols() - 1)
      _invalidate();

    if(valid()) {
      _set_dof();
      _increment_step();
      m_current.set_source(this->source(step()));
    }
    //std::cout << "At end, current value is : " << m_current.configdof().global_dof(m_strain_key).values().transpose() << "\n";
    //std::cout << "--FINISHED SEARCH " << _step()<< "--\n";
    return;
  }

}

