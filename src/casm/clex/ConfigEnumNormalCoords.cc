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
    PrimClex const &primclex,
    jsonParser const &_kwargs,
    Completer::EnumOption const &enum_opt,
    EnumeratorMap const *interface_map) {

    std::vector<ConfigEnumInput> in_configs = make_enumerator_input_configs(primclex, _kwargs, enum_opt, interface_map);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    DoFKey dof;
    try {
      if(!_kwargs.contains("dof")) {
        throw std::runtime_error("Field \"dof\" is required.\n");
      }
      from_json(dof, _kwargs["dof"]);
      DoF::traits(dof);
    }
    catch(std::exception &e) {
      throw std::runtime_error(std::string("Error parsing JSON arguments for ConfigEnumNormalCoords:") + e.what());
    }

    for(ConfigEnumInput const &config : in_configs) {
      Index result = run(primclex,
                         config,
                         dof,
                         filter_expr,
                         CASM::dry_run(_kwargs, enum_opt));
      if(result)
        return result;
    }

    return 0;
  }

  int ConfigEnumNormalCoords::run(PrimClex const &_primclex,
                                  ConfigEnumInput const &_in_config,
                                  DoFKey const &_dof,
                                  std::vector<std::string> const &_filter_expr,
                                  bool dry_run) {
    Configuration tconfig = _in_config.config();

    if(_in_config.sites().size() == 0) {
      tconfig.configdof().local_dof(_dof).values().setZero();
    }
    else {
      for(Index s : _in_config.sites()) {
        tconfig.configdof().local_dof(_dof).site_value(s).setZero();
      }
    }

    ConfigEnumInput config(tconfig, _in_config.sites());

    Eigen::MatrixXd norm_coords = collective_dof_normal_coords(config.sites().begin(),
                                                               config.sites().end(),
                                                               config.supercell().sym_info(),
                                                               _dof,
                                                               config.group()).transpose();

    std::cout << "norm_coords are \n" << norm_coords << "\n";

    auto constructor = [&](const ConfigEnumInput & _config) {
      return notstd::make_unique<ConfigEnumNormalCoords>(config,
                                                         _dof,
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

  ConfigEnumNormalCoords::ConfigEnumNormalCoords(const ConfigEnumInput &_init,
                                                 DoFKey const &_dof,
                                                 Eigen::Ref<const Eigen::MatrixXd> const &_coords):
    m_current(_init.config()),
    m_dof_key(_dof),
    m_sites(_init.sites()),
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

