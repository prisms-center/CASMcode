#include "casm/clex/SuperConfigEnum.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/clex/ConfigEnumEquivalents.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/container/Enumerator_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_SuperConfigEnum_interface() {
    return new CASM::EnumInterface<CASM::SuperConfigEnum>();
  }
}

namespace CASM {

  const std::string SuperConfigEnum::enumerator_name = "SuperConfigEnum";

  const std::string SuperConfigEnum::interface_help =

    "SuperConfigEnum: \n\n"

    "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
    "    Indicate supercells to enumerate super-configurations in using ScelEnum\n"
    "    input format, but the \"name\" option is not allowed. See 'ScelEnum'     \n"
    "    description for details. \n\n"

    "  subconfigs: string or JSON array\n"
    "    Indicate supercells to enumerate all occupational configurations in. May \n"
    "    be a JSON array of configuration names, or string specifying a           \n"
    "    configuration selection file anme. By default, all existing supercells   \n"
    "    are used. See 'ScelEnum' description for details. \n\n"

    "  primitive_only: bool (default=true)\n"
    "    If true, only the primitive form of a configuration is saved in the      \n"
    "    configuration list. Otherwise, both primitive and non-primitive          \n"
    "    configurations are saved.    \n\n"

    "  filter: string (optional, default=None)\n"
    "    A query command to use to filter which configurations are kept.          \n"
    "\n"
    "  Examples:\n"
    "    To enumerate super-configurations of listed sub-configurations:\n"
    "     casm enum --method SuperConfigEnum -i \n"
    "     '{ \n"
    "        \"supercells\": { \n"
    "          \"max\": 4, \n"
    "          \"unit_cell\": \"SCEL2_1_2_1_0_0_0\" \n"
    "        },\n"
    "        \"subconfigs\": [\n"
    "          \"SCEL1_1_1_1_0_0_0/0\",\n"
    "          \"SCEL2_1_2_1_0_0_0/0\",\n"
    "          \"SCEL2_1_2_1_0_0_0/1\"\n"
    "        ]\n"
    "      }' \n"
    "\n"
    "    To enumerate super-configurations of listed sub-configurations from a \n"
    "    selection file:\n"
    "     casm enum --method SuperConfigEnum -i \n"
    "     '{ \n"
    "        \"supercells\": { \n"
    "          \"max\": 4, \n"
    "          \"unit_cell\": \"SCEL2_1_2_1_0_0_0\" \n"
    "        }, \n"
    "        \"subconfigs\": \"selection_filename\"\n"
    "      }' \n"
    "\n";

  /// sub-routine for EnumInterface<SuperConfigEnum>::run,
  ///   generates primitive configurations from input
  void _generate_primitives(
    PrimClex &primclex,
    const Supercell &unit_cell,
    const jsonParser &kwargs,
    std::map<Configuration, std::string> &prim_subconfig) {

    auto check_is_supercell = [&](const Lattice & plat) {
      auto end = primclex.get_prim().factor_group().end();
      return is_supercell(unit_cell.get_real_super_lattice(),
                          plat,
                          primclex.get_prim().factor_group().begin(),
                          end,
                          primclex.crystallography_tol()).first != end;
    };

    if(!kwargs.contains("subconfigs")) {
      throw std::invalid_argument(
        "Error in SuperConfigEnum JSON input: 'subconfigs' is required");
    }
    Log &log = primclex.log();
    log << "Input sub-configurations:\n";
    if(kwargs["subconfigs"].is_array()) {
      const jsonParser &j = kwargs["subconfigs"];
      for(auto it = j.begin(); it != j.end(); ++it) {
        const Configuration &config = primclex.configuration(it->get<std::string>());
        log << "  " << config.name() << "\n";
        Configuration pconfig = config.primitive();
        if(!check_is_supercell(pconfig.ideal_lattice())) {
          primclex.err_log().error("Invalid subconfig");
          primclex.err_log() << "subconfig: " << config.name() << "\n";
          primclex.err_log() << "subconfig transf_mat: \n" << config.get_supercell().get_transf_mat() << "\n";
          primclex.err_log() << "prim transf_mat: \n" << pconfig.get_supercell().get_transf_mat() << "\n";
          primclex.err_log() << "unit cell: \n" << unit_cell.get_transf_mat() << "\n";
          throw std::invalid_argument("Error in SuperConfigEnum: subconfig does not fit in the unit cell");
        }
        prim_subconfig.insert(std::make_pair(pconfig, config.name()));
      }
    }
    else {
      ConstConfigSelection selection(primclex, kwargs["subconfigs"].get<fs::path>());
      for(auto it = selection.selected_config_begin(); it != selection.selected_config_end(); ++it) {
        log << "  " << it->name() << "\n";
        Configuration pconfig = it->primitive();
        if(!check_is_supercell(pconfig.ideal_lattice())) {
          primclex.err_log().error("Invalid subconfig");
          primclex.err_log() << "subconfig: " << it->name() << "\n";
          primclex.err_log() << "subconfig transf_mat: \n" << it->get_supercell().get_transf_mat() << "\n";
          primclex.err_log() << "prim transf_mat: \n" << pconfig.get_supercell().get_transf_mat() << "\n";
          primclex.err_log() << "unit cell: \n" << unit_cell.get_transf_mat() << "\n";
          throw std::invalid_argument("Error in SuperConfigEnum: subconfig does not fit in the unit cell");
        }
        prim_subconfig.insert(std::make_pair(pconfig, it->name()));
      }
    }
  }

  /// sub-routine for EnumInterface<SuperConfigEnum>::run,
  ///   generates equivalent configurations of all primitive configurations
  void _generate_equivalents(
    const PrimClex &primclex,
    Supercell &unit_cell,
    const std::map<Configuration, std::string> &prim_subconfig,
    std::vector<Configuration> &subconfig) {

    primclex.log() << "\nGenerating equivalents: \n";
    for(auto &_pair : prim_subconfig) {
      auto &pconfig = _pair.first;
      FillSupercell f(unit_cell, pconfig, primclex.crystallography_tol());
      ConfigEnumEquivalents e(f(pconfig));
      for(auto it = e.begin(); it != e.end(); ++it) {
        subconfig.push_back(*it);
        primclex.log() << "  Config: " << it->name() << "\n" << *it << "\n";
      }
    }
    primclex.log() << std::endl;
  }

  int SuperConfigEnum::run(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    // -- Disallow 'name' & --scelnames options --
    if(_kwargs.contains("name") || enum_opt.vm().count("scelnames")) {
      throw std::invalid_argument(
        "Error in SuperConfigEnum JSON input: 'name' is not allowed in the 'supercell' option");
    }

    // -- Option to include non-primitive --
    bool primitive_only = true;
    _kwargs.get_if(primitive_only, "primitive_only");

    // -- make <SupercellEnumerator<Lattice> & filter_expr from input
    std::unique_ptr<SupercellEnumerator<Lattice> > superlat_enum = make_enumerator_superlat_enum(primclex, _kwargs, enum_opt);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    // -- Unit cell --
    // the unit cell that will be the supercell of the subconfigs that are
    // tiled to form the super-configurations
    Supercell unit_cell(&primclex, *superlat_enum->begin());


    // -- get primitive of all subconfigs from input--
    std::map<Configuration, std::string> prim_subconfig;
    _generate_primitives(
      primclex,
      unit_cell,
      _kwargs,
      prim_subconfig);


    // -- generate equivalents of each primitive config --
    std::vector<Configuration> subconfig;
    _generate_equivalents(
      primclex,
      unit_cell,
      prim_subconfig,
      subconfig);

    // -- Enumerator construction --
    auto lambda = [&](Supercell & scel) {
      return notstd::make_unique<SuperConfigEnum>(scel, subconfig.begin(), subconfig.end());
    };

    int returncode = insert_configs_via_lattice_enum(
                       enumerator_name,
                       primclex,
                       superlat_enum->begin(),
                       superlat_enum->end(),
                       lambda,
                       filter_expr,
                       primitive_only);

    return returncode;
  }

  void SuperConfigEnum::_init() {

    // check that all sub-config have same supercell
    for(auto it = sub_config().begin(); it != sub_config().end(); ++it) {
      if(&it->get_supercell() != &(sub_config().begin()->get_supercell())) {
        throw std::runtime_error("Error constructing SuperConfigEnum: "
                                 "Sub-Configurations with different Supercells");
      }
    }
    m_sub_scel = &(m_sub_config.begin()->get_supercell());

    // construct PrimGrid
    m_prim_grid = notstd::make_cloneable<PrimGrid>(
                    _sub_supercell().get_real_super_lattice(),
                    _target_supercell().get_real_super_lattice()
                  );

    // initialize 'm_counter' to count over all possible sub-config on
    // each prim_grid site
    m_counter = Counter<Array<int> >(
                  Array<int>(prim_grid().size(), 0),
                  Array<int>(prim_grid().size(), sub_config().size() - 1),
                  Array<int>(prim_grid().size(), 1));

    // save indices for mapping occupants into super config
    // so that: m_current->occ(m_index_map[i][j]) = m_sub_scel[i].occ(j)
    // and same for all other site DoF
    m_index_map.resize(prim_grid().size());
    for(int i = 0; i < prim_grid().size(); ++i) {
      UnitCell ref = _sub_supercell().get_transf_mat().cast<Index>() * prim_grid().unitcell(i);
      for(int j = 0; j < _sub_supercell().num_sites(); ++j) {
        UnitCellCoord uccord = _sub_supercell().uccoord(j) + ref;
        Index linear_index = _target_supercell().find(uccord);
        m_index_map[i].push_back(linear_index);
      }
    }

    m_current = notstd::make_cloneable<Configuration>(_target_supercell());

    auto has_occ = [](const Configuration & c) {
      return c.has_occupation();
    };
    m_has_occ = std::any_of(m_sub_config.begin(), m_sub_config.end(), has_occ);

    auto has_disp = [](const Configuration & c) {
      return c.has_displacement();
    };
    m_has_disp = std::any_of(m_sub_config.begin(), m_sub_config.end(), has_disp);

    _initialize(&(*m_current));
    _fill(counter(), _current());

    // Make sure that current() satisfies requested conditions
    if(!_check_current()) {
      increment();
    }

    // set step to 0
    if(valid()) {
      _set_step(0);
    }
    _current().set_source(source(step()));
  }

  // **** Mutators ****
  // increment m_current and return a reference to it
  void SuperConfigEnum::increment() {

    bool is_valid_config {false};

    while(!is_valid_config && ++m_counter) {
      _fill(counter(), _current());
      is_valid_config = _check_current();
    }

    if(m_counter.valid()) {
      _increment_step();
    }
    else {
      _invalidate();
    }
  }

  /// Returns true if current() satisfies requested conditions
  bool SuperConfigEnum::_check_current() const {

    return true;
  }

  /// Fill DoF from sub_config into a Configuration
  ///
  /// \param counter_val The index of the sub_config on each PrimGrid site
  /// \param config The Configuration to set the occupation of
  ///
  void SuperConfigEnum::_fill(const Array<int> &counter_val,
                              Configuration &config) {

    // use m_sub_config, prim_grid(), and counter_val to set occupation in config
    if(m_has_occ) {
      m_current->init_occupation();
    }

    if(m_has_disp) {
      m_current->init_displacement();
    }

    ConfigDoF &to = config.configdof();
    for(Index i = 0; i < prim_grid().size(); ++i) {
      Configuration &from = _sub_config()[counter_val[i]];
      for(Index j = 0; j < _sub_supercell().num_sites(); ++j) {

        // copy site DoF
        if(from.has_occupation()) {
          to.occ(m_index_map[i][j]) = from.occ(j);
        }

        if(from.has_displacement()) {
          to.disp(m_index_map[i][j]) = from.disp(j);
        }
      }
    }
  }

}
