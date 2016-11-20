#include "casm/clex/SuperConfigEnum.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/clex/ConfigEnumEquivalents.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_SuperConfigEnum_interface() {
    return new CASM::EnumInterface<CASM::SuperConfigEnum>();
  }
}

namespace CASM {

  const std::string CASM_TMP::traits<CASM::SuperConfigEnum>::name = "SuperConfigEnum";

  const std::string CASM_TMP::traits<CASM::SuperConfigEnum>::help =

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
    "    otherwise, both primitive and non-primitive configurations are saved.    \n"
    "\n"
    "  filter: string (optional, default=None)\n"
    "    A query command to use to filter which configurations are kept.          \n"
    "\n"
    "  Examples:\n"
    "    To enumerate super-configurations of listed sub-configurations:\n"
    "     '{ \n"
    "        \"SuperConfigEnum\": {\n"
    "          \"supercells\": { \n"
    "            \"max\": 4, \n"
    "            \"unit_cell\": \"SCEL2_1_2_1_0_0_0\" \n"
    "          },\n"
    "          \"subconfigs\": [\n"
    "            \"SCEL1_1_1_1_0_0_0/0\",\n"
    "            \"SCEL2_1_2_1_0_0_0/0\",\n"
    "            \"SCEL2_1_2_1_0_0_0/1\"\n"
    "          ]\n"
    "        } \n"
    "      }' \n"
    "\n"
    "    To enumerate super-configurations of listed sub-configurations from a \n"
    "    selection file:\n"
    "     '{ \n"
    "        \"SuperConfigEnum\": {\n"
    "          \"supercells\": { \n"
    "            \"max\": 4, \n"
    "            \"unit_cell\": \"SCEL2_1_2_1_0_0_0\" \n"
    "          }, \n"
    "          \"subconfigs\": \"selection_filename\"\n"
    "        } \n"
    "      }' \n"
    "\n";

  int EnumInterface<SuperConfigEnum>::run(PrimClex &primclex, const jsonParser &_kwargs) const {
    if(_kwargs.contains("name")) {
      throw std::invalid_argument(
        "Error in SuperConfigEnum JSON input: 'name' is not allowed in the 'supercell' option");
    }

    // construct the superlattice enumerator
    jsonParser kwargs {_kwargs};
    if(kwargs.is_null()) {
      kwargs = jsonParser::object();
    }

    bool primitive_only = true;
    kwargs.get_if(primitive_only, "primitive_only");

    // default is use all existing Supercells
    jsonParser scel_input;
    scel_input["existing_only"] = true;
    kwargs.get_if(scel_input, "supercells");

    ScelEnumProps enum_props = make_scel_enum_props(primclex, scel_input);

    // the unit cell that will be the supercell of the subconfigs that are
    // tiled to form the super-configurations
    Supercell unit_cell(&primclex, enum_props.generating_matrix());

    // use the unit cell point group
    SupercellEnumerator<Lattice> superlat_enum(
      primclex.prim().lattice(),
      unit_cell.factor_group(),
      enum_props);

    auto check_is_supercell = [&](const Lattice & plat) {
      auto end = primclex.prim().factor_group().end();
      return is_supercell(unit_cell.real_super_lattice(),
                          plat,
                          primclex.prim().factor_group().begin(),
                          end,
                          primclex.crystallography_tol()).first != end;
    };

    // output initial info
    Log &log = primclex.log();
    Index Ninit = std::distance(primclex.config_begin(), primclex.config_end());
    log << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(name());

    // parse a filter expression
    std::vector<std::string> filter_expr;
    if(kwargs.contains("filter")) {
      filter_expr.push_back(kwargs["filter"].get<std::string>());
    };

    // -- get primitive of all subconfigs from input--
    std::map<Configuration, std::string> prim_subconfig;

    if(!kwargs.contains("subconfigs")) {
      throw std::invalid_argument(
        "Error in SuperConfigEnum JSON input: 'subconfigs' is required");
    }
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
          primclex.err_log() << "subconfig transf_mat: \n" << config.supercell().transf_mat() << "\n";
          primclex.err_log() << "prim transf_mat: \n" << pconfig.supercell().transf_mat() << "\n";
          primclex.err_log() << "unit cell: \n" << unit_cell.transf_mat() << "\n";
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
          primclex.err_log() << "subconfig transf_mat: \n" << it->supercell().transf_mat() << "\n";
          primclex.err_log() << "prim transf_mat: \n" << pconfig.supercell().transf_mat() << "\n";
          primclex.err_log() << "unit cell: \n" << unit_cell.transf_mat() << "\n";
          throw std::invalid_argument("Error in SuperConfigEnum: subconfig does not fit in the unit cell");
        }
        prim_subconfig.insert(std::make_pair(pconfig, it->name()));
      }
    }

    log << "\nGenerating equivalents: \n";
    std::vector<Configuration> subconfig;
    for(auto &_pair : prim_subconfig) {
      auto &pconfig = _pair.first;
      FillSupercell f(unit_cell, pconfig, primclex.crystallography_tol());
      ConfigEnumEquivalents e(f(pconfig));
      Index eq_count = 0;
      for(auto it = e.begin(); it != e.end(); ++it) {
        subconfig.push_back(*it);
        log << "  Config: " << it->name() << "\n" << *it << "\n";
      }
    }
    log << std::endl;

    // all supercells of the unit cell, unique by the 'unit_cell' factor group
    for(auto &superlat : superlat_enum) {

      Supercell target_scel(&primclex, superlat);

      Supercell &canon_scel = target_scel.canonical_form();
      log << "Enumerate configurations for " << canon_scel.name() << " ...  " << std::flush;

      // enumerate super-configurations
      SuperConfigEnum superconfig_enum(target_scel, subconfig.begin(), subconfig.end());
      Index num_before = canon_scel.config_list().size();

      if(kwargs.contains("filter")) {
        try {
          auto it = filter_begin(
                      superconfig_enum.begin(),
                      superconfig_enum.end(),
                      filter_expr,
                      primclex.settings().query_handler<Configuration>().dict());
          auto end = filter_end(superconfig_enum.end());
          for(; it != end; ++it) {
            it->insert(primitive_only);
          }
        }
        catch(std::exception &e) {
          primclex.err_log() << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
          return ERR_INVALID_ARG;
        }
      }
      else {
        auto it = superconfig_enum.begin();
        auto end = superconfig_enum.end();
        for(; it != end; ++it) {
          it->insert(primitive_only);
        }
      }
      log << (canon_scel.config_list().size() - num_before) << " configs." << std::endl;
    }
    log << "  DONE." << std::endl << std::endl;

    Index Nfinal = std::distance(primclex.config_begin(), primclex.config_end());

    log << "# new configurations: " << Nfinal - Ninit << "\n";
    log << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    std::cout << "Write SCEL..." << std::endl;
    primclex.print_supercells();
    log << "  DONE" << std::endl << std::endl;

    log << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    log << "  DONE" << std::endl;
    return 0;
  }



  void SuperConfigEnum::_init() {

    // check that all sub-config have same supercell
    for(auto it = sub_config().begin(); it != sub_config().end(); ++it) {
      if(&it->supercell() != &(sub_config().begin()->supercell())) {
        throw std::runtime_error("Error constructing SuperConfigEnum: "
                                 "Sub-Configurations with different Supercells");
      }
    }
    m_sub_scel = &(m_sub_config.begin()->supercell());

    // construct PrimGrid
    m_prim_grid = notstd::make_cloneable<PrimGrid>(
                    _sub_supercell().real_super_lattice(),
                    _target_supercell().real_super_lattice()
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
      UnitCell ref = _sub_supercell().transf_mat().cast<Index>() * prim_grid().unitcell(i);
      for(int j = 0; j < _sub_supercell().num_sites(); ++j) {
        UnitCellCoord uccord = _sub_supercell().uccoord(j) + ref;
        Index linear_index = _target_supercell().linear_index(uccord);
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
