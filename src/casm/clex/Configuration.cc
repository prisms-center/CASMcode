#include "casm/clex/Configuration.hh"

#include <sstream>
//#include "casm/misc/Time.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/clex/ConfigCompare.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigIOSelected.hh"
#include "casm/app/QueryHandler_impl.hh"

namespace CASM {

  const std::string QueryTraits<Configuration>::name = "Configuration";

  template class QueryHandler<Configuration>;

  namespace {
    typedef std::insert_iterator<std::map<std::string, std::shared_ptr<RuntimeLibrary> > > runtimelib_it_type;
    typedef std::insert_iterator<DataFormatterDictionary<Configuration> > config_dict_it_type;
  }

  template std::pair<config_dict_it_type, runtimelib_it_type> load_query_plugins(
    const ProjectSettings &set,
    config_dict_it_type dict_it,
    runtimelib_it_type lib_it);



  /// Construct a default Configuration
  Configuration::Configuration(Supercell &_supercell, const jsonParser &src, const ConfigDoF &_configdof)
    : id("none"), supercell(&_supercell), source_updated(true), multiplicity(-1),
      m_configdof(_configdof), prop_updated(true), m_selected(false) {
    set_source(src);
  }

  //*********************************************************************************
  /// Construct by reading from main data file (json)
  Configuration::Configuration(const jsonParser &json, Supercell &_supercell, Index _id)
    : supercell(&_supercell), source_updated(false), multiplicity(-1),
      m_configdof(_supercell.num_sites()), prop_updated(false) {

    std::stringstream ss;
    ss << _id;
    id = ss.str();

    read(json);
  }

  //********** MUTATORS  ***********

  void Configuration::set_id(Index _id) {
    std::stringstream ss;
    ss << _id;
    id = ss.str();

    prop_updated = true;
  }

  //*********************************************************************************
  void Configuration::set_source(const jsonParser &source) {
    if(source.is_null() || source.size() == 0) {
      m_source.put_array();
    }
    else if(!source.is_array()) {
      m_source.put_array();
      m_source.push_back(source);
    }
    else {
      m_source = source;
    }
    source_updated = true;
  }

  //*********************************************************************************
  void Configuration::push_back_source(const jsonParser &source) {

    if(source.is_null() || source.size() == 0) {
      return;
    }
    if(!source.is_array()) {

      // check if the new source is already listed, if it is do nothing
      for(int i = 0; i < m_source.size(); i++) {
        if(m_source[i] == source)
          return;
      }

      // else, add the new source
      m_source.push_back(source);

      source_updated = true;
    }
    else {

      // check all new sources, if already listed skip, if the any of the new sources is already listed, if it is do nothing

      for(int s = 0; s < source.size(); s++) {

        bool found = false;

        for(int i = 0; i < m_source.size(); i++) {
          if(m_source[i] == source[s]) {
            found = true;
            break;
          }
        }

        if(!found) {
          // else, add the new source
          m_source.push_back(source[s]);

          source_updated = true;
        }
      }
    }
  }

  //*********************************************************************************
  void Configuration::clear() {
    _invalidate_id();
    m_configdof.clear();
  }

  //*********************************************************************************
  void Configuration::init_occupation() {
    set_occupation(Array<int>(this->size(), 0));
  }

  //*********************************************************************************
  void Configuration::set_occupation(const Array<int> &new_occupation) {
    _invalidate_id();
    if(new_occupation.size() != this->size()) {
      default_err_log().error("Configuration::set_occupation size error");
      default_err_log() << "new_occupation.size(): " << new_occupation.size() << std::endl;
      default_err_log() << "Configuration size(): " << this->size() << std::endl;
      throw std::runtime_error("Error: Configuration::set_occupation with array of the wrong size");
    }
    m_configdof.set_occupation(new_occupation);
    return;
  }

  //*********************************************************************************
  void Configuration::set_occ(Index site_l, int val) {
    _invalidate_id();
    m_configdof.occ(site_l) = val;
  }

  //*********************************************************************************
  void Configuration::clear_occupation() {
    m_configdof.clear_occupation();
  }

  //*********************************************************************************
  void Configuration::init_displacement() {
    set_displacement(displacement_matrix_t::Zero(3, this->size()));
  }

  //*********************************************************************************
  void Configuration::set_displacement(const displacement_matrix_t &new_displacement) {
    _invalidate_id();
    if(new_displacement.cols() != this->size()) {
      default_err_log().error("Configuration::set_displacement size error");
      default_err_log() << "new_displacement.cols(): " << new_displacement.cols() << std::endl;
      default_err_log() << "Configuration size(): " << this->size() << std::endl;
      throw std::runtime_error("Error: Configuration::set_displacement with matrix of the wrong size");
    }
    m_configdof.set_displacement(new_displacement);
  }

  //*********************************************************************************
  void Configuration::set_disp(Index site_l, const Eigen::VectorXd &_disp) {
    _invalidate_id();
    m_configdof.disp(site_l) = _disp;
  }

  //*********************************************************************************
  void Configuration::clear_displacement() {
    _invalidate_id();
    m_configdof.clear_displacement();
  }

  //*********************************************************************************
  void Configuration::init_deformation() {
    set_deformation(Eigen::Matrix3d::Identity());
  }

  //*********************************************************************************

  void Configuration::set_deformation(const Eigen::Matrix3d &new_deformation) {
    _invalidate_id();
    m_configdof.set_deformation(new_deformation);
  }

  //*********************************************************************************

  void Configuration::clear_deformation() {
    _invalidate_id();
    m_configdof.clear_deformation();
  }

  //*******************************************************************************

  /// \brief Check if this is a primitive Configuration
  bool Configuration::is_primitive() const {
    return (find_translation() == get_supercell().translate_end());
  }

  //*******************************************************************************

  /// \brief Returns a PermuteIterator corresponding to the first non-zero pure
  /// translation that maps the Configuration onto itself.
  ///
  /// - If primitive, returns this->get_supercell().translate_end()
  PermuteIterator Configuration::find_translation() const {
    ConfigIsEquivalent f(*this, crystallography_tol());
    const Supercell &scel = get_supercell();
    auto begin = scel.translate_begin();
    auto end = scel.translate_end();
    if(++begin == end) {
      return end;
    }
    return std::find_if(begin, end, f);
  }

  //*******************************************************************************

  /// \brief Return a primitive Configuration
  ///
  /// - The result holds its own Supercell, so it must be considered only a
  ///   temporary Configuration
  Configuration Configuration::primitive() const {
    Configuration tconfig {*this};
    /*
    std::cout << "T: \n" << tconfig.get_supercell().get_transf_mat() << std::endl;
    std::cout << "L: \n" << tconfig.get_supercell().get_real_super_lattice().lat_column_mat() << std::endl;
    */

    std::unique_ptr<Supercell> next_scel;

    // check if config is primitive, and if not, obtain a translation that maps the config on itself
    while(true) {

      PermuteIterator result = tconfig.find_translation();
      if(result == tconfig.get_supercell().translate_end()) {
        break;
      }

      // replace one of the lattice vectors with the translation
      Lattice new_lat = replace_vector(
                          tconfig.ideal_lattice(),
                          result.sym_op().tau(),
                          crystallography_tol()).make_right_handed().get_reduced_cell();

      next_scel.reset(new Supercell(&get_primclex(), new_lat));
      /*
      std::cout << "T: \n" << next_scel->get_transf_mat() << std::endl;
      std::cout << "L: \n" << next_scel->get_real_super_lattice().lat_column_mat() << std::endl;
      */

      // create a sub configuration in the new supercell
      tconfig = sub_configuration(*next_scel, tconfig);
      //std::cout << "sub occ: \n" << tconfig.occupation() << std::endl;

      tconfig.m_supercell_ptr.reset(next_scel.release());
      tconfig.supercell = tconfig.m_supercell_ptr.get();

    }

    return tconfig;
  }

  //*******************************************************************************

  /// \brief Check if Configuration is in the canonical form
  bool Configuration::is_canonical() const {
    const Supercell &scel = get_supercell();
    ConfigIsEquivalent f(*this, crystallography_tol());
    return std::all_of(
             ++scel.permute_begin(),
             scel.permute_end(),
    [&](const PermuteIterator & p) {
      return f(p) || !f.is_less();
    });
  }

  //*******************************************************************************

  /// \brief Returns the operation that applied to *this returns the canonical form
  PermuteIterator Configuration::to_canonical() const {
    ConfigCompare f(*this, crystallography_tol());
    const Supercell &scel = get_supercell();
    return std::max_element(scel.permute_begin(), scel.permute_end(), f);
  }

  //*******************************************************************************

  /// \brief Returns the operation that applied to the canonical form returns *this
  PermuteIterator Configuration::from_canonical() const {
    return to_canonical().inverse();
  }

  //*******************************************************************************

  /// \brief Returns the canonical form Configuration in the same Supercell
  Configuration Configuration::canonical_form() const {
    return copy_apply(to_canonical(), (*this));
  }

  //*******************************************************************************

  /// \brief Returns the canonical form Configuration in the canonical Supercell
  ///
  /// - Will be a Supercell included in the PrimClex.get_supercell_list()
  Configuration Configuration::in_canonical_supercell() const {

    Supercell &canon_scel = get_supercell().canonical_form();

    FillSupercell f(canon_scel, *this, crystallography_tol());
    Configuration in_canon = f(*this);

    // only OK to use if canon_scel and this->get_supercell() are stored in
    //   primclex supercell list
    //Configuration in_canon = get_primclex().fill_supercell(canon_scel, *this);

    return in_canon.canonical_form();
  }

  //*******************************************************************************

  /// \brief Insert this configuration (in canonical form) in the canonical Supercell config list
  ///
  /// \param primitive_only If true, only the primitive Configuration is inserted.
  ///
  /// - By convention, the primitive canonical form of a configuration must
  ///   always be saved in the config list.
  /// - By default, both the primitive canonical Configuration and the equivalent
  ///   non-primitive Configuration in the canonical Supercell are saved
  /// - Optionally, this can insert just the primitive Configuration
  ///
  ConfigInsertResult Configuration::insert(bool primitive_only) const {

    ConfigInsertResult res;

    Configuration pconfig = this->primitive().in_canonical_supercell();
    Supercell &canon_scel = pconfig.get_supercell();
    Index config_index;

    res.insert_primitive = canon_scel.add_canon_config(pconfig, config_index);

    res.primitive_it = PrimClex::config_const_iterator(
                         &get_primclex(),
                         canon_scel.get_id(),
                         config_index);

    // if the primitive supercell is the same as the equivalent canonical supercell
    if(get_supercell().canonical_form() == pconfig.get_supercell()) {
      res.insert_canonical = res.insert_primitive;
      res.canonical_it = res.primitive_it;
    }
    else {
      if(primitive_only) {
        res.insert_canonical = false;
      }
      else {
        // primitive is returned as canonical form in canonical supercell
        Supercell &canon_scel = get_supercell().canonical_form();
        Index config_index;
        Supercell::permute_const_iterator permute_it;

        res.insert_canonical = canon_scel.add_config(this->in_canonical_supercell(), config_index, permute_it);

        res.canonical_it = PrimClex::config_const_iterator(
                             &get_primclex(),
                             canon_scel.get_id(),
                             config_index);
      }
    }
    return res;
  }

  //*******************************************************************************

  /// \brief Returns the subgroup of the Supercell factor group that leaves the
  ///        Configuration unchanged
  std::vector<PermuteIterator> Configuration::factor_group() const {
    std::vector<PermuteIterator> fg;
    ConfigIsEquivalent f(*this, crystallography_tol());
    const Supercell &scel = get_supercell();
    std::copy_if(scel.permute_begin(), scel.permute_end(), std::back_inserter(fg), f);
    return fg;
  }

  //*******************************************************************************

  /// \brief Returns the point group that leaves the Configuration unchanged
  SymGroup Configuration::point_group() const {
    SymGroup sym_group;
    sym_group.set_lattice(ideal_lattice());
    std::vector<PermuteIterator> config_factor_group;
    config_factor_group = factor_group();
    bool new_symop;
    for(int i = 0; i < config_factor_group.size(); i++) {
      new_symop = true;
      if(i > 0) {
        if(config_factor_group[i].factor_group_index() == config_factor_group[i - 1].factor_group_index())
          new_symop = false;
      }
      if(new_symop)
        sym_group.push_back(config_factor_group[i].sym_op());
    }
    return sym_group;
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  Configuration Configuration::fill_supercell(Supercell &scel, const SymOp &op) const {
    FillSupercell f(scel, op);
    return f(*this);

    // only OK to use if both supercells are stored in primclex supercell list:
    //return get_primclex().fill_supercell(scel, *this, op);
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  ///
  /// - Uses the first symop in g that such that scel is a supercell of op*(*this)
  Configuration Configuration::fill_supercell(Supercell &scel, const SymGroup &g) const {

    auto res = is_supercell(
                 scel.get_real_super_lattice(),
                 ideal_lattice(),
                 g.begin(),
                 g.end(),
                 crystallography_tol());

    if(res.first == g.end()) {

      std::cerr << "Requested supercell transformation matrix: \n"
                << scel.get_transf_mat() << "\n";
      std::cerr << "Requested motif Configuration: " <<
                name() << "\n";
      std::cerr << "Configuration transformation matrix: \n"
                << get_supercell().get_transf_mat() << "\n";

      throw std::runtime_error(
        "Error in 'Configuration::fill_supercell(const Supercell &scel, const SymGroup& g)'\n"
        "  The motif cannot be tiled onto the specified supercell."
      );
    }

    return fill_supercell(scel, *res.first);
  }

  //*********************************************************************************
  void Configuration::set_calc_properties(const jsonParser &calc) {
    prop_updated = true;
    calculated = calc;
  }

  //*********************************************************************************

  bool Configuration::read_calc_properties(jsonParser &parsed_props) const {
    //std::cout << "begin Configuration::read_calculated()" << std::endl;
    bool success = true;
    /// properties.calc.json: contains calculated properties
    ///   For default clex calctype only
    fs::path filepath = calc_properties_path();
    //std::cout << "filepath: " << filepath << std::endl;
    parsed_props = jsonParser();
    if(fs::exists(filepath)) {
      jsonParser json(filepath);

      //Record file timestamp
      parsed_props["data_timestamp"] = fs::last_write_time(filepath);

      std::vector<std::string> props = get_primclex().settings().properties();
      for(Index i = 0; i < props.size(); i++) {
        //std::cout << "checking for: " << props[i] << std::endl;
        if(json.contains(props[i])) {

          // normal by #prim cells for some properties
          if(props[i] == "energy" || props[i] == "relaxed_energy" || props[i] == "relaxed_magmom") {
            parsed_props[ props[i] ] = json[props[i]].get<double>() / get_supercell().volume();
          }
          else {
            parsed_props[props[i]] = json[props[i]];
          }
        }
        else
          success = false;
      }
      //Get relaxed magmom:
      if(json.contains("relaxed_magmom")) {
        parsed_props["relaxed_magmom"] = json["relaxed_magmom"].get<double>() / get_supercell().volume();
      }
      //Get RMS force:
      if(json.contains("relaxed_forces")) {
        if(json["relaxed_forces"].size()) {
          Eigen::MatrixXd forces;
          from_json(forces, json["relaxed_forces"]);
          parsed_props["rms_force"] = sqrt((forces.transpose() * forces).trace() / double(forces.rows()));
        }
        else {
          parsed_props["rms_force"] = 0.;
        }
      }
      //Get Magnetic moment per site:
      if(json.contains("relaxed_mag_basis")) {

        // get the number of each molecule type
        std::vector<int> num_each_molecule;
        std::vector<std::string> name_each_molecule;
        from_json(num_each_molecule, json["atoms_per_type"]);
        from_json(name_each_molecule, json["atom_type"]);

        // need to create an unsort_dict to put measured properties in the 'right' place
        auto struc_molecule_name = get_prim().get_struc_molecule_name();

        // Initialize container
        Eigen::VectorXd mag_each_molecule = Eigen::VectorXd::Constant(struc_molecule_name.size(), std::nan(""));

        Index i; // Index for which atom_type we're looking at
        Index j; // Index for which individual atom of that atom_type
        Index k = 0; // Global index over all atoms in the configuration
        double cum_molecule_mag; // Holder for the cumulative magmom for the current atom_type

        for(i = 0; i < num_each_molecule.size(); i++) {
          cum_molecule_mag = 0;
          for(j = 0; j < num_each_molecule[i]; j++) {
            cum_molecule_mag += json["relaxed_mag_basis"][k].get<double>();
            k += 1;
          }
          auto atom_idx = std::find(struc_molecule_name.begin(), struc_molecule_name.end(), name_each_molecule[i]) - struc_molecule_name.begin();
          if(atom_idx < struc_molecule_name.size()) {
            mag_each_molecule[atom_idx] = cum_molecule_mag / num_each_molecule[i];
          }
          parsed_props["relaxed_mag"] = mag_each_molecule;
        }
      }
    }
    else
      success = false;

    return success;
  }

  //********** ACCESSORS ***********

  const Lattice &Configuration::ideal_lattice()const {
    return get_supercell().get_real_super_lattice();
  }

  //*********************************************************************************

  std::string Configuration::get_id() const {
    return id;
  }

  //*********************************************************************************
  /// \brief Returns a Configuration name
  ///
  /// One of the following formats:
  /// - `$CANON_SCELNAME/$CANON_INDEX`
  ///   - For canonical forms in canonical supercells, whether primitive or not
  ///   - CANON_INDEX will be "none" if not in config list
  /// - `$PRIM_SCELNAME/$PRIM_CANON_INDEX.equiv.$FG_PERM.$TRANS_PERM`
  ///   - For primitive, but non-canonical configurations in a canonical supercell
  ///   - Primitive canonical form must exist already in config list or PRIM_CANON_INDEX will be "none"
  ///   - Applies PermuteIterator(FG_PERM, TRANS_PERM) to primitive canonical configuration
  /// - `$CANON_SCELNAME.$PRIM_FG_OP1/super.$PRIM_FG_OP2.$PRIM_SCELNAME/$PRIM_CANON_INDEX.equiv.$FG_PERM.$TRANS_PERM`
  ///   - If the supercell is non-canonical, or the configuration is non-primitive and non-canonical
  ///   - Primitive canonical form must exist already in config list or PRIM_CANON_INDEX will be "none"
  ///   - Applies PermuteIterator(FG_PERM, TRANS_PERM) to primitive canonical configuration
  ///   - Then applies prim Structure factor group op with index PRIM_FG_OP and
  ///     fills the supercell $CANON_SCELNAME.$PRIM_FG_OP1
  ///
  std::string Configuration::name() const {
    if(m_name.empty()) {
      _generate_name();
    }
    return m_name;
  }

  void Configuration::_generate_name() const {
    m_name = get_supercell().get_name() + "/" + get_id();
  }

  //*********************************************************************************
  std::string Configuration::calc_status() const {
    if(fs::exists(calc_status_path())) {
      jsonParser json(calc_status_path());
      if(json.contains("status"))
        return json["status"].get<std::string>();
    }
    return("not_submitted");
  }

  //*********************************************************************************
  std::string Configuration::failure_type() const {
    if(fs::exists(calc_status_path())) {
      jsonParser json(calc_status_path());
      if(json.contains("failure_type"))
        return json["failure_type"].get<std::string>();
    }
    return("none");
  }

  //*********************************************************************************
  const jsonParser &Configuration::source() const {
    return m_source;
  }

  //*********************************************************************************
  fs::path Configuration::get_path() const {
    return get_supercell().get_path() / get_id();
  }

  //*********************************************************************************
  ///Returns number of sites, NOT the number of primitives that fit in here
  Index Configuration::size() const {
    return get_supercell().num_sites();
  }

  //*********************************************************************************
  const Structure &Configuration::get_prim() const {
    return get_supercell().get_prim();
  }

  //*********************************************************************************
  //PrimClex &Configuration::get_primclex() {
  //return get_supercell().get_primclex();
  //}

  //*********************************************************************************
  PrimClex &Configuration::get_primclex() const {
    return get_supercell().get_primclex();
  }

  //*********************************************************************************
  Supercell &Configuration::get_supercell() const {
    return *supercell;
  }

  //*********************************************************************************
  double Configuration::crystallography_tol() const {
    return get_primclex().settings().crystallography_tol();
  }

  //*********************************************************************************
  UnitCellCoord Configuration::get_uccoord(Index site_l) const {
    return get_supercell().uccoord(site_l);
  }

  //*********************************************************************************
  int Configuration::get_b(Index site_l) const {
    return get_supercell().get_b(site_l);
  }

  //*********************************************************************************
  const Molecule &Configuration::get_mol(Index site_l) const {
    return get_prim().basis[ get_b(site_l) ].site_occupant()[ occ(site_l) ];
  }

  //*********************************************************************************
  const Properties &Configuration::calc_properties() const {
    return calculated;
  }

  //*********************************************************************************
  /*
    const DeltaProperties &Configuration::delta_properties() const {
      return delta;
    }
  */
  //*********************************************************************************

  const Properties &Configuration::generated_properties() const {
    return generated;
  }

  //*********************************************************************************

  /*const Correlation &Configuration::get_correlations() const {
    return correlations;
    }*/

  //*********************************************************************************

  /// Returns composition on each sublattice: sublat_comp[ prim basis site / sublattice][ molecule_type]
  ///   molucule_type is ordered as in the Prim structure's site_occupant list for that basis site (includes vacancies)
  ReturnArray< Array < double > > Configuration::get_sublattice_composition() const {

    // get the number of each molecule
    Array< Array < int > > sublat_num_each_molecule = get_sublat_num_each_molecule();
    Array< Array < double > > sublattice_composition(sublat_num_each_molecule.size());

    // divide by number of sites per sublattice ( supercell volume )
    for(Index i = 0; i < sublat_num_each_molecule.size(); i++) {
      sublattice_composition[i].resize(sublat_num_each_molecule[i].size());
      for(Index j = 0; j < sublat_num_each_molecule[i].size(); j++) {
        sublattice_composition[i][j] = (1.0 * sublat_num_each_molecule[i][j]) / get_supercell().volume();
      }
    }

    return sublattice_composition;
  }

  //*********************************************************************************
  /// Returns number of each molecule by sublattice:
  ///   sublat_num_each_molecule[ prim basis site / sublattice ][ molecule_type]
  ///   molucule_type is ordered as in the Prim structure's site_occupant list for that basis site
  ReturnArray< Array<int> > Configuration::get_sublat_num_each_molecule() const {

    Index i;

    // [basis_site][site_occupant_index]
    auto convert = get_index_converter(get_prim(), get_prim().get_struc_molecule());

    // create an array to count the number of each molecule
    Array< Array<int> > sublat_num_each_molecule;
    for(i = 0; i < get_prim().basis.size(); i++) {
      sublat_num_each_molecule.push_back(Array<int>(get_prim().basis[i].site_occupant().size(), 0));
    }

    // count the number of each molecule by sublattice
    for(i = 0; i < size(); i++) {
      sublat_num_each_molecule[ get_b(i) ][occ(i)]++;
    }

    return sublat_num_each_molecule;
  }

  //*********************************************************************************
  /// Returns composition, not counting vacancies
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's get_struc_molecule(), with [Va]=0.0
  ReturnArray<double> Configuration::get_composition() const {

    // get the number of each molecule type
    Array<int> num_each_molecule = get_num_each_molecule();

    /// get the total number of non-vacancy atoms
    int num_atoms = 0;

    // need to know which molecules are vacancies
    auto struc_molecule = get_prim().get_struc_molecule();

    Index i;
    for(i = 0; i < struc_molecule.size(); i++) {
      if(struc_molecule[i].is_vacancy()) {
        // set to zero, so the Va concentration is reported as 0.0
        num_each_molecule[i] = 0;
      }
      num_atoms += num_each_molecule[i];
    }

    // calculate the comp (not including vacancies) from the number of each molecule
    Array<double> comp;
    for(i = 0; i < num_each_molecule.size(); i++)
      comp.push_back((1.0 * num_each_molecule[i]) / double(num_atoms));

    return comp;

  }

  //*********************************************************************************
  /// Returns composition, including vacancies
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's get_struc_molecule()
  ReturnArray<double> Configuration::get_true_composition() const {

    Array<int> num_each_molecule = get_num_each_molecule();

    // calculate the true_comp (including vacancies) from the number of each molecule
    Array<double> comp;
    for(Index i = 0; i < num_each_molecule.size(); i++)
      comp.push_back((1.0 * num_each_molecule[i]) / size());

    return comp;
  }

  //*********************************************************************************
  /// Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::get_struc_molecule()
  ReturnArray<int> Configuration::get_num_each_molecule() const {
    return CASM::get_num_each_molecule(m_configdof, get_supercell());
  }

  //*********************************************************************************
  /// Returns parametric composition, as calculated using PrimClex::param_comp
  Eigen::VectorXd Configuration::get_param_composition() const {
    if(!get_primclex().has_composition_axes()) {
      std::cerr << "Error in Configuration::get_param_composition()" << std::endl;
      std::cerr << "  Composition axes are not set." << std::endl;
      exit(1);
    }

    return get_primclex().composition_axes().param_composition(get_num_each_component());
  }

  //*********************************************************************************
  /// Returns num_each_component[ component_type] per prim cell,
  ///   where 'component_type' is ordered as ParamComposition::get_components
  Eigen::VectorXd Configuration::get_num_each_component() const {

    // component order used for param_composition
    std::vector<std::string> v_components = get_primclex().composition_axes().components();

    // copy to CASM::Array
    std::vector<std::string> components;
    for(auto it = v_components.cbegin(); it != v_components.cend(); ++it) {
      components.push_back(*it);
    }

    // initialize
    Eigen::VectorXd num_each_component = Eigen::VectorXd::Zero(components.size());

    // [basis_site][site_occupant_index]
    auto convert = get_index_converter(get_prim(), components);

    // count the number of each component
    for(Index i = 0; i < size(); i++) {
      num_each_component[ convert[ get_b(i) ][occ(i)] ] += 1.0;
    }

    // normalize per prim cell
    for(Index i = 0; i < components.size(); i++) {
      num_each_component[i] /= get_supercell().volume();
    }

    return num_each_component;
  }



  //********* IO ************


  /// Writes the Configuration to a json object (the config list)
  ///   Uses PrimClex's current default settings to write the appropriate properties
  ///
  ///   'json' is a jsonParser JSON object (or will be set to a JSON object)
  ///
  ///   write_dof, source, selected: (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["dof"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["source"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["selected"]
  ///
  ///   write_properties: (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["generated"]
  jsonParser &Configuration::write(jsonParser &json) const {

    //std::cout << "begin Configuration::write()" << std::endl;

    const ProjectSettings &set = get_primclex().settings();
    std::string calc_string = "calctype." + set.default_clex().calctype;
    std::string ref_string = "ref." + set.default_clex().ref;

    /// write json object hierarchy if not existing
    jsonParser &json_scel = json["supercells"][get_supercell().get_name()];
    jsonParser &json_config = json_scel[get_id()];
    jsonParser &json_ref = json_config[calc_string][ref_string];
    jsonParser &json_prop = json_ref["properties"];

    json_config["selected"] = selected();

    if(!json_config.contains("dof")) {
      write_dof(json_config);
    }

    if(source_updated) {
      write_source(json_config);
    }

    if(prop_updated) {
      write_properties(json_prop);
    }

    //std::cout << "finish Configuration::write()" << std::endl;

    return json;
  }

  //*********************************************************************************

  void Configuration::write_pos() const {

    try {
      fs::create_directories(get_path());
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in Configuration::write_pos()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream file(get_pos_path());
    VaspIO::PrintPOSCAR p(*this);
    p.sort();
    p.print(file);
    return;
  }

  //*********************************************************************************

  void Configuration::print_occupation(std::ostream &stream) const {
    stream << occupation() << "\n";
    return;
  }

  //*********************************************************************************

  void Configuration::print_config_list(std::ostream &stream, int composition_flag) const {

    stream.width(10);
    stream.flags(std::ios::left);
    stream << id << " ";

    stream.width(10);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    stream << name() << " ";
    //Prints composition if comp_flag=1, true_composition if comp_flag=2
    // and the sublattice composition if comp_flag=3
    if(composition_flag == 1) {
      print_composition(stream);
    }
    else if(composition_flag == 2) {
      print_true_composition(stream);
    }
    else if(composition_flag == 3) {
      print_sublattice_composition(stream);
    }

    stream.width(8);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream << selected();

    stream << "\n";
  }

  //*********************************************************************************
  void Configuration::print_composition(std::ostream &stream) const {

    Array<double> comp = get_composition();
    auto mol_list = get_prim().get_struc_molecule();

    for(Index i = 0; i < mol_list.size(); i++) {
      if(mol_list[i].is_vacancy()) {
        continue;
      }
      stream.precision(6);
      stream.width(12);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream << comp[i] << " ";
    }

  }

  //*********************************************************************************
  void Configuration::print_true_composition(std::ostream &stream) const {

    Array<double> true_comp = get_true_composition();

    for(Index i = 0; i < true_comp.size(); i++) {
      stream.precision(6);
      stream.width(12);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream << true_comp[i] << " ";
    }

  }

  //*********************************************************************************
  void Configuration::print_sublattice_composition(std::ostream &stream) const {

    Array< Array<double> > sublattice_composition = get_sublattice_composition();

    for(Index i = 0; i < sublattice_composition.size(); i++) {
      for(Index j = 0; j < sublattice_composition[i].size(); j++) {
        stream.precision(6);
        stream.width(12);
        stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
        stream << sublattice_composition[i][j] << " ";
      }
    }

  }

  //*********************************************************************************

  /// Private members:

  /// Reads the Configuration from the config list
  ///   Uses PrimClex's current default settings to read in the appropriate properties
  ///
  /// This is private, because it is only called from the constructor:
  ///   Configuration(const Supercell &_supercell, Index _id)
  ///   It's called from the constructor because of the Supercell pointer
  ///
  ///   read dof, source, selected: (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["dof"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["source"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["selected"]
  ///
  ///   read properties: (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["generated"]
  ///
  void Configuration::read(const jsonParser &json) {

    //std::cout << "begin  Configuration::read()" << std::endl;

    const ProjectSettings &set = get_primclex().settings();

    std::string calc_string = "calctype." + set.default_clex().calctype;
    std::string ref_string = "ref." + set.default_clex().ref;

    // read dof
    if(!json.contains("supercells"))
      return;
    const jsonParser &json_scel = json["supercells"];
    if(!json_scel.contains(get_supercell().get_name()))
      return;
    if(!json_scel[get_supercell().get_name()].contains(get_id()))
      return;
    const jsonParser &json_config = json_scel[get_supercell().get_name()][get_id()];

    read_dof(json_config);


    // read properties: does not attempt to read in new calculation data
    if(!json_config.contains(calc_string))
      return;
    const jsonParser &json_calc = json_config[calc_string];
    if(!json_calc.contains(ref_string))
      return;
    const jsonParser &json_ref = json_calc[ref_string];
    if(!json_ref.contains("properties"))
      return;
    const jsonParser &json_prop = json_ref["properties"];

    read_properties(json_prop);

    //std::cout << "finish Configuration::read()" << std::endl;
  }

  //*********************************************************************************

  /// Read degree of freedom, source, and selected info
  ///
  ///   read dof, source, selected:  (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["dof"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["source"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["selected"]
  ///
  /// Tries to read dof from json["dof"]
  /// Tries to read source from json["source"]
  /// Tries to read selected from json["selected"]
  ///
  void Configuration::read_dof(const jsonParser &json) {

    /// json["dof"]: contains degree of freedom information
    if(!json.contains("dof")) {
      _invalidate_id();
      set_selected(false);
      return;
    }
    else {

      json.get_if(m_source, "source");
      json.get_else(m_selected, "selected", false);
      from_json(m_configdof, json["dof"]);
    }
  }

  //*********************************************************************************
  /// Read configuration properties
  /// - this does not automatically read new externally calculated properties
  ///
  ///   read properties: (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["generated"]
  ///
  /// Tries to read calculated from json["calc"]
  /// Tries to read generated from json["gen"]
  ///
  void Configuration::read_properties(const jsonParser &json) {
    if(json.contains("calc")) {
      from_json(calculated, json["calc"]);
    }

    if(json.contains("gen")) {
      from_json(generated, json["gen"]);
    }
  }

  //*********************************************************************************
  fs::path Configuration::get_pos_path() const {
    return get_primclex().dir().POS(name());
  }

  //*********************************************************************************
  fs::path Configuration::calc_dir() const {
    return get_primclex().dir().configuration_calc_dir(name(), get_primclex().settings().default_clex().calctype);
  }

  //*********************************************************************************
  fs::path Configuration::calc_properties_path() const {
    return get_primclex().dir().calculated_properties(name(), get_primclex().settings().default_clex().calctype);
  }

  //*********************************************************************************
  fs::path Configuration::calc_status_path() const {
    return get_primclex().dir().calc_status(name(), get_primclex().settings().default_clex().calctype);
  }

  //*********************************************************************************
  /// Write config.json file containing degree of freedom info
  ///
  ///   writes dof: (absolute path in config_list)
  ///     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["dof"]
  ///
  ///   adds: json["dof"]
  ///
  jsonParser &Configuration::write_dof(jsonParser &json) const {

    if(!json["dof"].is_obj()) {
      json["dof"].put_obj();
    }

    jsonParser &dof = json["dof"];

    if(occupation().size() == 0) {
      dof.erase("occupation");
    }
    else {
      dof["occupation"] = occupation();
    }

    if(displacement().size() == 0) {
      dof.erase("displacement");
    }
    else {
      dof["displacement"] = displacement();
    }

    if(!has_deformation()) {
      dof.erase("deformation");
    }
    else {
      dof["deformation"] = deformation();
    }

    return json;

  }

  //*********************************************************************************
  /// Write config.json file containing degree of freedom info
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID, adds: source
  ///
  jsonParser &Configuration::write_source(jsonParser &json) const {

    json["source"] = m_source;

    return json;

  }

  //*********************************************************************************
  /// Write POS file containing Structure
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID, adds: pos
  ///   If the configuration is completely vacant, json["pos"] = null
  ///
  jsonParser &Configuration::write_pos(jsonParser &json) const {

    // print POS to stringstream
    if(occupation() != get_supercell().vacant()) {
      std::stringstream ss;
      VaspIO::PrintPOSCAR p(*this);
      p.sort();
      p.print(ss);

      json["pos"] = ss.str();
    }
    else {
      json["pos"].put_null();
    }

    return json;

  }

  //*********************************************************************************
  /// Write param_composition.json file containing correlations
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID/CURR_CLEX/CURR_REF, adds: param_composition
  ///
  jsonParser &Configuration::write_param_composition(jsonParser &json) const {

    if(!get_primclex().has_composition_axes()) {
      json.erase("param_comp_formula");
      json.erase("param_composition");
      return json;
    }

    json["param_comp_formula"] = get_primclex().composition_axes().mol_formula();
    json["param_composition"] = get_param_composition();

    return json;

  }

  //*********************************************************************************
  /// Write properties.calc.json file containing calculated properties and param_composition
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID/CURR_CLEX/CURR_REF/properties, adds: calc
  ///
  jsonParser &Configuration::write_properties(jsonParser &json) const {

    if(calculated.size() == 0) {
      json.erase("calc");
    }
    else {
      json["calc"] = calculated;
    }

    if(generated.size() == 0) {
      json.erase("gen");
    }
    else {
      json["gen"] = generated;
    }

    return json;

  }

  //*********************************************************************************
  //--------------------------------------------------------------------------------------------------
  //Structure Factor
  Eigen::VectorXd Configuration::get_struct_fact_intensities() const {
    Eigen::VectorXd automatic_intensities(get_prim().get_struc_molecule().size());
    for(int i = 0; i < get_prim().get_struc_molecule().size(); i++)
      automatic_intensities(i) = i;
    return get_struct_fact_intensities(automatic_intensities);
  }

  Eigen::VectorXd Configuration::get_struct_fact_intensities(const Eigen::VectorXd &component_intensities) const {
    auto convert = get_index_converter(get_prim(), get_prim().get_struc_molecule());
    Eigen::VectorXd intensities(size());
    for(int i = 0; i < size(); i++) {
      intensities(i) = component_intensities(convert[get_b(i)][occ(i)]);
    }
    return intensities;
  }

  ///  Calculates the sublattice structure factors as:
  ///       intensities.segment<volume()>(i*volume()) * fourier_matrix = Q.column(i)
  ///       Q is arranged as: [Q1(k1) Q2(k1) ... Qn(k1)]
  ///                         [Q1(k2) Q2(k2) ... Qn(k2)]
  ///                         ...
  ///                         [Q1(kn) Q2(kn) ... Qn(kn)]
  ///  Q is called sublat_sf in the code
  void Configuration::calc_sublat_struct_fact(const Eigen::VectorXd &intensities) {
    //std::cout<<"Intensities"<<std::endl<<intensities<<std::endl;
    Eigen::MatrixXcd sublat_sf(supercell->basis_size(), supercell->fourier_matrix().cols());
    if(supercell->fourier_matrix().rows() == 0 || supercell->fourier_matrix().cols() == 0) {
      std::cerr << "ERROR in Configuration::calc_sublat_struct_fact. Did you "
                << "forget to initialize a fourier matrix in Supercell?"
                << " Quitting" << std::endl;
      exit(666);
    }
    //int num_kpoints = supercell->fourier_matrix().cols();
    int supercell_volume = supercell->volume();
    for(int i = 0; i < supercell->basis_size(); i++) {
      Eigen::VectorXd int_segment = intensities.segment(i * supercell_volume, supercell_volume);
      // std::cout<<"Intensity segment:"<<int_segment.transpose()<<std::endl;
      // std::cout<<"Fourier:"<<std::endl<<supercell->fourier_matrix()<<std::endl;
      sublat_sf.row(i) = int_segment.transpose() * supercell->fourier_matrix();
    }
    sublat_sf = sublat_sf / double(supercell->volume());
    // std::cout<<"Sublattice struct factor:"<<std::endl;
    // std::cout<<sublat_sf.transpose()<<std::endl;
    generated["sublat_struct_fact"] =  sublat_sf.transpose();
    prop_updated = true;
  }

  /// Structure factors are then calculated as S:
  /// S = (Q * m_phase_factor).diagonal().absolute_value()
  /// S is arranged as: [S(k1) S(k2) ... S(kn)]
  /// In the code: Q is called sublat_sf
  /// However, it would be useful to have a matrix that contained the coordinates of the k-points
  /// along with the intensities at those points. The matrix that is stored in generated is thus
  /// formatted as:
  ///   [k_x  k_y  k_z  S(k)]
  void Configuration::calc_struct_fact(const Eigen::VectorXd &intensities) {
    if(supercell->phase_factor().rows() == 0 || supercell->phase_factor().cols() == 0) {
      std::cerr << "ERROR in Configuration::calc_struct_fact. Did you "
                << "forget to initialize a phase-factor matrix in Supercell?"
                << " Quitting" << std::endl;
      exit(666);
    }
    calc_sublat_struct_fact(intensities);
    Eigen::MatrixXcd sublat_sf = generated["sublat_struct_fact"].get<Eigen::MatrixXcd>();
    //    std::cout<<"Multiplication matrix:"<<std::endl<<sublat_sf*supercell->phase_factor()<<std::endl;
    Eigen::VectorXcd raw_amplitudes = (sublat_sf * supercell->phase_factor()).diagonal();
    // std::cout.precision(4);
    // std::cout<<std::fixed;
    // std::cout<<"Raw amplitudes:"<<std::endl<<raw_amplitudes.transpose()<<std::endl;
    Eigen::MatrixXd sf_coords(raw_amplitudes.size(), 4);
    sf_coords.leftCols(3) = supercell->k_mesh();
    sf_coords.col(3) = raw_amplitudes.cwiseAbs().transpose() / double(supercell->basis_size());
    generated["struct_fact"] = sf_coords;
    prop_updated = true;
  }

  void Configuration::calc_struct_fact() {
    calc_struct_fact(get_struct_fact_intensities());
  }

  void Configuration::calc_sublat_struct_fact() {
    calc_sublat_struct_fact(get_struct_fact_intensities());
  }

  Eigen::MatrixXd Configuration::struct_fact() {
    if(!generated.contains("struct_fact"))
      calc_struct_fact();
    return generated["struct_fact"].get<Eigen::MatrixXd>();
  }

  Eigen::MatrixXcd Configuration::sublat_struct_fact() {
    if(!generated.contains("sublat_struct_fact"))
      calc_sublat_struct_fact();
    return generated["sublat_struct_fact"].get<Eigen::MatrixXcd>();
  }

  bool Configuration::is_equivalent(const Configuration &B) const {
    return this->canonical_form() == B.canonical_form();
  }

  bool Configuration::operator<(const Configuration &B) const {
    if(get_supercell() != B.get_supercell()) {
      return get_supercell() < B.get_supercell();
    }
    ConfigCompare f(*this, crystallography_tol());
    return f(B);
  }

  std::pair<std::string, Index> Configuration::split_name(std::string configname) {
    std::vector<std::string> splt_vec;
    boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);
    Index config_ind;
    if(splt_vec.size() != 2) {
      default_err_log().error("Parsing configuration name");
      default_err_log() << "configuration '" << configname << "' not valid." << std::endl;
      throw std::invalid_argument("Error in Configuration::split_name(const std::string &configname) const: Not valid");
    }

    try {
      config_ind = boost::lexical_cast<Index>(splt_vec[1]);
    }
    catch(boost::bad_lexical_cast &e) {
      default_err_log().error("Invalid config index");
      default_err_log() << "CRITICAL ERROR: In PrimClex::configuration(), malformed input:" << configname << "\n";
      throw e;
    }
    return std::make_pair(splt_vec[0], config_ind);
  }

  bool Configuration::_eq(const Configuration &B) const {
    if(get_supercell() != B.get_supercell()) {
      return false;
    }
    ConfigIsEquivalent f(*this, crystallography_tol());
    return f(B);
  }

  Configuration &apply(const PermuteIterator &it, Configuration &config) {
    apply(it, config.configdof());
    return config;
  }

  /// \brief Returns the sub-configuration that fills a particular Supercell
  ///
  /// \param sub_scel The Supercell of the sub-configuration
  /// \param super_config The super-configuration
  /// \param origin The UnitCell indicating the which unit cell in the
  ///        super-configuration is the origin in sub-configuration
  ///
  /// - Copies DoF from the super-configuration directly into the sub-configuration
  ///
  Configuration sub_configuration(
    Supercell &sub_scel,
    const Configuration &super_config,
    const UnitCell &origin) {

    //std::cout << "begin sub_configuration" << std::endl;
    if(&sub_scel.get_primclex() != &super_config.get_primclex()) {
      throw std::runtime_error(std::string("Error in 'sub_configuration:"
                                           " PrimClex of sub-Supercell and super-configuration are not the same"));
    }

    //std::cout << " here 1" << std::endl;
    Configuration sub_config {sub_scel};

    //std::cout << " here 2" << std::endl;
    // copy global dof
    if(super_config.has_deformation()) {
      sub_config.configdof().set_deformation(super_config.deformation());
    }


    //std::cout << " here 3" << std::endl;
    // initialize site dof
    if(super_config.has_occupation()) {
      sub_config.configdof().set_occupation(Array<int>(sub_config.size(), 0));
    }
    if(super_config.has_displacement()) {
      sub_config.configdof().set_displacement(
        ConfigDoF::displacement_matrix_t::Zero(3, sub_config.size()));
    }

    //std::cout << " here 4" << std::endl;
    // copy site dof
    for(Index i = 0; i < sub_config.size(); i++) {

      // unitcell of site i in sub_config
      UnitCellCoord unitcellcoord = sub_config.get_uccoord(i);

      // equivalent site in superconfig
      Index site_index = super_config.get_supercell().find(unitcellcoord + origin);

      // copy dof from superconfig to this:

      // occupation
      sub_config.configdof().occ(i) = super_config.occ(site_index);

      // displacement
      if(super_config.has_displacement()) {
        sub_config.configdof().disp(i) = super_config.disp(site_index);
      }

    }

    //std::cout << "finish sub_configuration" << std::endl;
    return sub_config;
  }


  namespace {
    std::vector<std::string> split(std::string s, char delim) {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep(&delim);
      tokenizer tok(s, sep);
      return std::vector<std::string>(tok.begin(), tok.end());
    }
  }

  /// \brief Make Configuration from name string
  ///
  /// Expects one of the following formats:
  /// - `$CANON_SCELNAME/$CANON_INDEX`
  ///   - For canonical forms, whether primitive or not
  ///   - Must exist already in config list
  /// - `$PRIM_SCELNAME/$PRIM_CANON_INDEX.equiv.$FG_PERM.$TRANS_PERM`
  ///   - For primitive, but non-canonical forms
  ///   - Primitive canonical form must exist already in config list
  ///   - Applies PermuteIterator(FG_PERM, TRANS_PERM) to primitive canonical configuration
  /// - `$CANON_SCELNAME.$PRIM_FG_OP1/super.$PRIM_FG_OP2.$PRIM_SCELNAME/$PRIM_CANON_INDEX.equiv.$FG_PERM.$TRANS_PERM`
  ///   - For non-primitive non-canonical forms
  ///   - Primitive canonical form must exist already in config list
  ///   - Supercell SCELNAME must exist already in supercell list
  ///   - Applies PermuteIterator(FG_PERM, TRANS_PERM) to primitive canonical configuration
  ///   - Then applies prim Structure factor group op with index PRIM_FG_OP and
  ///     fills the supercell SCELNAME
  ///
  Configuration make_configuration(PrimClex &primclex, std::string name) {

    // if $CANON_SCELNAME.$PRIM_FG_OP1/super.$PRIM_FG_OP2.$PRIMSCELNAME/$PRIM_CANON_INDEX.equiv.$FG_PERM.$TRANS_PERM
    auto pos = name.find("super");
    if(name.find("super") != std::string::npos) {
      std::string format = "$CANON_SCELNAME.$PRIM_FG_OP1/super.$PRIM_FG_OP2."
                           "$PRIM_SCELNAME/$PRIM_CANON_INDEX"
                           ".equiv.$FG_PERM.$TRANS_PERM";

      std::vector<std::string> tokens = split(name, '.');
      if(tokens.size() != 7) {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "tokens: " << tokens << std::endl;
        throw std::invalid_argument("Error in make_configuration: configuration name format error");
      }

      // prim equiv name
      Configuration prim_equiv = make_configuration(
                                   primclex,
                                   name.substr(pos + std::string("super").size() + 1));

      std::string scelname = name.substr(0, pos - 1);
      Index fg_op_index = boost::lexical_cast<Index>(tokens[1]);
      const auto &sym_op = primclex.get_prim().factor_group()[fg_op_index];

      if(sym_op.index() != fg_op_index) {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "read fg_op_index: " << fg_op_index << std::endl;
        primclex.err_log() << "primclex.get_prim().factor_group()[fg_op_index].index(): "
                           << primclex.get_prim().factor_group()[fg_op_index].index()
                           << std::endl << std::endl;
        throw std::runtime_error("Error in make_configuration: PRIM_FG_OP index mismatch");
      }

      FillSupercell f(
        primclex.get_supercell(scelname),
        sym_op);

      return f(prim_equiv);
    }

    // if $PRIM_SCELNAME/$PRIM_CANON_INDEX.equiv.$FG_PERM.$TRANS_PERM
    pos = name.find("equiv");
    if(pos != std::string::npos) {

      std::string format = "$PRIM_SCELNAME/$PRIM_CANON_INDEX"
                           ".equiv.$FG_PERM.$TRANS_PERM";

      //split $PRIM_SCELNAME/$PRIM_CANON_INDEX & $FG_PERM & $TRANS_PERM
      std::vector<std::string> tokens = split(name, '.');
      std::string primname = tokens[0];
      if(tokens.size() != 4) {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "tokens: " << tokens << std::endl;
        throw std::invalid_argument("Error in make_configuration: configuration name format error");
      }

      Configuration pconfig = primclex.configuration(primname);
      Index fg_index = boost::lexical_cast<Index>(tokens[2]);
      Index trans_index = boost::lexical_cast<Index>(tokens[3]);

      return apply(pconfig.get_supercell().permute_it(fg_index, trans_index), pconfig);
    }

    // if $CANON_SCELNAME/$CANON_INDEX
    return primclex.configuration(name);
  }

  /// \brief Returns correlations using 'clexulator'.
  Correlation correlations(const Configuration &config, Clexulator &clexulator) {
    return correlations(config.configdof(), config.get_supercell(), clexulator);
  }

  /// Returns parametric composition, as calculated using PrimClex::param_comp
  Eigen::VectorXd comp(const Configuration &config) {
    return config.get_param_composition();
  }

  /// \brief Returns the composition, as number of each species per unit cell
  Eigen::VectorXd comp_n(const Configuration &config) {
    return config.get_num_each_component();
  }

  /// \brief Returns the vacancy composition, as number per unit cell
  double n_vacancy(const Configuration &config) {
    if(config.get_primclex().vacancy_allowed()) {
      return comp_n(config)[config.get_primclex().vacancy_index()];
    }
    return 0.0;
  }

  /// \brief Returns the total number species per unit cell
  ///
  /// Equivalent to \code comp_n(config).sum() - n_vacancy(config) \endcode
  double n_species(const Configuration &config) {
    return comp_n(config).sum() - n_vacancy(config);
  }

  /// \brief Returns the composition as species fraction, with [Va] = 0.0, in the order of Structure::get_struc_molecule
  ///
  /// - Currently, this is really a Molecule fraction
  Eigen::VectorXd species_frac(const Configuration &config) {
    Eigen::VectorXd v = comp_n(config);
    if(config.get_primclex().vacancy_allowed()) {
      v(config.get_primclex().vacancy_index()) = 0.0;
    }
    return v / v.sum();
  }

  /// \brief Returns the composition as site fraction, in the order of Structure::get_struc_molecule
  Eigen::VectorXd site_frac(const Configuration &config) {
    return comp_n(config) / config.get_prim().basis.size();
  }

  /// \brief Returns the relaxed energy, normalized per unit cell
  double relaxed_energy(const Configuration &config) {
    return config.calc_properties()["relaxed_energy"].get<double>();
  }

  /// \brief Returns the relaxed energy, normalized per species
  double relaxed_energy_per_species(const Configuration &config) {
    return relaxed_energy(config) / n_species(config);
  }

  /// \brief Returns the reference energy, normalized per unit cell
  double reference_energy(const Configuration &config) {
    return reference_energy_per_species(config) * n_species(config);
  }

  /// \brief Returns the reference energy, normalized per species
  ///
  /// - Currently, this is per Molecule
  double reference_energy_per_species(const Configuration &config) {
    return config.get_primclex().chemical_reference()(config);
  }

  /// \brief Returns the formation energy, normalized per unit cell
  double formation_energy(const Configuration &config) {
    return relaxed_energy(config) - reference_energy(config);
  }

  /// \brief Returns the formation energy, normalized per species
  ///
  /// - Currently, this is really a Molecule fraction
  double formation_energy_per_species(const Configuration &config) {
    return formation_energy(config) / n_species(config);
  }

  /// \brief Returns the formation energy, normalized per unit cell
  double clex_formation_energy(const Configuration &config) {
    const auto &primclex = config.get_primclex();
    auto formation_energy = primclex.settings().clex("formation_energy");
    Clexulator clexulator = primclex.clexulator(formation_energy);
    const ECIContainer &eci = primclex.eci(formation_energy);

    if(eci.index().back() >= clexulator.corr_size()) {
      Log &err_log = default_err_log();
      err_log.error<Log::standard>("bset and eci mismatch");
      err_log << "using cluster expansion: 'formation_energy'" << std::endl;
      err_log << "basis set size: " << clexulator.corr_size() << std::endl;
      err_log << "max eci index: " << eci.index().back() << std::endl;
      throw std::runtime_error("Error: bset and eci mismatch");
    }

    return eci * correlations(config, clexulator);
  }

  /// \brief Returns the formation energy, normalized per unit cell
  double clex_formation_energy_per_species(const Configuration &config) {
    return clex_formation_energy(config) / n_species(config);
  }

  /// \brief Return true if all current properties have been been calculated for the configuration
  bool is_calculated(const Configuration &config) {
    const auto &set = config.get_primclex().settings();
    return std::all_of(set.properties().begin(),
                       set.properties().end(),
    [&](const std::string & key) {
      return config.calc_properties().contains(key);
    });
  }

  /// \brief Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)
  double rms_force(const Configuration &_config) {
    return _config.calc_properties()["rms_force"].get<double>();
  }

  /// \brief Cost function that describes the degree to which basis sites have relaxed
  double basis_deformation(const Configuration &_config) {
    return _config.calc_properties()["basis_deformation"].get<double>();
  }

  /// \brief Cost function that describes the degree to which lattice has relaxed
  double lattice_deformation(const Configuration &_config) {
    return _config.calc_properties()["lattice_deformation"].get<double>();
  }

  /// \brief Change in volume due to relaxation, expressed as the ratio V/V_0
  double volume_relaxation(const Configuration &_config) {
    return _config.calc_properties()["volume_relaxation"].get<double>();
  }

  /// \brief Returns the relaxed magnetic moment, normalized per unit cell
  double relaxed_magmom(const Configuration &_config) {
    return _config.calc_properties()["relaxed_magmom"].get<double>();
  }

  /// \brief Returns the relaxed magnetic moment, normalized per species
  double relaxed_magmom_per_species(const Configuration &_config) {
    return relaxed_magmom(_config) / n_species(_config);
  }

  /// \brief Returns the relaxed magnetic moment at each basis site
  Eigen::VectorXd relaxed_mag_basis(const Configuration &_config) {
    return _config.calc_properties()["relaxed_mag_basis"].get<Eigen::VectorXd>();
  }

  /// \brief Returns the relaxed magnetic moment for each molecule
  Eigen::VectorXd relaxed_mag(const Configuration &_config) {
    return _config.calc_properties()["relaxed_mag"].get<Eigen::VectorXd>();
  }

  /// \brief returns true if _config describes primitive cell of the configuration it describes
  bool is_primitive(const Configuration &_config) {
    return _config.is_primitive();
  }

  /// \brief returns true if _config no symmetry transformation applied to _config will increase its lexicographic order
  bool is_canonical(const Configuration &_config) {
    return _config.is_canonical();
  }

  bool has_relaxed_energy(const Configuration &_config) {
    return _config.calc_properties().contains("relaxed_energy");
  }

  bool has_reference_energy(const Configuration &_config) {
    return _config.get_primclex().has_composition_axes() &&
           _config.get_primclex().has_chemical_reference();
  }

  bool has_formation_energy(const Configuration &_config) {
    return has_relaxed_energy(_config) && has_reference_energy(_config);
  }

  bool has_rms_force(const Configuration &_config) {
    return _config.calc_properties().contains("rms_force");
  }

  bool has_basis_deformation(const Configuration &_config) {
    return _config.calc_properties().contains("basis_deformation");
  }

  bool has_lattice_deformation(const Configuration &_config) {
    return _config.calc_properties().contains("lattice_deformation");
  }

  bool has_volume_relaxation(const Configuration &_config) {
    return _config.calc_properties().contains("volume_relaxation");
  }

  bool has_relaxed_magmom(const Configuration &_config) {
    return _config.calc_properties().contains("relaxed_magmom");
  }

  bool has_relaxed_mag_basis(const Configuration &_config) {
    return _config.calc_properties().contains("relaxed_mag_basis");
  }

  /// \brief Constructor
  ///
  /// \param _scel Supercell to be filled
  /// \param _op SymOp that transforms the input motif before tiling into the
  ///        Supercell that is filled
  FillSupercell::FillSupercell(Supercell &_scel, const SymOp &_op) :
    m_scel(&_scel), m_op(&_op), m_motif_scel(nullptr) {}

  /// \brief Constructor
  ///
  /// \param _scel Supercell to be filled
  /// \param _motif Find the first SymOp that after application to _motif enables
  ///               tiling into _scel
  /// \param _tol tolerance
  ///
  FillSupercell::FillSupercell(Supercell &_scel, const Configuration &_motif, double _tol) :
    m_scel(&_scel), m_op(find_symop(_motif, _tol)), m_motif_scel(nullptr) {}

  Configuration FillSupercell::operator()(const Configuration &motif) const {

    if(&motif.get_supercell() != m_motif_scel) {
      _init(motif.get_supercell());
    }

    Configuration result(*m_scel);

    // ------- global dof ----------
    if(motif.has_deformation()) {
      result.set_deformation(m_op->matrix()*motif.deformation()*m_op->matrix().transpose());
    }

    // ------- site dof ----------
    Array<int> tscel_occ;
    ConfigDoF::displacement_matrix_t tscel_disp, motif_new_disp;

    // apply fg op
    if(motif.has_occupation()) {
      result.set_occupation(Array<int>(m_scel->num_sites(), 0));
    }
    //std::cout << "has_disp: " << motif.has_displacement() << std::endl;
    if(motif.has_displacement()) {
      result.init_displacement();
      //std::cout << "disp: " << motif.displacement() << std::endl;
      //std::cout << "mat: \n" << m_op->matrix() << std::endl;
      motif_new_disp = m_op->matrix() * motif.displacement();
      //std::cout << "new_disp: " << motif_new_disp << std::endl;

    }

    // copy transformed dof, as many times as necessary to fill the supercell
    for(Index s = 0; s < m_index_table.size(); ++s) {
      for(Index i = 0; i < m_index_table[s].size(); ++i) {
        Index scel_s = m_index_table[s][i];

        if(motif.has_occupation()) {
          result.configdof().occ(scel_s) = motif.occ(s);
        }
        if(motif.has_displacement()) {
          result.configdof().disp(scel_s) = motif_new_disp.col(s);
        }
      }
    }
    if(motif.has_displacement()) {
      //std::cout << "final disp: " << result.displacement() << std::endl;
    }
    return result;
  }

  /// \brief Find first SymOp in the prim factor group such that apply(op, motif)
  ///        can be used to fill the Supercell
  const SymOp *FillSupercell::find_symop(const Configuration &motif, double tol) {

    const Lattice &motif_lat = motif.get_supercell().get_real_super_lattice();
    const Lattice &scel_lat = m_scel->get_real_super_lattice();
    auto begin = m_scel->get_primclex().get_prim().factor_group().begin();
    auto end = m_scel->get_primclex().get_prim().factor_group().end();

    auto res = is_supercell(scel_lat, motif_lat, begin, end, tol);
    if(res.first == end) {

      std::cerr << "Requested supercell transformation matrix: \n"
                << m_scel->get_transf_mat() << "\n";
      std::cerr << "Requested motif Configuration: "
                << motif.name() << "\n";
      std::cerr << "Configuration transformation matrix: \n"
                << motif.get_supercell().get_transf_mat() << "\n";

      throw std::runtime_error(
        "Error in 'FillSupercell::find_symop':\n"
        "  The motif cannot be tiled onto the specified supercell."
      );
    }

    return &(*res.first);
  }

  void FillSupercell::_init(Supercell &_motif_scel) const {

    m_motif_scel = &_motif_scel;

    // ------- site dof ----------
    Lattice oriented_motif_lat = copy_apply(*m_op, m_motif_scel->get_real_super_lattice());

    // Create a PrimGrid linking the prim and the oriented motif each to the supercell
    // So we can tile the decoration of the motif config onto the supercell correctly
    PrimGrid prim_grid(oriented_motif_lat, m_scel->get_real_super_lattice());

    //std::cout << "m_op->matrix(): \n" << m_op->matrix() << std::endl;
    //std::cout << "m_op->tau(): \n" << m_op->tau() << std::endl;

    const Structure &prim = m_scel->get_prim();
    m_index_table.resize(m_motif_scel->num_sites());

    // for each site in motif
    for(Index s = 0 ; s < m_motif_scel->num_sites() ; s++) {

      //std::cout << "before: " << m_motif_scel->uccoord(s) << std::endl;

      // apply symmetry to re-orient and find unit cell coord
      UnitCellCoord oriented_uccoord = copy_apply(*m_op, m_motif_scel->uccoord(s), prim);
      //std::cout << "after: " << oriented_uccoord << std::endl;

      // for each unit cell of the oriented motif in the supercell, copy the occupation
      for(Index i = 0 ; i < prim_grid.size() ; i++) {

        Index prim_motif_tile_ind = m_scel->prim_grid().find(prim_grid.coord(i, PRIM));

        UnitCellCoord mc_uccoord =  m_scel->prim_grid().uccoord(prim_motif_tile_ind) + oriented_uccoord.unitcell();
        // b-index when doing UnitCellCoord addition is ambiguous; explicitly set it
        mc_uccoord.sublat() = oriented_uccoord.sublat();

        m_index_table[s].push_back(m_scel->find(mc_uccoord));
      }
    }

  }

  std::ostream &operator<<(std::ostream &sout, const Configuration &c) {

    sout << c.name() << "\n";
    if(c.has_deformation()) {
      sout << "Deformation:\n" << c.deformation() << std::endl;
    }
    for(Index i = 0; i < c.size(); ++i) {
      sout << "Linear index: " << i << "  UnitCellCoord: " << c.get_uccoord(i) << std::endl;
      if(c.has_occupation()) {
        sout << "  Occupation: " << c.occ(i) << "  (" << c.get_mol(i).name << ")\n";
      }
      if(c.has_displacement()) {
        sout << "  Displacement: " << c.disp(i).transpose() << "\n";
      }
    }
    return sout;
  }

}


