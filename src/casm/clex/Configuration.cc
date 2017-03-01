#include "casm/clex/Configuration.hh"

#include <sstream>

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/clex/ConfigCompare.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/app/QueryHandler_impl.hh"

namespace CASM {

  const std::string QueryTraits<Configuration>::name = "Configuration";

  template class QueryHandler<Configuration>;

  namespace {
    typedef std::insert_iterator<std::map<std::string, std::shared_ptr<RuntimeLibrary> > > runtimelib_it_type;
    typedef std::insert_iterator<DataFormatterDictionary<Configuration> > config_dict_it_type;

    std::vector<std::string> split(std::string s, char delim) {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep(&delim);
      tokenizer tok(s, sep);
      return std::vector<std::string>(tok.begin(), tok.end());
    }

    fs::path _calc_properties_path(const Configuration &config) {
      const PrimClex &primclex = config.primclex();
      return primclex.dir().calculated_properties(config.name(), primclex.settings().default_clex().calctype);
    }

    fs::path _calc_status_path(const Configuration &config) {
      const PrimClex &primclex = config.primclex();
      return primclex.dir().calc_status(config.name(), primclex.settings().default_clex().calctype);
    }

  }



  template std::pair<config_dict_it_type, runtimelib_it_type> load_query_plugins(
    const ProjectSettings &set,
    config_dict_it_type dict_it,
    runtimelib_it_type lib_it);



  /// Construct a default Configuration
  Configuration::Configuration(
    const Supercell &_supercell,
    const jsonParser &src,
    const ConfigDoF &_configdof) :
    m_supercell(&_supercell),
    m_source_updated(false),
    m_configdof(_configdof),
    m_dof_deps_updated(false),
    m_prop_updated(false) {

    set_source(src);
  }

  /// Construct a Configuration from JSON data
  Configuration::Configuration(
    const Supercell &_supercell,
    const std::string &_id,
    const jsonParser &_data) :
    m_source_updated(false),
    m_dof_deps_updated(false),
    m_prop_updated(false) {

    this->from_json(_data, _supercell, _id);

  }

  /// Construct a Configuration from JSON data
  Configuration::Configuration(const PrimClex &_primclex,
                               const std::string &_configname,
                               const jsonParser &_data) :
    Configuration(*_primclex.db<Supercell>().find(Configuration::split_name(_configname).first),
                  Configuration::split_name(_configname).second,
                  _data) {}

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
    m_source_updated = true;
    m_dof_deps_updated = true;
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

      m_source_updated = true;
      m_dof_deps_updated = true;
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

          m_source_updated = true;
          m_dof_deps_updated = true;
        }
      }
    }
  }

  //*********************************************************************************
  void Configuration::clear() {
    _modify_dof();
    m_configdof.clear();
  }

  //*********************************************************************************
  void Configuration::init_occupation() {
    _modify_dof();
    set_occupation(std::vector<int>(this->size(), 0));
  }

  //*********************************************************************************
  void Configuration::set_occupation(const std::vector<int> &new_occupation) {
    _modify_dof();
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
    _modify_dof();
    m_configdof.occ(site_l) = val;
  }

  //*********************************************************************************
  void Configuration::clear_occupation() {
    _modify_dof();
    m_configdof.clear_occupation();
  }

  //*********************************************************************************

  void Configuration::init_specie_id() {
    _modify_dof();
    m_configdof.specie_id().resize(this->size());
    for(Index i = 0; i < this->size(); ++i) {
      m_configdof.specie_id()[i].resize(mol(i).size(), 0);
    }
  }

  //*********************************************************************************

  std::vector<std::vector<Index> > &Configuration::specie_id() {
    _modify_dof();
    return m_configdof.specie_id();
  }

  //*********************************************************************************

  const std::vector<std::vector<Index> > &Configuration::specie_id() const {
    return m_configdof.specie_id();
  }

  //*********************************************************************************

  std::vector<Index> &Configuration::specie_id(Index site_l) {
    _modify_dof();
    return m_configdof.specie_id()[site_l];
  }

  //*********************************************************************************

  const std::vector<Index> &Configuration::specie_id(Index site_l) const {
    return m_configdof.specie_id()[site_l];
  }

  //*********************************************************************************

  void Configuration::clear_specie_id() {
    _modify_dof();
    m_configdof.clear_specie_id();
  }

  //*********************************************************************************
  void Configuration::init_displacement() {
    _modify_dof();
    set_displacement(displacement_matrix_t::Zero(3, this->size()));
  }

  //*********************************************************************************
  void Configuration::set_displacement(const displacement_matrix_t &new_displacement) {
    _modify_dof();
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
    _modify_dof();
    m_configdof.disp(site_l) = _disp;
  }

  //*********************************************************************************
  void Configuration::clear_displacement() {
    _modify_dof();
    m_configdof.clear_displacement();
  }

  //*********************************************************************************
  void Configuration::init_deformation() {
    _modify_dof();
    set_deformation(Eigen::Matrix3d::Identity());
  }

  //*********************************************************************************

  void Configuration::set_deformation(const Eigen::Matrix3d &new_deformation) {
    _modify_dof();
    m_configdof.set_deformation(new_deformation);
  }

  //*********************************************************************************

  void Configuration::clear_deformation() {
    _modify_dof();
    m_configdof.clear_deformation();
  }

  //*******************************************************************************

  /// \brief Check if this is a primitive Configuration
  bool Configuration::is_primitive() const {
    if(!cache().contains("is_primitive")) {
      bool result = (find_translation() == supercell().translate_end());
      cache_insert("is_primitive", result);
      return result;
    }
    return cache()["is_primitive"].get<bool>();
  }

  //*******************************************************************************

  /// \brief Returns a PermuteIterator corresponding to the first non-zero pure
  /// translation that maps the Configuration onto itself.
  ///
  /// - If primitive, returns this->supercell().translate_end()
  PermuteIterator Configuration::find_translation() const {
    ConfigIsEquivalent f(*this, crystallography_tol());
    const Supercell &scel = supercell();
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
    std::cout << "T: \n" << tconfig.supercell().transf_mat() << std::endl;
    std::cout << "L: \n" << tconfig.supercell().real_super_lattice().lat_column_mat() << std::endl;
    */

    std::unique_ptr<Supercell> next_scel;

    // check if config is primitive, and if not, obtain a translation that maps the config on itself
    while(true) {

      PermuteIterator result = tconfig.find_translation();
      if(result == tconfig.supercell().translate_end()) {
        break;
      }

      // replace one of the lattice vectors with the translation
      Lattice new_lat = replace_vector(
                          tconfig.ideal_lattice(),
                          result.sym_op().tau(),
                          crystallography_tol()).make_right_handed().reduced_cell();

      next_scel.reset(new Supercell(&primclex(), new_lat));
      /*
      std::cout << "T: \n" << next_scel->transf_mat() << std::endl;
      std::cout << "L: \n" << next_scel->real_super_lattice().lat_column_mat() << std::endl;
      */

      // create a sub configuration in the new supercell
      tconfig = sub_configuration(*next_scel, tconfig);
      //std::cout << "sub occ: \n" << tconfig.occupation() << std::endl;

      tconfig.m_supercell_ptr.reset(next_scel.release());
      tconfig.m_supercell = tconfig.m_supercell_ptr.get();

    }

    return tconfig;
  }

  //*******************************************************************************

  /// \brief Check if Configuration is in the canonical form
  bool Configuration::is_canonical() const {
    if(!cache().contains("to_canonical")) {
      const Supercell &scel = supercell();
      ConfigIsEquivalent f(*this, crystallography_tol());
      bool result = std::all_of(
                      ++scel.permute_begin(),
                      scel.permute_end(),
      [&](const PermuteIterator & p) {
        return f(p) || !f.is_less();
      });

      cache_insert("is_canonical", result);
      return result;
    }
    return cache()["is_calculated"].get<bool>();
  }

  //*******************************************************************************

  /// \brief Returns the operation that applied to *this returns the canonical form
  PermuteIterator Configuration::to_canonical() const {
    if(!cache().contains("to_canonical")) {
      ConfigCompare f(*this, crystallography_tol());
      const Supercell &scel = supercell();
      auto result = std::max_element(scel.permute_begin(), scel.permute_end(), f);
      cache_insert("to_canonical", result);
      return result;
    }
    else {
      return cache()["to_canonical"].get<PermuteIterator>(supercell());
    }
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
  /// - Will be a Supercell included in the PrimClex.supercell_list()
  Configuration Configuration::in_canonical_supercell() const {

    const Supercell &canon_scel = supercell().canonical_form();

    FillSupercell f(canon_scel, *this, crystallography_tol());
    Configuration in_canon = f(*this);

    // only OK to use if canon_scel and this->supercell() are stored in
    //   primclex supercell list
    //Configuration in_canon = primclex().fill_supercell(canon_scel, *this);

    return in_canon.canonical_form();
  }

  //*******************************************************************************

  /// \brief Insert this configuration (in primitive & canonical form) in the database
  ///
  /// \param primitive_only If true, only the primitive Configuration is inserted.
  ///
  /// - By convention, the primitive canonical form of a configuration must
  ///   always be saved in the config list.
  /// - If this is already known to be primitive & canonical, prefer to use
  ///   PrimClex::db<Configuration>.insert(config) directly.
  ConfigInsertResult Configuration::insert(bool primitive_only) const {

    ConfigInsertResult res;

    Configuration pconfig = this->primitive().in_canonical_supercell();
    std::tie(res.primitive_it, res.insert_primitive) =
      primclex().db<Configuration>().insert(pconfig);

    // if the primitive supercell is the same as the equivalent canonical supercell
    if(supercell().canonical_form() == pconfig.supercell()) {
      res.insert_canonical = res.insert_primitive;
      res.canonical_it = res.primitive_it;
    }
    else {
      if(primitive_only) {
        res.insert_canonical = false;
      }
      else {
        std::tie(res.canonical_it, res.insert_canonical) =
          primclex().db<Configuration>().insert(this->in_canonical_supercell());
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
    const Supercell &scel = supercell();
    std::copy_if(scel.permute_begin(), scel.permute_end(), std::back_inserter(fg), f);

    int mult = this->prim().factor_group().size() / fg.size();
    cache_insert("multiplicity", mult);

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
    cache_insert("point_group_name", sym_group.get_name());
    return sym_group;
  }

  //*******************************************************************************

  /// \brief Returns the point group that leaves the Configuration unchanged
  std::string Configuration::point_group_name() const  {
    if(!cache().contains("point_group_name")) {
      this->point_group();
    }
    return cache()["point_group_name"].get<std::string>();
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  Configuration Configuration::fill_supercell(Supercell &scel, const SymOp &op) const {
    FillSupercell f(scel, op);
    return f(*this);

    // only OK to use if both supercells are stored in primclex supercell list:
    //return primclex().fill_supercell(scel, *this, op);
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  ///
  /// - Uses the first symop in g that such that scel is a supercell of op*(*this)
  Configuration Configuration::fill_supercell(Supercell &scel, const SymGroup &g) const {

    auto res = is_supercell(
                 scel.real_super_lattice(),
                 ideal_lattice(),
                 g.begin(),
                 g.end(),
                 crystallography_tol());

    if(res.first == g.end()) {

      std::cerr << "Requested supercell transformation matrix: \n"
                << scel.transf_mat() << "\n";
      std::cerr << "Requested motif Configuration: " <<
                name() << "\n";
      std::cerr << "Configuration transformation matrix: \n"
                << supercell().transf_mat() << "\n";

      throw std::runtime_error(
        "Error in 'Configuration::fill_supercell(const Supercell &scel, const SymGroup& g)'\n"
        "  The motif cannot be tiled onto the specified supercell."
      );
    }

    return fill_supercell(scel, *res.first);
  }

  //*********************************************************************************
  void Configuration::set_calc_properties(const jsonParser &calc) {
    m_calculated = calc;
    m_prop_updated = true;
    m_dof_deps_updated = true;
  }

  //*********************************************************************************

  bool Configuration::read_calc_properties(jsonParser &parsed_props) const {
    //std::cout << "begin Configuration::read_calculated()" << std::endl;
    bool success = true;
    /// properties.calc.json: contains calculated properties
    ///   For default clex calctype only
    fs::path filepath = _calc_properties_path(*this);
    //std::cout << "filepath: " << filepath << std::endl;
    parsed_props = jsonParser();
    if(fs::exists(filepath)) {
      jsonParser json(filepath);

      //Record file timestamp
      parsed_props["data_timestamp"] = fs::last_write_time(filepath);

      std::vector<std::string> props = primclex().settings().properties();
      for(Index i = 0; i < props.size(); i++) {
        //std::cout << "checking for: " << props[i] << std::endl;
        if(json.contains(props[i])) {

          // normal by #prim cells for some properties
          if(props[i] == "energy" || props[i] == "relaxed_energy") {
            parsed_props[ props[i] ] = json[props[i]].get<double>() / supercell().volume();
          }
          else {
            parsed_props[props[i]] = json[props[i]];
          }
        }
        else
          success = false;
      }
      //Get RMS force:
      if(json.contains("relaxed_forces")) {
        if(json["relaxed_forces"].size()) {
          Eigen::MatrixXd forces;
          json["relaxed_forces"].get(forces);
          parsed_props["rms_force"] = sqrt((forces.transpose() * forces).trace() / double(forces.rows()));
        }
        else {
          parsed_props["rms_force"] = 0.;
        }
      }
    }
    else
      success = false;

    return success;
  }

  //********** ACCESSORS ***********

  const Lattice &Configuration::ideal_lattice()const {
    return supercell().real_super_lattice();
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
  std::string Configuration::generate_name() const {

    std::string result;
    // canonical forms in canonical supercells
    if(id() != "none") {
      result = supercell().name() + "/" + id();
    }
    else if(supercell().is_canonical() && is_canonical()) {
      result = supercell().name() + "/" + id();
    }
    else {
      result = supercell().name() + "/non_canonical_equivalent";
    }
    m_dof_deps_updated = true;
    return result;
  }

  //*********************************************************************************
  const jsonParser &Configuration::source() const {
    return m_source;
  }

  //*********************************************************************************
  ///Returns number of sites, NOT the number of primitives that fit in here
  Index Configuration::size() const {
    return supercell().num_sites();
  }

  //*********************************************************************************
  const Structure &Configuration::prim() const {
    return supercell().prim();
  }

  //*********************************************************************************
  //PrimClex &Configuration::primclex() {
  //return supercell().primclex();
  //}

  //*********************************************************************************
  const PrimClex &Configuration::primclex() const {
    return supercell().primclex();
  }

  //*********************************************************************************
  const Supercell &Configuration::supercell() const {
    return *m_supercell;
  }

  //*********************************************************************************
  double Configuration::crystallography_tol() const {
    return primclex().settings().crystallography_tol();
  }

  //*********************************************************************************
  UnitCellCoord Configuration::uccoord(Index site_l) const {
    return supercell().uccoord(site_l);
  }

  //*********************************************************************************
  Index Configuration::linear_index(const UnitCellCoord &bijk) const {
    return supercell().linear_index(bijk);
  }

  //*********************************************************************************
  int Configuration::sublat(Index site_l) const {
    return supercell().sublat(site_l);
  }

  //*********************************************************************************
  const Molecule &Configuration::mol(Index site_l) const {
    return prim().basis[ sublat(site_l) ].site_occupant()[ occ(site_l) ];
  }

  //*********************************************************************************
  const Properties &Configuration::calc_properties() const {
    return m_calculated;
  }

  //*********************************************************************************
  /// \brief Get symmetric multiplicity (i.e., size of configuration's factor_group)
  int Configuration::multiplicity() const {
    if(!cache().contains("multiplicity")) {
      this->factor_group();
    }
    return cache()["multiplicity"].get<int>();
  }

  //*********************************************************************************

  /// Returns composition on each sublattice: sublat_comp[ prim basis site / sublattice][ molecule_type]
  ///   molucule_type is ordered as in the Prim structure's site_occupant list for that basis site (includes vacancies)
  std::vector<Eigen::VectorXd> Configuration::sublattice_composition() const {

    // get the number of each molecule
    auto _sublat_num_each_molecule = sublat_num_each_molecule();
    std::vector<Eigen::VectorXd> sublattice_composition(_sublat_num_each_molecule.size());

    // divide by number of sites per sublattice ( supercell volume )
    for(Index i = 0; i < _sublat_num_each_molecule.size(); i++) {
      sublattice_composition[i] = Eigen::VectorXd::Zero(_sublat_num_each_molecule[i].size());
      for(Index j = 0; j < _sublat_num_each_molecule[i].size(); j++) {
        sublattice_composition[i][j] = (1.0 * _sublat_num_each_molecule[i][j]) / supercell().volume();
      }
    }

    return sublattice_composition;
  }

  //*********************************************************************************
  /// Returns number of each molecule by sublattice:
  ///   sublat_num_each_molecule[ prim basis site / sublattice ][ molecule_type]
  ///   molucule_type is ordered as in the Prim structure's site_occupant list for that basis site
  std::vector<Eigen::VectorXi> Configuration::sublat_num_each_molecule() const {

    Index i;

    // [basis_site][site_occupant_index]
    auto convert = index_converter(prim(), prim().struc_molecule());

    // create an array to count the number of each molecule
    std::vector<Eigen::VectorXi> sublat_num_each_molecule;
    for(i = 0; i < prim().basis.size(); i++) {
      sublat_num_each_molecule.push_back(Eigen::VectorXi::Zero(prim().basis[i].site_occupant().size()));
    }

    // count the number of each molecule by sublattice
    for(i = 0; i < size(); i++) {
      sublat_num_each_molecule[ sublat(i) ][occ(i)]++;
    }

    return sublat_num_each_molecule;
  }

  //*********************************************************************************
  /// Returns composition, not counting vacancies
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's struc_molecule(), with [Va]=0.0
  Eigen::VectorXd Configuration::composition() const {

    // get the number of each molecule type
    Eigen::VectorXi _num_each_molecule = num_each_molecule();

    /// get the total number of non-vacancy atoms
    int num_atoms = 0;

    // need to know which molecules are vacancies
    auto struc_molecule = prim().struc_molecule();

    Index i;
    for(i = 0; i < struc_molecule.size(); i++) {
      if(struc_molecule[i].is_vacancy()) {
        // set to zero, so the Va concentration is reported as 0.0
        _num_each_molecule[i] = 0;
      }
      num_atoms += _num_each_molecule[i];
    }

    // calculate the comp (not including vacancies) from the number of each molecule
    return _num_each_molecule.cast<double>() / double(num_atoms);

  }

  //*********************************************************************************
  /// Returns composition, including vacancies
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's struc_molecule()
  Eigen::VectorXd Configuration::true_composition() const {
    return num_each_molecule().cast<double>() / size();
  }

  //*********************************************************************************
  /// Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::struc_molecule()
  Eigen::VectorXi Configuration::num_each_molecule() const {
    return CASM::num_each_molecule(m_configdof, supercell());
  }

  //*********************************************************************************
  /// Returns parametric composition, as calculated using PrimClex::param_comp
  Eigen::VectorXd Configuration::param_composition() const {
    if(!primclex().has_composition_axes()) {
      std::cerr << "Error in Configuration::param_composition()" << std::endl;
      std::cerr << "  Composition axes are not set." << std::endl;
      exit(1);
    }

    return primclex().composition_axes().param_composition(num_each_component());
  }

  //*********************************************************************************
  /// Returns num_each_component[ component_type] per prim cell,
  ///   where 'component_type' is ordered as ParamComposition::components
  Eigen::VectorXd Configuration::num_each_component() const {

    // component order used for param_composition
    std::vector<std::string> components = primclex().composition_axes().components();

    // initialize
    Eigen::VectorXd num_each_component = Eigen::VectorXd::Zero(components.size());

    // [basis_site][site_occupant_index]
    auto convert = index_converter(prim(), components);

    // count the number of each component
    for(Index i = 0; i < size(); i++) {
      num_each_component[ convert[ sublat(i) ][occ(i)] ] += 1.0;
    }

    // normalize per prim cell
    for(Index i = 0; i < components.size(); i++) {
      num_each_component[i] /= supercell().volume();
    }

    return num_each_component;
  }



  //********* IO ************


  /// Writes the Configuration to a json object (the config list)
  ///   Uses PrimClex's current default settings to write the appropriate properties
  ///
  ///   "json" corresponds (entire file)["supercells"][scelname][id]
  ///   "configname" corresponds to "scelname/id"
  ///
  ///   write_dof, source: (absolute path in config_list)
  ///     json["dof"]
  ///     json["source"]
  ///
  ///   write_properties: (absolute path in config_list)
  ///     json["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]
  ///
  jsonParser &Configuration::to_json(jsonParser &json) const {

    //std::cout << "begin Configuration::to_json(jsonParser& json)" << std::endl;

    json.put_obj();

    jsonParser &dof = json["dof"];
    if(occupation().size()) {
      dof["occupation"] = occupation();
    }
    if(displacement().size()) {
      dof["displacement"] = displacement();
    }
    if(has_deformation()) {
      dof["deformation"] = deformation();
    }

    json["source"] = m_source;

    if(cache_updated()) {
      json["cache"] = cache();
    }

    if(m_prop_updated) {
      const ProjectSettings &set = primclex().settings();
      std::string calc_string = "calctype." + set.default_clex().calctype;
      std::string ref_string = "ref." + set.default_clex().ref;

      jsonParser &json_ref = json[calc_string][ref_string];
      jsonParser &json_prop = json_ref["properties"];

      if(m_calculated.size() == 0) {
        json.erase("calc");
      }
      else {
        json["calc"] = m_calculated;
      }
    }

    //std::cout << "finish Configuration::to_json(jsonParser& json)" << std::endl;
    return json;
  }

  //*********************************************************************************

  std::ostream &Configuration::write_pos(std::ostream &sout) const {
    VaspIO::PrintPOSCAR p(*this);
    p.sort();
    p.print(sout);
    return sout;
  }

  //*********************************************************************************

  void Configuration::write_pos() const {

    const auto &dir = primclex().dir();
    try {
      fs::create_directories(dir.configuration_dir(name()));
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in Configuration::write_pos()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream file(dir.POS(name()));
    write_pos(file);
    return;
  }

  //*********************************************************************************

  /// Private members:

  /// Reads the Configuration from JSON
  ///   Uses PrimClex's current default settings to read in the appropriate properties
  ///
  /// This is private, because it is only called from the constructor:
  ///   Configuration(const jsonParser& json, const PrimClex& _primclex, const std::string &configname)
  ///
  ///   "json" corresponds (entire file)["supercells"][scelname][id]
  ///   "configname" corresponds to "scelname/id"
  ///
  ///   read dof, source: (absolute path in config_list)
  ///     json["dof"]
  ///     json["source"]
  ///
  ///   read properties: (absolute path in config_list)
  ///     json["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]
  ///
  void Configuration::from_json(const jsonParser &json, const Supercell &scel, std::string _id) {

    //std::cout << "begin  Configuration::from_json()" << std::endl;

    m_supercell = &scel;

    this->clear_name();
    this->set_id(_id);

    json.get_if(m_source, "source");
    m_source_updated = false;

    CASM::from_json(m_configdof, json["dof"]);
    m_dof_deps_updated = false;

    CASM::from_json(cache(), json["cache"]);

    // read properties from 'json' input only: does not attempt to read in new
    // calculation data from the calc.properties.json file
    // - use read_calc_properties() to read new calc.properties.json files
    m_prop_updated = false;

    const ProjectSettings &set = primclex().settings();
    std::string calc_string = "calctype." + set.default_clex().calctype;
    auto calc_it = json.find(calc_string);
    if(calc_it == json.end()) {
      return;
    }

    std::string ref_string = "ref." + set.default_clex().ref;
    auto ref_it = calc_it->find(ref_string);
    if(ref_it == calc_it->end()) {
      return;
    }

    auto prop_it = ref_it->find("properties");
    if(prop_it == ref_it->end()) {
      return;
    }

    prop_it->get_if(m_calculated, "calc");

    //std::cout << "finish Configuration::from_json()" << std::endl;
  }

  void Configuration::from_json(const jsonParser &json, const PrimClex &primclex, std::string _configname) {
    auto name = Configuration::split_name(_configname);
    this->from_json(json, *primclex.db<Supercell>().find(name.first), name.second);
  }

  bool Configuration::is_equivalent(const Configuration &B) const {
    return this->canonical_form() == B.canonical_form();
  }

  bool Configuration::operator<(const Configuration &B) const {
    if(supercell() != B.supercell()) {
      return supercell() < B.supercell();
    }
    ConfigCompare f(*this, crystallography_tol());
    return f(B);
  }

  std::pair<std::string, std::string> Configuration::split_name(std::string configname) {
    std::vector<std::string> splt_vec;
    boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);
    if(splt_vec.size() != 2) {
      default_err_log().error("Parsing configuration name");
      default_err_log() << "configuration '" << configname << "' not valid." << std::endl;
      default_err_log() << "must have form: scelname/configid" << std::endl;
      throw std::invalid_argument("Error in Configuration::split_name(const std::string &configname) const: Not valid");
    }

    return std::make_pair(splt_vec[0], splt_vec[1]);
  }

  bool Configuration::_eq(const Configuration &B) const {
    if(supercell() != B.supercell()) {
      return false;
    }
    ConfigIsEquivalent f(*this, crystallography_tol());
    return f(B);
  }

  Configuration jsonConstructor<Configuration>::from_json(
    const jsonParser &json,
    const PrimClex &primclex,
    const std::string &configname) {

    return Configuration(primclex, configname, json);
  }

  Configuration jsonConstructor<Configuration>::from_json(
    const jsonParser &json,
    const Supercell &scel,
    const std::string &id) {

    return Configuration(scel, id, json);
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
    if(&sub_scel.primclex() != &super_config.primclex()) {
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
      sub_config.configdof().set_occupation(std::vector<int>(sub_config.size(), 0));
    }
    if(super_config.has_displacement()) {
      sub_config.configdof().set_displacement(
        ConfigDoF::displacement_matrix_t::Zero(3, sub_config.size()));
    }

    //std::cout << " here 4" << std::endl;
    // copy site dof
    for(Index i = 0; i < sub_config.size(); i++) {

      // unitcell of site i in sub_config
      UnitCellCoord unitcellcoord = sub_config.uccoord(i);

      // equivalent site in superconfig
      Index site_index = super_config.supercell().linear_index(unitcellcoord + origin);

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
      const auto &sym_op = primclex.prim().factor_group()[fg_op_index];

      if(sym_op.index() != fg_op_index) {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "read fg_op_index: " << fg_op_index << std::endl;
        primclex.err_log() << "primclex.prim().factor_group()[fg_op_index].index(): "
                           << primclex.prim().factor_group()[fg_op_index].index()
                           << std::endl << std::endl;
        throw std::runtime_error("Error in make_configuration: PRIM_FG_OP index mismatch");
      }

      FillSupercell f(
        *primclex.db<Supercell>().find(scelname),
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

      Configuration pconfig = *primclex.db<Configuration>().find(primname);
      Index fg_index = boost::lexical_cast<Index>(tokens[2]);
      Index trans_index = boost::lexical_cast<Index>(tokens[3]);

      return apply(pconfig.supercell().permute_it(fg_index, trans_index), pconfig);
    }

    // if $CANON_SCELNAME/$CANON_INDEX
    return *primclex.db<Configuration>().find(name);
  }

  /// \brief Returns correlations using 'clexulator'.
  Eigen::VectorXd correlations(const Configuration &config, Clexulator &clexulator) {
    return correlations(config.configdof(), config.supercell(), clexulator);
  }

  /// Returns parametric composition, as calculated using PrimClex::param_comp
  Eigen::VectorXd comp(const Configuration &config) {
    return config.param_composition();
  }

  /// \brief Returns the composition, as number of each species per unit cell
  Eigen::VectorXd comp_n(const Configuration &config) {
    return config.num_each_component();
  }

  /// \brief Returns the vacancy composition, as number per unit cell
  double n_vacancy(const Configuration &config) {
    if(config.primclex().vacancy_allowed()) {
      return comp_n(config)[config.primclex().vacancy_index()];
    }
    return 0.0;
  }

  /// \brief Returns the total number species per unit cell
  ///
  /// Equivalent to \code comp_n(config).sum() - n_vacancy(config) \endcode
  double n_species(const Configuration &config) {
    return comp_n(config).sum() - n_vacancy(config);
  }

  /// \brief Returns the composition as species fraction, with [Va] = 0.0, in the order of Structure::struc_molecule
  ///
  /// - Currently, this is really a Molecule fraction
  Eigen::VectorXd species_frac(const Configuration &config) {
    Eigen::VectorXd v = comp_n(config);
    if(config.primclex().vacancy_allowed()) {
      v(config.primclex().vacancy_index()) = 0.0;
    }
    return v / v.sum();
  }

  /// \brief Returns the composition as site fraction, in the order of Structure::struc_molecule
  Eigen::VectorXd site_frac(const Configuration &config) {
    return comp_n(config) / config.prim().basis.size();
  }

  /// \brief Status of calculation
  std::string calc_status(const Configuration &config) {
    fs::path p = _calc_status_path(config);
    if(fs::exists(p)) {
      jsonParser json(p);
      if(json.contains("status"))
        return json["status"].get<std::string>();
    }
    return("not_submitted");
  }

  // \brief Reason for calculation failure.
  std::string failure_type(const Configuration &config) {
    fs::path p = _calc_status_path(config);
    if(fs::exists(p)) {
      jsonParser json(p);
      if(json.contains("failure_type"))
        return json["failure_type"].get<std::string>();
    }
    return("none");
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
    return config.primclex().chemical_reference()(config);
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
    const auto &primclex = config.primclex();
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
    const auto &set = config.primclex().settings();
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
    return _config.primclex().has_composition_axes() &&
           _config.primclex().has_chemical_reference();
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

  /// \brief Constructor
  ///
  /// \param _scel Supercell to be filled
  /// \param _op SymOp that transforms the input motif before tiling into the
  ///        Supercell that is filled
  FillSupercell::FillSupercell(const Supercell &_scel, const SymOp &_op) :
    m_scel(&_scel), m_op(&_op), m_motif_scel(nullptr) {}

  /// \brief Constructor
  ///
  /// \param _scel Supercell to be filled
  /// \param _motif Find the first SymOp that after application to _motif enables
  ///               tiling into _scel
  /// \param _tol tolerance
  ///
  FillSupercell::FillSupercell(const Supercell &_scel, const Configuration &_motif, double _tol) :
    m_scel(&_scel), m_op(find_symop(_motif, _tol)), m_motif_scel(nullptr) {}

  Configuration FillSupercell::operator()(const Configuration &motif) const {

    if(&motif.supercell() != m_motif_scel) {
      _init(motif.supercell());
    }

    Configuration result(*m_scel);

    // ------- global dof ----------
    if(motif.has_deformation()) {
      result.set_deformation(m_op->matrix()*motif.deformation()*m_op->matrix().transpose());
    }

    // ------- site dof ----------
    std::vector<int> tscel_occ;
    ConfigDoF::displacement_matrix_t tscel_disp, motif_new_disp;

    // apply fg op
    if(motif.has_occupation()) {
      result.set_occupation(std::vector<int>(m_scel->num_sites(), 0));
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

    const Lattice &motif_lat = motif.supercell().real_super_lattice();
    const Lattice &scel_lat = m_scel->real_super_lattice();
    auto begin = m_scel->primclex().prim().factor_group().begin();
    auto end = m_scel->primclex().prim().factor_group().end();

    auto res = is_supercell(scel_lat, motif_lat, begin, end, tol);
    if(res.first == end) {

      std::cerr << "Requested supercell transformation matrix: \n"
                << m_scel->transf_mat() << "\n";
      std::cerr << "Requested motif Configuration: "
                << motif.name() << "\n";
      std::cerr << "Configuration transformation matrix: \n"
                << motif.supercell().transf_mat() << "\n";

      throw std::runtime_error(
        "Error in 'FillSupercell::find_symop':\n"
        "  The motif cannot be tiled onto the specified supercell."
      );
    }

    return &(*res.first);
  }

  void FillSupercell::_init(const Supercell &_motif_scel) const {

    m_motif_scel = &_motif_scel;

    // ------- site dof ----------
    Lattice oriented_motif_lat = copy_apply(*m_op, m_motif_scel->real_super_lattice());

    // Create a PrimGrid linking the prim and the oriented motif each to the supercell
    // So we can tile the decoration of the motif config onto the supercell correctly
    PrimGrid prim_grid(oriented_motif_lat, m_scel->real_super_lattice());

    //std::cout << "m_op->matrix(): \n" << m_op->matrix() << std::endl;
    //std::cout << "m_op->tau(): \n" << m_op->tau() << std::endl;

    const Structure &prim = m_scel->prim();
    m_index_table.resize(m_motif_scel->num_sites());

    // for each site in motif
    for(Index s = 0 ; s < m_motif_scel->num_sites() ; s++) {

      //std::cout << "before: " << m_motif_scel->uccoord(s) << std::endl;

      // apply symmetry to re-orient and find unit cell coord
      UnitCellCoord oriented_uccoord = copy_apply(*m_op, m_motif_scel->uccoord(s));
      //std::cout << "after: " << oriented_uccoord << std::endl;

      // for each unit cell of the oriented motif in the supercell, copy the occupation
      for(Index i = 0 ; i < prim_grid.size() ; i++) {

        Index prim_motif_tile_ind = m_scel->prim_grid().find(prim_grid.coord(i, PRIM));

        UnitCellCoord mc_uccoord(
          prim,
          oriented_uccoord.sublat(),
          m_scel->prim_grid().unitcell(prim_motif_tile_ind) + oriented_uccoord.unitcell());

        m_index_table[s].push_back(m_scel->linear_index(mc_uccoord));
      }
    }

  }

  std::ostream &operator<<(std::ostream &sout, const Configuration &c) {

    sout << c.name() << "\n";
    if(c.has_deformation()) {
      sout << "Deformation:\n" << c.deformation() << std::endl;
    }
    for(Index i = 0; i < c.size(); ++i) {
      sout << "Linear index: " << i << "  UnitCellCoord: " << c.uccoord(i) << std::endl;
      if(c.has_occupation()) {
        sout << "  Occupation: " << c.occ(i) << "  (" << c.mol(i).name << ")\n";
      }
      if(c.has_displacement()) {
        sout << "  Displacement: " << c.disp(i).transpose() << "\n";
      }
    }
    return sout;
  }

}


