#include "casm/clex/Configuration_impl.hh"

#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/casm_io/stream_io/container.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/database/DiffTransConfigDatabase.hh"

namespace CASM {

  template class HasPrimClex<Comparisons<Calculable<CRTPBase<Configuration> > > >;
  template class HasSupercell<Comparisons<Calculable<CRTPBase<Configuration> > > >;
  template class ConfigCanonicalForm<HasSupercell<Comparisons<Calculable<CRTPBase<Configuration> > > > >;

  namespace {
    std::vector<std::string> split(std::string s, char delim) {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep(&delim);
      tokenizer tok(s, sep);
      return std::vector<std::string>(tok.begin(), tok.end());
    }

  }


  /// Construct a default Configuration
  Configuration::Configuration(
    const Supercell &_supercell,
    const jsonParser &src,
    const ConfigDoF &_configdof) :
    m_supercell(&_supercell),
    m_configdof(_configdof) {

    set_source(src);
  }

  /// Construct a default Configuration that owns its Supercell
  Configuration::Configuration(
    const std::shared_ptr<Supercell> &_supercell,
    const jsonParser &source,
    const ConfigDoF &_configdof) :
    m_supercell(_supercell.get()),
    m_supercell_ptr(_supercell),
    m_configdof(_configdof) {

    set_source(source);
  }

  /// Construct a Configuration from JSON data
  Configuration::Configuration(
    const Supercell &_supercell,
    const std::string &_id,
    const jsonParser &_data) {
    if(_id == "none") {
      if(_data.contains("dof")) {
        ConfigDoF dof;
        CASM::from_json(dof, _data["dof"]);
        *this = Configuration(_supercell, _data, dof);
        return;
      }
      *this = Configuration(_supercell, _data);
      return;
    }
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

      // create a sub configuration in the new supercell
      tconfig = sub_configuration(*next_scel, tconfig);

      tconfig.m_supercell_ptr.reset(next_scel.release());
      tconfig.m_supercell = tconfig.m_supercell_ptr.get();

    }

    return tconfig;
  }

  //*******************************************************************************

  /// \brief Check if Configuration is an endpoint of an existing diff_trans_config
  bool Configuration::is_diff_trans_endpoint() const {
    auto it = primclex().db<Kinetics::DiffTransConfiguration>().scel_range(supercell().name()).begin();
    for(; it != primclex().db<Kinetics::DiffTransConfiguration>().scel_range(supercell().name()).end(); ++it) {
      if(is_equivalent(it->from_config()) || is_equivalent(it->to_config())) {
        return true;
      }
    }
    return false;
  }

  /// \brief tells which diff_trans this configuration is an endpoint of
  std::string Configuration::diff_trans_endpoint_of() const {
    std::set<std::string> collection;
    auto it = primclex().db<Kinetics::DiffTransConfiguration>().scel_range(supercell().name()).begin();
    for(; it != primclex().db<Kinetics::DiffTransConfiguration>().scel_range(supercell().name()).end(); ++it) {
      if(is_equivalent(it->from_config()) || is_equivalent(it->to_config())) {
        collection.insert(it->orbit_name());
      }
    }
    std::string result = "none";
    if(collection.size()) {
      result = "";
    }
    for(auto &n : collection) {
      result += n + ", ";
    }
    return result;
  }

  //*******************************************************************************

  /// \brief Returns the canonical form Configuration in the canonical Supercell
  ///
  /// - Canonical Supercell will be inserted in the PrimClex.db<Supercell>() if
  ///   necessary
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
  /// - Supercells are inserted in the Supercell database as necessary
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
    return invariant_subgroup();
  }

  //*******************************************************************************

  /// \brief Returns the subgroup of the Supercell factor group that leaves the
  ///        Configuration unchanged
  std::vector<PermuteIterator> Configuration::invariant_subgroup() const {
    std::vector<PermuteIterator> fg = ConfigurationBase::invariant_subgroup();
    int mult = this->prim().factor_group().size() / fg.size();
    cache_insert("multiplicity", mult);
    return fg;
  }

  //*******************************************************************************

  bool Configuration::is_canonical() const {
    if(!cache().contains("is_canonical")) {
      bool result = ConfigurationBase::is_canonical();
      cache_insert("is_canonical", result);
      return result;
    }
    return cache()["is_canonical"].get<bool>();
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
  Configuration Configuration::fill_supercell(const Supercell &scel, const SymOp &op) const {
    FillSupercell f(scel, op);
    return f(*this);

    // only OK to use if both supercells are stored in primclex supercell list:
    //return primclex().fill_supercell(scel, *this, op);
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  ///
  /// - Uses the first symop in g that such that scel is a supercell of op*(*this)
  Configuration Configuration::fill_supercell(const Supercell &scel, const SymGroup &g) const {

    auto res = is_supercell(
                 scel.lattice(),
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

  //********** ACCESSORS ***********

  const Lattice &Configuration::ideal_lattice()const {
    return supercell().lattice();
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
  std::string Configuration::generate_name_impl() const {

    // If 'id' is already known, just return configname
    if(id() != "none") {
      return supercell().name() + "/" + id();
    }

    const auto &db = primclex().db<Configuration>();

    // If primitive configuration in the canonical supercell:
    if(supercell().is_canonical() && is_primitive()) {

      auto prim_canon_config = canonical_form();
      auto find_it = db.search(prim_canon_config);

      // get configname (i.e. SCELV_A_B_C_D_E_F/I)
      std::string _name;
      if(find_it != db.end()) {
        return supercell().name() + "/" + id();
      }
      else {
        _name = supercell().name() + "/none";
      }

      auto op = from_canonical();
      if(op == supercell().permute_begin()) {
        // if canonical, return configname
        return _name;
      }
      else {
        // if non canonical, return configname.equiv.$FG_PERM.$TRANS_PERM
        return _name + ".equiv." + std::to_string(op.factor_group_index())
               + "." + std::to_string(op.translation_index());
      }
    }

    // If not primitive or not in the canonical supercell:

    // Find primitive canonical config:
    auto prim_canon_config = primitive().in_canonical_supercell();
    auto find_it = db.search(prim_canon_config);
    if(find_it != db.end()) {
      prim_canon_config = *find_it;
    }

    // Find the operation that takes the prim_canon_config and gives *this:
    //
    // op1: a prim_canon_config.supercell() PermuteIterator
    // op2: a prim().factor_group() SymOp
    // *this == copy_apply(op1, prim_canon_config).fill_supercell(supercell(), op2);

    auto f = equal_to();
    auto begin = prim_canon_config.supercell().permute_begin();
    auto end = prim_canon_config.supercell().permute_end();
    for(auto op1 = begin; op1 != end; ++op1) {
      auto test = copy_apply(op1, prim_canon_config);
      for(const auto &op2 : prim().factor_group()) {
        if(f(test.fill_supercell(supercell(), op2))) {
          return supercell().name() + "/super." + std::to_string(op2.index()) + "."
                 + prim_canon_config.name() + ".equiv." + std::to_string(op1.factor_group_index())
                 + "." + std::to_string(op1.translation_index());
        }
      }
    }

    throw std::runtime_error("Error in Configuration::generate_name_impl");
  }

  //*********************************************************************************
  ///Returns number of sites, NOT the number of primitives that fit in here
  Index Configuration::size() const {
    return supercell().num_sites();
  }

  //*********************************************************************************
  const Supercell &Configuration::supercell() const {
    return *m_supercell;
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
  /// \brief Get symmetric multiplicity, excluding translations
  ///
  /// - equal to prim.factor_group().size() / this->factor_group().size()
  int Configuration::multiplicity() const {
    if(!cache().contains("multiplicity")) {
      this->factor_group();
    }
    return cache()["multiplicity"].get<int>();
  }

  //*********************************************************************************

  Configuration &Configuration::apply_sym(const PermuteIterator &it) {
    configdof().apply_sym(it);
    return *this;
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
    auto convert = make_index_converter(prim(), prim().struc_molecule());

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
    auto convert = make_index_converter(prim(), components);

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
  ///
  ///   writes: dof, source, cache:
  ///     json["dof"]
  ///     json["source"]
  ///     json["cache"]
  ///
  jsonParser &Configuration::to_json(jsonParser &json) const {

    json.put_obj();

    CASM::to_json(m_configdof, json["dof"]);
    CASM::to_json(source(), json["source"]);

    json["cache"].put_obj();
    if(cache_updated()) {
      json["cache"] = cache();
    }

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

    m_supercell = &scel;

    this->clear_name();
    this->set_id(_id);

    auto source_it = json.find("source");
    if(source_it != json.end()) {
      set_source(*source_it);
    }
    CASM::from_json(m_configdof, json["dof"]);
    CASM::from_json(cache(), json["cache"]);

    // read properties from 'json' input only: does not attempt to read in new
    // calculation data from the calc.properties.json file
    // - use read_calc_properties() to read new calc.properties.json files

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

    auto calc_props_it = prop_it->find("calc");
    if(calc_props_it != prop_it->end()) {
      set_calc_properties(*calc_props_it);
    }

  }

  void Configuration::from_json(const jsonParser &json, const PrimClex &primclex, std::string _configname) {
    auto name = Configuration::split_name(_configname);
    this->from_json(json, *primclex.db<Supercell>().find(name.first), name.second);
  }

  bool Configuration::operator<(const Configuration &B) const {
    return less()(B);
  }

  ConfigCompare Configuration::less() const {
    return ConfigCompare(*this);
  }

  ConfigIsEquivalent Configuration::equal_to() const {
    return ConfigIsEquivalent(*this);
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

  bool Configuration::eq_impl(const Configuration &B) const {
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


    if(&sub_scel.primclex() != &super_config.primclex()) {
      throw std::runtime_error(std::string("Error in 'sub_configuration:"
                                           " PrimClex of sub-Supercell and super-configuration are not the same"));
    }


    Configuration sub_config {sub_scel};

    // copy global dof
    if(super_config.has_deformation()) {
      sub_config.configdof().set_deformation(super_config.deformation());
    }

    // initialize site dof
    if(super_config.has_occupation()) {
      sub_config.configdof().set_occupation(std::vector<int>(sub_config.size(), 0));
    }
    if(super_config.has_displacement()) {
      sub_config.configdof().set_displacement(
        ConfigDoF::displacement_matrix_t::Zero(3, sub_config.size()));
    }

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

  /// \brief Returns the relaxed energy, normalized per unit cell
  double relaxed_energy(const Configuration &config) {
    return config.calc_properties()["relaxed_energy"].get<double>() / config.supercell().volume();
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

  /// \brief Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)
  double rms_force(const Configuration &_config) {
    //Get RMS force:
    const jsonParser &props = _config.calc_properties();

    Eigen::MatrixXd forces;
    props["relaxed_forces"].get(forces);
    return sqrt((forces.transpose() * forces).trace() / double(forces.rows()));
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

  /// \brief returns true if no symmetry transformation applied to _config will increase its lexicographic order
  bool is_canonical(const Configuration &_config) {
    return _config.is_canonical();
  }

  /// \brief returns true if _config is an endpoint of an existing diff_trans_config in the database
  bool is_diff_trans_endpoint(const Configuration &_config) {
    return _config.is_diff_trans_endpoint();
  }

  /// \brief returns which diff_trans _config is an endpoint of
  std::string diff_trans_endpoint_of(const Configuration &_config) {
    return _config.diff_trans_endpoint_of();
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
    const jsonParser &props = _config.calc_properties();
    auto it = props.find("relaxed_forces");
    return it != props.end() && it->size();
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
    if(motif.has_displacement()) {
      result.init_displacement();

      motif_new_disp = m_op->matrix() * motif.displacement();

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
    return result;
  }

  /// \brief Find first SymOp in the prim factor group such that apply(op, motif)
  ///        can be used to fill the Supercell
  const SymOp *FillSupercell::find_symop(const Configuration &motif, double tol) {

    const Lattice &motif_lat = motif.supercell().lattice();
    const Lattice &scel_lat = m_scel->lattice();
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
    Lattice oriented_motif_lat = copy_apply(*m_op, m_motif_scel->lattice());

    // Create a PrimGrid linking the prim and the oriented motif each to the supercell
    // So we can tile the decoration of the motif config onto the supercell correctly
    PrimGrid prim_grid(oriented_motif_lat, m_scel->lattice());

    const Structure &prim = m_scel->prim();
    m_index_table.resize(m_motif_scel->num_sites());

    // for each site in motif
    for(Index s = 0 ; s < m_motif_scel->num_sites() ; s++) {

      // apply symmetry to re-orient and find unit cell coord
      UnitCellCoord oriented_uccoord = copy_apply(*m_op, m_motif_scel->uccoord(s));

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
        sout << "  Occupation: " << c.occ(i) << "  (" << c.mol(i).name() << ")\n";
      }
      if(c.has_displacement()) {
        sout << "  Displacement: " << c.disp(i).transpose() << "\n";
      }
    }
    return sout;
  }

}


