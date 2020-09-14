#include "casm/clex/Configuration_impl.hh"

#include <memory>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ClexParamPack.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/MappedPropertiesTools.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/Named_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

  template class HasPrimClex<Comparisons<Calculable<CRTPBase<Configuration> > > >;
  template class HasSupercell<Comparisons<Calculable<CRTPBase<Configuration> > > >;
  template class ConfigCanonicalForm<HasSupercell<Comparisons<Calculable<CRTPBase<Configuration> > > > >;

  namespace DB {
    template class Indexed<CRTPBase<Configuration> >;
    template class Named<CRTPBase<Configuration> >;
  }

  Configuration::Configuration(std::shared_ptr<Supercell> const &_supercell_ptr):
    m_supercell(_supercell_ptr.get()),
    m_supercell_ptr(_supercell_ptr),
    m_configdof(_supercell_ptr->zero_configdof(_supercell_ptr->prim().lattice().tol())) {}

  Configuration::Configuration(std::shared_ptr<Supercell> const &_supercell_ptr,
                               ConfigDoF const &_dof):
    m_supercell(_supercell_ptr.get()),
    m_supercell_ptr(_supercell_ptr),
    m_configdof(_dof) {}

  /// Construct a default Configuration
  Configuration::Configuration(
    const Supercell &_supercell,
    const jsonParser &src) :
    m_supercell(&_supercell),
    m_configdof(Configuration::zeros(_supercell).configdof()) {

    set_source(src);
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
    const std::shared_ptr<Supercell> &_supercell_ptr,
    const jsonParser &source) :
    m_supercell(_supercell_ptr.get()),
    m_supercell_ptr(_supercell_ptr),
    m_configdof(Configuration::zeros(*_supercell_ptr).configdof()) {

    set_source(source);
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
    const jsonParser &_data) :
    m_configdof(Configuration::zeros(_supercell).configdof()) {
    if(_id == "none") {
      if(_data.contains("dof")) {
        *this = Configuration(_supercell, _data, _data["dof"].get<ConfigDoF>(_supercell.prim()));
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

  Configuration Configuration::zeros(Supercell const &_scel) {
    return Configuration(_scel, jsonParser(), _scel.zero_configdof(_scel.prim().lattice().tol()));
  }

  //*********************************************************************************

  Configuration Configuration::zeros(Supercell const &_scel, double _tol) {
    return Configuration(_scel, jsonParser(), _scel.zero_configdof(_tol));
  }

  //*********************************************************************************

  Configuration Configuration::zeros(const std::shared_ptr<Supercell> &_supercell_ptr) {
    return Configuration {_supercell_ptr};
  }

  //*********************************************************************************

  Configuration Configuration::zeros(const std::shared_ptr<Supercell> &_supercell_ptr, double _tol) {
    return Configuration {_supercell_ptr, _supercell_ptr->zero_configdof(_tol)};
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

  void Configuration::set_occ(Index site_l, int val) {
    _modify_dof();
    m_configdof.occ(site_l) = val;
  }

  //*******************************************************************************

  /// \brief Check if this is a primitive Configuration
  bool Configuration::is_primitive() const {
    if(!cache().contains("is_primitive")) {
      bool result = (find_translation() == supercell().sym_info().translate_end());
      cache_insert("is_primitive", result);
      return result;
    }
    return cache()["is_primitive"].get<bool>();
  }

  //*******************************************************************************

  /// \brief Returns a PermuteIterator corresponding to the first non-zero pure
  /// translation that maps the Configuration onto itself.
  ///
  /// - If primitive, returns this->supercell().sym_info().translate_end()
  PermuteIterator Configuration::find_translation() const {
    ConfigIsEquivalent f(*this, crystallography_tol());
    const Supercell &scel = supercell();
    auto begin = scel.sym_info().translate_begin();
    auto end = scel.sym_info().translate_end();
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

    // check if config is primitive, and if not, obtain a translation that maps the config on itself
    while(true) {
      PermuteIterator result = tconfig.find_translation();

      if(result == tconfig.supercell().sym_info().translate_end()) {
        break;
      }

      // replace one of the lattice vectors with the translation
      Lattice new_lat = replace_vector(
                          tconfig.ideal_lattice(),
                          result.sym_op().tau(),
                          crystallography_tol()).make_right_handed().reduced_cell();

      // create a sub configuration in the new supercell
      tconfig = sub_configuration(std::make_shared<Supercell>(&primclex(), new_lat), tconfig);
    }

    return tconfig;
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
    if(fg.size() == 0) {
      default_err_log() << "Something went very wrong in invariant_subgroup returning group size 0" << std::endl;
    }
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
  std::vector<PermuteIterator> Configuration::point_group() const {
    std::vector<PermuteIterator> config_factor_group = factor_group();
    std::vector<PermuteIterator> result;
    for(int i = 0; i < config_factor_group.size(); i++) {
      if(i == 0 || config_factor_group[i].factor_group_index() != config_factor_group[i - 1].factor_group_index())
        result.push_back(config_factor_group[i]);
    }

    return result;
  }

  //*******************************************************************************

  /// \brief Returns the point group that leaves the Configuration unchanged
  std::string Configuration::point_group_name() const  {
    if(!cache().contains("point_group_name")) {
      cache_insert("point_group_name", make_sym_group(this->point_group(), this->supercell().sym_info().supercell_lattice()).get_name());
    }
    return cache()["point_group_name"].get<std::string>();
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with configuration
  Configuration Configuration::fill_supercell(const Supercell &scel) const {
    FillSupercell f(scel, prim().factor_group()[0]);
    return f(*this);
  }


  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  Configuration Configuration::fill_supercell(const Supercell &scel, const SymOp &op) const {
    FillSupercell f(scel, op);
    return f(*this);
  }

  //*******************************************************************************

  /// \brief Fills supercell 'scel' with reoriented configuration, op*(*this)
  ///
  /// - Uses the first symop in g that such that scel is a supercell of op*(*this)
  Configuration Configuration::fill_supercell(const Supercell &scel, const SymGroup &g) const {

    auto res = xtal::is_equivalent_superlattice(
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

  /// \brief Get operations that transform canonical primitive to this
  RefToCanonicalPrim::RefToCanonicalPrim(const Configuration &_config) :
    config(_config),
    // Find primitive canonical config:
    prim_canon_config(config.primitive().in_canonical_supercell()) {

    // Find the transformation of prim_canon_config and gives config:
    //
    // op1: a prim_canon_config.supercell() PermuteIterator
    // op2: a prim().point_group() SymOp
    // *this == copy_apply(op1, prim_canon_config).fill_supercell(supercell(), op2);

    // first find op2 == *res.first
    auto res = xtal::is_equivalent_superlattice(
                 config.ideal_lattice(),
                 prim_canon_config.ideal_lattice(),
                 config.prim().point_group().begin(),
                 config.prim().point_group().end(),
                 config.crystallography_tol());
    this->from_canonical_lat = *res.first;
    this->transf_mat = iround(res.second);

    // given op2, find op1
    auto f = config.equal_to();
    auto begin = prim_canon_config.supercell().sym_info().permute_begin();
    auto end = prim_canon_config.supercell().sym_info().permute_end();
    for(auto op1 = begin; op1 != end; ++op1) {
      auto test = copy_apply(op1, this->prim_canon_config);
      if(f(test.fill_supercell(config.supercell(), *res.first))) {
        this->from_canonical_config = op1;
        return;
      }
    }

    throw std::runtime_error("Error in RefToCanonicalPrim: could not find solution");
  }

  std::string RefToCanonicalPrim::name() const {
    return config.supercell().name() + "/super."
           + std::to_string(from_canonical_lat.index()) + "."
           + prim_canon_config.name() + ".equiv."
           + std::to_string(from_canonical_config.factor_group_index())
           + "." + std::to_string(from_canonical_config.translation_index());
  }

  //*********************************************************************************
  /// \brief Returns a Configuration name
  ///
  /// For configurations in supercells equivalent to the canonical supercell:
  ///   For canonical configurations:
  ///   - CANON_CONFIG_NAME = `$CANON_SCELNAME/$CONFIG_INDEX`
  ///   - The CANON_CONFIG is found in the config database by name.
  ///   For non-canonical configurations:
  ///   - NONCANON_CONFIG_NAME = `$CANON_CONFIG_NAME.equiv.$FG_PERM.$TRANS_PERM`
  ///   - The CANON_CONFIG is found in the config database,
  ///     then the FG_PERM-th factor_group permutation is applied,
  ///     follwed by TRANS_PERM-th translation permutation.
  /// For all other configurations:
  ///   - NONEQUIV_SCEL_CONFIG_NAME = `$SCEL_NAME/super.$PRIM_FG_OP2.`$NONCANON_CONFIG_NAME`
  ///   - SCEL_NAME may be for a canonical equivalent or non canonical equivalent supercell, in which case
  ///     it is represented by `CANON_SCEL_NAME.$PRIM_FG_OP1`
  ///   - The NONCANON_CONFIG is constructed, PRIM_FG_OP2 is applied, and then the SCEL is filled
  ///   - When generating the NONCANON_CONFIG_NAME, the primitive configuration is used
  ///
  std::string Configuration::generate_name_impl() const {
    /// if 'id' is known
    if(id() != "none") {
      return supercell().name() + "/" + id();
    }

    const auto &db = primclex().db<Configuration>();
    const Lattice &canon_scel_lat = supercell().canonical_form().lattice();
    bool is_canon_equiv_lat = xtal::is_equivalent(this->supercell().lattice(), canon_scel_lat);

    // If in the canonical equivalent supercell lattice:
    if(is_canon_equiv_lat) {

      // put Configuration in canonical supercell
      Configuration canon_scel_config = fill_supercell(supercell().canonical_form());

      // if canonical
      if(canon_scel_config.is_canonical()) {
        auto find_it = db.search(canon_scel_config);

        // if already in database
        if(find_it != db.end()) {
          // set id
          set_id(find_it->id());
          return supercell().name() + "/" + id();
        }
        // if not in database
        else {
          return supercell().name() + "/none";
        }
      }
      // if non-canonical
      else {
        // get canonical form and 'from_canonical' op
        Configuration canon_config = canon_scel_config.canonical_form();
        auto op = canon_scel_config.from_canonical();

        // get canonical form id if already in database, else 'none'
        auto find_it = db.search(canon_config);
        std::string canon_config_id = "none";
        if(find_it != db.end()) {
          canon_config_id = find_it->id();
        }

        // construct name
        return supercell().name() + "/" + canon_config_id
               + ".equiv." + std::to_string(op.factor_group_index())
               + "." + std::to_string(op.translation_index());
      }
    }

    RefToCanonicalPrim ref(*this);
    return ref.name();
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
    return prim().basis()[ sublat(site_l) ].occupant_dof()[ occ(site_l) ];
  }

  //*********************************************************************************
  /// \brief Get symmetric multiplicity, excluding translations
  ///
  /// - equal to prim.factor_group().size() / this->factor_group().size()
  int Configuration::multiplicity() const {
    if(!cache().contains("multiplicity")) {
      int result = this->prim().factor_group().size() / this->factor_group().size();
      cache_insert("multiplicity", result);
      return result;
    }
    return cache()["multiplicity"].get<int>();
  }

  //*********************************************************************************

  Configuration &Configuration::apply_sym(const PermuteIterator &op) {
    auto all_props = calc_properties_map();
    configdof().apply_sym(op);
    for(auto it = all_props.begin(); it != all_props.end(); ++it) {
      std::string calctype = it->first;
      set_calc_properties(copy_apply(op, it->second), calctype);
    }
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

    // create an array to count the number of each molecule
    std::vector<Eigen::VectorXi> sublat_num_each_molecule;
    for(i = 0; i < prim().basis().size(); i++) {
      sublat_num_each_molecule.push_back(Eigen::VectorXi::Zero(prim().basis()[i].occupant_dof().size()));
    }

    // count the number of each molecule by sublattice
    for(i = 0; i < size(); i++) {
      sublat_num_each_molecule[ sublat(i) ][occ(i)]++;
    }

    return sublat_num_each_molecule;
  }

  //*********************************************************************************
  /// Returns composition, not counting vacancies
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's xtal::struc_molecule_name(), with [Va]=0.0
  Eigen::VectorXd Configuration::composition() const {

    // get the number of each molecule type
    Eigen::VectorXi _num_each_molecule = num_each_molecule();

    /// get the total number of non-vacancy atoms
    int num_atoms = 0;

    // need to know which molecules are vacancies
    auto struc_mol = xtal::struc_molecule_name(prim());

    Index i;
    for(i = 0; i < struc_mol.size(); i++) {
      if(xtal::is_vacancy(struc_mol[i])) {
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
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's xtal::struc_molecule_name()
  Eigen::VectorXd Configuration::true_composition() const {
    return num_each_molecule().cast<double>() / size();
  }

  //*********************************************************************************
  /// Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::xtal::struc_molecule_name()
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
    CASM::from_json(m_configdof, json["dof"]);//, primclex().n_basis());
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
      set_calc_properties(calc_props_it->get<MappedProperties>(), set.default_clex().calctype);
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

  //*********************************************************************************

  std::string pos_string(Configuration const  &_config) {
    std::stringstream ss;
    VaspIO::PrintPOSCAR p(make_simple_structure(_config), _config.name());
    p.sort();
    p.print(ss);
    return ss.str();
  }

  //*********************************************************************************

  void write_pos(Configuration const &_config) {

    try {
      fs::create_directories(_config.primclex().dir().configuration_dir(_config.name()));
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in Configuration::write_pos()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream(_config.primclex().dir().POS(_config.name()))
        << pos_string(_config);
    return;
  }


  //*********************************************************************************

  std::string config_json_string(Configuration const  &_config) {
    std::stringstream ss;

    jsonParser tjson;
    to_json(make_simple_structure(_config), tjson);
    tjson.print(ss);
    return ss.str();
  }

  //*********************************************************************************

  void write_config_json(Configuration const &_config) {

    try {
      fs::create_directories(_config.primclex().dir().configuration_dir(_config.name()));
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in Configuration::write_pos()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream(_config.primclex().dir().config_json(_config.name()))
        << config_json_string(_config);
    return;
  }


  Configuration sub_configuration(
    std::shared_ptr<Supercell> sub_scel_ptr,
    const Configuration &super_config,
    const UnitCell &origin) {


    if(&sub_scel_ptr->primclex() != &super_config.primclex()) {
      throw std::runtime_error(std::string("Error in 'sub_configuration:"
                                           " PrimClex of sub-Supercell and super-configuration are not the same"));
    }

    Configuration sub_config {sub_scel_ptr};
    // copy global dof
    for(auto const &dof : super_config.configdof().global_dofs()) {
      sub_config.configdof().set_global_dof(dof.first, dof.second.values());
    }

    // copy site dof
    for(Index i = 0; i < sub_config.size(); i++) {

      // unitcell of site i in sub_config
      UnitCellCoord unitcellcoord = sub_config.uccoord(i);

      // equivalent site in superconfig
      Index site_index = super_config.supercell().linear_index(unitcellcoord + origin);

      // copy dof from superconfig to this:

      // occupation
      if(super_config.has_occupation()) {
        sub_config.configdof().occ(i) = super_config.occ(site_index);
      }

    }

    for(auto const &dof : super_config.configdof().local_dofs()) {
      LocalContinuousConfigDoFValues &res_ref(sub_config.configdof().local_dof(dof.first));
      for(Index i = 0; i < sub_config.size(); i++) {

        // unitcell of site i in sub_config
        UnitCellCoord unitcellcoord = sub_config.uccoord(i);

        // equivalent site in superconfig
        Index scel_site_index = super_config.supercell().linear_index(unitcellcoord + origin);

        // copy dof from superconfig to this:
        res_ref.site_value(i) = dof.second.site_value(scel_site_index);
      }
    }


    return sub_config;
  }

  namespace {

    /// \brief Make non-canonical Configuration (in canonical supercell) from name string
    ///
    /// Note: canonical config must be in database
    Configuration make_non_canon_configuration(const PrimClex &primclex, std::string name) {
      std::string format = "$CANON_CONFIGNAME.equiv.$FG_PERM.$TRANS_PERM";

      //split $CANON_CONFIGNAME & $FG_PERM & $TRANS_PERM
      std::vector<std::string> tokens;
      boost::split(tokens, name, boost::is_any_of("."), boost::token_compress_on);
      std::string canon_config_name = tokens[0];
      if(tokens.size() != 4) {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "tokens: " << tokens << std::endl;
        throw std::invalid_argument("Error in make_configuration: configuration name format error");
      }

      Configuration canon_config = *primclex.db<Configuration>().find(canon_config_name);
      Index fg_index = boost::lexical_cast<Index>(tokens[2]);
      Index trans_index = boost::lexical_cast<Index>(tokens[3]);

      return CASM::apply(canon_config.supercell().sym_info().permute_it(fg_index, trans_index), canon_config);
    }

    /// \brief Make general super Configuration from name string
    ///
    /// For non-primitive configurations, or configurations with supercells that are
    ///   not equivalent to the canonical supercell:
    ///   - NONEQUIV_SCEL_CONFIG_NAME = `$SCEL_NAME/super.$PRIM_FG_OP2.`$NONCANON_CONFIG_NAME`
    ///   - SCEL_NAME may be for a canonical equivalent or non canonical equivalent supercell, in which case
    ///     it is represented by `CANON_SCEL_NAME.$PRIM_FG_OP1`
    ///   - The NONCANON_CONFIG is constructed, PRIM_FG_OP2 is applied, and then the SCEL is filled
    ///
    Configuration make_super_configuration(const PrimClex &primclex, std::string name) {

      // expected format
      std::string format = "$SCEL_NAME/super.$PRIM_FG_OP2.$NONCANON_CONFIG_NAME";

      // tokenize name
      std::vector<std::string> tokens;
      boost::split(tokens, name, boost::is_any_of("./"), boost::token_compress_on);

      std::string scelname, non_canon_config_name;
      Index fg_op_index = 0;

      if(tokens[1] != "super" && tokens[2] != "super") {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "tokens: " << tokens << std::endl;

        throw std::invalid_argument("Error in make_configuration: configuration name format error");
      }

      try {
        // parse, if canonical equivalent lattice
        if(tokens[1] == "super") {
          scelname = tokens[0];
          fg_op_index = boost::lexical_cast<Index>(tokens[2]);
          non_canon_config_name = tokens[3] + '/' + tokens[4] + '.' + tokens[5] + '.' + tokens[6] + '.' + tokens[7];
        }
        // parse, if not canonical equivalent lattice
        else { // if(tokens[2] == "super")
          scelname = tokens[0] + '.' + tokens[1];
          fg_op_index = boost::lexical_cast<Index>(tokens[3]);
          non_canon_config_name = tokens[4] + '/' + tokens[5] + '.' + tokens[6] + '.' + tokens[7] + '.' + tokens[8];
        }
      }
      catch(...) {
        primclex.err_log().error("In make_configuration");
        primclex.err_log() << "expected format: " << format << "\n";
        primclex.err_log() << "name: " << name << std::endl;
        primclex.err_log() << "tokens: " << tokens << std::endl;

        throw std::invalid_argument("Error in make_configuration: configuration name format error");
      }

      // -- Generate primitive equivalent configuration

      // prim equiv name
      Configuration non_canon_config = make_non_canon_configuration(primclex, non_canon_config_name);

      const auto &sym_op = primclex.prim().factor_group()[fg_op_index];

      // canonical equivalent supercells can be found in the database
      if(tokens[1] == "super") {
        FillSupercell f(
          *primclex.db<Supercell>().find(scelname),
          sym_op);

        return f(non_canon_config);
      }
      // non canonical equivalent supercells must be constructed
      else {
        std::shared_ptr<Supercell> scel = make_shared_supercell(primclex, scelname);
        FillSupercell f(scel, sym_op);
        return f(non_canon_config);
      }
    }

  }

  /// \brief Make Configuration from name string
  ///
  /// For configurations in supercells equivalent to the canonical supercell:
  ///   For canonical configurations:
  ///   - CANON_CONFIG_NAME = `$CANON_SCELNAME/$CONFIG_INDEX`
  ///   - The CANON_CONFIG is found in the config database by name.
  ///   For non-canonical configurations:
  ///   - NONCANON_CONFIG_NAME = `$CANON_CONFIG_NAME.equiv.$FG_PERM.$TRANS_PERM`
  ///   - The CANON_CONFIG is found in the config database,
  ///     then the FG_PERM-th factor_group permutation is applied,
  ///     follwed by TRANS_PERM-th translation permutation.
  /// For all other configurations:
  ///   - NONEQUIV_SCEL_CONFIG_NAME = `$SCEL_NAME/super.$PRIM_FG_OP2.`$NONCANON_CONFIG_NAME`
  ///   - SCEL_NAME may be for a canonical equivalent or non canonical equivalent supercell, in which case
  ///     it is represented by `CANON_SCEL_NAME.$PRIM_FG_OP1`
  ///   - The NONCANON_CONFIG is constructed, PRIM_FG_OP2 is applied, and then the SCEL is filled
  ///
  Configuration make_configuration(const PrimClex &primclex, std::string name) {

    // if most general case:
    // format = $SCEL_NAME/super.$PRIM_FG_OP2.$NONCANON_CONFIG_NAME
    if(name.find("super") != std::string::npos) {
      return make_super_configuration(primclex, name);
    }

    // if non-canonical configuration in canonical equivalent supercell:
    // format = $CANON_CONFIG_NAME.equiv.$FG_PERM.$TRANS_PERM
    if(name.find("equiv") != std::string::npos) {
      return make_non_canon_configuration(primclex, name);
    }

    // if $CANON_SCELNAME/$CANON_INDEX
    return *primclex.db<Configuration>().find(name);
  }

  /// \brief Returns correlations using 'clexulator'.
  Eigen::VectorXd correlations(const Configuration &config, Clexulator const &clexulator) {
    return correlations(config.configdof(), config.supercell(), clexulator);
  }

  /// \brief Returns gradient correlations using 'clexulator', with respect to DoF 'dof_type'
  Eigen::MatrixXd gradcorrelations(const Configuration &config, Clexulator const &clexulator, DoFKey &key) {
    return gradcorrelations(config.configdof(), config.supercell(), clexulator, key);
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
    return comp_n(config) / config.prim().basis().size();
  }

  /// \brief Returns the relaxed energy, normalized per unit cell
  double relaxed_energy(const Configuration &config) {
    return config.calc_properties().scalar("relaxed_energy") / config.supercell().volume();
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
    Clexulator clexulator = primclex.clexulator(formation_energy.bset);
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

  /// \brief Cost function that describes the degree to which basis sites have relaxed
  double atomic_deformation(const Configuration &_config) {
    return _config.calc_properties().scalar("atomic_deformation");
  }

  /// \brief Cost function that describes the degree to which lattice has relaxed
  double lattice_deformation(const Configuration &_config) {
    return _config.calc_properties().scalar("lattice_deformation");
  }

  /// \brief Change in volume due to relaxation, expressed as the ratio V/V_0
  double volume_relaxation(const Configuration &_config) {
    StrainConverter tconvert("BIOT");
    return tconvert.rollup_E(_config.calc_properties().global.at("Ustrain")).determinant();
  }

  /// \brief Returns the relaxed magnetic moment, normalized per unit cell
  double relaxed_magmom(const Configuration &_config) {
    return _config.calc_properties().scalar("relaxed_magmom");
  }

  /// \brief Returns the relaxed magnetic moment, normalized per species
  double relaxed_magmom_per_species(const Configuration &_config) {
    return relaxed_magmom(_config) / n_species(_config);
  }

  /// \brief relaxed forces of configuration, determined from DFT (eV/Angstr.), as a 3xN matrix
  Eigen::MatrixXd relaxed_forces(const Configuration &_config) {
    return _config.calc_properties().site.at("relaxed_forces");
  }

  /// \brief Root-mean-square forces of relaxed configurations, determined from DFT (eV/Angstr.)
  double rms_force(const Configuration &_config) {
    //Get RMS force:
    Eigen::MatrixXd F = relaxed_forces(_config);

    return sqrt((F * F.transpose()).trace() / double(F.cols()));
  }


  /// \brief Returns an IntegralCluster representing the perturbation between the configs
  IntegralCluster config_diff(const Configuration &_config1, const Configuration &_config2) {
    if(_config1.supercell() != _config2.supercell()) {
      throw std::runtime_error("Misuse of basic config_diff: configs are not in same supercell");
    }
    std::vector<UnitCellCoord> uccoords;
    for(Index i = 0 ; i < _config1.occupation().size(); i++) {
      if(_config1.occ(i) != _config2.occ(i)) {
        uccoords.push_back(_config1.uccoord(i));
      }
    }
    IntegralCluster perturb(_config1.prim(), uccoords.begin(), uccoords.end());
    return perturb;
  }

  /// Returns a rotated/translated version of config 2 that leaves it closest to the occupation of config1
  Configuration closest_setting(const Configuration &_config1, const Configuration &_config2) {
    std::vector<PermuteIterator> non_fg_its;
    std::set_difference(_config2.supercell().sym_info().permute_begin(),
                        _config2.supercell().sym_info().permute_end(),
                        _config2.factor_group().begin(),
                        _config2.factor_group().end(),
                        std::back_inserter(non_fg_its));
    int min_size = config_diff(_config1, _config2).size();
    if(non_fg_its.size() == 0) {
      return _config2;
    }
    PermuteIterator best_it = *(non_fg_its.begin());
    for(auto it = non_fg_its.begin(); it != non_fg_its.end(); ++it) {
      int size = config_diff(_config1, copy_apply(*it, _config2)).size();
      if(size < min_size) {
        best_it = *it;
        min_size = size;
      }
    }
    if(min_size == config_diff(_config1, _config2).size()) {
      return _config2;
    }
    return copy_apply(best_it, _config2);

  }



  /// \brief Returns a Configuration with the sites in _clust clipped from _config and placed in _bg
  Configuration config_clip(const Configuration &_config, const Configuration &_bg, IntegralCluster &_clust) {
    if(_config.supercell() != _bg.supercell()) {
      throw std::runtime_error("Misuse of basic config_clip: configs are not in same supercell");
    }
    Configuration tmp = _bg;
    std::vector<Index> l_inds;
    std::vector<Index> l_values;
    for(auto &site : _clust) {
      l_inds.push_back(_config.linear_index(site));
      l_values.push_back(_config.occ(_config.linear_index(site)));
    }
    for(Index i = 0; i < l_inds.size(); i++) {
      tmp.set_occ(l_inds[i], l_values[i]);
    }
    return tmp;
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
    return _config.calc_properties().has_scalar("relaxed_energy");
  }

  bool has_reference_energy(const Configuration &_config) {
    return _config.primclex().has_composition_axes() &&
           _config.primclex().has_chemical_reference();
  }

  bool has_formation_energy(const Configuration &_config) {
    return has_relaxed_energy(_config) && has_reference_energy(_config);
  }

  bool has_rms_force(const Configuration &_config) {
    MappedProperties const &props = _config.calc_properties();
    auto it = props.site.find("relaxed_forces");
    return it != props.site.end() && it->second.cols();
  }

  bool has_atomic_deformation(const Configuration &_config) {
    return _config.calc_properties().has_scalar("atomic_deformation");
  }

  bool has_lattice_deformation(const Configuration &_config) {
    return _config.calc_properties().has_scalar("lattice_deformation");
  }

  bool has_volume_relaxation(const Configuration &_config) {
    return _config.calc_properties().global.count("Ustrain");
  }

  bool has_relaxed_magmom(const Configuration &_config) {
    return _config.calc_properties().has_scalar("relaxed_magmom");
  }

  bool has_relaxed_mag_basis(const Configuration &_config) {
    return _config.calc_properties().site.count("relaxed_mag_basis");
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

  /// \brief Constructor
  ///
  /// \param _scel Supercell to be filled
  /// \param _op SymOp that transforms the input motif before tiling into the
  ///        Supercell that is filled
  FillSupercell::FillSupercell(const std::shared_ptr<Supercell> &_scel, const SymOp &_op) :
    m_supercell_ptr(_scel), m_scel(m_supercell_ptr.get()), m_op(&_op), m_motif_scel(nullptr) {}

  /// \brief Constructor
  ///
  /// \param _scel Supercell to be filled
  /// \param _motif Find the first SymOp that after application to _motif enables
  ///               tiling into _scel
  /// \param _tol tolerance
  ///
  FillSupercell::FillSupercell(const std::shared_ptr<Supercell> &_scel, const Configuration &_motif, double _tol) :
    m_supercell_ptr(_scel), m_scel(m_supercell_ptr.get()), m_op(find_symop(_motif, _tol)), m_motif_scel(nullptr) {}

  Configuration FillSupercell::operator()(const Configuration &motif) const {

    if(&motif.supercell() != m_motif_scel) {
      _init(motif.supercell());
    }

    std::unique_ptr<Configuration> result;
    if(m_supercell_ptr) {
      result = notstd::make_unique<Configuration>(m_supercell_ptr);
    }
    else {
      result = notstd::make_unique<Configuration>(*m_scel);
    }

    // We reorien the starting configuration (motif) by m_op in two steps.
    // In the first step, we transform the DoFs of motif by m_op, without permuting their site indices
    // Then, we utilize m_index_table to decorate the new cell with these transformed DoFs.
    // m_index_table is built in FillSupercel::_init() by explicitly transforming each UnitCellCoordinate of
    // the starting supercell (*m_motif_scel) and finding it in the resultant supercell (*m_supercell_ptr)
    // This ensures that sublattice permutation and microscopic translations are handled appropriately
    ConfigDoF trans_motif(motif.configdof());
    trans_motif.apply_sym_no_permute(*m_op);

    // copy transformed dof, as many times as necessary to fill the supercell
    for(auto const &dof : trans_motif.global_dofs())
      result->configdof().set_global_dof(dof.first, dof.second.values());

    for(Index s = 0; s < m_index_table.size(); ++s) {
      for(Index i = 0; i < m_index_table[s].size(); ++i) {
        Index scel_s = m_index_table[s][i];

        if(trans_motif.has_occupation()) {
          result->configdof().occ(scel_s) = trans_motif.occ(s);
        }
      }
    }

    for(auto const &dof : trans_motif.local_dofs()) {
      LocalContinuousConfigDoFValues &res_ref(result->configdof().local_dof(dof.first));
      for(Index s = 0; s < m_index_table.size(); ++s) {
        for(Index i = 0; i < m_index_table[s].size(); ++i) {
          Index scel_s = m_index_table[s][i];
          res_ref.site_value(scel_s) = dof.second.site_value(s);
        }
      }
    }
    return *result;
  }

  /// \brief Find first SymOp in the prim factor group such that apply(op, motif)
  ///        can be used to fill the Supercell
  const SymOp *FillSupercell::find_symop(const Configuration &motif, double tol) {

    const Lattice &motif_lat = motif.supercell().lattice();
    const Lattice &scel_lat = m_scel->lattice();
    auto begin = m_scel->prim().factor_group().begin();
    auto end = m_scel->prim().factor_group().end();

    auto res = xtal::is_equivalent_superlattice(scel_lat, motif_lat, begin, end, tol);
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
    Lattice oriented_motif_lat = sym::copy_apply(*m_op, m_motif_scel->lattice());

    auto oriended_motif_lattice_points = xtal::make_lattice_points(oriented_motif_lat, m_scel->lattice(), TOL);

    const Structure &prim = m_scel->prim();
    m_index_table.resize(m_motif_scel->num_sites());

    // for each site in motif
    for(Index s = 0 ; s < m_motif_scel->num_sites() ; s++) {

      // apply symmetry to re-orient and find unit cell coord
      UnitCellCoord oriented_uccoord = sym::copy_apply(*m_op, m_motif_scel->uccoord(s), prim.lattice(), prim.basis_permutation_symrep_ID());

      // for each unit cell of the oriented motif in the supercell, copy the occupation
      for(const UnitCell &oriented_motif_uc : oriended_motif_lattice_points) {
        UnitCell oriented_motif_uc_relative_to_prim = oriented_motif_uc.reset_tiling_unit(oriented_motif_lat, prim.lattice());

        Index prim_motif_tile_ind = m_scel->sym_info().unitcell_index_converter()(oriented_motif_uc_relative_to_prim);

        UnitCellCoord mc_uccoord(
          oriented_uccoord.sublattice(),
          m_scel->sym_info().unitcell_index_converter()(prim_motif_tile_ind) + oriented_uccoord.unitcell());

        m_index_table[s].push_back(m_scel->linear_index(mc_uccoord));
      }
    }
  }

  std::ostream &operator<<(std::ostream &sout, const Configuration &c) {
    sout << c.name() << "\n";
    //if(c.has_deformation()) {
    //sout << "Deformation:\n" << c.deformation() << std::endl;
    //}

    for(Index i = 0; i < c.size(); ++i) {
      sout << "Linear index: " << i << "  UnitCellCoord: " << c.uccoord(i) << std::endl;
      if(c.has_occupation()) {
        sout << "  Occupation: " << c.occ(i) << "  (" << c.mol(i).name() << ")\n";
      }
      //if(c.has_displacement()) {
      //sout << "  Displacement: " << c.disp(i).transpose() << "\n";
      //}
    }

    return sout;
  }

  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell &scel, Clexulator const &clexulator) {

    //Size of the supercell will be used for normalizing correlations to a per primitive cell value
    int scel_vol = scel.volume();

    Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

    //Inform Clexulator of the bitstring

    //Holds contribution to global correlations from a particular neighborhood
    Eigen::VectorXd tcorr = correlations;
    //std::vector<double> corr(clexulator.corr_size(), 0.0);

    for(int v = 0; v < scel_vol; v++) {
      //Fill up contributions
      clexulator.calc_global_corr_contribution(configdof,
                                               scel.nlist().sites(v).data(),
                                               end_ptr(scel.nlist().sites(v)),
                                               tcorr.data(),
                                               end_ptr(tcorr));

      correlations += tcorr;
    }

    correlations /= (double) scel_vol;

    return correlations;
  }

  /// \brief Returns gradient correlations using 'clexulator', with respect to DoF 'dof_type'
  Eigen::MatrixXd gradcorrelations(const ConfigDoF &configdof, const Supercell &scel, Clexulator const &clexulator, DoFKey &key) {
    ClexParamKey paramkey;
    ClexParamKey corr_key(clexulator.param_pack().key("corr"));
    ClexParamKey dof_key;
    if(key == "occ") {
      paramkey = clexulator.param_pack().key("diff/corr/" + key + "_site_func");
      dof_key = clexulator.param_pack().key("occ_site_func");
    }
    else {
      paramkey = clexulator.param_pack().key("diff/corr/" + key + "_var");
      dof_key = clexulator.param_pack().key(key + "_var");
    }

    std::string em_corr, em_dof;
    em_corr = clexulator.param_pack().eval_mode(corr_key);
    em_dof = clexulator.param_pack().eval_mode(dof_key);

    // this const_cast is not great...
    // but it seems like the only place passing const Clexulator is a problem and it is not
    // actually changing clexulator before/after this function
    const_cast<Clexulator &>(clexulator).param_pack().set_eval_mode(corr_key, "DIFF");
    const_cast<Clexulator &>(clexulator).param_pack().set_eval_mode(dof_key, "DIFF");

    Eigen::MatrixXd gcorr;
    Index scel_vol = scel.volume();
    if(DoF::BasicTraits(key).global()) {
      Eigen::MatrixXd gcorr_func = configdof.global_dof(key).values();
      gcorr.setZero(gcorr_func.size(), clexulator.corr_size());
      //Holds contribution to global correlations from a particular neighborhood

      //std::vector<double> corr(clexulator.corr_size(), 0.0);
      for(int v = 0; v < scel_vol; v++) {

        //Fill up contributions
        clexulator.calc_global_corr_contribution(configdof,
                                                 scel.nlist().sites(v).data(),
                                                 end_ptr(scel.nlist().sites(v)));

        for(Index c = 0; c < clexulator.corr_size(); ++c)
          gcorr.col(c) += clexulator.param_pack().read(paramkey(c));


      }
    }
    else {

      Eigen::MatrixXd gcorr_func;
      gcorr.setZero(configdof.local_dof(key).values().size(), clexulator.corr_size());
      //Holds contribution to global correlations from a particular neighborhood
      Index l;
      for(int v = 0; v < scel_vol; v++) {
        //Fill up contributions
        clexulator.calc_global_corr_contribution(configdof,
                                                 scel.nlist().sites(v).data(),
                                                 end_ptr(scel.nlist().sites(v)));

        for(Index c = 0; c < clexulator.corr_size(); ++c) {
          gcorr_func = clexulator.param_pack().read(paramkey(c));

          for(Index n = 0; n < scel.nlist().sites(v).size(); ++n) {
            l = scel.nlist().sites(v)[n];
            //for(Index i=0; i<gcorr_func.cols(); ++i){
            gcorr.block(l * gcorr_func.rows(), c, gcorr_func.rows(), 1) += gcorr_func.col(n);
            //std::cout << "Block: (" << l * gcorr_func.rows() << ", " << c << ", " << gcorr_func.rows() << ", " << 1 << ") += " << gcorr_func.col(n).transpose() << "\n";
            //}
          }
        }
      }
    }
    const_cast<Clexulator &>(clexulator).param_pack().set_eval_mode(corr_key, em_corr);
    const_cast<Clexulator &>(clexulator).param_pack().set_eval_mode(dof_key, em_dof);


    return gcorr;
  }


  /// \brief Returns num_each_molecule(molecule_type), where 'molecule_type' is ordered as Structure::xtal::struc_molecule_name()
  Eigen::VectorXi num_each_molecule(const ConfigDoF &configdof, const Supercell &scel) {

    auto mol_names = xtal::struc_molecule_name(scel.prim());
    // [basis_site][site_occupant_index]
    auto convert = make_index_converter(scel.prim(), mol_names);

    // create an array to count the number of each molecule
    Eigen::VectorXi num_each_molecule = Eigen::VectorXi::Zero(mol_names.size());

    // count the number of each molecule
    for(Index i = 0; i < configdof.size(); i++) {
      num_each_molecule(convert[ scel.sublat(i) ][ configdof.occ(i)])++;
    }

    return num_each_molecule;
  }

  /// \brief Returns comp_n, the number of each molecule per primitive cell, ordered as Structure::xtal::struc_molecule_name()
  Eigen::VectorXd comp_n(const ConfigDoF &configdof, const Supercell &scel) {
    return num_each_molecule(configdof, scel).cast<double>() / scel.volume();
  }

}
