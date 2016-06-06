#include "casm/clex/Configuration.hh"

#include <sstream>
//#include "casm/misc/Time.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/symmetry/SymBasisPermute.hh"

namespace CASM {

  /// Construct a default Configuration
  Configuration::Configuration(Supercell &_supercell, const jsonParser &src, const ConfigDoF &_configdof)
    : m_id("none"), m_supercell(&_supercell), m_source_updated(true), m_multiplicity(-1), m_dof_updated(true),
      m_configdof(_configdof), m_prop_updated(true), m_selected(false) {
    set_source(src);
  }

  //*********************************************************************************
  /// Construct by reading from main data file (json)
  Configuration::Configuration(const jsonParser &json, Supercell &_supercell, Index _id)
    : m_supercell(&_supercell), m_source_updated(false), m_multiplicity(-1), m_dof_updated(false),
      m_configdof(_supercell.num_sites()), m_prop_updated(false) {

    std::stringstream ss;
    ss << _id;
    m_id = ss.str();

    read(json);
  }


  //********** MUTATORS  ***********

  void Configuration::set_id(Index _id) {
    std::stringstream ss;
    ss << _id;
    m_id = ss.str();

    m_dof_updated = true;
    m_prop_updated = true;
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
    m_source_updated = true;
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
        }
      }
    }
  }

  //*********************************************************************************
  void Configuration::set_occupation(const Array<int> &new_occupation) {
    m_dof_updated = true;
    m_configdof.set_occupation(new_occupation);
    return;
  }

  //*********************************************************************************
  void Configuration::set_occ(Index site_l, int val) {
    //if( i > occupation.size() || i < 0)
    //{
    //    std::cout << "ERROR in Configuration::set_occ(). i: " << i << " occupation.size(): "<< occupation.size() << "  val: " << val << std::endl;
    //    exit(1);
    //}
    //std::cout << "Configuration::set_occ(). i: " << i << " occupation.size(): "<< occupation.size() << "  val: " << val << std::endl;
    m_dof_updated = true;
    m_configdof.occ(site_l) = val;
  }

  //*********************************************************************************

  void Configuration::set_displacement(const displacement_matrix_t &new_displacement) {
    m_dof_updated = true;
    m_configdof.set_displacement(new_displacement);
  }

  //*********************************************************************************

  void Configuration::set_deformation(const Eigen::Matrix3d &new_deformation) {
    m_dof_updated = true;
    m_configdof.set_deformation(new_deformation);
  }

  //*********************************************************************************

  Configuration Configuration::canonical_form(PermuteIterator it_begin, PermuteIterator it_end, PermuteIterator &it_canon, double tol) const {
    Configuration tconfig(*this);
    tconfig.m_configdof = m_configdof.canonical_form(it_begin, it_end, it_canon, tol);
    return tconfig;
  }

  //*********************************************************************************
  void Configuration::set_calc_properties(const jsonParser &calc) {
    m_prop_updated = true;
    m_calculated = calc;
    //delta = calculated - reference;
  }

  //*********************************************************************************

  bool Configuration::read_calc_properties(jsonParser &parsed_props) const {
    //std::cout << "begin Configuration::read_calculated()" << std::endl;
    bool success = true;
    /// properties.calc.json: contains calculated properties
    ///   Currently only loading those properties that have references
    fs::path filepath = calc_properties_path();
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
        Eigen::MatrixXd forces;
        from_json(forces, json["relaxed_forces"]);
        parsed_props["rms_force"] = sqrt((forces.transpose() * forces).trace() / double(forces.rows()));
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

  std::string Configuration::id() const {
    return m_id;
  }

  //*********************************************************************************
  std::string Configuration::name() const {
    return supercell().name() + "/" + id();
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
  fs::path Configuration::path() const {
    return supercell().path() / id();
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
  Supercell &Configuration::supercell() {
    return *m_supercell;
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

  const Properties &Configuration::generated_properties() const {
    return m_generated;
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
    auto convert = index_converter(prim(), prim().get_struc_molecule());

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
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's get_struc_molecule(), with [Va]=0.0
  Eigen::VectorXd Configuration::composition() const {

    // get the number of each molecule type
    Eigen::VectorXi _num_each_molecule = num_each_molecule();

    /// get the total number of non-vacancy atoms
    int num_atoms = 0;

    // need to know which molecules are vacancies
    auto struc_molecule = prim().get_struc_molecule();

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
  ///    composition[ molecule_type ]: molecule_type ordered as prim structure's get_struc_molecule()
  Eigen::VectorXd Configuration::true_composition() const {
    return num_each_molecule().cast<double>() / size();
  }

  //*********************************************************************************
  /// Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::get_struc_molecule()
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
  ///   where 'component_type' is ordered as ParamComposition::get_components
  Eigen::VectorXd Configuration::num_each_component() const {

    // component order used for param_composition
    std::vector<std::string> v_components = primclex().composition_axes().components();

    // copy to CASM::Array
    std::vector<std::string> components;
    for(auto it = v_components.cbegin(); it != v_components.cend(); ++it) {
      components.push_back(*it);
    }

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


  /// Writes the Configuration to a json object
  ///   Uses PrimClex's current settings to write the appropriate
  ///   Properties, DeltaProperties and Correlations files
  ///
  ///   'json' is a jsonParser JSON object (or will be set to a JSON object)
  ///   Configuration data is saved in several object, we write the *'d objects:
  ///
  ///   *config.json:             json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["dof"]
  ///   *corr.json:               json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CLEX"]["corr"]
  ///    properties.calc.json:    casmroot/supercells/SCEL_NAME/CONFIG_ID/CURR_CALCTYPE/properties.calc.json
  ///   *param_composition.json:  json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["param_composition"]
  ///   *properties.ref.json:     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["ref"]
  ///   *properties.calc.json:    json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]     // contains param_comp
  ///   *properties.delta.json:   json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["delta"]
  ///   *properties.generated.json:   json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["generated"]
  jsonParser &Configuration::write(jsonParser &json) const {

    //std::cout << "begin Configuration::write()" << std::endl;

    const ProjectSettings &set = primclex().settings();
    std::string calc_string = "calctype." + set.calctype();
    std::string ref_string = "ref." + set.ref();

    /// write json object hierarchy if not existing
    jsonParser &json_scel = json["supercells"][supercell().name()];
    jsonParser &json_config = json_scel[id()];
    //jsonParser &json_bset = json_config[primclex().curr_bset()];
    jsonParser &json_ref = json_config[calc_string][ref_string];
    jsonParser &json_prop = json_ref["properties"];

    json_config["selected"] = selected();

    if(!json_config.contains("dof")) {
      write_dof(json_config);
    }

    if(!json_config.contains("pos")) {
      //write_pos(json_config);
    }

    if(m_source_updated) {
      write_source(json_config);
    }

    if(!json_ref.contains("param_composition") || m_prop_updated) {
      //write_param_composition(json_ref);
    }

    if(m_prop_updated) {
      write_properties(json_prop);
    }

    //std::cout << "finish Configuration::write()" << std::endl;

    return json;
  }

  //*********************************************************************************

  void Configuration::write_pos() const {

    try {
      fs::create_directories(path());
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in Configuration::write_pos()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream file(pos_path());
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
    stream << id() << " ";

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

    Eigen::VectorXd comp = composition();
    auto mol_list = prim().get_struc_molecule();

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

    Eigen::VectorXd true_comp = true_composition();

    for(Index i = 0; i < true_comp.size(); i++) {
      stream.precision(6);
      stream.width(12);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream << true_comp[i] << " ";
    }

  }

  //*********************************************************************************
  void Configuration::print_sublattice_composition(std::ostream &stream) const {

    std::vector<Eigen::VectorXd> sublattice_comp = sublattice_composition();

    for(Index i = 0; i < sublattice_comp.size(); i++) {
      for(Index j = 0; j < sublattice_comp[i].size(); j++) {
        stream.precision(6);
        stream.width(12);
        stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
        stream << sublattice_comp[i][j] << " ";
      }
    }

  }

  //*********************************************************************************

  /// Private members:

  /// Reads the Configuration from the expected casm directory
  ///   Uses PrimClex's current settings to read in the appropriate
  ///   Properties, DeltaProperties and Correlations files if they exist
  ///
  /// This is private, because it is only called from the constructor:
  ///   Configuration(const Supercell &_supercell, Index _id)
  ///   It's called from the constructor because of the Supercell pointer
  ///
  ///   Configuration data is saved in several object, we write the *'d objects:
  ///
  ///   *config.json:             json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["dof"]
  ///   *corr.json:               json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CLEX"]["corr"]
  ///    properties.calc.json:    son["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["properties"]["calc"]
  ///    param_composition.json:  json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["param_composition"]
  ///   *properties.ref.json:     json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["ref"]
  ///   *properties.calc.json:    json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["calc"]
  ///    properties.delta.json:   json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["delta"]
  ///   *properties.generated.json:   json["supercells"]["SCEL_NAME"]["CONFIG_ID"]["CURR_CALCTYPE"]["CURR_REF"]["properties"]["generated"]
  ///
  void Configuration::read(const jsonParser &json) {

    //std::cout << "begin  Configuration::read()" << std::endl;

    const ProjectSettings &set = primclex().settings();

    std::string calc_string = "calctype." + set.calctype();
    std::string ref_string = "ref." + set.ref();

    // read dof
    if(!json.contains("supercells"))
      return;
    const jsonParser &json_scel = json["supercells"];
    if(!json_scel.contains(supercell().name()))
      return;
    if(!json_scel[supercell().name()].contains(id()))
      return;
    const jsonParser &json_config = json_scel[supercell().name()][id()];

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

  /// Read source and degree of freedom info
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID
  ///
  void Configuration::read_dof(const jsonParser &json) {

    /// json["dof"]: contains degree of freedom information
    if(!json.contains("dof")) {
      m_id = "none";
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
  ///   location: json = casmroot/supercells/SCEL_NAME/CONFIG_ID/CURR_CALCTYPE/CURR_REF/properties
  ///
  /// Tries to read reference from json["ref"], if it doesn't exist, tries to generate references
  /// Tries to read calculated from json["calc"]
  /// Tries to read generated from json["gen"]
  ///
  /// Calculates delta = calculated - reference
  ///
  void Configuration::read_properties(const jsonParser &json) {
    /*
        if(!json.contains("ref")) {
          generate_reference();
        }
        else {
          from_json(reference, json["ref"]);
        }
    */
    if(json.contains("calc")) {
      from_json(m_calculated, json["calc"]);
    }

    if(json.contains("gen")) {
      from_json(m_generated, json["gen"]);
    }

    //    delta = calculated - reference;

  }

  //*********************************************************************************
  fs::path Configuration::pos_path() const {
    return primclex().dir().POS(name());
  }

  //*********************************************************************************
  fs::path Configuration::calc_properties_path() const {
    return primclex().dir().calculated_properties(name(), primclex().settings().calctype());
  }

  //*********************************************************************************
  fs::path Configuration::calc_status_path() const {
    return primclex().dir().calc_status(name(), primclex().settings().calctype());
  }

  //*********************************************************************************
  /// Write config.json file containing degree of freedom info
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID, adds: dof
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

    if(!is_strained()) {
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
    if(occupation() != supercell().vacant()) {
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

    if(!primclex().has_composition_axes()) {
      json.erase("param_comp_formula");
      json.erase("param_composition");
      return json;
    }

    json["param_comp_formula"] = primclex().composition_axes().mol_formula();
    json["param_composition"] = param_composition();

    return json;

  }

  //*********************************************************************************
  /// Write properties.calc.json file containing calculated properties and param_composition
  ///   location: json = supercells/SCEL_NAME/CONFIG_ID/CURR_CLEX/CURR_REF/properties, adds: calc
  ///
  jsonParser &Configuration::write_properties(jsonParser &json) const {

    if(!primclex().has_composition_axes()) {

      //json.erase("calc");
      json.erase("ref");
      json.erase("delta");
      json.erase("gen");

      //return json;
    }

    if(m_calculated.size() == 0) {
      json.erase("calc");
    }
    else {
      json["calc"] = m_calculated;
    }
    /*
        if(reference.size() == 0) {
          json.erase("ref");
        }
        else {
          json["ref"] = reference;
        }

        if(delta.size() == 0) {
          json.erase("delta");
        }
        else {
          json["delta"] = delta;
        }
    */
    if(m_generated.size() == 0) {
      json.erase("gen");
    }
    else {
      json["gen"] = m_generated;
    }

    return json;

  }

  /// \brief Returns correlations using 'clexulator'.
  ///
  /// This still assumes that the PrimClex and Supercell are set up for this, so make sure you've called:
  /// - primclex.populate_basis_tables(basis_type);
  /// - primclex.populate_cluster_basis_function_tables();
  /// - primclex.generate_full_nlist();
  /// - primclex.populate_clex_supercell_nlists();
  ///
  /// In the future that shouldn't be necessary
  ///
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

  /// \brief Returns the composition as species fraction, with [Va] = 0.0, in the order of Structure::get_struc_molecule
  ///
  /// - Currently, this is really a Molecule fraction
  Eigen::VectorXd species_frac(const Configuration &config) {
    Eigen::VectorXd v = comp_n(config);
    if(config.primclex().vacancy_allowed()) {
      v(config.primclex().vacancy_index()) = 0.0;
    }
    return v / v.sum();
  }

  /// \brief Returns the composition as site fraction, in the order of Structure::get_struc_molecule
  Eigen::VectorXd site_frac(const Configuration &config) {
    return comp_n(config) / config.prim().basis.size();
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
    Clexulator clexulator = config.primclex().global_clexulator();
    return config.primclex().global_eci("formation_energy") * correlations(config, clexulator);
  }

  /// \brief Returns the formation energy, normalized per unit cell
  double clex_formation_energy_per_species(const Configuration &config) {
    return clex_formation_energy(config) / n_species(config);
  }

  /// \brief Return true if all current properties have been been calculated for the configuration
  bool is_calculated(const Configuration &config) {
    return std::all_of(config.primclex().settings().properties().begin(),
                       config.primclex().settings().properties().end(),
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
    return _config.is_primitive(_config.supercell().permute_begin());
  }

  /// \brief returns true if _config no symmetry transformation applied to _config will increase its lexicographic order
  bool is_canonical(const Configuration &_config) {
    return _config.is_canonical(_config.supercell().permute_begin(), _config.supercell().permute_end());
  }

  bool has_relaxed_energy(const Configuration &_config) {
    return _config.calc_properties().contains("relaxed_energy");
  }

  bool has_reference_energy(const Configuration &_config) {
    return _config.primclex().has_chemical_reference();
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


  /// \brief Application results in filling supercell 'scel' with reoriented motif, op*motif
  ///
  /// Currently only applies to occupation
  Configuration &apply(const ConfigTransform &f, Configuration &motif) {

    Configuration result(f.scel);
    const PrimClex &primclex = motif.primclex();
    const Structure &prim = motif.prim();

    Lattice oriented_motif_lat = copy_apply(f.op, motif.supercell().real_super_lattice());

    // Create a PrimGrid linking the prim and the oriented motif each to the supercell
    // So we can tile the decoration of the motif config onto the supercell correctly
    PrimGrid prim_grid(oriented_motif_lat, f.scel.real_super_lattice());

    // For each site in the motif, re-orient and then translate by all possible
    // translations from prim_grid, index in the mc_prim_prim, and use index to
    // assign occupation

    // std::vector of occupations in the MonteCarlo cell
    Array<int> tscel_occ(motif.size()*prim_grid.size());

    // for each site in motif
    for(Index s = 0 ; s < motif.size() ; s++) {

      // apply symmetry to re-orient and find unit cell coord
      UnitCellCoord oriented_uccoord = motif.uccoord(s).copy_apply(f.op);

      // for each unit cell of the oriented motif in the supercell, copy the occupation
      for(Index i = 0 ; i < prim_grid.size() ; i++) {

        Index prim_motif_tile_ind = f.scel.prim_grid().find(prim_grid.coord(i, PRIM));

        UnitCellCoord mc_uccoord(
          prim,
          oriented_uccoord.sublat(),
          f.scel.prim_grid().unitcell(prim_motif_tile_ind) + oriented_uccoord.unitcell()
        );

        Index occ_ind = f.scel.find(mc_uccoord);

        tscel_occ[occ_ind] = motif.occ(s);
      }
    }

    result.set_occupation(tscel_occ);

    return motif = result;

  }

}


