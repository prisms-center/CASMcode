#include <exception>

#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/casm_io/Log.hh"

#include "casm/casm_io/json_io/container.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/basis_set/DoFIsEquivalent.hh"
#include "casm/basis_set/DoFIsEquivalent_impl.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {


  Site::Site(const Lattice &init_home) :
    Coordinate(init_home),
    m_label(-1),
    m_type_ID(-1),
    m_site_occupant(MoleculeOccupant(DoFType::occupation(), "s" /*variable name*/, std::vector<Molecule>()).clone()) {
    //site_occupant.set_value(0);
  }

  //****************************************************

  Site::Site(const Coordinate &init_pos, const std::string &occ_name) :
    Coordinate(init_pos),
    m_label(-1),
    m_type_ID(-1),
    m_site_occupant(MoleculeOccupant(DoFType::occupation(), "s"/* variable name*/, std::vector<Molecule>()).clone()) {

    std::vector<Molecule> tocc;
    tocc.push_back(Molecule::make_atom(occ_name));
    m_site_occupant->set_domain(tocc);
    m_site_occupant->set_value(0);

    return;
  }

  /// \brief Construct site with initial position and the allowed Molecule
  Site::Site(const Coordinate &init_pos, std::initializer_list<Molecule> site_occ) :
    Coordinate(init_pos),
    m_label(-1),
    m_type_ID(-1),
    m_site_occupant(MoleculeOccupant(DoFType::occupation(), "s"/* variable name*/, std::vector<Molecule>(site_occ)).clone()) {

  }

  Site::~Site() {}

  //****************************************************

  const MoleculeOccupant &Site::site_occupant() const {
    return *m_site_occupant;
  }

  //****************************************************

  Index Site::dof_size() const {
    return m_dof_map.size();
  }

  //****************************************************

  DoFSet const &Site::dof(std::string const &_dof_type) const {
    auto it = m_dof_map.find(_dof_type);
    if(it != m_dof_map.end())
      return (it->second);
    else
      throw std::runtime_error(std::string("In Structure::dof(), this structure does not contain any global DoF's of type ") + _dof_type);
  }

  //****************************************************
  bool Site::has_dof(std::string const &_dof_type) const {
    return site_occupant().size() > 1 || m_dof_map.find(_dof_type) != m_dof_map.end();
  }

  //****************************************************

  std::vector<std::string> Site::dof_types() const {
    std::vector<std::string> result;
    //if(site_occupant().size() > 1)
    //result.push_back(site_occupant().type_name());
    for(auto it = m_dof_map.begin(); it != m_dof_map.end(); ++it)
      result.push_back(it->first);
    return result;
  }

  //****************************************************

  bool Site::time_reversal_active() const {
    for(auto it = m_dof_map.begin(); it != m_dof_map.end(); ++it)
      if(DoF::traits(it->first).time_reversal_active())
        return true;

    return false;
  }

  //****************************************************

  bool Site::is_vacant()const {
    return site_occupant().is_specified()
           && site_occupant().occ().is_vacancy();
  }

  //****************************************************

  Index Site::label() const {
    return m_label;
  };

  //****************************************************

  std::string Site::occ_name() const {
    if(!site_occupant().is_specified())
      return "?";
    else
      return site_occupant().occ().name();
  }

  //****************************************************

  const Molecule &Site::occ() const {
    return site_occupant().occ();
  }

  //****************************************************
  Site &Site::_apply_sym_attributes(const SymOp &op) {
    for(Index i = 0; i < site_occupant().size(); i++)
      (*m_site_occupant)[i].apply_sym(op);

    auto it = m_dof_map.begin();
    for(; it != m_dof_map.end(); ++it)
      it->second = copy_apply(op, it->second);

    m_type_ID = -1;
    return *this;
  }
  //****************************************************
  /**
   *
   */
  //****************************************************

  Site &Site::apply_sym(const SymOp &op) {
    Coordinate::apply_sym(op);
    _apply_sym_attributes(op);
    return *this;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  Site &Site::apply_sym_no_trans(const SymOp &op) {
    Coordinate::apply_sym_no_trans(op);
    _apply_sym_attributes(op);

    return *this;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  Site &Site::operator+=(const Coordinate &translation) {
    Coordinate::operator += (translation);
    return *this;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  Site &Site::operator-=(const Coordinate &translation) {
    Coordinate::operator -= (translation);
    return *this;
  }

  //*******************************************************************************************
  /**
   *
   */
  //*******************************************************************************************

  bool Site::compare(const Coordinate &test_coord) const {
    return (min_dist(test_coord) < lattice().tol());
  }

  //*******************************************************************************************
  /**
   *
   */
  //*******************************************************************************************

  bool Site::compare(const Site &test_site) const {
    return (compare_type(test_site) && min_dist(test_site) < lattice().tol());
  }

  //*******************************************************************************************
  /**
   *
   *
   */
  //*******************************************************************************************

  bool Site::compare(const Site &test_site, const Coordinate &shift) const {

    return (compare_type(test_site)) && (min_dist(test_site + shift) < lattice().tol());

  }

  //*******************************************************************************************
  /**
   *
   */
  //*******************************************************************************************

  bool Site::compare_type(const Site &test_site) const {
    assert(((site_occupant().size() <= 1 || test_site.site_occupant().size() <= 1)
            || (site_occupant().is_specified() == test_site.site_occupant().is_specified()))
           && "In Site::compare_type() comparing initialized occupant to uninitialized occupant!  This isn't a good idea!");

    return (_type_ID() == test_site._type_ID()) && site_occupant().value() == test_site.site_occupant().value();
  }

  //*******************************************************************************************

  bool Site::operator==(const Site &test_site) const {
    return (compare_type(test_site) && Coordinate::operator==(test_site));
  }

  //*******************************************************************************************

  bool Site::almost_equal(const Site &test_site) const {
    return (compare_type(test_site) && dist(test_site) < lattice().tol());
  }

  //****************************************************

  bool Site::contains(const std::string &name) const {
    for(Index i = 0; i < site_occupant().size(); i++)
      if(site_occupant()[i].contains(name)) {
        return true;
      }
    return false;
  }

  //****************************************************

  bool Site::contains(const std::string &name, int &index) const {
    for(Index i = 0; i < site_occupant().size(); i++)
      if(site_occupant()[i].contains(name)) {
        index = i;
        return true;
      }
    return false;
  }

  //****************************************************

  void Site::set_allowed_occupants(std::vector<Molecule> const &_occ_domain) {
    m_site_occupant->set_domain(_occ_domain);
    m_type_ID = -1;
  }

  //****************************************************

  void Site::set_occ_value(int new_val) {
    m_site_occupant->set_value(new_val);
  }

  //****************************************************

  void Site::set_occ(const Molecule &new_occ) {
    m_site_occupant->set_current_state(new_occ);
  }

  //****************************************************

  void Site::set_local_dofs(std::vector<DoFSet> const &_new_dofs) {
    m_dof_map.clear();
    for(DoFSet const &_dof : _new_dofs)
      m_dof_map.emplace(std::make_pair(_dof.type_name(), _dof));
  }


  //****************************************************

  std::vector<std::string> Site::allowed_occupants() const {
    std::vector<std::string> occ_list;
    for(Index i = 0; i < site_occupant().size(); i++) {
      occ_list.push_back(site_occupant()[i].name());
    }
    return occ_list;
  }

  //****************************************************

  void Site::set_basis_ind(Index new_ind) {
    Coordinate::set_basis_ind(new_ind);
    m_site_occupant->set_ID(new_ind);
    for(auto &dof : m_dof_map)
      dof.second.set_ID(new_ind);
  }

  //****************************************************

  void Site::set_label(Index new_ind) {
    if(new_ind == m_label)
      return;
    m_label = new_ind;

    m_type_ID = -1;
    return;
  }

  //****************************************************
  //   read site, including all possible occupants
  void Site::read(std::istream &stream, bool SD_is_on) {
    set_label(-1);

    char ch;

    AtomPosition::sd_type SD_flag;

    Coordinate::read(stream, COORD_MODE::CHECK());
    if(SD_is_on) {
      for(int i = 0; i < 3; i++) {
        stream >> ch;
        if(ch == 'T') {
          SD_flag[i] = true;
        }
        else if(ch == 'F') {
          SD_flag[i] = false;
        }
      }
    }

    std::vector<Molecule> tocc;

    ch = stream.peek();
    //    while(ch != '\n' && !stream.eof()) {
    while(ch != '\n' && ch != ':' && !stream.eof()) {
      //      while((ch < 'a' || ch > 'z') && (ch < 'A' || ch > 'Z') && ch != '\n' && !stream.eof()) {
      while((ch < 'a' || ch > 'z') && (ch < 'A' || ch > 'Z') && ch != '\n' && ch != ':' && !stream.eof()) {
        stream.ignore();
        ch = stream.peek();
      }
      if(ch != '\n' && ch != ':' && !stream.eof()) {
        //Need to change this part for real molecules
        std::string mol_name;
        stream >> mol_name;
        tocc.push_back(Molecule::make_atom(mol_name));
      }
      ch = stream.peek();
    }

    if(ch == ':') {
      stream.ignore();
      stream.ignore();

      std::string mol_name;
      stream >> mol_name;

      if(tocc.size()) {
        m_site_occupant->set_domain(tocc);
        Index index = tocc.size();
        for(Index i = 0; i < tocc.size(); i++)
          if(tocc[i].name() == mol_name) {
            index = i;
            break;
          }
        if(index == tocc.size()) {
          default_err_log() << "ERROR in Site::read(). Occupying molecule not listed in possible occupants" << std::endl;
          default_err_log() << "  occupying molecule name: " << mol_name << "  index: " << index << std::endl;
          default_err_log() << "  possible occupants: ";
          for(Index i = 0; i < tocc.size(); i++)
            default_err_log() << tocc[i].name() << " ";
          default_err_log() << " " << std::endl;
          exit(1);
        }
        else
          m_site_occupant->set_value(index);

      }
      else {
        default_err_log() << "WARNING: Trying to read Site info, but no valid input was received." << std::endl;
      }
      m_type_ID = -1;
      return;
    }

    if(tocc.size()) {
      m_site_occupant->set_domain(tocc);
    }
    else {
      default_err_log() << "WARNING: Trying to read Site info, but no valid input was received." << std::endl;
    }
    stream.ignore(1000, '\n');

    m_type_ID = -1;
    return;
  }

  //****************************************************
  // read site, using 'elem' as site occupant domain
  void Site::read(std::istream &stream, std::string &elem, bool SD_is_on) {
    char ch;

    set_label(-1);

    AtomPosition::sd_type SD_flag;

    Coordinate::read(stream, COORD_MODE::CHECK());
    if(SD_is_on) {
      for(int i = 0; i < 3; i++) {
        stream >> ch;
        if(ch == 'T') {
          SD_flag[i] = true;
        }
        else if(ch == 'F') {
          SD_flag[i] = false;
        }
      }
    }

    std::vector<Molecule> tocc;

    tocc.push_back(Molecule::make_atom(elem, SD_flag));

    if(tocc.size()) {
      m_site_occupant->set_domain(tocc);
    }
    else {
      default_err_log() << "WARNING: Trying to read Site info, but no valid input was received." << std::endl;
    }
    stream.ignore(1000, '\n');

    m_type_ID = -1;
    return;
  }


  //****************************************************
  /**	Print coordinate of site with name of all possible
  *		occupying molecule
   */
  //****************************************************

  void Site::print(std::ostream &stream, Eigen::IOFormat format) const {
    //site_occupant().occ().print(stream, *this, SD_is_on);
    Coordinate::print(stream, 0, format);
    stream << " ";
    site_occupant().print(stream);
    stream << std::flush;

    return;
  }

  //****************************************************
  /**	Print coordinate of site with name of occupying molecule
   */
  //****************************************************

  void Site::print_occ(std::ostream &stream, Eigen::IOFormat format) const {
    //site_occupant().occ().print(stream, *this, SD_is_on);
    Site::print(stream, format);
    stream << " :: ";
    site_occupant().print_occ(stream);
    stream << std::flush;

    return;
  }

  //****************************************************
  /**	Print each AtomPosition in current molecule,
   *		with name of occupying atom
   */
  //****************************************************

  void Site::print_mol(std::ostream &stream,
                       int spaces,
                       char delim,
                       bool SD_is_on) const {

    Eigen::Matrix3d c2f = Eigen::Matrix3d::Identity();
    // change this to use FormatFlag
    if(COORD_MODE::IS_FRAC())
      c2f = home().inv_lat_column_mat();
    site_occupant().occ().print(stream,
                                cart(),
                                c2f,
                                spaces,
                                delim,
                                SD_is_on);
    return;
  }

  //****************************************************

  std::ostream &operator<< (std::ostream &stream, const Site &site) {
    site.print(stream);
    return stream;
  }

  //****************************************************

  Site jsonConstructor<Site>::from_json(const jsonParser &json, Lattice const &_home, COORD_TYPE coordtype) {
    Site result(_home);
    CASM::from_json(result, json, _home, coordtype);
    return result;
  }

  //****************************************************

  jsonParser &to_json(const Site &site, jsonParser &json, COORD_TYPE coordtype) {
    json.put_obj();

    // class Site : public Coordinate
    if(coordtype == FRAC)
      to_json_array(site.frac(), json["coordinate"]);
    else
      to_json_array(site.cart(), json["coordinate"]);

    // MoleculeOccupant site_occupant;
    Eigen::Matrix3d c2f = Eigen::Matrix3d::Identity();
    // change this to use FormatFlag
    if(coordtype == FRAC)
      c2f = site.home().inv_lat_column_mat();
    CASM::to_json(site.site_occupant(), json["site_occupant"], c2f);


    // Index m_label
    if(valid_index(site.label()))
      json["label"] = site.label();


    return json;

  }

  //****************************************************

  void from_json(Site &site, const jsonParser &json, Lattice const &_home, COORD_TYPE coordtype) {
    site.set_lattice(_home, coordtype);
    if(coordtype == FRAC)
      site.frac() = json["coordinate"].get<Eigen::Vector3d>();
    else
      site.cart() = json["coordinate"].get<Eigen::Vector3d>();

    // MoleculeOccupant -- allowed occupants at site
    //MoleculeOccupant t_occ(DoF::traits("occ"));
    //if(coordtype == FRAC)
    //CASM::from_json(t_occ, json["site_occupant"], _home.inv_lat_column_mat());
    //else
    //CASM::from_json(t_occ, json["site_occupant"], Eigen::Matrix3d::Identity());
    //if(t_occ.size())
    //site.set_allowed_occupants(t_occ.domain());
    //else
    //site.set_allowed_occupants({Molecule::make_unknown()});

    // Index m_label -- must be greater than zero
    Index _label = -1;
    if(json.contains("label")) {
      CASM::from_json(_label, json["label"]);
      if(!valid_index(_label))
        throw std::runtime_error("JSON specification of site has {\"label\" : " + std::to_string(_label) + "}, but \"label\" must be greater than 0.\n");
    }
    site.set_label(_label);

    // Local continuous dofs
    std::vector<DoFSet> _dof_vec;
    if(json.contains("dof")) {
      std::map<std::string, DoFSet> _dof_map;
      auto it = json["dof"].begin(), end_it = json["dof"].end();
      for(; it != end_it; ++it) {
        if(_dof_map.count(it.name()))
          throw std::runtime_error("Error parsing global field \"dof\" from JSON. DoF type " + it.name() + " cannot be repeated.");

        try {
          _dof_map.emplace(std::make_pair(it.name(), it->get<DoFSet>(DoF::traits(it.name()))));
        }
        catch(std::exception &e) {
          throw std::runtime_error("Error parsing global field \"dof\" from JSON. Failure for DoF type " + it.name() + ": " + e.what());
        }

      }
      for(auto const &_dof : _dof_map)
        _dof_vec.push_back(_dof.second);
    }
    site.set_local_dofs(_dof_vec);


  }

  //*******************************************************************************************

  bool Site::_compare_type_no_ID(const Site &_other) const {
    //compare domain but not value
    if(!(label() == _other.label() && OccupantDoFIsEquivalent<Molecule>(site_occupant())(_other.site_occupant())))
      return false;

    if(m_dof_map.size() != _other.m_dof_map.size())
      return false;

    auto it1 = m_dof_map.begin(), it2 = _other.m_dof_map.begin();
    for(; it1 != m_dof_map.end(); ++it1, ++it2)
      if(!DoFIsEquivalent(it1->second)(it2->second))
        return false;


    return true;
  }

  //*******************************************************************************************

  Index Site::_type_ID() const {
    if(!valid_index(m_type_ID)) {
      for(m_type_ID = 0; m_type_ID < _type_prototypes().size(); m_type_ID++) {
        if(_compare_type_no_ID(_type_prototypes()[m_type_ID]))
          return m_type_ID;
      }
      //std::cout << "NEW TYPE PROTOTYPE!\n";
      //print_occ(std::cout);
      //std::cout << " : TYPE_ID-> " << m_type_ID << "\n";
      _type_prototypes().push_back(*this);
    }
    return m_type_ID;
  }

  //****************************************************

  Site operator*(const SymOp &LHS, const Site &RHS) {
    return Site(RHS).apply_sym(LHS);
  }

  //****************************************************

  Site operator+(const Site &LHS, const Coordinate &RHS) {
    return Site(LHS) += RHS;
  }

  //****************************************************

  Site operator+(const Coordinate &LHS, const Site &RHS) {
    return Site(RHS) += LHS;
  }



};
