#include "casm/clex/ScelEnum_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/app/enum.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ScelEnum_interface() {
    return new CASM::EnumInterface<CASM::ScelEnum>();
  }
}

namespace CASM {

  const std::string ScelEnumByName::enumerator_name = "ScelEnumByName";

  /// \brief Construct with PrimClex a initializer_list of Supercell names
  ///
  /// \param primclex A PrimClex for which to enumerate Supercells
  /// \param scelnames A list of names of Supercells to enumerate
  ///
  ScelEnumByName::ScelEnumByName(
    const PrimClex &primclex,
    std::initializer_list<std::string> scelnames) :
    RandomAccessEnumeratorBase<Supercell>(scelnames.size()),
    m_primclex(&primclex) {

    m_scelptr.reserve(scelnames.size());
    for(auto it = scelnames.begin(); it != scelnames.end(); ++it) {
      m_scelptr.push_back(&*m_primclex->db<Supercell>().find(*it));
    }
    _init();
  }

  /// \brief Construct with PrimClex and JSON array containing supercell names
  ///
  /// \param primclex A PrimClex for which to enumerate Supercells
  /// \param input A JSON array of names of Supercells to enumerate
  ///
  ScelEnumByName::ScelEnumByName(
    const PrimClex &primclex,
    const jsonParser &input) :
    RandomAccessEnumeratorBase<Supercell>(input.size()),
    m_primclex(&primclex) {

    m_scelptr.reserve(input.size());
    for(auto it = input.begin(); it != input.end(); ++it) {
      m_scelptr.push_back(&*m_primclex->db<Supercell>().find(it->get<std::string>()));
    }
    _init();
  }

  std::string ScelEnumByName::name() const {
    return ScelEnumByName::enumerator_name;
  }

  /// Random access implementation
  const Supercell *ScelEnumByName::at_step(step_type n) {
    return m_scelptr[n];
  }

  void ScelEnumByName::_init() {

    this->_set_size(m_scelptr.size());
    if(this->size() > 0) {
      this->_initialize(m_scelptr[0]);
    }
    else {
      this->_invalidate();
    }
  }


  const std::string ScelEnumByProps::enumerator_name = "ScelEnumByProps";

  /// \brief Construct with PrimClex and ScelEnumProps settings
  ///
  /// \param primclex A PrimClex for which to enumerate Supercells
  /// \param enum_props Specifies which Supercells to enumerate
  /// \param existing_only Skip Supercells specified by enum_props but not already
  ///        existing in the Supercell list
  ///
  ScelEnumByProps::ScelEnumByProps(const PrimClex &primclex, const ScelEnumProps &enum_props, bool existing_only) :
    m_primclex(&primclex),
    m_existing_only(existing_only) {
    m_lattice_enum.reset(new SupercellEnumerator<Lattice>(
                           m_primclex->prim().lattice(),
                           m_primclex->prim().point_group(),
                           enum_props
                         ));

    m_lat_it = m_lattice_enum->begin();
    m_lat_end = m_lattice_enum->end();

    while(!_include(*m_lat_it) && m_lat_it != m_lat_end) {
      ++m_lat_it;
    }

    if(m_lat_it != m_lat_end) {
      Supercell scel(m_primclex, *m_lat_it);
      this->_initialize(&*scel.insert().first);
    }
    else {
      this->_invalidate();
    }
  }

  namespace {

    bool _get_else(const jsonParser &json, std::string key, bool default_value) {
      bool tmp;
      json.get_else<bool>(tmp, key, default_value);
      return tmp;
    }
  }

  /// \brief Construct with PrimClex and ScelEnumProps JSON settings
  ///
  /// \param primclex A PrimClex for which to enumerate Supercells
  /// \param input JSON used to make an ScelEnumProps object specifying which
  ///        Supercells to enumerate, via ::make_scel_enum_props
  ///
  /// - The JSON input is also checked for the boolean property "existing_only",
  ///   which if true indicates that Supercell not already included in the
  ///   supercell list should be skipped
  ///
  ScelEnumByProps::ScelEnumByProps(const PrimClex &primclex, const jsonParser &input) :
    ScelEnumByProps(primclex, make_scel_enum_props(primclex, input), _get_else(input, "existing_only", false)) {}

  std::string ScelEnumByProps::name() const {
    return ScelEnumByProps::enumerator_name;
  }

  /// Implements increment over supercells
  void ScelEnumByProps::increment() {
    ++m_lat_it;

    while(!_include(*m_lat_it) && m_lat_it != m_lat_end) {
      ++m_lat_it;
    }

    if(m_lat_it != m_lat_end) {
      Supercell scel(m_primclex, *m_lat_it);

      this->_set_current_ptr(&*scel.insert().first);
      this->_increment_step();
    }
    else {
      this->_invalidate();
    }
  }

  /// Check for existing supercells
  bool ScelEnumByProps::_include(const Lattice &lat) const {
    if(m_existing_only) {
      std::string name = canonical_scelname(m_primclex->prim(), *m_lat_it);
      return m_primclex->db<Supercell>().find(name) != m_primclex->db<Supercell>().end();
    }
    return true;
  }

  std::string ScelEnum::name() const {
    return ScelEnum::enumerator_name;
  }

  /// \relates ::ScelEnumT
  const std::string ScelEnum::enumerator_name = "ScelEnum";

  /// \relates ::ScelEnumT
  std::string ScelEnum::interface_help() {
    return

      "ScelEnum: \n\n"

      "  min: int, >0 (default=1)\n"
      "    The minimum volume supercell to enumerate. The volume is measured\n"
      "    relative the unit cell being used to generate supercells.\n"
      "\n"
      "  max: int, >= min (default=max existing scel_size)\n"
      "    The maximum volume supercell to enumerate. The volume is measured\n"
      "    relative the unit cell being used to generate supercells.\n"
      "\n"
      "  existing_only: bool (default=true)\n"
      "    If true, only existing supercells are used. This is useful when it\n"
      "    is used as input to a Configuration enumeration method.\n"
      "\n"
      "  dirs: string (default=\"abc\")\n"
      "    This option may be used to restrict the supercell enumeration to 1, \n"
      "    2 or 3 of the lattice vectors, to get 1-, 2-, or 3-dimensional      \n"
      "    supercells. By specifying combinations of 'a', 'b', and 'c', you    \n"
      "    determine which of the unit cell lattice vectors you want to        \n"
      "    enumerate over. For example, to enumerate 1-dimensional supercells  \n"
      "    along the 'c' use \"dirs\":\"c\". If you want 2-dimensional        \n"
      "    supercells along the 'a' and 'c' lattice vectors, specify           \n"
      "    \"dirs\":\"ac\". \n"
      "\n"
      "  unit_cell: 3x3 matrix of int, or string (default=identity matrix)     \n"
      "    This option may be used to specify the unit cell. It may be         \n"
      "    specified using a 3x3 matrix of int, representing the transformation\n"
      "    matrix, T, such that U = P*T, where P are the primitive lattice     \n"
      "    and U are the unit cell lattice vectors. For example, a unit cell   \n"
      "    that whose lattice vectors are (2*a+b, b, c) (with respect to the   \n"
      "    the primitive cell vectors) could be specified using:\n"
      "\n"
      "      \"unit_cell\" : [\n"
      "        [2, 0, 0],\n"
      "        [1, 1, 0],\n"
      "        [0, 0, 1]\n"
      "       ]\n"
      "\n"
      "    Or it may be specified by  \n"
      "    the name of the existing supercell to use as the unit cell, for     \n"
      "    example: \n"
      "\n"
      "      \"unit_cell\" : \"SCEL2_1_1_2_0_0_0\"\n"
      "\n"
      "  names: JSON array of string (optional)\n"
      "    As an alternative to the above options, an array of existing supercell\n"
      "    names to explicitly indicate which supercells to act on. If this is \n"
      "    included, other properties are ignored. This is useful as an input \n"
      "    to other enumeration methods, such as ConfigEnumAllOccupations when \n"
      "    supercells have already been enumerated. \n"
      "\n"
      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n"
      "\n"
      "Examples:\n"
      "\n"
      "    To enumerate supercells up to and including size 4:\n"
      "      casm enum --method ScelEnum -i '{\"max\": 4}' \n"
      "\n"
      "    To enumerate 2d supercells up to and including size 4:\n"
      "      casm enum --method ScelEnum -i '{\"max\": 4, \"dirs\": \"ab\"}' \n"
      "\n"
      "    If the prim is primitive FCC, two dimensional supercells of the \n"
      "    conventional FCC unit cell up to and including 4x the unit cell volume\n"
      "    could be enumerated using:\n"
      "\n"
      "     casm enum --method ScelEnum -i \n"
      "     '{\n"
      "        \"min\": 1,\n"
      "        \"max\": 4,\n"
      "        \"dirs\": \"ab\",\n"
      "        \"unit_cell\" : [\n"
      "          [-1,  1,  1],\n"
      "          [ 1, -1,  1],\n"
      "          [ 1,  1, -1]\n"
      "        ]\n"
      "      }'\n"
      "\n";
  }

  /// \relates ::ScelEnumT
  int ScelEnum::run(
    const PrimClex &primclex,
    const jsonParser &kwargs,
    const Completer::EnumOption &enum_opt,
    EnumeratorMap const *interface_map) {

    Log &log = primclex.log();
    log.begin(enumerator_name);

    Index list_size = primclex.db<Supercell>().size();

    bool verbose = true;

    jsonParser input {kwargs};
    // check supercell shortcuts
    if(enum_opt.vm().count("min")) {
      input["min"] = enum_opt.min_volume();
    }
    if(enum_opt.vm().count("max")) {
      input["max"] = enum_opt.max_volume();
    }

    ScelEnum scel_enum(primclex, input);

    bool dry_run = CASM::dry_run(kwargs, enum_opt);
    std::string dry_run_msg = CASM::dry_run_msg(dry_run);
    for(auto &scel : scel_enum) {
      if(verbose) {
        if(primclex.db<Supercell>().size() != list_size) {
          log << dry_run_msg << "  Generated: " << scel.name() << "\n";
        }
        else {
          log << dry_run_msg << "  Generated: " << scel.name() << " (already existed)\n";
        }
      }
      list_size = primclex.db<Supercell>().size();
    }
    log << dry_run_msg << "  DONE." << std::endl << std::endl;

    if(!dry_run) {
      log << "Write supercells..." << std::endl;
      primclex.db<Supercell>().commit();
      log << "  DONE" << std::endl << std::endl;
    }
    return 0;
  }

  /// \brief Construct with PrimClex and JSON settings
  ///
  /// \see EnumInterface<ScelEnum>::run
  ScelEnum::ScelEnum(const PrimClex &primclex, const jsonParser &input) {
    if(input.contains("names")) {
      m_enum.ptr.reset(new ScelEnumByName(primclex, input["names"]));
    }
    else {
      m_enum.ptr.reset(new ScelEnumByProps(primclex, input));
    }

    m_it = m_enum.begin();
    m_end = m_enum.end();

    if(m_it != m_end) {
      this->_initialize(&(*m_it));
      this->_set_step(0);
    }
    else {
      this->_invalidate();
    }
  }


  /// Implements increment over all occupations
  void ScelEnum::increment() {
    ++m_it;
    if(m_it != m_end) {
      this->_set_current_ptr(&(*m_it));
      this->_increment_step();
    }
    else {
      this->_invalidate();
    }
  }

}
