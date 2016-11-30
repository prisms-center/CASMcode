#ifndef CASM_ScelEnum_impl
#define CASM_ScelEnum_impl

#include "casm/clex/ScelEnum.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  template<bool IsConst>
  const std::string CASM_TMP::traits<ScelEnumByNameT<IsConst> >::name = "ScelEnumByName";

  /// \brief Construct with PrimClex and ScelEnumProps settings
  template<bool IsConst>
  ScelEnumByNameT<IsConst>::ScelEnumByNameT(
    PrimClex &primclex,
    std::initializer_list<std::string> scelnames) :
    RandomAccessEnumeratorBase<Supercell, IsConst>(scelnames.size()),
    m_primclex(&primclex) {

    m_scelptr.reserve(scelnames.size());
    for(auto it = scelnames.begin(); it != scelnames.end(); ++it) {
      m_scelptr.push_back(&m_primclex->get_supercell(*it));
    }
    _init();
  }

  /// \brief Construct with PrimClex and ScelEnumProps settings
  template<bool IsConst>
  template<typename ScelNameIterator>
  ScelEnumByNameT<IsConst>::ScelEnumByNameT(
    PrimClex &primclex,
    ScelNameIterator begin,
    ScelNameIterator end) :
    RandomAccessEnumeratorBase<Supercell, IsConst>(),
    m_primclex(&primclex) {

    for(auto it = begin; it != end; ++it) {
      m_scelptr.push_back(&m_primclex->get_supercell(*it));
    }
    _init();
  }

  /// \brief Construct with PrimClex and array of supercell names
  template<bool IsConst>
  ScelEnumByNameT<IsConst>::ScelEnumByNameT(
    PrimClex &primclex,
    const jsonParser &input) :
    RandomAccessEnumeratorBase<Supercell, IsConst>(input.size()),
    m_primclex(&primclex) {

    m_scelptr.reserve(input.size());
    for(auto it = input.begin(); it != input.end(); ++it) {
      m_scelptr.push_back(&m_primclex->get_supercell(it->get<std::string>()));
    }
    _init();
  }

  /// Random access implementation
  template<bool IsConst>
  Supercell *ScelEnumByNameT<IsConst>::at_step(step_type n) {
    return m_scelptr[n];
  }

  template<bool IsConst>
  void ScelEnumByNameT<IsConst>::_init() {

    this->_set_size(m_scelptr.size());
    if(this->size() > 0) {
      this->_initialize(m_scelptr[0]);
    }
    else {
      this->_invalidate();
    }
  }


  template<bool IsConst>
  const std::string CASM_TMP::traits<ScelEnumByPropsT<IsConst> >::name = "ScelEnumByProps";

  /// \brief Construct with PrimClex and ScelEnumProps settings
  template<bool IsConst>
  ScelEnumByPropsT<IsConst>::ScelEnumByPropsT(PrimClex &primclex, const ScelEnumProps &enum_props, bool existing_only) :
    m_primclex(&primclex),
    m_existing_only(existing_only) {

    m_lattice_enum.reset(new SupercellEnumerator<Lattice>(
                           m_primclex->get_prim().lattice(),
                           m_primclex->get_prim().factor_group(),
                           enum_props
                         ));

    m_lat_it = m_lattice_enum->begin();
    m_lat_end = m_lattice_enum->end();

    while(!_include(*m_lat_it) && m_lat_it != m_lat_end) {
      std::cout << "here 0" << std::endl;
      ++m_lat_it;
    }

    if(m_lat_it != m_lat_end) {
      this->_initialize(&m_primclex->get_supercell(m_primclex->add_supercell(*m_lat_it)));
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
  template<bool IsConst>
  ScelEnumByPropsT<IsConst>::ScelEnumByPropsT(PrimClex &primclex, const jsonParser &input) :
    ScelEnumByPropsT(primclex, make_scel_enum_props(primclex, input), _get_else(input, "existing_only", false)) {}

  /// Implements increment over supercells
  template<bool IsConst>
  void ScelEnumByPropsT<IsConst>::increment() {
    ++m_lat_it;

    while(!_include(*m_lat_it) && m_lat_it != m_lat_end) {
      std::cout << "here 1" << std::endl;
      ++m_lat_it;
    }

    if(m_lat_it != m_lat_end) {
      this->_set_current_ptr(&m_primclex->get_supercell(m_primclex->add_supercell(*m_lat_it)));
      this->_increment_step();
    }
    else {
      this->_invalidate();
    }
  }

  /// Check for existing supercells
  template<bool IsConst>
  bool ScelEnumByPropsT<IsConst>::_include(const Lattice &lat) const {
    if(m_existing_only) {
      Lattice canon_lat = canonical_equivalent_lattice(
                            *m_lat_it,
                            m_primclex->get_prim().point_group(),
                            m_primclex->crystallography_tol());
      Supercell tmp(m_primclex, canon_lat);
      return m_primclex->contains_supercell(tmp);
    }
    return true;
  }

  template<bool IsConst>
  const std::string CASM_TMP::traits<ScelEnumT<IsConst> >::name = "ScelEnum";

  template<bool IsConst>
  const std::string CASM_TMP::traits<ScelEnumT<IsConst> >::help =

    "ScelEnum: \n\n"

    "  min: int, >0 (default=1)\n"
    "    The minimum volume supercell to enumerate. The volume is measured\n"
    "    relative the unit cell being used to generate supercells.\n"
    "\n"
    "  max: int, >= min (default=max existing scel_size)\n"
    "    The maximum volume supercell to enumerate. The volume is measured\n"
    "    relative the unit cell being used to generate supercells.\n"
    "\n"
    "  existing_only: bool (default=false)\n"
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
    "  name: JSON array of string (optional)\n"
    "    As an alternative to the above options, an array of existing supercell\n"
    "    names to explicitly indicate which supercells to act on. If this is \n"
    "    included, other properties are ignored. This is useful as an input \n"
    "    to other enumeration methods, such as ConfigEnumAllOccupations when \n"
    "    supercells have already been enumerated. \n"
    "\n"
    "Examples:\n"
    "\n"
    "    To enumerate supercells up to and including size 4:\n"
    "      casm enum --method SuperConfigEnum -i '{\"max\": 4}' \n"
    "\n"
    "    To enumerate 2d supercells up to and including size 4:\n"
    "      casm enum --method SuperConfigEnum -i '{\"max\": 4, \"dirs\": \"ab\"}' \n"
    "\n"
    "    If the prim is primitive FCC, two dimensional supercells of the \n"
    "    conventional FCC unit cell up to and including 4x the unit cell volume\n"
    "    could be enumerated using:\n"
    "\n"
    "     casm enum --method SuperConfigEnum -i \n"
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


  template<bool IsConst>
  int EnumInterface<ScelEnumT<IsConst> >::run(
    PrimClex &primclex,
    const jsonParser &kwargs,
    const Completer::EnumOption &enum_opt) const {

    Log &log = primclex.log();
    log.begin(name());

    auto &supercell_list = primclex.get_supercell_list();
    Index list_size = supercell_list.size();

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
    for(auto &scel : scel_enum) {
      if(verbose) {
        if(supercell_list.size() != list_size) {
          log << "  Generated: " << scel.get_name() << "\n";
        }
        else {
          log << "  Generated: " << scel.get_name() << " (already existed)\n";
        }
      }
      list_size = supercell_list.size();
    }
    log << "  DONE." << std::endl << std::endl;

    log << "Write SCEL..." << std::endl;
    primclex.print_supercells();
    log << "  DONE" << std::endl << std::endl;

    log << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    log << "  DONE" << std::endl;

    return 0;
  }

  /// \brief Construct with PrimClex and JSON settings
  ///
  /// \see EnumInterface<ScelEnumT<IsConst> >::run
  template<bool IsConst>
  ScelEnumT<IsConst>::ScelEnumT(PrimClex &primclex, const jsonParser &input) {

    if(input.contains("name")) {
      m_enum.ptr.reset(new ScelEnumByName(primclex, input["name"]));
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
  template<bool IsConst>
  void ScelEnumT<IsConst>::increment() {
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

#endif
