#ifndef SYMGROUPREPID_HH
#define SYMGROUPREPID_HH

#include "casm/casm_io/jsonParser.hh"
#include <iostream>

namespace CASM {

  /** \ingroup SymGroup
   *  @{
   */

  /// \brief Type-safe ID object for communicating and accessing Symmetry representation info

  /// SymGroupRepIDs can be returned by routines that creat a SymGroupRep inside the MasterSymGroup
  /// and then be shared among all objects that transform according to that representation.
  /// The representation can be accessed from the MasterSymGroup or from any SymOpRepresentation that
  /// belongs to the MasterSymGroup
  class SymGroupRepID {
  public:
    /// \brief Construct from group index (i.e., MasterSymGroup::group_index()) and representation index
    /// This constructor is typically only used by SymGroup class
    /// throws if _group_index or _rep_index are out-of-bounds (less than zero, for signed Index)
    SymGroupRepID(Index _group_index, Index _rep_index) :
      m_group_index(_group_index),
      m_rep_index(_rep_index) {
      if(!(valid_index(_group_index) && valid_index(_rep_index))) {
        throw std::runtime_error(std::string("SymGroupRepID initialized with out-of-bounds values!\n"));
      }
    }

    /// \brief Default constructor initializes group_index and rep_index as out-of-bounds
    /// SymGroupRepID::empty() returns true for default-constructed object
    SymGroupRepID() : SymGroupRepID(-1, -1, true) {}

    /// \brief Static function to construct an ID for identity representations.
    /// @param dim is dimension of vector space
    static SymGroupRepID identity(Index dim) {
      return SymGroupRepID(-2, dim, true);
    }

    /// \brief Index of master group in which the corresponding SymGroupRep is stored
    /// Used internally to MasterSymGroup to verify provenance of the SymGroupRepID
    Index group_index() const {
      return m_group_index;
    }

    /// \brief Index of SymGroupRep within the master group
    /// Used internally to MasterSymGroup to access the correct representation
    Index rep_index() const {
      return m_rep_index;
    }

    /// \brief Returns true if SymGroupRepID corresponds to an Identity representation
    bool is_identity() const {
      return group_index() == Index(-2);
    }

    /// \brief Returns true if SymGroupRepID has not been initialized with valid group_index or rep_index
    bool empty() const {
      return group_index() == Index(-1) || rep_index() == Index(-1);
    }

    /// \brief Output internal state to JSON
    jsonParser const &from_json(jsonParser const &json) {
      json.get_else(m_group_index, "group_index", Index(-1));
      json.get_else(m_rep_index, "group_index", Index(-1));
      return json;
    }

  private:
    /// Private constructor skips bounds checks -- Used for identity() and default construction
    SymGroupRepID(Index _group_index, Index _rep_index, bool override) :
      m_group_index(_group_index),
      m_rep_index(_rep_index) {

    }

    Index m_group_index;
    Index m_rep_index;
  };

  /// Compares true if group_index() and rep_index() are equal
  inline
  bool operator==(SymGroupRepID const &a, SymGroupRepID const &b) {
    return a.group_index() == b.group_index() && a.rep_index() == b.rep_index();
  }

  /// Compares false if group_index() and rep_index() are equal
  inline
  bool operator!=(SymGroupRepID const &a, SymGroupRepID const &b) {
    return !(a == b);
  }

  /// Less-than comparison for use in STL containers (std::set, std::map, etc)
  inline
  bool operator<(SymGroupRepID const &a, SymGroupRepID const &b) {
    return (a.group_index() < b.group_index()) ||
           (a.group_index() == b.group_index() && a.rep_index() < b.rep_index());
  }

  inline
  jsonParser &to_json(SymGroupRepID const &_id, jsonParser &json) {
    json["group_index"] = _id.group_index();
    json["rep_index"] = _id.rep_index();
    return json;
  }

  inline
  jsonParser const &from_json(SymGroupRepID &_id, jsonParser const &json) {
    return _id.from_json(json);
  }

  inline
  std::ostream &operator<<(std::ostream &out, SymGroupRepID const &_id) {
    out << "{group_index = " << _id.group_index() << ", rep_index = " << _id.rep_index() << "}";
    return out;
  }

  /** @} */
}
#endif
