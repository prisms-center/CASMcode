#ifndef SYMGROUPREP_HH
#define SYMGROUPREP_HH

#include <iostream>
#include <string>
#include <iomanip>

#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/container/multivector.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
  class SymGroupRep;
  class MasterSymGroup;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class SymGroupRepHandle;

  /** \ingroup SymGroup
   *  @{
   */

  ///\brief SymGroupRep is an alternative representation of a SymGroup for something other than real space.
  /// There is a one-to-one correspondence of SymOps in some SymGroup with the SymOpRepresentations in SymGroupRep
  /// SymGroupRep does not know or care about the specifics of what the SymOpRepresentations describe
  /// or how they are implemented
  class SymGroupRep : public std::vector<SymOpRepresentation *> {
  public:
    friend class SymGroup;

    typedef SymGroupRepHandle RemoteHandle;

    /// SymGroupRep::NO_HOME is safeguard to prevent usage of this constructor when full constructor should be used
    enum NullInitializer {NO_HOME};

    // No default construction allowed
    SymGroupRep() = delete;


    SymGroupRep(const SymGroup &_head, SymGroupRepID _rep_ID = SymGroupRepID()) :
      m_rep_ID(_rep_ID), m_master_group(&_head.master_group()) {
      if(has_valid_master())
        resize(master_group().size(), nullptr);
      else {
        std::cerr << "WARNING: SymGroupRep initialized without valid master_group!  Bad things may happen!\n";
        assert(0);
      }
    }

    /// \brief Standalone constructor
    /// Use this constructor when MasterSymGroup is unknown or doesn't exist
    /// Ex: SymGroupRep my_rep(SymGroupRep::NO_HOME, 48);
    SymGroupRep(SymGroupRep::NullInitializer init, Index _size) :
      std::vector<SymOpRepresentation * >(_size, nullptr),
      m_master_group(nullptr) { }

    /// \brief explicit copy constructor
    /// Necessary because SymGroupRep directly manages heap-allocated SymOpRepresentations
    SymGroupRep(const SymGroupRep &RHS);

    /// \brief explicit destructor
    /// Necessary because SymGroupRep directly manages heap-allocated SymOpRepresentations
    ~SymGroupRep();

    /// \brief Reference to MasterSymGroup for which this SymGroupRep is a group representation
    /// If a MasterSymGroup contains this SymGroupRep, it must be the master_group() of this SymGroupRep
    /// A SymGroupRep that is not contained in a MasterSymGroup may (and often does) still have a valid master_group()
    const MasterSymGroup &master_group() const {
      assert(m_master_group);
      return *m_master_group;
    }

    /// \brief Returns true if this SymGroupRep has valid pointer to a MasterSymGroup
    bool has_valid_master() const {
      return m_master_group != nullptr;
    }

    /// \brief Assign ma
    void set_master_group(const MasterSymGroup &master, const SymGroupRepID &_rep_ID);

    /// \brief Sets the representation of operation at entry 'op_index' of this group representation
    /// Throws if this group representation already contains entry at 'op_index'
    /// \param op_index Index of group operation to be set
    /// \param new_rep New representation to be recorded in this SymGrouprep for operation at entry 'op_index' of master_group()
    void set_rep(Index op_index, const SymOpRepresentation &new_rep);

    /// \brief Sets the representation of operation at entry 'base_rep.index()' of this group representation
    /// Equivalent to this->set_rep(base_rep.index(),new_rep)
    /// Throws if this group representation already contains entry at 'base_rep.index()'
    /// \param base_rep Fully initialized SymOpRepresentation that has valid 'base_rep.index()' value
    /// \param new_rep New representation to be recorded in this SymGroupRep for operation at entry 'base_rep.endex()' of master_group()
    void set_rep(const SymOpRepresentation &base_rep, const SymOpRepresentation &new_rep);

    /// \brief Sets the representation of operation at entry 'base_rep.index()' of this group representation
    /// Throws if this group representation already contains entry at 'base_rep.index()'
    /// If base_rep.index() is invalid, it will attempt to find base_rep in the master_group() to resolve index
    /// \param base_rep SymOp corresponding to an operation in master_group()
    /// \param new_rep New representation to be recorded in this SymGroupRep for operation at entry 'base_rep.endex()' of master_group()
    void set_rep(const SymOp &base_rep, const SymOpRepresentation &new_rep);

    /// \brief explicit assignment operator
    /// Necessary because SymGroupRep directly manages heap-allocated SymOpRepresentations
    SymGroupRep &operator=(const SymGroupRep &RHS);

    /// \brief Returns pointer to unmanaged copy of this SymGroupRep
    SymGroupRep *copy()const {
      return new SymGroupRep(*this);
    }

    /// \brief pointer to MatrixXd corresponding to SymOpRepresentation at entry 'i' of this SymGroupRep
    /// Returns nullptr if entry 'i' is uninitialized or if entry 'i' cannot be represented as a Matrix
    Eigen::MatrixXd const *MatrixXd(Index i) const;

    /// \brief pointer to MatrixXd corresponding to SymOpRepresentation at entry '_op.index()' of this SymGroupRep
    /// Returns nullptr if entry '_op.index()' is uninitialized or if entry '_op.index()' cannot be represented as a Matrix
    Eigen::MatrixXd const *MatrixXd(const SymOpRepresentation &_op) const;

    /// \brief pointer to Permutation corresponding to SymOpRepresentation at entry 'i' of this SymGroupRep
    /// Returns nullptr if entry 'i' is uninitialized or if entry 'i' cannot be represented as a Permutation
    Permutation const *permutation(Index i) const;

    /// \brief pointer to Permutation corresponding to SymOpRepresentation at entry '_op.index()' of this SymGroupRep
    /// Returns nullptr if entry '_op.index()' is uninitialized or if entry '_op.index()' cannot be represented as a Permutation
    Permutation const *permutation(const SymOpRepresentation &_op) const;

    /// \brief Returns SymGroupRepID of this SymGroupRep in its master_group. If this is not contained in a master_group, symrep_ID().empty()==true
    SymGroupRepID symrep_ID() const {
      return m_rep_ID;
    }

    jsonParser &to_json(jsonParser &json) const;

    // If 'm_master_group' is not nullptr, should be initialized accordingly
    void from_json(const jsonParser &json);
  private:

    void clear();

    /// rep_ID is unique identifier of a specific SymGroupRep instantiation
    mutable SymGroupRepID m_rep_ID;

    /// Pointer to the home_group that generated this SymGroupRep
    MasterSymGroup const *m_master_group;

    using std::vector<SymOpRepresentation *>::push_back;
    using std::vector<SymOpRepresentation *>::resize;

    /// Pointer version of constructor is private for internal construction of master-less representations
    SymGroupRep(const MasterSymGroup *_home) :  m_master_group(_home) { }


  };

  jsonParser &to_json(const SymGroupRep &rep, jsonParser &json);

  // If 'm_home_group' is not nullptr, should be initialized accordingly
  void from_json(SymGroupRep &rep, const jsonParser &json);

  /// \brief Make a copy of representation on vector space 'V' that is transformed into a representation on vector space 'W'
  /// 'trans_mat' is the unitary matrix that isomorphically maps 'V'->'W' (i.e., [w = trans_mat * v] and [v = trans_mat.transpose() * w] )
  /// If the original representation to be transformed is just a temporary standalone SymGroupRep, be sure to delete it before falling out of scope
  SymGroupRep coord_transformed_copy(SymGroupRep const &_rep, const Eigen::MatrixXd &trans_mat);



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// SymGroupRepHandle is used to provide easy remote access to a SymGroup representation
  /// The head group may be a subgroup of the MasterSymGroup where the SymGroupRep is stored.
  class SymGroupRepHandle { // <-- typedefed as SymGroupRep::RemoteHandle
    SymGroupRep const *m_group_rep;
    std::vector<Index> m_subgroup_op_inds;
  public:
    SymGroupRepHandle():
      m_group_rep(nullptr) {}

    SymGroupRepHandle(const SymGroup &head_group, SymGroupRepID symrep_ID) :
      m_group_rep(nullptr), m_subgroup_op_inds(head_group.op_indices()) {
      if(head_group.size() == 0 || !head_group[0].has_valid_master()) {
        std::cerr << "CRITICAL ERROR: Requested representation of an improperly initialized SymGroup.\n"
                  << "                Exiting...\n";
        assert(0);
      }
      assert(!symrep_ID.empty());
      m_group_rep = &(head_group[0].master_group().representation(symrep_ID));
      assert(m_group_rep);
    }

    /// Size of the associate SymGroup
    Index size() const {
      return m_subgroup_op_inds.size();
    }

    SymGroupRepID symrep_ID() const {
      return rep_ptr()->symrep_ID();
    }

    /// Matrix dimension of representation
    Index dim() const {
      return (m_group_rep->MatrixXd(m_subgroup_op_inds[0]))->cols();
    }

    bool empty() const {
      return m_group_rep == nullptr;
    }

    SymGroupRep const *rep_ptr() const {
      return m_group_rep;
    }

    SymGroupRep const *operator->() const {
      return rep_ptr();
    }

    SymOpRepresentation const *operator[](Index i) const {
      return m_group_rep->at(m_subgroup_op_inds[i]);
    }

    bool operator==(const SymGroupRepHandle &RHS) const {
      return m_group_rep == RHS.m_group_rep && m_subgroup_op_inds == RHS.m_subgroup_op_inds;
    }

    void set_rep(const SymGroup &head_group, SymGroupRepID symrep_ID) {
      m_subgroup_op_inds = head_group.op_indices();
      assert(head_group.size() && !symrep_ID.empty() && head_group[0].has_valid_master());
      m_group_rep = &(head_group[0].master_group().representation(symrep_ID));
      assert(m_group_rep);

    }
  };

  /** @} */
}

#endif
