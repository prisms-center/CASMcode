#ifndef SYMOPREPRESENTATION_HH
#define SYMOPREPRESENTATION_HH

#include <iostream>
#include <cmath>

#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"
#include "casm/symmetry/SymGroupRepID.hh"

namespace CASM {
  class MasterSymGroup;
  class Permutation;
  class UnitCellCoord;
  class SymBasisPermute;

  /** \ingroup Symmetry
   *
   *  @{
  */

  ///\brief SymOpRepresentation is the base class for anything describes a symmetry operation
  class SymOpRepresentation {
  public:
    SymOpRepresentation() : m_master_group(nullptr), m_op_index(-1) {}

    SymOpRepresentation(const MasterSymGroup &_master_group, SymGroupRepID _rep_ID, Index _op_index) :
      SymOpRepresentation(&_master_group, _rep_ID, _op_index) {}

    /// SymOpRepresentation specifies how a symmetry operation acts on a certain type of object.
    /// It is a virtual class, since we don't need to know the specifics of its behavior at higher levels of abstraction
    virtual ~SymOpRepresentation() {} // = 0;

    /// \brief Make copy of Derived object through Base-class interface
    virtual SymOpRepresentation *copy() const = 0;

    /// \brief Calculates character of the representation (if well-defined)
    virtual double character() const {
      return NAN;
    }

    virtual Permutation  const *get_permutation() const {
      return nullptr;
    }

    virtual Eigen::MatrixXd const *get_MatrixXd() const {
      return nullptr;
    }

    virtual SymBasisPermute const *get_ucc_permutation() const {
      return nullptr;
    }


    /// get pointer to matrix representation corresponding to rep_ID
    Eigen::MatrixXd const *get_matrix_rep(SymGroupRepID _rep_ID) const;

    /// get pointer to permutation representation corresponding to _rep_ID
    Permutation const *get_permutation_rep(SymGroupRepID _rep_ID) const;

    /// get pointer to BasisPermute representation corresponding to _rep_ID
    SymBasisPermute const *get_basis_permute_rep(SymGroupRepID _rep_ID) const;

    /// get array of pointers to matrix representations for representations corresponding to _rep_IDs
    Array<Eigen::MatrixXd const * > get_matrix_reps(Array<SymGroupRepID> _rep_IDs) const;

    /// set representation for SymOp corresponding to _rep_ID
    void set_rep(SymGroupRepID _rep_ID, const SymOpRepresentation &op_rep) const;

    /// Change m_master_group and determine op_index
    void set_identifiers(const MasterSymGroup &new_group, SymGroupRepID new_rep_ID);

    /// Set m_master_group, _rep_ID, and op_index
    void set_identifiers(const MasterSymGroup &new_group, SymGroupRepID new_rep_ID, Index new_op_index);

    /// const access of head group
    const MasterSymGroup &master_group() const {
      assert(m_master_group);
      return *m_master_group;
    }

    /// \brief check if this representation is registered with a MasterSymGroup
    bool has_valid_master() const {
      return m_master_group != nullptr;
    }

    void invalidate_index() {
      m_op_index = -1;
    }

    /// Index of this operation within the master_group
    Index index()const {
      return m_op_index;
    }

    /// ID of representation that this operation belongs to within the master_group
    SymGroupRepID rep_ID()const {
      return m_rep_ID;
    }

    /// Get the operation index of the inverse of this operation, using the master_group's multiplication table
    Index ind_inverse()const;

    /// Get the operation index of operation that is the product of this operation and @param RHS,
    /// using the master_group's multiplication table
    Index ind_prod(const SymOpRepresentation &RHS)const;

    virtual jsonParser &to_json(jsonParser &json) const = 0;
    virtual void from_json(const jsonParser &json) = 0;

  protected:
    /// Protected constructor to allow internal construction of masterless symops
    SymOpRepresentation(MasterSymGroup const *_master_group_ptr, SymGroupRepID _rep_ID, Index _op_index) :
      m_master_group(_master_group_ptr), m_rep_ID(_rep_ID), m_op_index(_op_index) {}

    /// Pointer to the MasterSymGroup where prototype of this SymOp lives
    MasterSymGroup const *m_master_group;

    /// ID of this representation within the master_group.  Default is uninitialized.
    SymGroupRepID m_rep_ID;

    ///Index into MasterSymGroup that specifies the operation
    Index m_op_index;
  };

  struct SymRepIndexCompare {

    SymRepIndexCompare() {}

    bool operator()(const SymOpRepresentation &A, const SymOpRepresentation &B) const {
      return A.index() < B.index();
    }
  };

  jsonParser &to_json(const SymOpRepresentation *rep, jsonParser &json);
  /// This allocates a new object to 'rep'.
  void from_json(SymOpRepresentation *rep, const jsonParser &json);

  /** @}*/

}
#endif
