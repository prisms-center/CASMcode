#ifndef SYMGROUPREP_HH
#define SYMGROUPREP_HH

#include <iostream>
#include <string>
#include <iomanip>

#include "casm/container/Array.hh"
#include "casm/container/multivector.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"


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
  class SymGroupRep : public Array<SymOpRepresentation *> {
  public:
    typedef SymGroupRepHandle RemoteHandle;
    enum NullInitializer {NO_HOME};

    /// Use this constructor when MasterSymGroup is unknown or doesn't exist
    /// You must promise that you know what you're doing

    SymGroupRep(SymGroupRep::NullInitializer init, Index _size) :
      Array<SymOpRepresentation * >(_size, nullptr),
      m_master_group(nullptr) { }

    SymGroupRep(const SymGroup &_head, SymGroupRepID _rep_ID = SymGroupRepID()) :
      m_rep_ID(_rep_ID), m_master_group(&_head.master_group()) {
      if(has_valid_master())
        resize(master_group().size(), nullptr);
      else {
        std::cerr << "WARNING: SymGroupRep initialized without valid master_group!  Bad things may happen!\n";
        assert(0);
      }
    }

    SymGroupRep(const SymGroupRep &RHS);

    ~SymGroupRep();

    //void push_back(const SymOpRepresentation &new_op);

    /// SPECIAL METHOD: const, but it modifies symop representation 'i',
    /// after performing appropriate checks. Do not refer to this as an example
    /// of how to write mutators.
    /// Critical Failure if symop representation 'i' is already initialized
    void set_rep(const SymOpRepresentation &base_rep, const SymOpRepresentation &new_rep) const;
    void set_rep(const SymOp &base_rep, const SymOpRepresentation &new_rep) const;
    void set_rep(Index op_index, const SymOpRepresentation &new_rep) const;

    const MasterSymGroup &master_group() const {
      assert(m_master_group);
      return *m_master_group;
    }

    bool has_valid_master() const {
      return m_master_group != nullptr;
    }
    void set_master_group(const MasterSymGroup &master, const SymGroupRepID &_rep_ID);

    /// Adds copy of this representation its home_group
    SymGroupRepID add_copy_to_master() const;

    SymGroupRep &operator=(const SymGroupRep &RHS);

    void clear();

    SymGroupRep *copy()const {
      return new SymGroupRep(*this);
    }

    //John G 050513

    void print_permutation(std::ostream &stream) const;
    void print_MatrixXd(std::ostream &stream) const;
    void print_MatrixXd(std::ostream &stream, const SymGroup &subgroup) const;

    //\John G

    // block_shape_matrix is sum of squares of each (i,j) matrix element over all operations in SymGroupRep
    // It reveals the block_diagonalization of the symgrouprep
    Eigen::MatrixXd block_shape_matrix() const;
    Eigen::MatrixXd block_shape_matrix(const SymGroup &subgroup) const;

    // counts number of blocks in the block_shape_matrix
    Index num_blocks() const;
    Index num_blocks(const SymGroup &subgroup) const;

    Eigen::MatrixXd const *get_MatrixXd(Index i) const;
    Eigen::MatrixXd const *get_MatrixXd(const SymOpRepresentation &) const;
    Permutation const *get_permutation(Index i) const;
    Permutation const *get_permutation(const SymOpRepresentation &) const;

    // adds a copy of '_pushed' temporary workaround for 'masterless' SymGroupReps; to be deprecated
    // will throw if valid master exists.
    void push_back_copy(const SymOpRepresentation &_pushed);

    SymGroupRepID get_ID() const {
      return m_rep_ID;
    }

    std::vector<Eigen::MatrixXd> irreducible_wedges(const SymGroup &head_group, std::vector<Index> &multiplicities)const;
    multivector<Eigen::VectorXd>::X<3> calc_special_total_directions(const SymGroup &subgroup)const;
    ReturnArray<Array< Eigen::MatrixXd> > calc_special_subspaces(const SymGroup &subgroup)const;

    ReturnArray<Index> num_each_irrep() const;
    ReturnArray<Index> num_each_irrep(const SymGroup &sub_group, bool verbose = false) const;
    ReturnArray<Index> num_each_real_irrep(const SymGroup &subgroup, bool verbose = false) const;

    ReturnArray<SymGroupRepID> get_irrep_IDs(const SymGroup &subgroup) const;

    bool is_irrep() const;
    bool is_irrep(const SymGroup &head_group) const;

    /// Make a copy of representation on vector space 'V' that is transformed into a representation on vector space 'W'
    /// 'trans_mat' is the unitary matrix that isomorphically maps 'V'->'W' (i.e., [w = trans_mat * v] and [v = trans_mat.transpose() * w] )
    /// If the original representation to be transformed is just a temporary standalone SymGroupRep, be sure to delete it before falling out of scope
    SymGroupRep coord_transformed_copy(const Eigen::MatrixXd &trans_mat) const;

    /// Make copy of (*this) that is transformed so that axes are oriented along high-symmetry direction
    /// and confined to subspaces that transform as irreps.
    SymGroupRep symmetry_adapted_copy(const SymGroup &head_group) const {
      return coord_transformed_copy(get_irrep_trans_mat(head_group));
    }

    /// Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
    /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
    /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
    Eigen::MatrixXd get_irrep_trans_mat(const SymGroup &head_group) const;
    /// Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
    /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
    /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
    /// Also populate 'subspaces' with lists of columns that form irreps
    Eigen::MatrixXd get_irrep_trans_mat(const SymGroup &head_group, std::vector<std::vector<Index> > &subspaces) const;

    Eigen::MatrixXd get_irrep_trans_mat_blind(const SymGroup &head_group) const;
    ReturnArray<Eigen::MatrixXd> get_projection_operators() const;

    jsonParser &to_json(jsonParser &json) const;

    // If 'm_master_group' is not nullptr, should be initialized accordingly
    void from_json(const jsonParser &json);
  private:

    /// rep_ID is unique identifier of a specific SymGroupRep instantiation
    mutable SymGroupRepID m_rep_ID;

    /// Pointer to the home_group that generated this SymGroupRep
    MasterSymGroup const *m_master_group;
    // don't use default constructor
    SymGroupRep() {
      std::cerr << "Cannot perform default construction of SymGroupRep.\nExiting...\n";
      exit(1);
    }

    using Array<SymOpRepresentation *>::push_back;
    using Array<SymOpRepresentation *>::resize;

    /// Pointer version of constructor is private for internal construction of master-less representations
    SymGroupRep(const MasterSymGroup *_home) :  m_master_group(_home) { }

    void calc_new_irreps(int max_iter = 1000) const;
    void calc_new_irreps(const SymGroup &sub_group, int max_iter = 1000) const;

    ReturnArray<Array<Eigen::VectorXd> > _calc_special_irrep_directions(const SymGroup &subgroup)const;

    /// Find a new coordinate system oriented along high-symmetry directions in vector space 'V' as determined by
    /// the subset of SymOpRepresentations specified by 'subgroup'.
    /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
    /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
    Eigen::MatrixXd _symmetrized_irrep_trans_mat(const SymGroup &subgroup) const;

  };

  SymGroupRep subset_permutation_rep(const SymGroupRep &permute_rep, const Array<Index>::X2 &subsets);
  SymGroupRep permuted_direct_sum_rep(const SymGroupRep &permute_rep, const Array<SymGroupRep const *> &sum_reps);
  SymGroupRep kron_rep(const SymGroupRep &LHS, const SymGroupRep &RHS);

  jsonParser &to_json(const SymGroupRep &rep, jsonParser &json);

  // If 'm_home_group' is not nullptr, should be initialized accordingly
  void from_json(SymGroupRep &rep, const jsonParser &json);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// SymGroupRepHandle is used to provide easy remote access to a SymGroup representation
  /// The head group may be a subgroup of the MasterSymGroup where the SymGroupRep is stored.
  class SymGroupRepHandle { // <-- typedefed as SymGroupRep::RemoteHandle
    SymGroupRep const *m_group_rep;
    Array<Index> m_subgroup_op_inds;
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

    /// Matrix dimension of representation
    Index dim() const {
      return (m_group_rep->get_MatrixXd(m_subgroup_op_inds[0]))->cols();
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

    const SymOp &sym_op(Index i) const {
      return (m_group_rep->master_group())[m_subgroup_op_inds[i]];
    }

    Index ind_inverse(Index i) const {
      return m_subgroup_op_inds.find((m_group_rep->master_group()).ind_inverse(m_subgroup_op_inds[i]));
    }

    Index ind_prod(Index i, Index j) const {
      return m_subgroup_op_inds.find((m_group_rep->master_group()).ind_prod(m_subgroup_op_inds[i], m_subgroup_op_inds[j]));
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
