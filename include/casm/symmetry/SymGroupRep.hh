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


  /// \brief Find irrep decomposition of _rep and add them to irrep listing in master_group
  void calc_new_irreps(SymGroupRep const &_rep, int max_iter = 1000);

  /// \brief Find irrep decomposition of _rep and add them to irrep listing in sub_group
  void calc_new_irreps(SymGroupRep const &_rep, const SymGroup &sub_group, int max_iter = 1000);

  /// \brief Assuming that _rep is an irrep of head_group, find high-symmetry directions
  /// throws if _rep is not an irrep
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// to a unique subgroup of 'head_group' (i.e., no other direction in the space, except the negative of that direction, is invariant to that subgroup)
  /// result[i] is an orbit of symmetrically equivalent directions, result[i][j] is an individual direction. Direction vectors are normalized to unit length.
  /// The total set of all directions is guaranteed to span the space.
  std::vector<std::vector<Eigen::VectorXd> > special_irrep_directions(SymGroupRep const &_rep, const SymGroup &head_group);

  /// Find a new coordinate system oriented along high-symmetry directions in vector space 'V' as determined by
  /// the subset of SymOpRepresentations specified by 'subgroup'.
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd symmetrized_irrep_trans_mat(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Matrix B such that B(i,j) is sum of squares of each (i,j) matrix element over all SymOpRepresentations in _rep
  /// It reveals the block_diagonalization of SymGroupRep _rep
  Eigen::MatrixXd block_shape_matrix(SymGroupRep const &_rep);

  /// \brief Matrix B such that B(i,j) is sum of squares of each (i,j) matrix element over SymOpRepresentations in _rep corresponding to operations of head_group
  /// It reveals the block_diagonalization of SymGroupRep _rep for the subset of operations contained in head_group
  Eigen::MatrixXd block_shape_matrix(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief counts number of blocks in the block_shape_matrix of _rep
  /// Reveals number of invariant subspaces that comprise the vector space supporting _rep
  Index num_blocks(SymGroupRep const &_rep);

  /// \brief counts number of blocks in the block_shape_matrix of _rep
  /// Reveals number of invariant subspaces (with respect to head_group) that comprise the vector space supporting _rep
  Index num_blocks(SymGroupRep const &_rep, const SymGroup &head_group);


  /// \brief finds high-symmetry directions within vector space supporting _rep, wrt symmetry of head_group
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// a subgroup of 'head_group'. These are constructed by finding irreps of _rep and then calling special_irrep_directions on each
  /// result[i] is the set of special directions belonging to the i'th irrep constituting _rep
  /// result[i][j] s an orbit of symmetrically equivalent directions, and result[i][j][k] is an individual direction.
  /// Direction vectors are normalized to unit length. The total set of all directions is guaranteed to span the space.
  multivector<Eigen::VectorXd>::X<3> special_total_directions(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief finds high-symmetry directions within vector space supporting _rep wrt symmetry of head_group using the small_subgroups of head_group
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// a subgroup of 'head_group'. These are constructed by finding irreps of _rep and then calling special_irrep_directions on each
  /// result[i] is the set of special directions belonging to the i'th irrep constituting _rep
  /// result[i][j] s an orbit of symmetrically equivalent directions, and result[i][j][k] is an individual direction.
  /// Direction vectors are normalized to unit length. This routine is EXPERIMENTAL, USE AT YOUR OWN RISK
  multivector< Eigen::VectorXd >::X<3> calc_special_total_directions_experimental(SymGroupRep const &_rep, const SymGroup &head_group, double vector_norm_compare_tolerance = 0.001);

  /// \brief
  bool is_new_direction(const multivector<Eigen::VectorXd>::X<3> &special_directions, const Eigen::VectorXd &test_direction, double vector_norm_compare_tolerance);

  /// \brief generates the orbit of equivalent directions to direction by applying the operations of head_group
  /// \result A set of directions that are equivalent to direction under the application of head_group operations
  std::vector<Eigen::VectorXd> generate_special_direction_orbit(Eigen::VectorXd direction, const SymGroupRep &_rep, const SymGroup &head_group, double vector_norm_compare_tolerance);

  /// \brief finds high-symmetry subspaces within vector space supporting _rep, wrt symmetry of head_group
  /// High-symmetry subspaces are closed under the action of a nontrivial subgroup of head_group, without spanning
  /// the entire vector space supporting _rep
  /// \result Set of matrices MxN matrices (M>N) in the vector space on which '_rep' is defined, such that the column space of the
  /// matrix is invariant (i.e., closed) under the action of a nontrivial subgroup of head_group
  /// result[i] is an orbit of special subspaces that are equivalent by symmetry
  /// result[i][j] contains the spanning vectors of a specific special subspace
  /// Columns of each subspace matrix are orthogonal and normalized to unit length
  std::vector<std::vector< Eigen::MatrixXd> > special_subspaces(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Returns number of each irrep comprising _rep, listed in same order as the character table of master_group
  std::vector<Index> num_each_irrep(SymGroupRep const &_rep);

  /// \brief Returns number of each irrep of head_group comprising _rep, listed in same order as the character table of head_group
  std::vector<Index> num_each_irrep(SymGroupRep const &_rep, const SymGroup &head_group, bool verbose = false);

  /// \brief Returns number of each irrep of head_group comprising _rep, listed in same order as the character table of head_group
  /// _rep is assumed to be real, and irreps that are otherwise complex conjugates of each other are combined into a single real 'irrep'
  std::vector<Index> num_each_real_irrep(SymGroupRep const &_rep, const SymGroup &head_group, bool verbose = false);

  /// \brief Find irrep decomposition of _rep wrt group head_group and returns it as a list of SymGroupRepIDs corresponding to representtions of master_group
  std::vector<SymGroupRepID> get_irrep_IDs(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Returns true if _rep is irreducible wrt master_group (does not use character table information)
  bool is_irrep(SymGroupRep const &_rep);

  /// \brief Returns true if _rep is irreducible wrt head_group (does not use character table information)
  bool is_irrep(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Make a copy of representation on vector space 'V' that is transformed into a representation on vector space 'W'
  /// 'trans_mat' is the unitary matrix that isomorphically maps 'V'->'W' (i.e., [w = trans_mat * v] and [v = trans_mat.transpose() * w] )
  /// If the original representation to be transformed is just a temporary standalone SymGroupRep, be sure to delete it before falling out of scope
  SymGroupRep coord_transformed_copy(SymGroupRep const &_rep, const Eigen::MatrixXd &trans_mat);

  /// Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd get_irrep_trans_mat_old(SymGroupRep const &_rep, const SymGroup &head_group);

  /// Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// Also populate 'subspaces' with lists of columns that form irreps
  std::pair<Eigen::MatrixXd, std::vector<Index>> get_irrep_trans_mat_and_dims_old(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// This routine additionally orients the resulting basis vectors along high-symmetry directions of the
  /// vector space on which they are defined.
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result Transformation matrix with the ROWS comprising the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd get_irrep_trans_mat(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result Pair, with first element being the transformation matrix with the ROWS comprising
  /// the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// The second element is the dimension of irreducible subspaces, ordered identically to the rows of the transformation matrix
  std::pair<Eigen::MatrixXd, std::vector<Index>> get_irrep_trans_mat_and_dims(SymGroupRep const &_rep,
                                                                              const SymGroup &head_group,
                                                                              std::function<Eigen::MatrixXd(const SymGroupRep &,
                                                                                  const SymGroup &head_group)> symmetrizer_func);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result Pair, with first element being the transformation matrix with the ROWS comprising
  /// the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// The second element is the dimension of irreducible subspaces, ordered identically to the rows of the transformation matrix
  std::pair<Eigen::MatrixXd, std::vector<Index>> get_irrep_trans_mat_and_dims(SymGroupRep const &_rep,
                                                                              const SymGroup &head_group,
                                                                              std::function<Eigen::MatrixXd(const SymGroupRep &,
                                                                                  const SymGroup &head_group)> symmetrizer_func,
                                                                              Eigen::Ref<const Eigen::MatrixXd> const &_subspace);


  /// \brief Make copy of (*this) that is transformed so that axes are oriented along high-symmetry direction
  /// and confined to subspaces that transform as irreps.
  inline
  SymGroupRep symmetry_adapted_copy(SymGroupRep const &_rep, const SymGroup &head_group) {
    return coord_transformed_copy(_rep, get_irrep_trans_mat(_rep, head_group));
  }

  /// \brief Finds the project operators that project an arbitrary vector in the space supporting _rep onto irrep basis vectors
  /// \result The column space of result[i] spans the i'th irrep of master_group (if it is contained in _rep)
  std::vector<Eigen::MatrixXd> get_projection_operators(SymGroupRep const &_rep);


  SymGroupRep subset_permutation_rep(const SymGroupRep &permute_rep, const std::vector<std::vector<Index>> &subsets);
  SymGroupRep permuted_direct_sum_rep(const SymGroupRep &permute_rep, const std::vector<SymGroupRep const *> &sum_reps);
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

    const SymOp &sym_op(Index i) const {
      return (m_group_rep->master_group())[m_subgroup_op_inds[i]];
    }

    Index ind_inverse(Index i) const {
      return find_index(m_subgroup_op_inds, (m_group_rep->master_group()).ind_inverse(m_subgroup_op_inds[i]));
    }

    Index ind_prod(Index i, Index j) const {
      return find_index(m_subgroup_op_inds, (m_group_rep->master_group()).ind_prod(m_subgroup_op_inds[i], m_subgroup_op_inds[j]));
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
