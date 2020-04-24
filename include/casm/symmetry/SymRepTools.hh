#ifndef CASM_SymRepTools
#define CASM_SymRepTools

#include "casm/container/multivector.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"
namespace CASM {
  class SymGroupRep;
  class SymGroupRepID;
  class SymGroup;

  namespace SymRepTools {
    struct IrrepWedge {
      IrrepWedge() {}
      IrrepWedge(Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                 std::vector<Index> _mult) :
        axes(_axes),
        mult(_mult) {}
      Eigen::MatrixXd axes;
      std::vector<Index> mult;
    };

    class SubWedge {
    public:
      SubWedge(std::vector<IrrepWedge> const &_iwedges) :
        m_iwedges(_iwedges),
        m_trans_mat(_subwedge_to_trans_mat(m_iwedges)) {

      }

      std::vector<IrrepWedge> const &irrep_wedges() const {
        return m_iwedges;
      }

      Eigen::MatrixXd const &trans_mat() const {
        return m_trans_mat;
      }

    private:
      std::vector<IrrepWedge> m_iwedges;
      Eigen::MatrixXd m_trans_mat;

      static Eigen::MatrixXd _subwedge_to_trans_mat(std::vector<IrrepWedge> const &_iwedges);
    };


    std::vector<IrrepWedge> irreducible_wedges(const SymGroup &head_group, SymGroupRepID id);
    std::vector<SubWedge> symrep_subwedges(SymGroup const &_group, SymGroupRepID id);
  }

  bool rep_check(SymGroupRep const &_rep, SymGroup const &head_group, bool verbose);

  /// \brief Find irrep decomposition of _rep and add them to irrep listing in sub_group
  void calc_new_irreps(SymGroupRep const &_rep, const SymGroup &sub_group, int max_iter = 1000);

  /// \brief Assuming that _rep is an irrep of head_group, find high-symmetry directions
  /// throws if _rep is not an irrep
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// to a unique subgroup of 'head_group' (i.e., no other direction in the space, except the negative of that direction, is invariant to that subgroup)
  /// result[i] is an orbit of symmetrically equivalent directions, result[i][j] is an individual direction. Direction vectors are normalized to unit length.
  /// The total set of all directions is guaranteed to span the space.
  multivector<Eigen::VectorXd>::X<2> special_irrep_directions(SymGroupRep const &_rep,
                                                              SymGroup const &head_group,
                                                              double vec_compare_tol);


  /// \brief Assuming that _rep is an irrep of head_group, find high-symmetry directions
  /// throws if _rep is not an irrep
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// to a unique subgroup of 'head_group' (i.e., no other direction in the space, except the negative of that direction, is invariant to that subgroup)
  /// result[i] is an orbit of symmetrically equivalent directions, result[i][j] is an individual direction. Direction vectors are normalized to unit length.
  /// The total set of all directions is guaranteed to span the space.
  /// @param all_subgroups denotes whether all subgroups of head_group should be used for symmetry analysis (if true), or only cyclic subgroups (if false)
  multivector<Eigen::VectorXd>::X<2> special_irrep_directions(SymGroupRep const &_rep,
                                                              SymGroup const &head_group,
                                                              Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                                              double vec_compare_tol,
                                                              bool all_subgroups = false);


  /// Find a new coordinate system oriented along high-symmetry directions in vector space 'V' as determined by
  /// the subset of SymOpRepresentations specified by 'subgroup'.
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_symmetrizer(SymGroupRep const &_rep,
                                    const SymGroup &head_group,
                                    double vec_compare_tol);

  /// Find a new coordinate system oriented along high-symmetry directions in vector space spanned by
  /// matrix argument '_subspace', using as reference symmetry the subset of SymOpRepresentations specified by 'subgroup'.
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_symmetrizer(SymGroupRep const &_rep,
                                    const SymGroup &head_group,
                                    Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                    double vec_compare_tol);


  /// \brief Matrix B such that B(i,j) is sum of squares of each (i,j) matrix element over SymOpRepresentations in _rep corresponding to operations of head_group
  /// It reveals the block_diagonalization of SymGroupRep _rep for the subset of operations contained in head_group
  Eigen::MatrixXd block_shape_matrix(SymGroupRep const &_rep, const SymGroup &head_group);

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

  /// \brief finds high-symmetry subspaces within vector space supporting _rep, wrt symmetry of head_group
  /// High-symmetry subspaces are closed under the action of a nontrivial subgroup of head_group, without spanning
  /// the entire vector space supporting _rep
  /// \result Set of matrices MxN matrices (M>N) in the vector space on which '_rep' is defined, such that the column space of the
  /// matrix is invariant (i.e., closed) under the action of a nontrivial subgroup of head_group
  /// result[i] is an orbit of special subspaces that are equivalent by symmetry
  /// result[i][j] contains the spanning vectors of a specific special subspace
  /// Columns of each subspace matrix are orthogonal and normalized to unit length
  std::vector<std::vector< Eigen::MatrixXd> > special_subspaces(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Returns number of each irrep of head_group comprising _rep, listed in same order as the character table of head_group
  std::vector<Index> num_each_irrep(SymGroupRep const &_rep, const SymGroup &head_group, bool verbose = false);

  /// \brief Returns number of each irrep of head_group comprising _rep, listed in same order as the character table of head_group
  /// _rep is assumed to be real, and irreps that are otherwise complex conjugates of each other are combined into a single real 'irrep'
  std::vector<Index> num_each_real_irrep(SymGroupRep const &_rep, const SymGroup &head_group, bool verbose = false);

  /// \brief Find irrep decomposition of _rep wrt group head_group and returns it as a list of SymGroupRepIDs corresponding to representtions of master_group
  std::vector<SymGroupRepID> irrep_IDs(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Returns true if _rep is irreducible wrt head_group (does not use character table information)
  bool is_irrep(SymGroupRep const &_rep, const SymGroup &head_group);

  /// Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_trans_mat_old(SymGroupRep const &_rep, const SymGroup &head_group);

  /// Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// Also populate 'subspaces' with lists of columns that form irreps
  std::pair<Eigen::MatrixXd, std::vector<Index>> irrep_trans_mat_and_dims_old(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// This routine additionally orients the resulting basis vectors along high-symmetry directions of the
  /// vector space on which they are defined.
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result Transformation matrix with the ROWS comprising the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_trans_mat(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result Pair, with first element being the transformation matrix with the ROWS comprising
  /// the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// The second element is the dimension of irreducible subspaces, ordered identically to the rows of the transformation matrix
  std::pair<Eigen::MatrixXd, std::vector<Index>> irrep_trans_mat_and_dims(SymGroupRep const &_rep,
                                                                          const SymGroup &head_group,
                                                                          std::function<Eigen::MatrixXd(Eigen::Ref<const Eigen::MatrixXd> const &f_subspace)> symmetrizer_func);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result Pair, with first element being the transformation matrix with the ROWS comprising
  /// the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// The second element is the dimension of irreducible subspaces, ordered identically to the rows of the transformation matrix
  std::pair<Eigen::MatrixXd, std::vector<Index>> irrep_trans_mat_and_dims(SymGroupRep const &_rep,
                                                                          const SymGroup &head_group,
                                                                          std::function<Eigen::MatrixXd(Eigen::Ref<const Eigen::MatrixXd> const &f_subspace)> symmetrizer_func,
                                                                          Eigen::Ref<const Eigen::MatrixXd> const &_subspace);


  /// \brief Make copy of (*this) that is transformed so that axes are oriented along high-symmetry direction
  /// and confined to subspaces that transform as irreps.
  SymGroupRep symmetry_adapted_copy(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Finds the project operators that project an arbitrary vector in the space supporting _rep onto irrep basis vectors
  /// \result The column space of result[i] spans the i'th irrep of master_group (if it is contained in _rep)
  std::vector<Eigen::MatrixXd> projection_operators(SymGroupRep const &_rep, SymGroup const &head_group);


  //----- Operations involving representations without consideration of 'head_group' ----

  SymGroupRep subset_permutation_rep(const SymGroupRep &permute_rep, const std::vector<std::vector<Index>> &subsets);
  SymGroupRep permuted_direct_sum_rep(const SymGroupRep &permute_rep, const std::vector<SymGroupRep const *> &sum_reps);
  SymGroupRep kron_rep(const SymGroupRep &LHS, const SymGroupRep &RHS);


}
#endif
