#ifndef CASM_SymRepTools
#define CASM_SymRepTools

#include "casm/container/multivector.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
  class SymGroupRep;
  class SymGroupRepID;
  class SymGroup;

  namespace SymRepTools {

    using Symmetrizer = std::pair<Eigen::MatrixXcd, multivector<Eigen::VectorXcd>::X<2> >;
    using SymmetrizerFunction = std::function<Symmetrizer(Eigen::Ref<const Eigen::MatrixXcd> const &f_subspace)>;

    struct IrrepInfo {
      template<typename Derived>
      static IrrepInfo make_dummy(Eigen::MatrixBase<Derived> const &_trans_mat) {
        Eigen::VectorXcd tchar(1);
        tchar[0] = std::complex<double>(double(_trans_mat.rows()), 0.);
        return IrrepInfo(_trans_mat. template cast<std::complex<double> >(),
                         tchar);
      }

      IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters);

      // Dimension of irreducible vector space
      Index irrep_dim() const {
        return trans_mat.cols();
      }

      // Dimension of initial vector space
      Index vector_dim() const {
        return trans_mat.rows();
      }

      /// true if any character has non-zero imaginary component, false otherwise
      bool complex;

      /// true if irrep is real but was created as direct sum of two complex irreps
      /// in this case, the 'irrep' is reducible, but this is the most-reduced
      /// representation that can still have real basis vectors
      bool pseudo_irrep;
      Index index;

      /// vector_dim() x irrep_dim() matrix describes of the 'big' vector space that is the irrep
      Eigen::MatrixXcd trans_mat;

      /// complex characters of each operation in the group for this irrep
      Eigen::VectorXcd characters;

      /// High-symmetry directions in the 'big' vector space that correspond to this irrep
      std::vector<std::vector<Eigen::VectorXd> > directions;
    };


    struct IrrepWedge {

      IrrepWedge(IrrepInfo _irrep_info,
                 Eigen::MatrixXd _axes) :
        irrep_info(std::move(_irrep_info)),
        axes(std::move(_axes)) {}

      IrrepInfo irrep_info;

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


    std::vector<IrrepWedge> irrep_wedges(SymGroup const &head_group,
                                         SymGroupRepID id,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

    std::vector<IrrepWedge> irrep_wedges(SymGroupRep const &_rep,
                                         SymGroup const &head_group,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

    std::vector<SubWedge> symrep_subwedges(SymGroup const &_group,
                                           SymGroupRepID id);

    std::vector<SubWedge> symrep_subwedges(SymGroup const &_group,
                                           SymGroupRepID id,
                                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

    std::vector<SubWedge> symrep_subwedges(SymGroupRep const &_rep,
                                           SymGroup const &head_group);

    std::vector<SubWedge> symrep_subwedges(SymGroupRep const &_rep,
                                           SymGroup const &head_group,
                                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

  }


  struct VectorSpaceSymReport {
    /// \brief Matrix representation for each operation in the group
    std::vector<Eigen::MatrixXd> symgroup_rep;

    /// \brief A list of all irreducible representation that make up the big representation
    std::vector<SymRepTools::IrrepInfo> irreps;

    /// \brief Irreducible wedge in the vector
    std::vector<SymRepTools::SubWedge> irreducible_wedge;

    /// \brief Symmetry-oriented subspace of the vector space (columns are the basis vectors)
    Eigen::MatrixXd symmetry_adapted_dof_subspace;
  };

  VectorSpaceSymReport vector_space_sym_report(SymGroupRep const &_rep,
                                               SymGroup const &head_group,
                                               Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                               bool calc_wedges = false);

  /// \brief Assuming that _rep is an irrep of head_group, find high-symmetry directions
  /// throws if _rep is not an irrep
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// to a unique subgroup of 'head_group' (i.e., no other direction in the space, except the negative of that direction, is invariant to that subgroup)
  /// result[i] is an orbit of symmetrically equivalent directions, result[i][j] is an individual direction. Direction vectors are normalized to unit length.
  /// The total set of all directions is guaranteed to span the space.
  multivector<Eigen::VectorXcd>::X<2> special_irrep_directions(SymGroupRep const &_rep,
                                                               SymGroup const &head_group,
                                                               double vec_compare_tol);


  /// \brief Assuming that _rep is an irrep of head_group, find high-symmetry directions
  /// throws if _rep is not an irrep
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// to a unique subgroup of 'head_group' (i.e., no other direction in the space, except the negative of that direction, is invariant to that subgroup)
  /// result[i] is an orbit of symmetrically equivalent directions, result[i][j] is an individual direction. Direction vectors are normalized to unit length.
  /// The total set of all directions is guaranteed to span the space.
  /// @param all_subgroups denotes whether all subgroups of head_group should be used for symmetry analysis (if true), or only cyclic subgroups (if false)
  multivector<Eigen::VectorXcd>::X<2> special_irrep_directions(SymGroupRep const &_rep,
                                                               SymGroup const &head_group,
                                                               Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                                               double vec_compare_tol,
                                                               bool all_subgroups = false);


  /// Find a new coordinate system oriented along high-symmetry directions in vector space 'V' as determined by
  /// the subset of SymOpRepresentations specified by 'subgroup'.
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  SymRepTools::Symmetrizer irrep_symmetrizer_and_directions(SymGroupRep const &_rep,
                                                            const SymGroup &head_group,
                                                            double vec_compare_tol);

  /// Find a new coordinate system oriented along high-symmetry directions in vector space spanned by
  /// matrix argument '_subspace', using as reference symmetry the subset of SymOpRepresentations specified by 'subgroup'.
  /// The ROWS of trans_mat are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  SymRepTools::Symmetrizer irrep_symmetrizer_and_directions(SymGroupRep const &_rep,
                                                            const SymGroup &head_group,
                                                            Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                                            double vec_compare_tol);

  Eigen::MatrixXcd irrep_symmetrizer_from_directions(multivector<Eigen::VectorXcd>::X<2> const &special_directions,
                                                     Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                                     double vec_compare_tol);



  /// \brief Matrix B such that B(i,j) is sum of squares of each (i,j) matrix element over SymOpRepresentations in _rep corresponding to operations of head_group
  /// It reveals the block_diagonalization of SymGroupRep _rep for the subset of operations contained in head_group
  Eigen::MatrixXd block_shape_matrix(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief counts number of blocks in the block_shape_matrix of _rep
  /// Reveals number of invariant subspaces (with respect to head_group) that comprise the vector space supporting _rep
  Index num_blocks(SymGroupRep const &_rep, const SymGroup &head_group);

  Eigen::MatrixXd full_trans_mat(std::vector<SymRepTools::IrrepInfo> const &irreps);

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
  std::pair<Eigen::MatrixXd, std::vector<Eigen::VectorXcd> > irrep_trans_mat_and_characters(SymGroupRep const &_rep,
      const SymGroup &head_group,
      SymRepTools::SymmetrizerFunction symmetrizer_func);

  /// \brief Finds the transformation matrix that block-diagonalizes this representation of head_group into irrep blocks
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// \param head_group The group with respect to which irreps are determined, which may be a subset of all operations in this representation
  /// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with identity first (if present)
  /// repeated irreps are sequential, and are distinguished by IrrepInfo::index
  //OLD the new basis vectors in terms of the old such that
  //OLD new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  //OLD The second element is the dimension of irreducible subspaces, ordered identically to the rows of the transformation matrix
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          bool allow_complex);

  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                                          bool allow_complex);

  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                          bool allow_complex);


  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                          Eigen::MatrixXd subspace,
                                                          bool allow_complex);


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
