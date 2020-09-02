#ifndef CASM_SymRepTools
#define CASM_SymRepTools

#include <set>
#include "casm/container/multivector.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"

namespace CASM {
  class SymGroupRep;
  class SymGroupRepID;
  class SymGroup;

  namespace SymRepTools {

    /// A Symmetrizer comprises a MatrixXcd, which defines a rotation of a generically-oriented irreducible subspace
    /// that aligns its components along high-symmetry directions, and a vector of orbits of high-symmetry directions
    using Symmetrizer = std::pair<Eigen::MatrixXcd, multivector<Eigen::VectorXcd>::X<2> >;

    /// A SymmetrizerFunction takes a matrix whose columns specify an irreducible subspace, and returns a Symmetrizer
    /// for symmetrizing that subspace
    using SymmetrizerFunction = std::function<Symmetrizer(Eigen::Ref<const Eigen::MatrixXcd> const &f_subspace)>;


    struct IrrepInfo {

      /// \brief Named constructor to initialize an IrrepInfo with a user-specified transformtion matrix and character vector of [(dim,0)]
      /// where 'dim' is dimension of irrep
      template<typename Derived>
      static IrrepInfo make_dummy(Eigen::MatrixBase<Derived> const &_trans_mat) {
        Eigen::VectorXcd tchar(1);
        tchar[0] = std::complex<double>(double(_trans_mat.rows()), 0.);
        return IrrepInfo(_trans_mat. template cast<std::complex<double> >(),
                         tchar);
      }

      /// \brief Construct an IrrepInfo with transformation matrix and vector of irreducible characters
      IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters);

      /// Dimension of irreducible vector space (less than or equal to vector_dim())
      Index irrep_dim() const {
        return trans_mat.rows();
      }

      // Dimension of initial vector space (greater than or equal to irrep_dim())
      Index vector_dim() const {
        return trans_mat.cols();
      }

      /// true if any character has non-zero imaginary component, false otherwise
      bool complex;

      /// true if irrep is real but was created as direct sum of two complex irreps
      /// in this case, the 'irrep' is reducible, but this is the most-reduced
      /// representation that can still have real basis vectors
      bool pseudo_irrep;

      /// sequentially-assigned index used to distinguish between identical irreps
      /// irreps are identical if they have the same character vectors
      Index index;

      /// irrep_dim() x vector_dim() matrix that transforms a vector from the initial
      /// vector space into a vector in the irreducible vector space
      Eigen::MatrixXcd trans_mat;

      /// vector containing complex character of each group operation's action on the irreducible vector space
      Eigen::VectorXcd characters;

      /// Vectors in the initial vector space that correspond to high-symmetry directions in the
      /// irreducible vector space. directions[i] is the i'th orbit of equivalent high-symmetry directions
      /// and directions[i].size() is the symmetric multiplicity of a direction in that orbit
      std::vector<std::vector<Eigen::VectorXd> > directions;
    };


    ///\brief An irreducible wedge in an irreducible vector space
    struct IrrepWedge {

      IrrepWedge(IrrepInfo _irrep_info,
                 Eigen::MatrixXd _axes) :
        irrep_info(std::move(_irrep_info)),
        axes(std::move(_axes)) {}

      /// The description of the associated irreducible vector space
      IrrepInfo irrep_info;

      /// columns of 'axes' are high-symmetry direction that form the edges
      /// of the irreducible wedge. The number of columns should be equal to the
      /// dimension of the irreducible subspace
      Eigen::MatrixXd axes;

      /// Symmetric multiplicity of each column of 'axes'
      std::vector<Index> mult;
    };

    ///\brief SubWedge is a vector of IrrepWedge, one from each irreducible subspace
    /// Together, the IrrepWedges that comprise the Subwedge span the entire space
    /// However, it is likely that the orbit of equivalent SubWedges does not *fill* the entire space
    class SubWedge {
    public:
      SubWedge(std::vector<IrrepWedge> const &_iwedges) :
        m_iwedges(_iwedges),
        m_trans_mat(_subwedge_to_trans_mat(m_iwedges)) {
      }

      /// IrrepWedges comprising the Subwedge
      std::vector<IrrepWedge> const &irrep_wedges() const {
        return m_iwedges;
      }

      /// Transformation matrix to convert from a vector in terms of the SubWedge axes
      /// to a vector in the original vector space
      Eigen::MatrixXd const &trans_mat() const {
        return m_trans_mat;
      }

    private:
      std::vector<IrrepWedge> m_iwedges;
      Eigen::MatrixXd m_trans_mat;

      static Eigen::MatrixXd _subwedge_to_trans_mat(std::vector<IrrepWedge> const &_iwedges);
    };


    /// Find IrrepWedges for a group-represented vector space, knowing only the group and the representation ID
    /// @param _subspace matrix whose columns span a subspace of the full vector space (pass Identity matrix to decompose full space)
    std::vector<IrrepWedge> irrep_wedges(SymGroup const &head_group,
                                         SymGroupRepID id,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

    /// Find IrrepWedges for a group-represented vector space, knowing only the group and representation
    /// @param _subspace matrix whose columns span a subspace of the full vector space (pass Identity matrix to decompose full space)
    std::vector<IrrepWedge> irrep_wedges(SymGroupRep const &_rep,
                                         SymGroup const &head_group,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

    /// Find full irreducible wedge of a group-represented vector space, as a vector of SubWedges, knowing only the SymGroup and representation ID
    std::vector<SubWedge> symrep_subwedges(SymGroup const &_group,
                                           SymGroupRepID id);

    /// Find full irreducible wedge of a group-represented vector space, as a vector of SubWedges, knowing only the SymGroup and representation ID
    /// @param _subspace matrix whose columns span a subspace of the full vector space
    std::vector<SubWedge> symrep_subwedges(SymGroup const &_group,
                                           SymGroupRepID id,
                                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

    /// Find full irreducible wedge of a group-represented vector space, as a vector of SubWedges, knowing only the SymGroup and representation
    std::vector<SubWedge> symrep_subwedges(SymGroupRep const &_rep,
                                           SymGroup const &head_group);

    /// Find full irreducible wedge of a group-represented vector space, as a vector of SubWedges, knowing only the SymGroup and representation
    /// @param _subspace matrix whose columns span a subspace of the full vector space
    std::vector<SubWedge> symrep_subwedges(SymGroupRep const &_rep,
                                           SymGroup const &head_group,
                                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace);

  }

  ///\brief Summary of data associated with the action of a symmetry group on a vector space
  struct VectorSpaceSymReport {

    /// \brief Matrix representation for each operation in the group -- defines action of group on vector space
    std::vector<Eigen::MatrixXd> symgroup_rep;

    /// \brief A list of all irreducible representation that make up the full representation
    std::vector<SymRepTools::IrrepInfo> irreps;

    /// \brief Irreducible wedge in the vector space
    /// encoded as a vector of symmetrically distinct SubWedges
    std::vector<SymRepTools::SubWedge> irreducible_wedge;

    /// \brief Symmetry-oriented subspace of the vector space (columns are the basis vectors)
    Eigen::MatrixXd symmetry_adapted_dof_subspace;

    /// \brief Names given to individual axes in initial (un-adapted) vector space, corresponding to rows of symmetry_adapted_dof_subspace
    std::vector<std::string> axis_glossary;
  };

  ///\brief Construct the VectorSpaceSymReport for
  /// @param _rep matrix representation of head_group, this defines group action on the underlying vector space
  /// @param head_group group for which the sym report is to be generated
  /// @param _subspace matrix such that _subspace.rows()==_rep.dim() and whose columns specify a subspace of underlying vector space
  /// @param calc_wedges if true, 'irreducible_wedge' of returned object is initialized, if false, 'irreducible_wedge' is empty
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
  /// @param vec_compare_tol tolerance for elementwise floating-point comparisons of vectors
  multivector<Eigen::VectorXcd>::X<2> special_irrep_directions(SymGroupRep const &_rep,
                                                               SymGroup const &head_group,
                                                               double vec_compare_tol);


  /// \brief Assuming that _rep is an irrep of head_group, find high-symmetry directions
  /// throws if _rep is not an irrep
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// to a unique subgroup of 'head_group' (i.e., no other direction in the space, except the negative of that direction, is invariant to that subgroup)
  /// result[i] is an orbit of symmetrically equivalent directions, result[i][j] is an individual direction. Direction vectors are normalized to unit length.
  /// The total set of all directions is guaranteed to span the space.
  /// @param _subspace matrix such that _subspace.rows()==_rep.dim() and whose columns specify a subspace of underlying vector space
  /// @param vec_compare_tol tolerance for elementwise floating-point comparisons of vectors
  /// @param all_subgroups denotes whether all subgroups of head_group should be used for symmetry analysis (if true), or only cyclic subgroups (if false)
  multivector<Eigen::VectorXcd>::X<2> special_irrep_directions(SymGroupRep const &_rep,
                                                               SymGroup const &head_group,
                                                               Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                                               double vec_compare_tol,
                                                               bool all_subgroups = false);


  /// Find a new coordinate system oriented along high-symmetry directions in underlying vector space
  /// @param _rep matrix representation specifying group action on underlying vector space
  /// @param head_group subgroup of operations in _rep to be used for symmetry analysis
  /// @param vec_compare_tol tolerance for elementwise floating-point comparisons of vectors
  /// The ROWS of result.first are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = result.first * old_symrep_matrix * result.first.transpose();
  SymRepTools::Symmetrizer irrep_symmetrizer(SymGroupRep const &_rep,
                                             const SymGroup &head_group,
                                             double vec_compare_tol);

  /// Find a new coordinate system oriented along high-symmetry directions in vector space spanned by columns of @param _subspace
  /// @param _rep matrix representation specifying group action on underlying vector space
  /// @param head_group subgroup of operations in _rep to be used for symmetry analysis
  /// @param _subspace matrix such that _subspace.rows()==_rep.dim() and whose columns specify a subspace of underlying vector space
  /// @param vec_compare_tol tolerance for elementwise floating-point comparisons of vectors
  /// The ROWS of result.first are the new basis vectors in terms of the old such that
  /// new_symrep_matrix = result.first * old_symrep_matrix * result.first.transpose();
  SymRepTools::Symmetrizer irrep_symmetrizer(SymGroupRep const &_rep,
                                             const SymGroup &head_group,
                                             Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                             double vec_compare_tol);


  /// \brief counts number of nonzero blocks in matrix representation of head_group as specified by  _rep
  /// Reveals number of invariant subspaces (with respect to head_group) that comprise the vector space supporting _rep
  Index num_blocks(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Assumes that irreps are real, and concatenates their individual trans_mats to form larger trans_mat
  Eigen::MatrixXd full_trans_mat(std::vector<SymRepTools::IrrepInfo> const &irreps);

  /// \brief finds high-symmetry directions within vector space supporting _rep, wrt symmetry of head_group
  /// \result Set of directions in the vector space on which '_rep' is defined, such that each direction is invariant
  /// a subgroup of 'head_group'. These are constructed by finding irreps of _rep and then calling special_irrep_directions on each
  /// result[i] is the set of special directions belonging to the i'th irrep constituting _rep
  /// result[i][j] is an orbit of symmetrically equivalent directions, and result[i][j][k] is an individual direction.
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

  /// \brief Find irrep decomposition of _rep wrt group head_group and returns it as a list of SymGroupRepIDs corresponding to representtions of master_group
  std::vector<SymGroupRepID> irrep_IDs(SymGroupRep const &_rep, const SymGroup &head_group);

  /// \brief Returns true if _rep is irreducible wrt head_group (does not use character table information)
  bool is_irrep(SymGroupRep const &_rep, const SymGroup &head_group);

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
  /// @param symmetrizer_func, a function object that takes an unsymmetrized irrep subspace and returns its Symmetrizer object
  /// \result Pair, with result.first being the transformation matrix with the ROWS comprising
  /// the new basis vectors in terms of the old such that
  /// new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  /// result.second contains the complex character vectors of the irreducible subspaces, in same order as they appear in result.first
  std::pair<Eigen::MatrixXd, std::vector<Eigen::VectorXcd> > irrep_trans_mat_and_characters(SymGroupRep const &_rep,
      const SymGroup &head_group,
      SymRepTools::SymmetrizerFunction symmetrizer_func);

  /// \brief Finds irreducible subspaces that comprise an underlying subspace
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// @param _rep matrix representation of head_group, this defines group action on the underlying vector space
  /// @param head_group group for which the irreps are to be found
  /// @param allow_complex if true, irreps may be complex-valued, if false, complex irreps would be combined to form real representations
  /// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with identity first (if present)
  /// repeated irreps are sequential, and are distinguished by IrrepInfo::index
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          bool allow_complex);

  /// \brief Finds irreducible subspaces that comprise an underlying subspace
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// @param _rep matrix representation of head_group, this defines group action on the underlying vector space
  /// @param head_group group for which the irreps are to be found
  /// @param allow_complex if true, irreps may be complex-valued, if false, complex irreps would be combined to form real representations
  /// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with identity first (if present)
  /// repeated irreps are sequential, and are distinguished by IrrepInfo::index
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                                          bool allow_complex);

  /// \brief Finds irreducible subspaces that comprise an underlying subspace
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// @param _rep matrix representation of head_group, this defines group action on the underlying vector space
  /// @param head_group group for which the irreps are to be found
  /// @param allow_complex if true, irreps may be complex-valued, if false, complex irreps would be combined to form real representations
  /// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with identity first (if present)
  /// repeated irreps are sequential, and are distinguished by IrrepInfo::index
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                          bool allow_complex);


  /// \brief Finds irreducible subspaces that comprise an underlying subspace
  /// It does not rely on the character table, but instead utilizes a brute-force method
  /// @param _rep matrix representation of head_group, this defines group action on the underlying vector space
  /// @param head_group group for which the irreps are to be found
  /// @param allow_complex if true, irreps may be complex-valued, if false, complex irreps would be combined to form real representations
  /// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with identity first (if present)
  /// repeated irreps are sequential, and are distinguished by IrrepInfo::index
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                          Eigen::MatrixXd subspace,
                                                          bool allow_complex);

  //----- Operations involving representations without consideration of 'head_group' ----

  /// Given a permutation representation that permutes indices, calculate a new permutation representation of lower or equal dimension
  /// that describes the permutation of ordered subsets of indices. The total number of indices contained in @param subsets must be <= permute_rep.dim()
  /// and each Index in @param subsets must be unique an on the range [0,permute_rep.dim())
  SymGroupRep subset_permutation_rep(const SymGroupRep &permute_rep, const std::vector<std::set<Index>> &subsets);

  /// Given a permutation of positions, and a list of matrix representations (one per position), build a large matrix representation
  /// that has blocks specified by the smaller matrices in @param sum_reps and whose positions within the large matrix correspond to
  /// the permutation from @param permute_rep
  SymGroupRep permuted_direct_sum_rep(const SymGroupRep &permute_rep, const std::vector<SymGroupRep const *> &sum_reps);

  /// Build a large matrix representation by forming kronecker products of a matrix from @param LHS with a matrix from @param RHS
  SymGroupRep kron_rep(const SymGroupRep &LHS, const SymGroupRep &RHS);


}
#endif
