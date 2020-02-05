#ifndef CASM_DoFSet
#define CASM_DoFSet

#include <vector>
#include <string>
#include "casm/basis_set/DoF.hh"



namespace CASM {

  template<typename T> struct jsonConstructor;


  class AnisoValTraits;

  namespace xtal {
    struct SymOp;
  }

  class SymGroup;

  // All identifying information for a vector of continuous independent variables (Degrees of Freedom / DoFs)
  // DoFSets are associated with a specific DoF 'type', which has a predefined 'standard' coordinate system
  // ex: displacement -> 3-vector (x,y,z) -> displacement components (relative to fixed laboratory frame)
  //     strain -> 6-vector (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy) -> tensor elements
  // DoFSets have a typename, which specifies the type, and a set of basis vectors, which are denoted relative
  // to the DoF type's standard axes. This allows the DoFSet components to be specified by the user,
  // including the ability to only allow DoF values within a subspace of the standard valeus.
  // DoFSet records (at least) the DoF typename, the names of the vector components, the axes of the
  // vector components (relative to a set of standard axes), and the SymGroupRepID of the DoFSet (which
  // is a key for retrieving the SymGroupRep that encodes how the DoFSet transforms with symmetry.
  class DoFSet;

  /// \brief struct for consolidating the SymGroupRepID and coordinates axes for continuous vector DoF
  ///  This is the minimal data object required for performing most crystallographic symmetry analysis
  struct DoFSetInfo {

    /// \brief construct DoFSetInfo with SymGroupRepID and _basis matrix
    /// \param _basis {columns are vectors in stard vector space that specify the alignment of DoFSet components}
    /// A particulr DoFSet column-vector value, 'v', represented in the standard coordinate space is given by
    /// v_standard = _basis * v
    DoFSetInfo(SymGroupRepID _id, Eigen::Ref<const Eigen::MatrixXd> const &_basis) :
      m_symrep_ID(_id),
      m_basis(_basis) {
      if(basis().cols() > 0 && basis().rows() > 0)
        m_inv_basis = basis().transpose().colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(dim(), dim())).transpose();
    }

    /// \brief returns SymGroupRepID of associated DoFSet
    SymGroupRepID const &symrep_ID() const {
      return m_symrep_ID;
    }

    /// \brief sets SymGroupRepID
    void set_symrep_ID(SymGroupRepID _id) {
      m_symrep_ID = _id;
    }

    /// \brief returns const reference to DoFSet coordinate axes
    Eigen::MatrixXd const &basis() const {
      return m_basis;
    }

    /// \brief Sets the DoFSet coordinate axes and caches the pseudoinverse (which is used for symmetry analysis, among other things)
    Eigen::MatrixXd set_basis(Eigen::Ref<const Eigen::MatrixXd> const &_basis) {
      m_basis = _basis;
      if(basis().cols() > 0 && basis().rows() > 0)
        m_inv_basis = basis().transpose().colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(dim(), dim())).transpose();
      return m_basis;
    }

    /// \brief pseudoinverse of DoFSet coordinate axes, useful for transforming from standard basis to user-specified basis
    Eigen::MatrixXd const  &inv_basis() const {
      return m_inv_basis;
    }

    /// \brief Dimension of the DoFSet, equivalent to basis().cols()
    Index dim() const {
      return basis().cols();
    }

  private:
    SymGroupRepID m_symrep_ID;
    Eigen::MatrixXd m_basis;
    Eigen::MatrixXd m_inv_basis;
  };


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// \briefDoFSet specifies all identifying information for a vector of continuous independent variables (Degrees of Freedom / DoFs)
  /// DoFSets are associated with a specific DoF 'type', which has a predefined 'standard' coordinate system
  /// ex: displacement -> 3-vector (x,y,z) -> displacement components (relative to fixed laboratory frame)
  ///     strain -> 6-vector (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy) -> tensor elements
  /// DoFSets have a typename, which specifies the type, and a set of basis vectors, which are denoted relative
  /// to the DoF type's standard axes. This allows the DoFSet components to be specified by the user,
  /// including the ability to only allow DoF values within a subspace of the standard valeus.
  /// DoFSet records (at least) the DoF typename, the names of the vector components, the axes of the
  /// vector components (relative to a set of standard axes), and the SymGroupRepID of the DoFSet (which
  /// is a key for retrieving the SymGroupRep that encodes how the DoFSet transforms with symmetry.
  class DoFSet {
  public:
    using BasicTraits = AnisoValTraits;

    using TypeFunc =  std::function<notstd::cloneable_ptr<BasicTraits>()>;

    using Container = std::vector<ContinuousDoF>;

    using const_iterator = std::vector<ContinuousDoF>::const_iterator;

    /// \brief Factor function for DoFSet having standard dimension and coordinate basis for given type.
    /// Example usage:
    ///    \code
    ///    DoFSet my_disp_dof = DoFSet::make_default(DoF::traits("disp"));
    ///    \endcode
    static DoFSet make_default(DoFSet::BasicTraits const &_type);

    /// \brief Construct with BasicTraits object, in order to catch typos near point of origin
    /// DoFSet is initialized to zero dimension. Use static function DoFSet::make_default() to construct
    /// A viable DoFSet having standard parameters
    /// Example usage:
    ///    \code
    ///    DoFSet my_disp_dof(DoF::traits("disp"))
    ///    \endcode
    DoFSet(BasicTraits const &_type);

    /// \brief Returns number of components in this DoFSet
    Index size() const {
      return m_components.size();
    }

    /// \brief Returns type_name of DoFSet, which should be a standardized DoF type (e.g., "disp", "magspin", "GLstrain")
    std::string const &type_name() const {
      return m_traits.name();
    }

    /// \brief Returns traits object for the DoF type of this DoFSet
    BasicTraits const &traits() const {
      return m_traits;
    }

    /// \brief Set identifying index (ID) of this DoFSet
    /// ID is context dependent but is used to different distinct DoFSets that are otherwise equivalent
    /// (i.e., have same type_name(), components, etc)
    /// Example, two sites in the primitive cell may have the same DoFs, but they have different indices within
    /// the primitive basis. Their corresponding DoFSets can be assigned IDs in order to distinguish the DoFSets from each other
    void set_ID(Index _ID) {
      for(auto &c : m_components)
        c.set_ID(_ID);
    }

    /// \brief Returns reference to DoFSetInfo (which contains SymGroupRepID and basis vectors of coordinate system
    DoFSetInfo const &info() const {
      return m_info;
    }

    /// \brief Return i'th component of DoFSet
    /// DoFSet components are represented as ContinuousDoF, which, in addition to a type_name also has a variable name
    ContinuousDoF const &operator[](Index i) const {
      return m_components[i];
    }

    /// \brief Iterator pointing to first component
    const_iterator begin() const {
      return m_components.cbegin();
    }

    /// \brief Iterator pointing one past last component
    const_iterator end() const {
      return m_components.cend();
    }

    /// \brief const iterator to first component
    const_iterator cbegin() const {
      return m_components.cbegin();
    }

    /// \brief const iterator pointing one past last component
    const_iterator cend() const {
      return m_components.cend();
    }

    /// \brief Returns true if DoFSet is inactive (e.g., takes zero values) when specified occupant is present
    bool is_excluded_occ(std::string const &_occ_name) const {
      return m_excluded_occs.count(_occ_name);
    }

    Index dim() const {
      return m_info.dim();
    }

    /// \brief Matrix that relates DoFSet variables to a conventional coordiante system
    /// columns of coordinate_space() matrix are directions in conventional coordinate system
    /// so that  conventional_coord = DoFSet.coordinate_space()*DoFSet.values()
    /// coordinate_space() matrix has dimensions (N x size()), where N >= size()
    Eigen::MatrixXd const &basis() const {
      return m_info.basis();
    }

    /// \brief SymGroupRepID of this DoFSet, necessary for applying symmetry transformations
    SymGroupRepID const &symrep_ID() const {
      return m_info.symrep_ID();
    }

    /// \brief Allocates an empty symmetry representation in \param _group and records its SymGroupRepID
    /// This representation becomes THE representation for this DoFSet, to be initialized accordingly
    void allocate_symrep(SymGroup const &_group) const;

    /// \brief returns true of \param rhs has identical components and basis to this DoFSet
    bool identical(DoFSet const &rhs) const;

    /// \brief Equivalent to m_basis=trans_mat*m_basis. Invalidates SymGroupRepID
    void transform_basis(Eigen::Ref<const Eigen::MatrixXd> const &trans_mat);

    /// \brief Convenience function, when DoFSet IDs are being set, en masse
    /// Each ID before_IDs[i] will get changed to corresponding ID after_IDs[i].
    /// Determines if this DoFSet has an ID contained in before_IDs -- If so, changes it to corresponding
    /// ID in after_IDs and returns true; otherwise returns false.
    bool update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs);

    /// \brief Parse DoFSet basis vectors, components, and excluded_occs from JSON input
    void from_json(jsonParser const &json);

    /// \brief Write DoFSet basis vectors, components, and excluded_occs to JSON
    jsonParser &to_json(jsonParser &json) const;


  private:

    //AnisovalTraits, defines the intrinsic space of a type of DoF
    BasicTraits m_traits;
    std::vector<ContinuousDoF> m_components;
    mutable DoFSetInfo m_info;

    std::set<std::string> m_excluded_occs;

  };

  //********************************************************************

  inline
  bool operator==(DoFSet const &A, DoFSet const &B) {
    return A.identical(B);
  }

  //********************************************************************

  inline
  bool operator!=(DoFSet const &A, DoFSet const &B) {
    return !A.identical(B);
  }

  //********************************************************************

  /// \brief Apply SymOp to a DoFSet
  DoFSet &apply(const xtal::SymOp &op, DoFSet &_dof);

  //********************************************************************

  /// \brief Copy and apply SymOp to a DoFSet
  DoFSet copy_apply(const xtal::SymOp &op, const DoFSet &_dof);

  //********************************************************************
  template<>
  struct jsonConstructor<DoFSet> {
    static DoFSet from_json(const jsonParser &json, DoFSet::BasicTraits const &_type);
  };

  //********************************************************************
  inline
  void from_json(DoFSet &_dof, jsonParser &json) {
    return _dof.from_json(json);
  }

  //********************************************************************
  inline
  jsonParser &to_json(DoFSet const &_dof, jsonParser &json) {
    return _dof.to_json(json);
  }



}

#endif
