#ifndef CASM_DoFSet
#define CASM_DoFSet

#include <vector>
#include "casm/basis_set/DoF.hh"


namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<typename T> struct jsonConstructor;

  class SymOp;
  class SymGroup;

  class DoFSet {
  public:
    using BasicTraits = DoFType::BasicTraits;

    using TypeFunc =  std::function<notstd::cloneable_ptr<BasicTraits>()>;

    using Container = std::vector<ContinuousDoF>;

    using const_iterator = std::vector<ContinuousDoF>::const_iterator;

    DoFSet(BasicTraits const &_type) :
      m_type_name(_type.type_name()),
      m_info(SymGroupRepID(), Eigen::MatrixXd::Zero(_type.dim(), 0)) {}

    Index size() const {
      return m_components.size();
    }

    std::string const &type_name() const {
      return m_type_name;
    }

    DoFType::Traits const &traits() const;

    void set_ID(Index _ID) {
      for(auto &c : m_components)
        c.set_ID(_ID);
    }

    DoFSetInfo const &info() const {
      return m_info;
    }

    ContinuousDoF const &operator[](Index i) const {
      return m_components[i];
    }

    const_iterator begin() const {
      return m_components.cbegin();
    }

    const_iterator end() const {
      return m_components.cend();
    }

    const_iterator cbegin() const {
      return m_components.cbegin();
    }

    const_iterator cend() const {
      return m_components.cend();
    }

    bool is_excluded_occ(std::string const &_occ_name) const {
      return m_excluded_occs.count(_occ_name);
    }

    /// \brief Matrix that relates DoFSet variables to a conventional coordiante system

    /// columns of coordinate_space() matrix are directions in conventional coordinate system
    /// so that  conventional_coord = DoFSet.coordinate_space()*DoFSet.values()
    /// coordinate_space() matrix has dimensions (N x size()), where N >= size()
    Eigen::MatrixXd const &basis() const {
      return m_info.basis();
    }

    SymGroupRepID const &symrep_ID() const {
      return m_info.symrep_ID();
    }

    void allocate_symrep(SymGroup const &_group) const;

    bool identical(DoFSet const &rhs) const;

    /// \brief Equivalent to m_basis=trans_mat*m_basis. Invalidates SymGroupRepID
    void transform_basis(Eigen::Ref<const Eigen::MatrixXd> const &trans_mat);

    bool update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs);

    void from_json(jsonParser const &json);

    jsonParser &to_json(jsonParser &json) const;


  private:
    std::string m_type_name;
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
  DoFSet &apply(const SymOp &op, DoFSet &_dof);

  //********************************************************************

  /// \brief Copy and apply SymOp to a DoFSet
  DoFSet copy_apply(const SymOp &op, const DoFSet &_dof);

  //********************************************************************
  template<>
  struct jsonConstructor<DoFSet> {
    static DoFSet from_json(const jsonParser &json, DoF::BasicTraits const &_type);
  };

  //********************************************************************
  inline
  jsonParser &to_json(DoFSet const &_dof, jsonParser &json) {
    return _dof.to_json(json);
  }



}

#endif
