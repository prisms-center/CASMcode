#include "casm/basis_set/DoFSet.hh"

#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {

  DoFSet::DoFSet(AnisoValTraits const &_type,
                 std::vector<std::string> components,
                 DoFSetInfo _info,
                 std::unordered_set<std::string> _excluded_occupants)
    : m_traits(_type), m_info(_info), m_excluded_occs(std::move(_excluded_occupants)) {

    if(components.empty()) {
      components = traits().standard_var_names();
      m_info.set_basis(Eigen::MatrixXd::Identity(traits().dim(), traits().dim()));
    }


    if(m_info.basis().rows() != traits().dim()) {
      throw std::runtime_error("Cannot construct DoFSet of type " + type_name() + ", number of rows in basis matrix is "
                               + std::to_string(m_info.basis().rows()) + " but must be "
                               + std::to_string(traits().dim()));
    }
    if(m_info.basis().cols() != components.size()) {
      throw std::runtime_error("Cannot construct DoFSet of type " + type_name() + ", number of columns in basis matrix is "
                               + std::to_string(m_info.basis().cols()) + " but number of variable names is "
                               + std::to_string(components.size()) + ". These numbers must be equal.");
    }

    for(std::string const &name : components)
      m_components.push_back(ContinuousDoF(traits(), name,
                                           -1, // ID
                                           -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()));
  }

  //********************************************************************

  bool DoFSet::identical(DoFSet const &rhs) const {
    if(!std::equal(begin(), end(), rhs.begin(), compare_no_value)) {
      return false;
    }

    return almost_equal(basis(), rhs.basis());
  }

  //********************************************************************

  void DoFSet::transform_basis(Eigen::Ref<const Eigen::MatrixXd> const &trans_mat) {
    m_info.set_basis(trans_mat * basis());
    m_info.set_symrep_ID(SymGroupRepID());
  }

  //********************************************************************

  void DoFSet::lock_IDs() {
    for(ContinuousDoF &dof : this->m_components) {
      dof.lock_ID();
    }
    return;
  }

  void DoFSet::set_sequential_IDs() {
    for(int i = 0; i < this->m_components.size(); ++i) {
      m_components[i].set_ID(i);
    }
    return;
  }

  bool DoFSet::update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs) {

    Index ID_ind;
    bool is_updated(false);
    for(Index i = 0; i < m_components.size(); i++) {
      // IMPORTANT: Do before_IDs.find(), NOT m_components().find() (if such a thing existed)
      ID_ind = find_index(before_IDs, m_components[i].ID());
      // Only set ID if DoF doesn't have an ID lock
      if(ID_ind < after_IDs.size() && !m_components[i].is_locked()) {
        m_components[i].set_ID(after_IDs[ID_ind]);
        // The new ID only changes the formula if the corresponding coeff is nonzero
        // if(!almost_zero(m_coeffs[i]))
        is_updated = true;
      }
    }
    return is_updated;
  }

  DoFSet DoFSet::make_default(DoFSet::BasicTraits const &_type) {
    DoFSet result(_type,
                  _type.standard_var_names(),
    {SymGroupRepID(), Eigen::MatrixXd::Identity(_type.dim(), _type.dim())});

    if(_type.global()) {
      Index i = 0;
      for(auto &dof : result.m_components) {
        dof.set_ID(i++);
        dof.lock_ID();
      }
    }
    return result;
  }
  //********************************************************************

} // namespace CASM
