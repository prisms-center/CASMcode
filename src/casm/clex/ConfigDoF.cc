#include "casm/basis_set/Adapter.hh"
#include "casm/clex/ConfigDoFValues.hh"
#include "casm/clex/ConfigDoF_impl.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymPermutation.hh"

namespace CASM {

/// Number of sites in the ConfigDoF
Index ConfigDoF::size() const { return occupation().size(); }

/// Integer volume of ConfigDoF
Index ConfigDoF::n_vol() const { return m_occupation.n_vol(); }

/// Number of sublattices in ConfigDoF
Index ConfigDoF::n_sublat() const { return m_occupation.n_sublat(); }

/// Tolerance for comparison of continuous DoF values
double ConfigDoF::tol() const { return m_tol; }

/// Set all DoF values to zero
void ConfigDoF::setZero() {
  m_occupation.setZero();
  for (auto &local_dof : m_local_dofs) {
    local_dof.second.setZero();
  }
  for (auto &global_dof : m_global_dofs) {
    global_dof.second.setZero();
  }
}

/// Set occupation values
///
/// \throws std::runtime_error ("Size mismatch in ConfigDoF::set_occupation...")
/// if input vector size does not match current size
void ConfigDoF::set_occupation(
    Eigen::Ref<const Eigen::VectorXi> const &_occupation) {
  if (occupation().size() != _occupation.size()) {
    std::stringstream msg;
    msg << "Size mismatch in ConfigDoF::set_occupation(): "
        << "Expected size=" << occupation().size()
        << ", received size=" << _occupation.size();
    throw std::runtime_error(msg.str());
  }
  m_occupation.set_values(_occupation);
}

bool ConfigDoF::has_occupation() const {
  return size() != 0 && occupation().size() == size();
}

std::map<DoFKey, GlobalContinuousConfigDoFValues> const &
ConfigDoF::global_dofs() const {
  return m_global_dofs;
}

GlobalContinuousConfigDoFValues const &ConfigDoF::global_dof(
    DoFKey const &_key) const {
  auto it = m_global_dofs.find(_key);
  if (it == m_global_dofs.end())
    throw std::runtime_error(
        "Attempting to access uninitialized ConfigDoF value for '" + _key +
        "'");
  return it->second;
}

GlobalContinuousConfigDoFValues &ConfigDoF::global_dof(DoFKey const &_key) {
  auto it = m_global_dofs.find(_key);
  if (it == m_global_dofs.end())
    throw std::runtime_error(
        "Attempting to access uninitialized ConfigDoF value for '" + _key +
        "'");
  return it->second;
}

bool ConfigDoF::has_global_dof(DoFKey const &_key) const {
  return global_dofs().count(_key);
}

/// Set global continuous DoF values
///
/// \throws std::runtime_error ("Attempting to assign global ConfigDoF
/// values...") if "_key" does not exist in global dofs
///
/// \throws std::runtime_error ("Size mismatch in ConfigDoF::set_global_dof...")
/// if input vector size does not match current size
void ConfigDoF::set_global_dof(DoFKey const &_key,
                               Eigen::Ref<const Eigen::VectorXd> const &_val) {
  auto it = m_global_dofs.find(_key);
  if (it == m_global_dofs.end()) {
    throw std::runtime_error(
        "Attempting to assign global ConfigDoF values of type " + _key +
        " but no such value type has been allocated.\n");
  }
  if (it->second.values().size() != _val.size()) {
    std::stringstream msg;
    msg << "Size mismatch in ConfigDoF::set_global_dof(): "
        << "For key=\"" << _key
        << "\", expected size=" << it->second.values().size()
        << ", received size=" << _val.size();
    throw std::runtime_error(msg.str());
  }

  it->second.set_values(_val);
}

std::map<DoFKey, LocalContinuousConfigDoFValues> const &ConfigDoF::local_dofs()
    const {
  return m_local_dofs;
}

LocalContinuousConfigDoFValues const &ConfigDoF::local_dof(
    DoFKey const &_key) const {
  auto it = m_local_dofs.find(_key);
  if (it == m_local_dofs.end())
    throw std::runtime_error(
        "Attempting to access uninitialized ConfigDoF value for '" + _key +
        "'");
  return it->second;
}

LocalContinuousConfigDoFValues &ConfigDoF::local_dof(DoFKey const &_key) {
  auto it = m_local_dofs.find(_key);
  if (it == m_local_dofs.end())
    throw std::runtime_error(
        "Attempting to access uninitialized ConfigDoF value for '" + _key +
        "'");
  return it->second;
}

bool ConfigDoF::has_local_dof(DoFKey const &_key) const {
  return local_dofs().count(_key);
}

/// Set local continuous DoF values
///
/// \throws std::runtime_error ("Attempting to assign local ConfigDoF
/// values...") if "_key" does not  exist in local dofs
//
/// \throws std::runtime_error ("Size mismatch in ConfigDoF::set_local_dof...")
/// if input matrix size does not match current matrix size
void ConfigDoF::set_local_dof(DoFKey const &_key,
                              Eigen::Ref<const Eigen::MatrixXd> const &_val) {
  auto it = m_local_dofs.find(_key);
  if (it == m_local_dofs.end()) {
    throw std::runtime_error(
        "Attempting to assign local ConfigDoF values of type " + _key +
        " but no such value type has been allocated.\n");
  }

  if (it->second.values().rows() != _val.rows() ||
      it->second.values().cols() != _val.cols()) {
    std::stringstream msg;
    msg << "Size mismatch in ConfigDoF::set_local_dof(): "
        << "For key=\"" << _key
        << "\", expected rows=" << it->second.values().rows()
        << ", received rows=" << _val.rows()
        << ", expected cols=" << it->second.values().cols()
        << ", received cols=" << _val.cols();
    throw std::runtime_error(msg.str());
  }

  it->second.set_values(_val);
}

/// Update DoF values using the effect of symmetry, including permutation
/// among sites
///
/// The effect is: *this = permute_iterator * (*this)
ConfigDoF &ConfigDoF::apply_sym(PermuteIterator const &it) {
  for (auto &dof : m_global_dofs) {
    dof.second.set_values(*(it.global_dof_rep(dof.first).MatrixXd()) *
                          dof.second.values());
  }

  Permutation tperm(it.combined_permute());
  if (occupation().size()) {
    if (it.sym_info().has_aniso_occs()) {
      Index l = 0;
      for (Index b = 0; b < n_sublat(); ++b) {
        for (Index n = 0; n < n_vol(); ++n, ++l) {
          occ(l) = (*(it.occ_rep(b).permutation()))[occ(l)];
        }
      }
    }
    set_occupation(tperm * occupation());
  }

  for (auto &dof : m_local_dofs) {
    LocalContinuousConfigDoFValues tmp = dof.second;

    for (Index b = 0; b < tmp.n_sublat(); ++b) {
      Index rows = it.local_dof_rep(dof.first, b).MatrixXd()->rows();
      tmp.sublat(b).topRows(rows) =
          *(it.local_dof_rep(dof.first, b).MatrixXd()) *
          dof.second.sublat(b).topRows(rows);
    }
    for (Index l = 0; l < size(); ++l) {
      dof.second.site_value(l) = tmp.site_value(tperm[l]);
    }
  }

  return *this;
}

/// Update DoF values using only the effect of symmetry on the value at each
/// site, without permutation among sites
ConfigDoF &ConfigDoF::apply_sym_no_permute(SymOp const &_op) {
  for (auto &dof : m_global_dofs) {
    dof.second.set_values(
        *(_op.representation(dof.second.info().symrep_ID()).MatrixXd()) *
        dof.second.values());
  }

  if (occupation().size()) {
    Index l = 0;
    for (Index b = 0; b < n_sublat(); ++b) {
      if (!m_occupation.symrep_IDs()[b].is_identity()) {
        SymPermutation const &permrep(
            *_op.get_permutation_rep(m_occupation.symrep_IDs()[b]));
        l = b * n_vol();
        for (Index n = 0; n < n_vol(); ++n, ++l) {
          occ(l) = (*permrep.permutation())[occ(l)];
        }
      }
    }
  }

  for (auto &dof : m_local_dofs) {
    LocalContinuousConfigDoFValues tmp = dof.second;

    for (Index b = 0; b < tmp.n_sublat(); ++b)
      dof.second.sublat(b) =
          *(_op.representation(dof.second.info()[b].symrep_ID()).MatrixXd()) *
          dof.second.sublat(b);
  }

  return *this;
}

void ConfigDoF::swap(ConfigDoF &RHS) {
  std::swap(m_occupation, RHS.m_occupation);
  std::swap(m_local_dofs, RHS.m_local_dofs);
  std::swap(m_global_dofs, RHS.m_global_dofs);
  std::swap(m_tol, RHS.m_tol);
}

void swap(ConfigDoF &A, ConfigDoF &B) { A.swap(B); }

}  // namespace CASM
