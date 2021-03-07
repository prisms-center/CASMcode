#include "casm/clex/ConfigEnumStrain.hh"

#include <algorithm>

#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/symmetry/IrrepWedge.hh"

namespace CASM {

namespace {

// Increment past any values outside ellipsoid defined by shape_factor
// matrix (if trim_corners==true)
void advance_past_corners(EigenCounter<Eigen::VectorXd> &counter,
                          Eigen::MatrixXd shape_factor) {
  while (counter.valid() &&
         double(counter().transpose() * shape_factor * counter()) > 1.0 + TOL) {
    ++counter;
  }
}

}  // namespace

ConfigEnumStrain::ConfigEnumStrain(ConfigEnumInput const &initial_state,
                                   ConfigEnumStrainParams const &params)
    : ConfigEnumStrain(initial_state, params.wedges, params.min_val,
                       params.max_val, params.inc_val, params.dof,
                       params.auto_range, params.trim_corners) {}

ConfigEnumStrain::ConfigEnumStrain(
    ConfigEnumInput const &initial_state,
    std::vector<SymRepTools_v2::SubWedge> const &wedges,
    Eigen::VectorXd min_val, Eigen::VectorXd max_val, Eigen::VectorXd inc_val,
    DoFKey const &strain_key, bool auto_range, bool trim_corners)
    : m_strain_key(strain_key),
      m_trim_corners(trim_corners),
      m_current(initial_state.configuration()),
      m_equiv_ind(0),
      m_wedges(wedges),
      m_shape_factor(
          Eigen::MatrixXd::Identity(min_val.size(), min_val.size())) {
  // Check there are wedges to counter over
  if (m_wedges.size() == 0) return;

  // Check for dimensional mismatches
  if (min_val.size() != m_wedges[0].trans_mat().cols() ||
      min_val.size() != inc_val.size() || min_val.size() != max_val.size()) {
    std::stringstream msg;
    msg << "Error in ConfigEnumStrain: dimension mismatch";
    throw std::runtime_error(msg.str());
  }

  // Check for very small inc_val that would cause problems. Currently it is
  // hard-coded that there should not be >1e4 points along a dimension
  for (int d = 0; d < inc_val.size(); ++d) {
    if (inc_val(d) == 0. || (max_val(d) - min_val(d)) / inc_val(d) > 1e4) {
      std::stringstream msg;
      msg << "Error in ConfigEnumStrain: Increment along dimension " << d + 1
          << " is too small: " << inc_val(d);
      throw std::runtime_error(msg.str());
    }
  }

  // If auto_range==true, set min_val[i]=-max_val[i], for strain components, i,
  // with multiplicity==1. Typically, set to true if using symmetry adapted axes
  if (auto_range) {
    int d = 0;
    for (auto const &irrep_wedge : m_wedges[0].irrep_wedges()) {
      for (int i = 0; i < irrep_wedge.axes.cols(); ++i) {
        if (irrep_wedge.mult[i] == 1) {
          min_val(d) = -max_val(d);
        }
        d++;
      }
    }
  }

  // If trim_corners==true, set shape_factor matrix to exclude grid points
  // that lie outside an ellipsoid defined by most extreme min_val/max_val
  // Excluded if:
  //     m_counter().transpose() * m_shape_factor * m_counter() > 1.0 + TOL
  if (m_trim_corners) {
    int dim = min_val.size();
    for (int d = 0; d < dim; ++d) {
      double mag = max(abs(min_val(d)), abs(max_val(d)));
      m_shape_factor(d, d) = 1. / (mag * mag);
    }
  }

  // Find first valid config
  m_counter = EigenCounter<Eigen::VectorXd>(min_val, max_val, inc_val);
  m_counter.reset();
  if (m_trim_corners) advance_past_corners(m_counter, m_shape_factor);

  reset_properties(m_current);
  m_current.configdof().set_global_dof(
      m_strain_key, m_wedges[m_equiv_ind].trans_mat() * m_counter());
  this->_initialize(&m_current);

  if (!m_counter.valid()) {
    this->_invalidate();
  }

  m_current.set_source(this->source(step()));
}

// Implements _increment
void ConfigEnumStrain::increment() {
  ++m_counter;
  if (m_trim_corners) advance_past_corners(m_counter, m_shape_factor);

  // Move to next part of wedge if necessary
  if (!m_counter.valid() && m_equiv_ind + 1 < m_wedges.size()) {
    m_counter.reset();
    ++m_equiv_ind;

    if (m_trim_corners) advance_past_corners(m_counter, m_shape_factor);
  }

  if (m_counter.valid()) {
    m_current.configdof().set_global_dof(
        m_strain_key, m_wedges[m_equiv_ind].trans_mat() * m_counter());

    _increment_step();
  } else {
    _invalidate();
  }

  m_current.set_source(this->source(step()));
  return;
}

std::string ConfigEnumStrain::name() const { return enumerator_name; }

const std::string ConfigEnumStrain::enumerator_name = "ConfigEnumStrain";

Index ConfigEnumStrain::subwedge_index() const { return m_equiv_ind; }

Eigen::VectorXd ConfigEnumStrain::normal_coordinate() const {
  return m_counter();
}

}  // namespace CASM
