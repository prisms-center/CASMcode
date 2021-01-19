#include "casm/clex/ConfigEnumStrain.hh"

#include <algorithm>

#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {

ConfigEnumStrain::ConfigEnumStrain(ConfigEnumInput const &initial_state,
                                   ConfigEnumStrainParams const &params)
    : ConfigEnumStrain(initial_state, params.wedges, params.min_val,
                       params.max_val, params.inc_val, params.dof,
                       params.auto_range, params.trim_corners) {}

ConfigEnumStrain::ConfigEnumStrain(
    ConfigEnumInput const &initial_state,
    std::vector<SymRepTools::SubWedge> const &wedges, Eigen::VectorXd min_val,
    Eigen::VectorXd max_val, Eigen::VectorXd inc_val, DoFKey const &strain_key,
    bool auto_range, bool trim_corners)
    : m_strain_key(strain_key),
      m_trim_corners(trim_corners),
      m_current(initial_state.configuration()),
      m_equiv_ind(0),
      m_wedges(wedges),
      m_shape_factor(
          Eigen::MatrixXd::Identity(min_val.size(), min_val.size())) {
  // Condition range arrays and build shape_factor matrix
  Index nc = 0;
  if (m_wedges.size() == 0) return;
  Index nsub = m_wedges[0].irrep_wedges().size();
  for (Index s = 0; s < nsub; s++) {
    double abs_mag = 0.;
    for (Index i = 0; i < m_wedges[0].irrep_wedges()[s].axes.cols(); i++)
      abs_mag = max<double>(abs_mag, max(abs(max_val(nc)), abs(min_val(nc))));

    for (Index i = 0; i < m_wedges[0].irrep_wedges()[s].axes.cols();
         i++, nc++) {
      if (m_wedges[0].irrep_wedges()[s].mult[i] == 1 && auto_range) {
        min_val(nc) = -max_val(nc);
      }

      if (abs_mag > TOL) m_shape_factor(nc, nc) /= abs_mag * abs_mag;

      if (inc_val(nc) == 0. ||
          (max_val(nc) - min_val(nc)) / inc_val(nc) > 1e4) {
        max_val(nc) = min_val(nc) + TOL;
        inc_val(nc) = 10 * TOL;
      }
    }
  }

  // Find first valid config
  m_counter = EigenCounter<Eigen::VectorXd>(min_val, max_val, inc_val);

  m_counter.reset();

  // Increment past any invalid values, including those that are outside
  // specified ellipsoid (if trim_corners==true)

  while (
      m_counter.valid() &&
      (trim_corners &&
       double(m_counter().transpose() *
              m_wedges[m_equiv_ind].trans_mat().transpose() * m_shape_factor *
              m_wedges[m_equiv_ind].trans_mat() * m_counter()) > 1.0 + TOL)) {
    ++m_counter;
  }

  reset_properties(m_current);
  m_current.configdof().set_global_dof(
      m_strain_key, m_wedges[m_equiv_ind].trans_mat() * m_counter());
  this->_initialize(&m_current);

  if (!m_counter.valid()) {
    this->_invalidate();
  }

  // auto &log = CASM::log();
  // log.subsection().begin_section<Log::debug>();
  // log.indent() << "normal coordinates: " << m_counter().transpose() <<
  // std::endl; log.indent() << "dof value: " <<
  // (m_wedges[m_equiv_ind].trans_mat() * m_counter()).transpose() << std::endl;
  // log.indent() << std::endl;
  // log.end_section();

  m_current.set_source(this->source(step()));
}

// Implements _increment
void ConfigEnumStrain::increment() {
  // Increment past any invalid values, including those that are outside
  // specified ellipsoid (if trim_corners==true)
  while (
      ++m_counter &&
      (m_trim_corners &&
       double(m_counter().transpose() *
              m_wedges[m_equiv_ind].trans_mat().transpose() * m_shape_factor *
              m_wedges[m_equiv_ind].trans_mat() * m_counter()) > 1.0 + TOL)) {
  }

  // move to next part of wedge if necessary
  if (!m_counter.valid() && m_equiv_ind + 1 < m_wedges.size()) {
    m_counter.reset();
    ++m_equiv_ind;

    // Increment past any invalid values, including those that are outside
    // specified ellipsoid (if trim_corners==true)
    // this time it's for the new wedge

    while (
        m_counter &&
        (m_trim_corners &&
         double(m_counter().transpose() *
                m_wedges[m_equiv_ind].trans_mat().transpose() * m_shape_factor *
                m_wedges[m_equiv_ind].trans_mat() * m_counter()) > 1.0 + TOL)) {
      ++m_counter;
    }
  }

  if (m_counter.valid()) {
    // auto &log = CASM::log();
    // log.subsection().begin_section<Log::debug>();
    // log.indent() << "normal coordinates: " << m_counter().transpose() <<
    // std::endl; log.indent() << "dof value: " <<
    // (m_wedges[m_equiv_ind].trans_mat() * m_counter()).transpose() <<
    // std::endl; log.indent() << std::endl; log.end_section();

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
