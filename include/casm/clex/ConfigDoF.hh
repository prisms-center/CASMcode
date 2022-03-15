#ifndef ConfigDoF_HH
#define ConfigDoF_HH

#include <vector>

#include "casm/clex/ConfigDoFValues.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

class PermuteIterator;
class SymGroupRepID;
class SymOp;

/// A container class for the different degrees of freedom values a
/// Configuration might have
///
/// \ingroup Configuration
///
/// ConfigDoF has:
/// - The values of the discrete site DoF within the supercell ("occupation",
/// Eigen::VectorXi)
///
///   Example: Occupation values, accessed via `Eigen::VectorXi const
///   &ConfigDoF::occupation()`, with integer value corresponding to which
///   Molecule in the the a Site::occupant_dof vector is occupying a particular
///   site:
///
///       [<- sublattice 0 "occ" values -> | <- sublattice 1 "occ" values -> |
///       ... ]
///
/// - The values of the continuous site DoF within the supercell ("local_dofs",
///   std::map<DoFKey, LocalContinuousConfigDoFValues>)
///
///   Example: Displacement values, with prim DoF basis equal to the standard
///   basis (dx, dy, dz), accessed via `Eigen::MatrixXd const
///   &ConfigDoF::local_dofs("disp").values()`:
///
///       [<- sublattice 0 dx values -> | <- sublattice 1 dx values -> | ... ]
///       [<- sublattice 0 dy values -> | <- sublattice 1 dy values -> | ... ]
///       [<- sublattice 0 dz values -> | <- sublattice 1 dz values -> | ... ]
///
///   Example: Displacement values, with non-standard prim DoF basis:
///
///       "basis" : [ {
///           "coordinate": [c0x, c0y, c0z],
///           "occupants": [...],
///           "dofs": {
///             "disp" : {
///               "axis_names" : ["dxy", "dz"],
///               "axes" : [[1.0, 1.0, 0.0],
///                       [0.0, 0.0, 1.0]]}}
///         },
///         {
///           "coordinate": [c1x, c1y, c1z],
///           "occupants": [...],
///           "dofs": {
///             "disp" : {
///               "axis_names" : ["d\bar{x}y", "dz"],
///               "axes" : [[-1.0, 1.0, 0.0],
///                       [0.0, 0.0, 1.0]]}}
//          },
///         ...
///       }
///
///       [<- sublat 0 dxy values -> | <- sublat 1 d\bar{x}y values ->| ... ]
///       [<- sublat 0 dz values  -> | <- sublat 1 dz values ->       | ... ]
///
///       Note that the values matrix has only two rows, this is the maximum
///       site basis dimension.
///
///   Example: Displacement values, with varying prim DoF basis:
///
///       "basis" : [ {
///           "coordinate": [c0x, c0y, c0z],
///           "occupants": [...],
///           "dofs": {
///             "disp" : {}
///         },
///         {
///           "coordinate": [c1x, c1y, c1z],
///           "occupants": [...],
///           "dofs": {
///             "disp" : {
///               "axis_names" : ["dxy", "dz"],
///               "axes" : [[-1.0, 1.0, 0.0],
///                       [0.0, 0.0, 1.0]]}}
//          },
///         ...
///       }
///
///       [<- sublattice 0 dx values -> | <- sublattice 1 dxy values ->| ... ]
///       [<- sublattice 0 dy values -> | <- sublattice 1 dz values -> | ... ]
///       [<- sublattice 0 dz values -> | <- 0.0 values ->             | ... ]
///
///       Note that the values matrix has three rows, this is the maximum
///       site basis dimension, but for sublattices with lower site basis
///       dimension it is padded with fixed zeros.
///
/// - The values of the continuous global DoF ("global_dofs",
///   std::map<DoFKey, GlobalContinuousConfigDoFValues>).
///
///   Example: GLstrain values, with prim DoF basis equal to the standard basis,
///   accessed via `Eigen::VectorXd const
///   &ConfigDoF::global_dofs("GLstrain").values()`:
///
///       [e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy]
///
/// Note: Continuous DoF values stored in memory in ConfigDoF are coordinates in
/// the prim DoF basis (the "axes" given in "prim.json" which set
/// xtal::SiteDoFSet::basis() or xtal::DoFSet::basis()), but when saved to file
/// (i.e. `.casm/config/config_list.json`) they are saved as coordinates in the
/// standard DoF basis (with axes meaning as described by
/// AnisoValTraits::standard_var_names()).
///
class ConfigDoF {
 public:
  // Can treat as a Eigen::VectorXd
  // typedef displacement_matrix_t::ColXpr displacement_t;
  // typedef displacement_matrix_t::ConstColXpr const_displacement_t;

  /// Initialize with number of sites, and dimensionality of global and local
  /// DoFs GlobalInfoContainerType is an iterable container having value_type
  /// std::pair<DoFKey,ContinuousDoFInfo> LocalInfoContainerType is an iterable
  /// container having value_type
  /// std::pair<DoFKey,std::vector<ContinuousDoFInfo>  > OccInfoContainerType is
  /// an iterable container having value_type SymGroupRepID
  template <typename GlobalInfoContainerType, typename LocalInfoContainerType>
  ConfigDoF(Index _N_sublat, Index _N_vol,
            GlobalInfoContainerType const &global_dof_info,
            LocalInfoContainerType const &local_dof_info,
            std::vector<SymGroupRepID> const &occ_symrep_IDs, double _tol);

  ConfigDoF(ConfigDoF const &RHS);

  ConfigDoF &operator=(ConfigDoF const &RHS);

  /// Number of sites in the ConfigDoF
  Index size() const;

  /// Integer volume of ConfigDoF
  Index n_vol() const;

  /// Number of sublattices in ConfigDoF
  Index n_sublat() const;

  /// Tolerance for comparison of continuous DoF values
  double tol() const;

  /// Set all DoF values to zero
  void setZero();

  /// Reference occupation value on site i
  int &occ(Index i) {  // keep inline
    return m_dof_values.occupation(i);
  }

  /// Const reference to occupation value on site i
  const int &occ(Index i) const {  // keep inline
    return m_dof_values.occupation(i);
  }

  /// Set occupation values
  void set_occupation(Eigen::Ref<const Eigen::VectorXi> const &_occupation);

  /// Const reference occupation values
  Eigen::VectorXi const &occupation() const {  // keep inline
    return m_dof_values.occupation;
  }

  bool has_occupation() const;

  std::map<DoFKey, GlobalContinuousConfigDoFValues> const &global_dofs() const;

  GlobalContinuousConfigDoFValues const &global_dof(DoFKey const &_key) const;

  GlobalContinuousConfigDoFValues &global_dof(DoFKey const &_key);

  bool has_global_dof(DoFKey const &_key) const;

  /// Set global continuous DoF values
  void set_global_dof(DoFKey const &_key,
                      Eigen::Ref<const Eigen::VectorXd> const &_val);

  std::map<DoFKey, LocalContinuousConfigDoFValues> const &local_dofs() const;

  LocalContinuousConfigDoFValues const &local_dof(DoFKey const &_key) const;

  LocalContinuousConfigDoFValues &local_dof(DoFKey const &_key);

  bool has_local_dof(DoFKey const &_key) const;

  /// Set local continuous DoF values
  void set_local_dof(DoFKey const &_key,
                     Eigen::Ref<const Eigen::MatrixXd> const &_val);

  /// Update DoF values using the effect of symmetry, including permutation
  /// among sites
  ConfigDoF &apply_sym(PermuteIterator const &it);

  /// Update DoF values using only the effect of symmetry on the value at each
  /// site, without permutation among sites
  ConfigDoF &apply_sym_no_permute(SymOp const &_op);

  /// Access underlying data structure. *Do not resize DoF value containers*.
  clexulator::ConfigDoFValues &values() { return m_dof_values; }

  /// Access underlying data structure.
  clexulator::ConfigDoFValues const &values() const { return m_dof_values; }

 private:
  // populate m_local_dofs and m_global_dofs to reference this object's
  // m_dof_values
  void _make_continuous_dof_values(ConfigDoF const &RHS);


  Index m_N_sublat;

  Index m_N_vol;

  // this contains raw data structures, without all the methods that
  // apply symmetry and converting to/from standard values
  clexulator::ConfigDoFValues m_dof_values;

  std::vector<SymGroupRepID> m_occ_symrep_IDs;

  std::map<std::string, GlobalContinuousConfigDoFValues> m_global_dofs;

  std::map<std::string, LocalContinuousConfigDoFValues> m_local_dofs;

  /// Tolerance used for transformation to canonical form -- used also for
  /// continuous DoF comparisons, since comparisons are only meaningful to
  /// within the tolerance used for finding the canonical form.
  mutable double m_tol;
};


inline void reset_properties(ConfigDoF &_dof) { return; }

}  // namespace CASM

#endif
