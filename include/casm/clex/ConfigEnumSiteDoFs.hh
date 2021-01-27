#ifndef CASM_ConfigEnumSiteDoFs
#define CASM_ConfigEnumSiteDoFs

#include "casm/clex/Configuration.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"

namespace CASM {

class ConfigEnumInput;

/// Parameters controlling ConfigEnumSiteDoFs
///
/// ConfigEnumSiteDoFs is a method to generate Configurations that represent a
/// sampling of a combination of subspaces of site (local, continuous) DoF.
///
struct ConfigEnumSiteDoFsParams {
  /// Name of site degree of freedom which will be sampled
  ///
  /// DoFKey is a typedef for std::string
  DoFKey dof;

  /// Axes defining normal coordinates used to sample the site DoF space:
  ///
  ///     site_DoF_values = axes * normal_coordinates
  ///
  /// Note:
  /// - Number of rows equals the total dimensionality of the space
  ///   - If the dimensionality of the site DoF is the same on every site, then
  ///   the total
  ///     dimensionality is (number of sites in the configuration / supercell) x
  ///     (dimension of the site DoF value vector representation).
  ///   - DoF values are unrolled in row-major order, such that values from a
  ///   particular site are
  ///     listed in contiguous rows.
  /// - Each column represents a normal mode (collective excitation) which will
  /// be applied to the
  ///   initial enumeration state in linear combinations to construct
  ///   configurations
  /// - The axes may be rank deficient, indicating the subspace of the total
  /// site DoF that will
  ///   be sampled.
  ///   - The number of columns equals the dimensionality of the subspace being
  ///   sampled
  ///
  ///   Example:
  ///     [[1,  1,  1]
  ///      [1, -1, -1]
  ///      [1,  0,  0]
  ///      [1, -1,  1]
  ///      [1,  1, -1]
  ///      [1,  0,  0]]
  ///
  /// NOTE:
  /// - This "axes" is the transpose of the JSON input value used with `casm
  /// enum -m ConfigEnumSiteDoFs`
  /// - Use VectorSpaceSymReport to generate symmetry adapted axes and an
  /// "axis_glossary" which
  ///   describes which site DoF correspond to each row for a particular choice
  ///   of DoFSpace (choice of dof type and ConfigEnumInput).
  Eigen::MatrixXd axes;

  // The following specify min / increment / max amplitude of normal modes
  // (columns of axes) applied to the initial enumeration state

  /// Normal mode amplitude minimum values.
  /// Dimension must be equal to number of columns of "axes".
  ///
  ///     min_val[i]: the minimum amplitude for axes.col(i)
  ///
  Eigen::VectorXd min_val;

  /// Normal mode amplitude maximum values.
  /// Dimension must be equal to number of columns of "axes".
  ///
  ///     max_val[i]: the maximum amplitude for axes.col(i)
  ///
  Eigen::VectorXd max_val;

  /// Normal mode amplitude increment values.
  /// Dimension must be equal to number of columns of "axes".
  ///
  ///     inv_val[i]: the amplitude spacing for axes.col(i)
  ///
  Eigen::VectorXd inc_val;

  // Even if "axes" are rank deficient, the site DoF space defined by axes may
  // quickly become very high dimensional (number of sites x mean site DoF
  // dimension), so rather than sample the entire space, ConfigEnumSiteDoFs
  // perturbs the input configuration by applying a linear combination of normal
  // modes (collective excitations represented by columns of axes).

  // The parameters "min_nonzero" and "max_nonzero" specifies how many normal
  // mode amplitudes should be nonzero (inclusive range [min_nonzero,
  // max_nonzero]). The method generates all n choose k (n=site DoF space
  // dimension, k iterates through [min_nonzer, max_nonzero]) combinations of
  // normal modes in that range, and for each combination applies all the k
  // chosen normal modes with amplitudes specified by "min" / "increment" /
  // "max". Note that this may quickly become very large, depending on n, k, and
  // the range specified by "min" / "increment" / "max".;

  /// Minimum number of number modes included in the linear combination applied
  /// to the initial state
  Index min_nonzero;

  /// Maximum number of number modes included in the linear combination applied
  /// to the initial state
  Index max_nonzero;

  /// A boolean flag which lets you turn on or off homogeneous modes
  /// They are included by default and this flag should be used to turn them off
  bool exclude_homogeneous_modes;
};

/// Enumerate site (continuous local) DoFs
///
/// Note: See ConfigEnumSiteDoFsParams for method and parameter details
///
class ConfigEnumSiteDoFs : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  /// See `ConfigEnumSiteDoFsParams` for method and parameter details
  ConfigEnumSiteDoFs(ConfigEnumInput const &_in_config,
                     ConfigEnumSiteDoFsParams const &params);

  /// See `ConfigEnumSiteDoFsParams` for method and parameter details
  ConfigEnumSiteDoFs(ConfigEnumInput const &_init, DoFKey const &_dof,
                     Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                     Eigen::Ref<const Eigen::VectorXd> const &min_val,
                     Eigen::Ref<const Eigen::VectorXd> const &max_val,
                     Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                     bool exclude_homogeneous_modes, Index _min_nonzero,
                     Index _max_nonzero);

  std::string name() const override;

  static const std::string enumerator_name;

  /// Implements increment over all strain states
  void increment() override;

 private:
  void _set_dof();

  bool _increment_combo();

  bool _check_sparsity() const;

  bool _check_current() const;

  // -- Unique -------------------
  notstd::cloneable_ptr<Configuration> m_current;

  LocalContinuousConfigDoFValues *m_dof_vals;

  DoFKey m_dof_key;

  Index m_min_nonzero;

  Index m_max_nonzero;

  Eigen::MatrixXd m_axes;

  Eigen::VectorXd m_min;

  Eigen::VectorXd m_max;

  Eigen::VectorXd m_inc;

  bool m_unit_length;

  std::vector<Index> m_sites;

  std::vector<Index> m_dof_dims;

  bool m_subset_mode;

  bool m_exclude_homogeneous_modes;

  Index m_combo_index;

  std::vector<Index> m_combo;

  EigenCounter<Eigen::VectorXd> m_counter;
};

/// Checks if all the dof values at each site are the same
bool are_all_dof_vals_same(const std::vector<Eigen::VectorXd> &dof_vals);

/// Constructs a homogeneous mode space of the dof space represented in terms of
/// user defined axes. Returns a column vectors matrix
Eigen::MatrixXd make_homogeneous_mode_space(
    const std::vector<DoFSetInfo> &dof_info);

/// Returns true if the homogeneous modes are mixed within two irreps.
/// This will be used to issue a warning to the user to look at stuff more
/// carefully
bool are_homogeneous_modes_mixed_in_irreps(
    const Eigen::MatrixXd &axes, const Eigen::MatrixXd &homogeneous_mode_space);
}  // namespace CASM

#endif
