#ifndef CASM_DoFSpace
#define CASM_DoFSpace

#include "casm/crystallography/DoFDecl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {
  struct VectorSpaceSymReport;

  class DoFSpaceSymReport;

  ///\brief Compositiong of enumeration environment (as ConfigEnumInput), DoF type (as DoFKey),
  // and vector subpsace. This simplifies generation of a VectorSpaceSymReport for a range of inputs
  struct DoFSpace {

    /// Definition of Configuration region for which enumeration/analysis is to be performed
    ConfigEnumInput config_region;

    /// Type of degree of freedom that is under consideration (e.g., "disp", "Hstrain", "occ")
    DoFKey dof_key;

    /// Columns of dof_subspace span a subspace of allowed DoF values in configuration
    /// Number of rows must equal dimension of intrinsic DoF vector space (e.g., 3*N for displacement)
    /// Number columns must be less than or equal to number of rows
    Eigen::MatrixXd dof_subspace;

    /// Construct with ConfigEnuminput, DoFKey (which specifies DoF type), and subspace matrix
    /// Constructing with Supercell or Configuration is allowed via implicit conversiont to ConfigEnumInput.
    /// By default, _dof_subspace matrix (i.e., an empty matrix) results in dof_subspace being initialized to
    /// Identity matrix having same dimension as complete DoF vector space of the provided ConfigEnumInput.
    /// Columns of _dof_subspace are vectors in the DoF vector space of the provided ConfigEnumInput, and so
    /// the number of rows of _dof_subspace is determined by the ConfigEnumInput and DoF type.
    /// The rank of _dof_subspace must be equal to its number of columns (i.e., columns are linearly independent).
    DoFSpace(ConfigEnumInput _config_region,
             DoFKey _dof_key,
             Eigen::MatrixXd _dof_subspace = Eigen::MatrixXd());

  };

  Index get_dof_space_dimension(DoFKey dof_key,
                                Configuration const &configuration,
                                std::set<Index> const &sites);

  std::vector<std::string> make_axis_glossary(DoFKey dof_key,
                                              Configuration const &configuration,
                                              std::set<Index> const &sites);

  template<typename PermuteIteratorIt>
  VectorSpaceSymReport vector_space_sym_report(DoFSpace const &dof_space,
                                               PermuteIteratorIt group_begin,
                                               PermuteIteratorIt group_end,
                                               bool calc_wedges = true);

}

#endif
