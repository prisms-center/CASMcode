#ifndef CASM_ConfigEnumStrain
#define CASM_ConfigEnumStrain

#include "casm/clex/Configuration.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  namespace SymRepTools {
    class SubWedge;
  }

  class ConfigEnumInput;

  struct ConfigEnumStrainParams {

    /// Type of strain
    DoFKey dof;

    /// The positive coordinate space defined by the SubWedge axes are a symmetrically unique portion
    /// of the strain space. It may require multiple SubWedge to fill strain space, hence the vector.
    ///
    /// Note:
    /// - Symmetry adapted wedges can be obtained from VectorSpaceSymReport::irreducible_wedge
    /// - The "SubWedge::make_dummy" static function can be used to create a fake SubWedge to
    ///   generate strain states directly on custom axes
    std::vector<SymRepTools::SubWedge> wedges;

    // The options "min_val", "inc_val", "max_val" defines a grid of points to sample along each
    // wedge axis.
    // Must have same dimension as the number of subwedge axes columns.

    Eigen::VectorXd min_val;
    Eigen::VectorXd max_val;
    Eigen::VectorXd inc_val;

    // If true, set min_val[i]=-max_val[i], for strain components, i, with multiplicity==1
    // Typically, set to true if using symmetry adapted axes
    bool auto_range = true;

    /// If true, fit an ellipsoid inside the grid defined by min_val/max_val and exclude grid points
    /// that lie outside the ellipsoid
    bool trim_corners = true;

  };


  /// Enumerate strained Configurations
  ///
  /// \ingroup ConfigEnumGroup
  ///
  class ConfigEnumStrain : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumStrain(ConfigEnumInput const &initial_state,
                     ConfigEnumStrainParams const &params);

    ConfigEnumStrain(ConfigEnumInput const &initial_state,
                     std::vector<SymRepTools::SubWedge> const &wedges,
                     Eigen::VectorXd min_val,
                     Eigen::VectorXd max_val,
                     Eigen::VectorXd inc_val,
                     DoFKey const &strain_key,
                     bool auto_range,
                     bool trim_corners);


    std::string name() const override;

    static const std::string enumerator_name;

  private:

    /// Implements increment over all strain states
    void increment() override;


    // -- Unique -------------------
    DoFKey m_strain_key;

    bool m_trim_corners;

    Configuration m_current;

    // counts over strain grid
    EigenCounter<Eigen::VectorXd> m_counter;

    // counts over transformation matrices
    Index m_equiv_ind;

    std::vector<SymRepTools::SubWedge> m_wedges;

    //set of non-equivalent transformation matrices matrices that, along with m_counter define irreducible space
    //std::vector<Eigen::MatrixXd> m_trans_mats;
    PermuteIterator m_perm_begin, m_perm_end;
    Eigen::MatrixXd m_shape_factor;

    const PermuteIterator &_perm_begin() {
      return m_perm_begin;
    }
    const PermuteIterator &_perm_end() {
      return m_perm_end;
    }

  };

}

#endif
