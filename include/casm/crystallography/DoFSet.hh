#ifndef CASM_XTAL_DOFSET
#define CASM_XTAL_DOFSET

#include "casm/crystallography/AnisoValTraits.hh"
#include <unordered_set>

namespace CASM {
  namespace xtal {
    /**
     * SiteDoFSet specifies all identifying information for a vector of continuous independent variables (Degrees of Freedom / DoFs)
     * DoFSets are associated with a specific DoF 'type', which has a predefined 'standard' coordinate system
     * ex: displacement -> 3-vector (x,y,z) -> displacement components (relative to fixed laboratory frame)
     *     strain -> 6-vector (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy) -> tensor elements
     * DoFSets have a typename, which specifies the type, and a set of basis vectors, which are denoted relative
     * to the DoF type's standard axes. This allows the DoFSet components to be specified by the user,
     * including the ability to only allow DoF values within a subspace of the standard values.
     * DoFSet records  the DoF typename, the names of the vector components, and the axes of the
     * vector components (relative to a set of standard axes)
     */
    class DoFSet {
    public:
      using BasicTraits = AnisoValTraits;

      /// \brief Returns type_name of DoFSet, which should be a standardized DoF type (e.g., "disp", "magspin", "GLstrain")
      const std::string &type_name() const {
        return m_traits.name();
      }

      /// \brief Returns traits object for the DoF type of this DoFSet
      BasicTraits const &traits() const {
        return m_traits;
      }

      /// Returns the number of dimension of the DoF, corresponding to the number of axes in the vector space
      Index dim() const {
        return this->basis().cols();
      }

      /// Matrix that relates DoFSet variables to a conventional coordiante system
      Eigen::MatrixXd const &basis() const {
        return m_basis;
      }

      /// \brief returns true of \param rhs has identical components and basis to this DoFSet
      bool is_identical(DoFSet const &rhs) const;

    private:
      /// AnisoValTraits. Describes the type of DoF, and can convert Cartesian symmetry representations into the appropriate representation
      BasicTraits m_traits;

      /// Names for each axis of the basis, for example "x". "y", "z" for displacement
      std::vector<std::string> m_component_names;

      /// The basis defines the space of the DoF, which should be a linear combination of the AnisoValTraits conventional coordinates.
      /// For example, you may want to define displacements that only happen along a particular direction
      Eigen::MatrixXd m_basis;
      Eigen::MatrixXd m_basis_inverse;
    };

    /**
     * Identical to xtal::DoFSet, but also keeps track of a list of molecule names that the DoFSet does not
     * apply to. For example, don't apply displacements to a vacancy.
     */

    class SiteDoFSet : public DoFSet {
    public:
      /// \brief Returns true if DoFSet is inactive (e.g., takes zero values) when specified occupant is present
      bool is_excluded_occ(std::string const &_occ_name) const {
        return m_excluded_occs.count(_occ_name);
      }

    private:
      std::unordered_set<std::string> m_excluded_occs;
    };

  } // namespace xtal
} // namespace CASM

#endif
