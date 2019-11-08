#include <set>
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace xtal {
    // This is a helper function that isn't part of the Niggli interface, and is
    // defined in Niggli.cc
    std::set<Eigen::Matrix3d, StandardOrientationCompare>
    _niggli_set(const Lattice &in_lat, double compare_tol, bool keep_handedness);

    namespace canonical {
      /// Return canonical equivalent lattice, and 'to_canonical' SymOp
      ///
      /// The 'to_canonical' SymOp, to_canonical_op, satisfies:
      ///   canonical_lat == niggli(to_canonical_op.matrix() * ref_lat.lat_column_mat())
      /// to within the specified tolerance.
      ///   where ref_lat = niggli(in_lat, compare_tol, false)
      ///
      std::pair<Lattice, Index>
      _equivalent_lattice_and_index(const Lattice &in_lat, const SymOpVector &point_grp, double compare_tol) {

        auto lat_set = _niggli_set(in_lat, compare_tol, true);
        if(lat_set.empty())
          throw std::runtime_error(
            "In _canonical_equivalent_lattice(), did not find any niggli representations of provided lattice.");
        bool first = true;
        Eigen::Matrix3d most_canonical_lat_mat, trans_lat_mat;
        Index to_canonical_ix = 0;

        // The niggli cell is based purely on the metric tensor, and so is invariant cartesian rotation
        // (i.e., the niggli cell uniquely orders lattice parameters and angles but does not uniquely specify orientation)
        // We find all niggli representations for 'in_lat', and then loop over point group operations.
        // For each point group operation, reorient all the niggli representations, keeping track of
        // which orientation is 'most' canonical, according to casm standard orientation
        Index ix = 0;
        for(auto it = point_grp.begin(); it != point_grp.end(); ++it, ++ix) {
          // Skip operations that change the handedness of the lattice
          if(get_matrix(*it).determinant() <= 0.0) {
            continue;
          }

          for(Eigen::MatrixXd const &lat_mat : lat_set) {
            trans_lat_mat = get_matrix(*it) * lat_mat;
            assert(is_niggli(lat_mat, compare_tol) && "Result of 'niggli()' is not a Niggli cell");

            if(first || standard_orientation_compare(most_canonical_lat_mat, trans_lat_mat, compare_tol)) {
              first = false;
              most_canonical_lat_mat = trans_lat_mat;
              to_canonical_ix = ix;
            }
          }
        }

        return std::make_pair(Lattice(most_canonical_lat_mat, in_lat.tol()), to_canonical_ix);
      }

      Lattice equivalent(const Lattice &in_lat, const SymOpVector &point_grp, double compare_tol) {
        return _equivalent_lattice_and_index(in_lat, point_grp, compare_tol).first;
      }

      Lattice equivalent(const Lattice &lat, SymOpVector const &g) {
        return canonical::equivalent(lat, g, lat.tol());
      }

      Lattice equivalent(const Lattice &lat) {
        return canonical::equivalent(lat, make_point_group(lat));
      }

      bool check(const Lattice &lat, SymOpVector const &g) {
        return almost_equal(lat.lat_column_mat(), canonical::equivalent(lat, g).lat_column_mat(), lat.tol());
      }

      bool check(const Lattice &lat) {
        return canonical::check(lat, make_point_group(lat));
      }

      Index operation_index(const Lattice &in_lat, const SymOpVector &point_grp, double compare_tol) {
        return _equivalent_lattice_and_index(in_lat, point_grp, compare_tol).second;
      }

      Index operation_index(const Lattice &lat, SymOpVector const &g) {
        return canonical::operation_index(lat, g, lat.tol());
      }

    } // namespace canonical

  } // namespace xtal
} // namespace CASM
