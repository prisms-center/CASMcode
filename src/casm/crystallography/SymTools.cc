#include "casm/crystallography/SymTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"

namespace {
  using CASM::xtal::Lattice;

  //TODO
  //I removed this from the library, because it seems like something that exists solely for generating the point group
  //If you need this routine somewhere else, go right ahead and put it back in the xtal namespace, and uncomment its
  //declaration in SymTools.hh
  /// Relax the vectors of the given lattice such that it obeys the symmetry of the given group,
  /// where the symmetry operations are given in fractional representations
  Lattice symmetrized_with_fractional(const Lattice &lattice, const std::vector<Eigen::Matrix3i> &fractional_point_group) {
    Eigen::Matrix3d tLat2(Eigen::Matrix3d::Zero());

    for(const auto &frac_mat : fractional_point_group) {
      tLat2 += frac_mat.cast<double>().transpose() * lattice.lat_column_mat().transpose() * lattice.lat_column_mat() * frac_mat.cast<double>();
    }
    tLat2 /= double(fractional_point_group.size());

    // tLat2 has the symmetrized lengths and angles -- it is equal to L.transpose()*L, where L=lat_column_mat()
    // we will find the sqrt of tLat2 and then reorient it so that it matches the original lattice
    Eigen::Matrix3d tMat;

    Eigen::JacobiSVD<Eigen::Matrix3d> tSVD(tLat2, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // tMat is sqrt of symlat.transpose()*symlat
    tMat = tSVD.matrixU() * tSVD.singularValues().cwiseSqrt().asDiagonal() * tSVD.matrixV().transpose();

    // if starting lattice were perfect, we would have origlat=rotation*tMat
    tSVD.compute(tMat * lattice.inv_lat_column_mat());

    tLat2 = tSVD.matrixV() * tSVD.matrixU().transpose() * tMat;

    return Lattice(tLat2, lattice.tol());
  }
}

namespace CASM {
  namespace xtal {
    /// \brief Construct the subgroup that leaves a lattice unchanged
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, std::vector<SymOp> const &super_grp) {
      std::vector<Index> result;
      invariant_subgroup_indices(lat, super_grp, std::back_inserter(result));
      return result;
    }

    Lattice symmetrize(const Lattice &lattice, const std::vector<SymOp> &enforced_group) {
      Eigen::Matrix3d tLat2(Eigen::Matrix3d::Zero());
      Eigen::Matrix3d frac_mat;
      for(Index ng = 0; ng < enforced_group.size(); ng++) {
        frac_mat = iround(lattice.inv_lat_column_mat() * get_matrix(enforced_group[ng]) * lattice.lat_column_mat()).cast<double>();
        tLat2 += frac_mat.transpose() * lattice.lat_column_mat().transpose() * lattice.lat_column_mat() * frac_mat;
      }
      tLat2 /= double(enforced_group.size());

      // tLat2 has the symmetrized lengths and angles -- it is equal to L.transpose()*L, where L=lat_column_mat()
      // we will find the sqrt of tLat2 and then reorient it so that it matches the original lattice
      Eigen::Matrix3d tMat;

      Eigen::JacobiSVD<Eigen::Matrix3d> tSVD(tMat, Eigen::ComputeFullU | Eigen::ComputeFullV);

      // tMat is sqrt of symlat.transpose()*symlat
      tMat = tSVD.matrixU() * tSVD.singularValues().cwiseSqrt().asDiagonal() * tSVD.matrixV().transpose();

      // if starting lattice were perfect, we would have origlat=rotation*tMat
      tSVD.compute(tMat * lattice.inv_lat_column_mat());

      tLat2 = tSVD.matrixV() * tSVD.matrixU().transpose() * tMat;

      return Lattice(tLat2, lattice.tol());
    }

    Lattice symmetrize(const Lattice &lattice, double point_group_tolerance) {
      auto point_group = make_point_group(lattice, point_group_tolerance);
      return symmetrize(lattice, point_group);
    }

    std::vector<SymOp> make_point_group(Lattice const &_lat) {
      return make_point_group(_lat, _lat.tol());
    }

    std::vector<SymOp> make_point_group(Lattice const &_lat, double tol) {
      std::vector<Eigen::Matrix3i> frac_point_group;
      frac_point_group.reserve(48);

      //Enumerate all possible matrices with elements equal to -1, 0, or 1
      //These represent operations that reorder lattice vectors or replace one
      //or more lattice vectors with a face or body diagonal.

      //For this algorithm to work, lattice needs to be in reduced form.
      Lattice tlat_reduced(_lat.reduced_cell());
      tlat_reduced.set_tol(tol);
      IsPointGroupOp is_equiv(tlat_reduced);
      for(Eigen::Matrix3i const &mat : unimodular_matrices()) {
        if(is_equiv(mat)) {
          frac_point_group.push_back(mat);
        }
      }

      // Find group closure using fractional form (skip for groups sizes that are very unlikely to be unconverged)
      if(frac_point_group.size() != 48 && frac_point_group.size() != 24) {
        Eigen::Matrix3i t_op;
        for(Index i = 0; i < frac_point_group.size(); ++i) {
          t_op = inverse(frac_point_group[i]);
          if(!contains(frac_point_group, t_op)) {
            frac_point_group.push_back(t_op);
          }
          for(Index j = i; j < frac_point_group.size(); ++j) {
            t_op = frac_point_group[i] * frac_point_group[j];
            if(!contains(frac_point_group, t_op)) {
              frac_point_group.push_back(t_op);
            }
          }
        }
      }

      std::vector<SymOp> result;
      result.reserve(frac_point_group.size());
      Eigen::Matrix3d t_cart, t_diff;
      Lattice symlat = symmetrized_with_fractional(tlat_reduced, frac_point_group);
      //std::cout << "Symmetrized lattice is\n " << symlat.lat_column_mat() << "\n";
      for(Eigen::Matrix3i const &frac : frac_point_group) {
        t_cart = symlat.lat_column_mat() * frac.cast<double>() * symlat.inv_lat_column_mat();
        t_diff = t_cart * _lat.lat_column_mat() - _lat.lat_column_mat();
        result.push_back(SymOp::point_operation(t_cart));
      }
      return result;
    }

    Lattice &apply(const SymOp &op, Lattice &lat) {
      return lat = Lattice(get_matrix(op) * lat.lat_column_mat(), lat.tol());
    }

    Lattice copy_apply(const SymOp &op, const Lattice &lat) {
      return Lattice(get_matrix(op) * lat.lat_column_mat(), lat.tol());
    }

  } // namespace xtal
} // namespace CASM

