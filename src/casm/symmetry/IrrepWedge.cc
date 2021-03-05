#include "casm/symmetry/IrrepWedge.hh"

#include "casm/container/Counter.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/symmetry/IrrepDecompositionImpl.hh"
#include "casm/symmetry/VectorSymCompare_v2.hh"

namespace CASM {

namespace SymRepTools_v2 {

namespace IrrepWedgeImpl {

/// \param _rep The dimension of _rep should match the dimension of the irrep
static IrrepWedge _wedge_from_pseudo_irrep(IrrepInfo const &irrep,
                                           MatrixRep const &_rep,
                                           GroupIndices const &head_group) {
  Eigen::MatrixXd t_axes = irrep.trans_mat.transpose().real();
  Eigen::MatrixXd axes = vector_space_prepare(t_axes, TOL);
  Eigen::VectorXd v = axes.col(0);
  Eigen::VectorXd vbest;

  axes.setZero(irrep.vector_dim(), irrep.irrep_dim());

  axes.col(0) = v;

  for (Index i = 1; i < axes.cols(); ++i) {
    double bestproj = -1;
    for (Index element_index : head_group) {
      v = _rep[element_index] * axes.col(0);
      // std::cout << "v: " << v.transpose() << std::endl;
      bool skip_op = false;
      for (Index j = 0; j < i; ++j) {
        if (almost_equal(v, axes.col(j))) {
          skip_op = true;
          break;
        }
      }
      if (skip_op) continue;

      // std::cout << "bproj: " << bestproj << "  proj: " <<
      // (v.transpose()*axes).transpose() << std::endl;
      if (bestproj < (v.transpose() * axes).sum()) {
        bestproj = (v.transpose() * axes).sum();
        vbest = v;
      }
    }
    axes.col(i) = vbest;
  }
  return IrrepWedge(irrep, axes);
}

}  // namespace IrrepWedgeImpl

IrrepWedge::IrrepWedge(IrrepInfo _irrep_info, Eigen::MatrixXd _axes)
    : irrep_info(std::move(_irrep_info)), axes(std::move(_axes)) {}

/// Construct a "dummy" IrrepWedge with user specified axes
///
/// Note:
/// - Sets the multiplicity for each IrrepWedge axis to 1. This means each
///   direction in the columns of axes are unique and have no orbits
/// - This is a way to directly specify a custom enumeration grid for methods
///   that are written to enumerate on a vector of SubWedge spaces (i.e.
///   ConfigEnumStrain).
IrrepWedge make_dummy_irrep_wedge(Eigen::MatrixXd const &axes) {
  IrrepWedge irrep_wedge(make_dummy_irrep_info(axes), axes);
  irrep_wedge.mult.reserve(axes.cols());
  for (Index i = 0; i < axes.cols(); ++i) {
    irrep_wedge.mult.push_back(1);
  }
  return irrep_wedge;
}

SubWedge::SubWedge(std::vector<IrrepWedge> const &_iwedges)
    : m_iwedges(_iwedges), m_trans_mat(_subwedge_to_trans_mat(m_iwedges)) {}

Eigen::MatrixXd SubWedge::_subwedge_to_trans_mat(
    std::vector<IrrepWedge> const &_iwedges) {
  if (_iwedges.empty()) return Eigen::MatrixXd();
  Eigen::MatrixXd result(_iwedges[0].axes.rows(), _iwedges[0].axes.rows());
  Index i = 0;
  for (Index w = 0; w < _iwedges.size(); ++w) {
    result.block(0, i, _iwedges[w].axes.rows(), _iwedges[w].axes.cols()) =
        _iwedges[w].axes;
    i += _iwedges[w].axes.cols();
  }
  if (i < result.cols()) result.conservativeResize(Eigen::NoChange, i);
  return result;
}

/// Makes a "dummy" SubWedge from a single "dummy" IrrepWedge with given axes
///
/// Note:
/// - Sets the multiplicity for each IrrepWedge axis to 1. This means each
///   direction in the columns of axes are unique and have no orbits
/// - This is a way to directly specify a custom enumeration grid for methods
///   that are written to enumerate on a vector of SubWedge spaces (i.e.
///   ConfigEnumStrain).
SubWedge make_dummy_subwedge(Eigen::MatrixXd const &axes) {
  return SubWedge({make_dummy_irrep_wedge(axes)});
}

/// Make IrrepWedges from an IrrepDecomposition
///
/// \result A vector of IrrepWedge. The IrrepWedge axes have number of cols ==
///     irrep dimension and number of rows equal to the full space dimension.
std::vector<IrrepWedge> make_irrep_wedges(
    IrrepDecomposition const &irrep_decomposition) {
  std::vector<IrrepInfo> const &irreps = irrep_decomposition.irreps;
  MatrixRep const &fullspace_rep = irrep_decomposition.fullspace_rep;
  GroupIndices const &head_group = irrep_decomposition.head_group;

  std::vector<IrrepWedge> wedges;
  wedges.reserve(irreps.size());
  double best_proj, tproj;

  // std::cout << "subspace: \n" << _subspace << std::endl;
  // std::cout << "irrep decomposition size: " << irreps.size() << std::endl;

  for (IrrepInfo const &irrep : irreps) {
    // 1D irreps directions can have positive and negative directions, but we
    // only want to include one. only 1D irreps can have singly-degenerate
    // directions and two singly-degenerate directions indicate the same vector
    // duplicated in positive and negative direction (because they are not
    // equivalent by symmetry) If irrep.directions[0] is singly degenerate
    // (orbits size == 1) then irrepdim is 1 and we only need one direction to
    // define wedge
    wedges.emplace_back(
        irrep, Eigen::MatrixXd::Zero(irrep.vector_dim(), irrep.irrep_dim()));
    // std::cout << "Irrep characters: \n" << irrep.characters << std::endl;
    // std::cout << "Irrep directions: " << irrep.directions.size() <<
    // std::endl;
    if (irrep.directions.empty()) {
      wedges.back() = IrrepWedgeImpl::_wedge_from_pseudo_irrep(
          irrep, fullspace_rep, head_group);
      continue;
    }

    // std::cout << "Irrep direction orbit" << 0 << " : " <<
    // irrep.directions[0].size() << std::endl; std::cout << "Irrep direction: "
    // << irrep.directions[0][0].transpose() << std::endl;
    wedges.back().axes.col(0) = irrep.directions[0][0];
    wedges.back().mult.push_back(irrep.directions[0].size());
    for (Index i = 1; i < irrep.irrep_dim(); i++) {
      // std::cout << "Irrep direction orbit" << i << " : " <<
      // irrep.directions[i].size() << std::endl; std::cout << "Irrep direction:
      // " << irrep.directions[i][0].transpose() << std::endl;
      Index j_best = 0;
      best_proj =
          (wedges.back().axes.transpose() * irrep.directions[i][0]).sum();
      for (Index j = 1; j < irrep.directions[i].size(); j++) {
        tproj = (wedges.back().axes.transpose() * irrep.directions[i][j]).sum();
        if (tproj > best_proj) {
          best_proj = tproj;
          j_best = j;
        }
      }

      wedges.back().axes.col(i) = irrep.directions[i][j_best];
      wedges.back().mult.push_back(irrep.directions[i].size());
    }
    // std::cout << "New irrep wedge: \n" << wedges.back().axes.transpose() <<
    // std::endl;
  }
  return wedges;
}

/// \brief Find full irreducible wedge of a group-represented vector space, as
/// a vector of SubWedges, from an IrrepDecomposition
std::vector<SubWedge> make_symrep_subwedges(
    IrrepDecomposition const &irrep_decomposition) {
  auto irrep_wedge_compare = [](const IrrepWedge &a,
                                const IrrepWedge &b) -> bool {
    return Eigen::almost_equal(a.axes, b.axes);
  };

  auto tot_wedge_compare = [irrep_wedge_compare](
                               const std::vector<IrrepWedge> &a,
                               const std::vector<IrrepWedge> &b) -> bool {
    if (a.size() != b.size()) return false;
    for (auto ita = a.begin(), itb = b.begin(); ita != a.end(); ++ita, ++itb)
      if (!irrep_wedge_compare(*ita, *itb)) return false;
    // std::cout << "Equal: \n" << SubWedge(a).trans_mat() << ",\n" <<
    // SubWedge(b).trans_mat() << "\n\n";
    return true;
  };

  std::vector<IrrepWedge> init_wedges = make_irrep_wedges(irrep_decomposition);
  MatrixRep const &fullspace_rep = irrep_decomposition.fullspace_rep;
  GroupIndices const &head_group = irrep_decomposition.head_group;

  std::vector<SubWedge> result;

  // irrep_wedge_orbits[w] is orbit of wedges[w]
  multivector<IrrepWedge>::X<2> irrep_wedge_orbits;
  irrep_wedge_orbits.reserve(init_wedges.size());
  // max_equiv[w] is irrep_wedge_orbits[w].size()-1
  std::vector<Index> max_equiv;
  max_equiv.reserve(init_wedges.size());
  // std::cout << "irreducible wedges for group of order " << head_group.size()
  // << std::endl;
  Index imax = 0;
  multivector<Index>::X<2> subgroups;
  for (IrrepWedge const &wedge : init_wedges) {
    // std::cout << "Working wedge with axes: \n" << wedge.axes.transpose() <<
    // std::endl;
    irrep_wedge_orbits.push_back({wedge});

    // Start getting orbit of wedges[w]
    subgroups.push_back({});
    for (Index element_index : head_group) {
      IrrepWedge test_wedge{wedge};
      test_wedge.axes = fullspace_rep[element_index] * wedge.axes;
      Index o = 0;
      for (; o < irrep_wedge_orbits.back().size(); ++o) {
        if (irrep_wedge_compare(irrep_wedge_orbits.back()[o], test_wedge)) {
          if (o == 0) {
            subgroups.back().push_back(element_index);
          }
          break;
        }
      }
      if (o < irrep_wedge_orbits.back().size()) continue;
      irrep_wedge_orbits.back().push_back(test_wedge);
    }
    // std::cout << "wedge mult: " << irrep_wedge_orbits.back().size();
    // std::cout << "; subgroups[" << subgroups.size() << "]: " <<
    // subgroups.back() << std::endl; std::cout << "N equiv wedges found: " <<
    // irrep_wedge_orbits.back().size() << std::endl;
    max_equiv.push_back(irrep_wedge_orbits.back().size() - 1);
    if (max_equiv.back() > max_equiv[imax]) imax = max_equiv.size() - 1;
  }
  max_equiv[imax] = 0;

  // Counter over combinations of equivalent wedges
  Counter<std::vector<Index>> wcount(std::vector<Index>(init_wedges.size(), 0),
                                     max_equiv,
                                     std::vector<Index>(init_wedges.size(), 1));

  // std::cout << "max_equiv: " << max_equiv << std::endl;
  // std::cout << "init wcount: " << wcount() << std::endl;
  multivector<IrrepWedge>::X<3> tot_wedge_orbits;
  // std::cout << "Starting slow bit!\n";
  for (; wcount; ++wcount) {
    std::vector<IrrepWedge> twedge = init_wedges;
    for (Index i = 0; i < init_wedges.size(); i++)
      twedge[i].axes = irrep_wedge_orbits[i][wcount[i]].axes;

    if (contains_if(
            tot_wedge_orbits,
            [tot_wedge_compare, &twedge](
                const multivector<IrrepWedge>::X<2> &wedge_orbit) -> bool {
              return contains(wedge_orbit, twedge, tot_wedge_compare);
            }))
      continue;

    tot_wedge_orbits.push_back({twedge});
    result.emplace_back(twedge);
    for (Index p : subgroups[imax]) {
      for (Index i = 0; i < twedge.size(); i++)
        twedge[i].axes =
            fullspace_rep[p] * result.back().irrep_wedges()[i].axes;
      if (!contains(tot_wedge_orbits.back(), twedge, tot_wedge_compare)) {
        tot_wedge_orbits.back().push_back(twedge);
      }
    }
  }
  // std::cout << "Num subwedges: " << result.size() << std::endl;
  return result;
}

}  // namespace SymRepTools_v2

}  // namespace CASM
