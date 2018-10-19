#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/misc/algorithm.hh"
#include "casm/container/Counter.hh"

namespace CASM {
  namespace SymRepTools {

    Eigen::MatrixXd SubWedge::_subwedge_to_trans_mat(std::vector<IrrepWedge> const &_iwedges) {

      if(_iwedges.empty())
        return Eigen::MatrixXd();
      Eigen::MatrixXd result(_iwedges[0].axes.rows(), _iwedges[0].axes.rows());
      Index i = 0;
      for(Index w = 0; w < _iwedges.size(); ++w) {
        result.block(0, i, _iwedges[w].axes.rows(), _iwedges[w].axes.cols()) = _iwedges[w].axes;
        i += _iwedges[w].axes.cols();
      }
      if(i < result.cols())
        result.conservativeResize(Eigen::NoChange, i);
      return result;
    }

    std::vector<IrrepWedge> irreducible_wedges(const SymGroup &head_group, SymGroupRepID id) {
      SymGroupRep const &srep(head_group.master_group().representation(id));
      Index dim = srep.MatrixXd(head_group[0])->cols();

      //Handle for strain symrep
      SymGroupRep::RemoteHandle trep(head_group, id);

      multivector<Eigen::VectorXd>::X<3> sdirs = srep.calc_special_total_directions(head_group);
      std::vector<IrrepWedge> result(sdirs.size());
      double best_proj, tproj;

      for(Index s = 0; s < sdirs.size(); s++) {
        result[s].axes = Eigen::MatrixXd::Zero(dim, sdirs[s].size());
        result[s].axes.col(0) = sdirs[s][0][0];
        result[s].mult.push_back(sdirs[s][0].size());
        for(Index i = 1; i < sdirs[s].size(); i++) {
          Index j_best = 0;
          best_proj = (result[s].axes.transpose() * sdirs[s][i][0]).sum();
          for(Index j = 1; j < sdirs[s][i].size(); j++) {
            tproj = (result[s].axes.transpose() * sdirs[s][i][j]).sum();
            if(tproj > best_proj) {
              best_proj = tproj;
              j_best = j;
            }
          }

          result[s].axes.col(i) = sdirs[s][i][j_best];
        }
      }
      return result;
    }




    //*******************************************************************************************

    std::vector<SubWedge > symrep_subwedges(SymGroup const &head_group, SymGroupRepID id) {
      auto irrep_wedge_compare = [](const IrrepWedge & a, const IrrepWedge & b)->bool {
        return Eigen::almost_equal(a.axes, b.axes);
      };

      auto tot_wedge_compare = [irrep_wedge_compare](const std::vector<IrrepWedge> &a, const std::vector<IrrepWedge> &b)->bool {
        for(auto ita = a.begin(), itb = b.begin(); ita != a.end(); ++ita, ++itb)
          if(!irrep_wedge_compare(*ita, *itb))
            return false;
        //std::cout << "Equal: \n" << SubWedge(a).trans_mat() << ",\n" << SubWedge(b).trans_mat() << "\n\n";
        return true;
      };


      std::vector<SubWedge> result;
      SymGroupRep const &srep(head_group.master_group().representation(id));
      if(!srep[0]->MatrixXd())
        throw std::runtime_error("In symrep_subwedges, SymGroupRep does not describe matrix representation");
      Index dim = srep[0]->MatrixXd()->cols();

      //for(SymOp const &op : head_group) {
      //std::cout << "OP " << op.index() << ":\n" << *(srep[op.index()]->MatrixXd()) << "\n\n";
      //}

      std::vector<IrrepWedge> init_wedges = irreducible_wedges(head_group, id);
      //for(Index i=0; i<init_wedges.size(); ++i){
      //std::cout << "IrrepWedge #" << i+1 << ":\n" << init_wedges[i].axes.transpose() << "\n\n";
      //}


      //Handle for strain symrep
      SymGroupRep::RemoteHandle trep(head_group, id);
      //irrep_wedge_orbits[w] is orbit of wedges[w]
      multivector<IrrepWedge>::X<2> irrep_wedge_orbits;
      irrep_wedge_orbits.reserve(init_wedges.size());
      //max_equiv[w] is irrep_wedge_orbits[w].size()-1
      std::vector<Index> max_equiv;
      max_equiv.reserve(init_wedges.size());

      for(IrrepWedge const &wedge : init_wedges) {
        irrep_wedge_orbits.push_back({wedge});

        //Start getting orbit of wedges[w]
        for(Index p = 0; p < trep.size(); p++) {
          IrrepWedge test_wedge((*(trep[p]->MatrixXd()))*wedge.axes,
                                wedge.mult);

          if(contains(irrep_wedge_orbits.back(), test_wedge, irrep_wedge_compare))
            continue;
          irrep_wedge_orbits.back().push_back(test_wedge);
        }

        max_equiv.push_back(irrep_wedge_orbits.back().size() - 1);
        max_equiv[0] = 0;
      }



      //Counter over combinations of equivalent wedges

      Counter<std::vector<Index> > wcount(std::vector<Index>(init_wedges.size(), 0), max_equiv, std::vector<Index>(init_wedges.size(), 1));
      multivector<IrrepWedge>::X<3> tot_wedge_orbits;
      for(; wcount; ++wcount) {
        std::vector<IrrepWedge> twedge = init_wedges;
        for(Index i = 1; i < init_wedges.size(); i++)
          twedge[i].axes = irrep_wedge_orbits[i][wcount[i]].axes;


        if(contains_if(tot_wedge_orbits,
        [tot_wedge_compare, &twedge](const multivector<IrrepWedge>::X<2> &wedge_orbit)->bool {
        return contains(wedge_orbit,
                        twedge,
                        tot_wedge_compare);
        }))
        continue;

        tot_wedge_orbits.push_back({twedge});
        result.push_back({twedge});
        for(Index p = 0; p < trep.size(); p++) {
          for(Index i = 0; i < twedge.size(); i++)
            twedge[i].axes = (*(trep[p]->MatrixXd())) * result.back().irrep_wedges()[i].axes;
          if(!contains(tot_wedge_orbits.back(), twedge, tot_wedge_compare)) {
            //std::cout << "Adding subwedge!\n";
            tot_wedge_orbits.back().push_back(twedge);
          }
        }
      }

      return result;
    }

  }
}
