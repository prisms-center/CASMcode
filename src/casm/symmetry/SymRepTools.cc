#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/misc/algorithm.hh"
#include "casm/container/Counter.hh"

namespace CASM {
  namespace SymRepTools {

    std::pair<std::vector<Index>, std::vector<Eigen::MatrixXd> > symrep_subwedges(SymGroup const &_group, SymGroupRepID id) {
      auto eigen_compare = [](const Eigen::MatrixXd & a, const Eigen::MatrixXd & b)->bool {
        // required because Eigen::almost_equal takes 3 args (function pointers don't know about default args)
        return Eigen::almost_equal(a, b);
      };

      SymGroupRep const &srep(_group.master_group().representation(id));
      if(!srep[0]->get_MatrixXd())
        throw std::runtime_error("In symrep_subwedges, SymGroupRep does not describe matrix representation");
      Index dim = srep[0]->get_MatrixXd()->cols();

      std::vector<Eigen::MatrixXd> result_wedges;
      //Eigen::MatrixXd axes=m_strain_calc.sop_transf_mat();
      std::vector<Index> result_dims;
      std::vector<Eigen::MatrixXd> wedges = srep.irreducible_wedges(_group, result_dims);

      //Handle for strain symrep
      SymGroupRep::RemoteHandle trep(_group, id);
      //wedge_orbits[w] is orbit of wedges[w]
      multivector<Eigen::MatrixXd>::X<2> wedge_orbits(wedges.size());
      //max_equiv[w] is wedge_orbits[w].size()-1
      std::vector<Index> max_equiv(wedges.size());
      for(Index w = 0; w < wedges.size(); w++) {
        //don't bother with symmetry if we aren't perturbing the current wedge
        //if(absmags[w] < TOL) {
        //wedge_orbits[w].push_back(0.0 * wedges[w]);
        //continue;
        //}

        //Start getting orbit of wedges[w]
        for(Index p = 0; p < trep.size(); p++) {
          Eigen::MatrixXd twedge((*(trep[p]->get_MatrixXd()))*wedges[w]);
          if(contains(wedge_orbits[w], twedge, eigen_compare))
            continue;
          wedge_orbits[w].push_back(twedge);
        }
        //Finish getting orbit of wedge[w]
        max_equiv[w] = wedge_orbits[w].size() - 1;
      }

      //Counter over combinations of equivalent wedges
      Counter<std::vector<Index> > wcount(std::vector<Index>(wedges.size(), 0), max_equiv, std::vector<Index>(wedges.size(), 1));
      multivector<Eigen::MatrixXd>::X<2> trans_mat_orbits;
      for(; wcount; ++wcount) {
        Eigen::MatrixXd ttrans(dim, dim);
        Index l = 0;
        for(Index i = 0; i < wedges.size(); i++) {
          for(Index j = 0; j < wedges[i].cols(); j++)
            ttrans.col(l++) = wedge_orbits[i][wcount[i]].col(j);
        }

        if(contains_if(trans_mat_orbits,
        [&ttrans, eigen_compare](const std::vector<Eigen::MatrixXd> &_orbit)->bool {
        return contains(_orbit, ttrans, eigen_compare);
        })
          ) continue;
        trans_mat_orbits.push_back(std::vector<Eigen::MatrixXd>(1, ttrans));
        result_wedges.push_back(ttrans);
        for(Index p = 0; p < trep.size(); p++) {
          Eigen::MatrixXd symtrans((*(trep[p]->get_MatrixXd()))*ttrans);
          if(!contains(trans_mat_orbits.back(), symtrans, eigen_compare))
            trans_mat_orbits.back().push_back(symtrans);
        }
      }




    }

  }
}
