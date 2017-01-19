#include "casm/clex/ConfigEnumStrain.hh"
#include <algorithm>
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumStrain_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumStrain>();
  }
}

namespace CASM {

  const std::string ConfigEnumStrain::enumerator_name = "ConfigEnumStrain";

  const std::string ConfigEnumStrain::interface_help =
    "ConfigEnumStrain: \n\n"

    "  ... include help documentation here ... \n\n";

  int ConfigEnumStrain::run(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {
    throw std::runtime_error("EnumInterface<Strain>::run is not implemented");
  }

  ConfigEnumStrain::ConfigEnumStrain(Supercell &_scel,
                                     const Configuration &_init,
                                     const std::vector<Index> &linear_partitions,
                                     const std::vector<double> &magnitudes,
                                     std::string _mode) :
    m_current(_init),
    m_equiv_ind(0),
    m_strain_calc(_mode),
    m_perm_begin(_scel.permute_begin()),
    m_perm_end(_scel.permute_end()),
    m_shape_factor(Eigen::MatrixXd::Identity(m_strain_calc.dim(), m_strain_calc.dim())) {

    auto eigen_compare = [](const Eigen::MatrixXd & a, const Eigen::MatrixXd & b)->bool {
      // required because Eigen::almost_equal takes 3 args (function pointers don't know about default args)
      return Eigen::almost_equal(a, b);
    };

    Index sdim(m_strain_calc.dim());
    // Force magnitudes to be positive
    std::vector<double> absmags;
    std::transform(magnitudes.cbegin(), magnitudes.cend(), std::back_inserter(absmags), std::abs<double>);

    m_strain_calc.set_symmetrized_sop(_scel.get_primclex().get_prim().point_group());

    //Eigen::MatrixXd axes=m_strain_calc.sop_transf_mat();
    std::vector<Index> mult;
    std::vector<Eigen::MatrixXd> wedges = m_strain_calc.irreducible_wedges(_scel.get_primclex().get_prim().point_group(), mult);
    Eigen::VectorXd init(sdim), final(sdim), inc(sdim);
    Index num_sub = wedges.size();

    Index nc = 0;
    for(Index s = 0; s < num_sub; s++) {
      double wedgevol = sqrt((wedges[s].transpose() * wedges[s]).determinant());
      Index N = round(pow(linear_partitions[s], wedges[s].cols()));
      //double density = double(N) / pow(magnitudes[s], wedges[s].cols());
      std::cout << "wedgevol: " << wedgevol << ", N: " << N;// << ", density: " << density;
      N = max(1, (int) ceil(pow(wedgevol * double(N), 1.0 / double(wedges[s].cols())) - TOL));
      std::cout << ", linearN: " << N << "\n";
      std::cout << "mult.size() is: " << mult.size() << "  and mult is ";
      for(auto m : mult)
        std::cout  << m << "  ";
      std::cout << "\n";

      for(Index i = 0; i < wedges[s].cols(); i++, nc++) {
        if(mult[nc] == 1 && linear_partitions[s] > 1) {
          init(nc) = -absmags[s];
          //inc(nc)=2*absmags[s]/double(subspace_partitions[s]);
        }
        else {
          init(nc) = 0.0;
          //inc(nc)=absmags[s]/double(subspace_partitions[s]);
        }

        if(absmags[s] > TOL)
          m_shape_factor(nc, nc) /= absmags[s] * absmags[s];

        if(linear_partitions[s] < 2) {
          final(nc) = init(nc);
          inc(nc) = 10 * TOL;
        }
        else {
          final(nc) = absmags[s] + TOL;
          inc(nc) = absmags[s] / double(N);
        }
      }
    }

    //Handle for strain symrep
    SymGroupRep::RemoteHandle trep(_scel.get_prim().factor_group(), m_strain_calc.symrep_ID());
    //wedge_orbits[w] is orbit of wedges[w]
    multivector<Eigen::MatrixXd>::X<2> wedge_orbits(wedges.size());
    //max_equiv[w] is wedge_orbits[w].size()-1
    std::vector<Index> max_equiv(wedges.size());
    for(Index w = 0; w < wedges.size(); w++) {
      //don't bother with symmetry if we aren't perturbing the current wedge
      if(absmags[w] < TOL) {
        wedge_orbits[w].push_back(0.0 * wedges[w]);
        continue;
      }

      //Start getting orbit of wedges[w]
      for(Index p = 0; p < trep.size(); p++) {
        Eigen::MatrixXd twedge(m_strain_calc.sop_transf_mat() * (*(trep[p]->get_MatrixXd()))*m_strain_calc.sop_transf_mat().transpose()*wedges[w]);
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
      Eigen::MatrixXd ttrans(sdim, sdim);
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
      m_trans_mats.push_back(ttrans);
      for(Index p = 0; p < trep.size(); p++) {
        Eigen::MatrixXd symtrans(m_strain_calc.sop_transf_mat() * (*(trep[p]->get_MatrixXd()))*m_strain_calc.sop_transf_mat().transpose()*ttrans);
        if(!contains(trans_mat_orbits.back(), symtrans, eigen_compare))
          trans_mat_orbits.back().push_back(symtrans);
      }
    }



    m_counter = EigenCounter<Eigen::VectorXd>(init, final, inc);

    std::cout << "Project matrices are \n";
    for(Index i = 0; i < m_trans_mats.size(); i++)
      std::cout << i << " of " << m_trans_mats.size() << "\n" << m_trans_mats[i] << "\n\n";
    m_shape_factor = m_strain_calc.sop_transf_mat().transpose() * m_shape_factor * m_strain_calc.sop_transf_mat();

    m_counter.reset();
    while(m_counter.valid() && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL) {
      ++m_counter;
    }

    reset_properties(m_current);
    this->_initialize(&m_current);

    if(!m_counter.valid()) {
      this->_invalidate();
    }
    _current().set_source(this->source(step()));
  }

  // Implements _increment
  void ConfigEnumStrain::increment() {
    //bool is_valid_config(false);
    //std::cout << "Incrementing...\n";

    while(++m_counter && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL) {
      //just burning throught the count
    }

    // move to next part of wedge if necessary
    if(!m_counter.valid() && m_equiv_ind + 1 < m_trans_mats.size()) {
      m_counter.reset();
      ++m_equiv_ind;
      std::cout << "INCREMENTING m_equiv_ind to " << m_equiv_ind << "\n";
    }

    while(m_counter && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL) {
      //just burning throught the count
      ++m_counter;
    }

    if(m_counter.valid()) {
      _current().set_deformation(m_strain_calc.unrolled_strain_metric_to_F(m_trans_mats[m_equiv_ind] * m_counter()));
      std::cout << "Counter is " << m_counter().transpose() << "\n\n";
      //std::cout << "strain vector is \n" << m_trans_mats[m_equiv_ind]*m_counter() << "\n\n";
      //std::cout << "DEFORMATION IS\n" << _current().deformation() << "\n\n";
      //is_valid_config = current().is_canonical(_perm_begin(), _perm_end());
      //std::cout << "counter() is: " << m_counter() << ";  is_valid_config: " << is_valid_config
      //<< ";  is_valid_counter: " << m_counter.valid() << "\n";
      _increment_step();
    }
    else {
      //std::cout << "REACHED END OF THE LINE!\n";
      _invalidate();
    }
    _current().set_source(this->source(step()));
    //std::cout << "--FINISHED SEARCH " << _step()<< "--\n";
    return;
  }

}

