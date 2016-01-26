#include <algorithm>
#include "casm/clex/Supercell.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

  template<typename ConfigType>
  ConfigEnumStrain<ConfigType>::ConfigEnumStrain(Supercell &_scel,
                                                 const value_type &_init,
                                                 const std::vector<Index> &linear_partitions,
                                                 const std::vector<double> &magnitudes,
                                                 std::string _mode) :
    ConfigEnum<ConfigType>(_init, _init, -1),
    //m_counter(_init_vec, _final_vec, Eigen::VectorXd::Constant(_init_vec.size(), _inc)),
    m_strain_calc(_mode),
    m_proj(m_strain_calc.dim(), m_strain_calc.dim()),
    m_perm_begin(_scel.permute_begin()),
    m_perm_end(_scel.permute_end()),
    m_shape_factor(Eigen::VectorXd::Ones(m_strain_calc.dim())) {

    // Force magnitudes to be positive
    std::vector<double> absmags;
    std::transform(magnitudes.cbegin(), magnitudes.cend(), std::back_inserter(absmags), std::abs<double>);

    m_strain_calc.set_symmetrized_sop(_scel.get_primclex().get_prim().point_group());

    //Eigen::MatrixXd axes=m_strain_calc.sop_transf_mat();
    std::vector<Index> mult;
    std::vector<Eigen::MatrixXd> wedges = m_strain_calc.irreducible_wedges(_scel.get_primclex().get_prim().point_group(), mult);
    Eigen::VectorXd init(m_proj.cols()), final(m_proj.cols()), inc(m_proj.cols());
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
      for( auto m : mult)
        std::cout  << m << "  ";
      std::cout << "\n";

      for(Index w = 0; w < wedges[s].cols(); w++, nc++) {
        m_proj.col(nc) = wedges[s].col(w);
        if(mult[nc] == 1 && linear_partitions[s] > 1) {
          init(nc) = -absmags[s];
          //inc(nc)=2*absmags[s]/double(subspace_partitions[s]);
        }
        else {
          init(nc) = 0.0;
          //inc(nc)=absmags[s]/double(subspace_partitions[s]);
        }

        if(absmags[s] > TOL)
          m_shape_factor[nc] /= absmags[s] * absmags[s];

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
    m_counter = EigenCounter<Eigen::VectorXd>(init, final, inc);
    _source() = "strain_enumeration";

    std::cout << "Project matrix is \n" << m_proj << "\n";
    //m_shape_factor=m_proj.transpose()*m_shape_factor*m_proj;

    m_counter.reset();
    while(m_counter.valid() && double(m_counter().transpose()*m_shape_factor.asDiagonal()*m_counter()) > 1.0 + TOL) {
      ++m_counter;
    }

    if(!m_counter.valid()) {
      std::cout << "COUNTER IS INVALID\n";
      _step() = -1;
    }
    else {
      _step() = 0;
    }
  }
  //*******************************************************************************************
  // **** Mutators ****
  // increment m_current and return a reference to it
  template<typename ConfigType>
  const typename ConfigEnumStrain<ConfigType>::value_type &ConfigEnumStrain<ConfigType>::increment() {
    //bool is_valid_config(false);
    //std::cout << "Incrementing...\n";
    //while(!is_valid_config && ++m_counter) {
    while(++m_counter && double(m_counter().transpose()*m_shape_factor.asDiagonal()*m_counter()) > 1.0 + TOL) {
      //just burning throught the count
    }

    if(m_counter.valid()) {
      _current().set_deformation(m_strain_calc.unrolled_strain_metric_to_F(m_proj * m_counter()));
      std::cout << "Counter is " << m_counter().transpose() << "\n\n";
      //std::cout << "strain vector is \n" << m_proj*m_counter() << "\n\n";
      //std::cout << "DEFORMATION IS\n" << _current().deformation() << "\n\n";
      //is_valid_config = current().is_canonical(_perm_begin(), _perm_end());
      //std::cout << "counter() is: " << m_counter() << ";  is_valid_config: " << is_valid_config
      //<< ";  is_valid_counter: " << m_counter.valid() << "\n";
      _step()++;
    }
    else {
      //std::cout << "REACHED END OF THE LINE!\n";
      _step() = -1;
    }
    //std::cout << "--FINISHED SEARCH " << _step()<< "--\n";
    _current().set_source(source());
    return current();
  }

  //*******************************************************************************************
  // set m_current to correct value at specified step and return a reference to it
  template<typename ConfigType>
  const typename ConfigEnumStrain<ConfigType>::value_type &ConfigEnumStrain<ConfigType>::goto_step(step_type _step) {
    std::cerr << "CRITICAL ERROR: Class ConfigEnumStrain does not implement a goto_step() method. \n"
              << "                You may be using a ConfigEnumIterator in an unsafe way!\n"
              << "                Exiting...\n";
    assert(0);
    exit(1);
    return current();
  };
}

