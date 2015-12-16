#include <algorithm>
#include "casm/clex/Supercell.hh"

namespace CASM {

  template<typename ConfigType>
  ConfigEnumStrain<ConfigType>::ConfigEnumStrain(Supercell &_scel,
                                                 const value_type &_init,
                                                 const std::vector<Index> & linear_partitions,
                                                 const std::vector<double> &magnitudes,
                                                 std::string _mode) :
    ConfigEnum<ConfigType>(_init,_init,-1),
    //m_counter(_init_vec, _final_vec, Eigen::VectorXd::Constant(_init_vec.size(), _inc)),
    //m_proj(proj_mat),
    m_strain_calc(_mode),
    m_perm_begin(_scel.permute_begin()),
    m_perm_end(_scel.permute_end()) {

    std::vector<double> absmags;
    std::transform(magnitudes.cbegin(),magnitudes.cend(),std::back_inserter(absmags),std::abs<double>);
                   
    m_strain_calc.set_symmetrized_sop(_scel.get_primclex().get_prim().point_group());

    //Eigen::MatrixXd axes=m_strain_calc.sop_transf_mat();
    std::vector<Index> mult, subspaces;
    m_proj=m_strain_calc.irreducible_wedge(_scel.get_primclex().get_prim().point_group(),mult,subspaces);
    Eigen::VectorXd init(m_proj.cols()),final(m_proj.cols()),inc(m_proj.cols());
    Index num_sub=subspaces.back()+1;

    Index nc=0;
    for(Index s=0; s<num_sub; s++){
      for(;nc<subspaces.size() && subspaces[nc]==s; nc++){
        if(mult[nc]==1){
          init(nc)=-absmags[s];
          //inc(nc)=2*absmags[s]/double(subspace_partitions[s]);
        }
        else{
          init(nc)=0.0;
          //inc(nc)=absmags[s]/double(subspace_partitions[s]);
        }

        if(linear_partitions[s]<2){
          final(nc)=init(nc);
          inc(nc)=10*TOL;
        }
        else{
          final(nc)=absmags[s]+TOL;
          inc(nc)=absmags[s]/double(linear_partitions[s]);
        }
      }
    }
    m_counter=EigenCounter<Eigen::VectorXd>(init,final,inc);
    _source() = "strain_enumeration";

    std::cout << "Project matrix is \n" << m_proj << "\n";

    m_counter.reset();
    if(!m_counter.valid()){
      std::cout << "COUNTER IS INVALID\n";
      _step() = -1;
    }
    else{
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
    if(++m_counter) {
      _current().set_deformation(m_strain_calc.unrolled_strain_metric_to_F(m_proj*m_counter()));
      std::cout << "Counter is " << m_counter().transpose() << "\n\n";
      //std::cout << "strain vector is \n" << m_proj*m_counter() << "\n\n";
      //std::cout << "DEFORMATION IS\n" << _current().deformation() << "\n\n";
      //is_valid_config = current().is_canonical(_perm_begin(), _perm_end());
      //std::cout << "counter() is: " << m_counter() << ";  is_valid_config: " << is_valid_config
      //<< ";  is_valid_counter: " << m_counter.valid() << "\n";
    }

    if(m_counter.valid())
      _step()++;
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

