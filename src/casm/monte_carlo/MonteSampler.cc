#include "casm/monte_carlo/MonteSampler.hh"

#include "casm/monte_carlo/MonteCarlo.hh"

namespace CASM {
  
  // ---- MonteSampler Definitions ---------------------------------
  
  /// \brief Construct sampler that does not need to converge
  MonteSampler::MonteSampler(const std::string &print_name, 
                             double data_confidence, 
                             size_type data_initsize):
    m_must_converge(false),
    m_conf(data_confidence),
    m_data(data_initsize),
    m_name(print_name) {
  
  }

  /// \brief Construct sampler that must converge
  MonteSampler::MonteSampler(const std::string &print_name, 
                             double data_prec, 
                             double data_confidence, 
                             size_type data_initsize):
    m_must_converge(true),
    m_prec(data_prec),
    m_conf(data_confidence),
    m_data(data_initsize),
    m_name(print_name) {
    
  }
  
  
  // ---- ScalarMonteSampler Definitions ---------------------------------
  
  /// \brief Construct sampler that does not need to converge
  ///
  /// \param _property_name Name of scalar property to sample, ex: "formation_energy"
  /// \param print_name Name to be printed, ex: "formation_energy"
  /// \param data_initsize For constructing MCData object
  ///
  ScalarMonteSampler::ScalarMonteSampler(std::string _property_name, 
                     std::string print_name, 
                     double data_confidence, 
                     size_type data_initsize) :
    MonteSampler(print_name, data_confidence, data_initsize),
    m_property_name(_property_name) {}
  
  /// \brief Construct sampler that must converge
  ///
  /// \param _property_name Name of scalar property to sample, ex: "formation_energy"
  /// \param print_name Name to be printed, ex: "formation_energy"
  /// \param data_prec Required precision level
  /// \param data_confidence Required confidence level
  /// \param data_initsize For constructing MCData object
  ///
  ScalarMonteSampler::ScalarMonteSampler(std::string _property_name, 
                     std::string print_name, 
                     double data_prec, 
                     double data_confidence, 
                     size_type data_initsize) :
    MonteSampler(print_name, data_prec, data_confidence, data_initsize),
    m_property_name(_property_name) {}
  
  
  /// \brief Sample data from a MonteCarlo calculation
  void ScalarMonteSampler::sample(const MonteCarlo &mc) {
    data().push_back(mc.scalar_property(m_property_name));
  }
  
  /// \brief Clone this object
  std::unique_ptr<MonteSampler> ScalarMonteSampler::clone() const {
    return std::unique_ptr<MonteSampler>(new ScalarMonteSampler(*this));
  }
  
  
  
  // ---- VectorMonteSampler Definitions ---------------------------------
  
  /// \brief Construct sampler that does not need to converge
  ///
  /// \param _property_name Name of vector property to sample, ex: "corr"
  /// \param _index Index of individual element of vector property to sample, ex: 10
  /// \param print_name Name to be printed, ex: "corr(10)"
  /// \param data_initsize For constructing MCData object
  ///
  VectorMonteSampler::VectorMonteSampler(std::string _property_name,
                                         size_type _index, 
                                         std::string print_name, 
                                         double data_confidence, 
                                         size_type data_initsize) :
    MonteSampler(print_name, data_confidence, data_initsize),
    m_property_name(_property_name),
    m_index(_index) {}
    
  /// \brief Construct sampler that must converge
  ///
  /// \param _property_name Name of vector property to sample, ex: "corr"
  /// \param _index Index of individual element of vector property to sample, ex: 10
  /// \param print_name Name to be printed, ex: "corr(10)"
  /// \param data_prec Required precision level
  /// \param data_confidence Required confidence level
  /// \param data_initsize For constructing MCData object
  ///
  VectorMonteSampler::VectorMonteSampler(std::string _property_name,
                                         size_type _index, 
                                         std::string print_name, 
                                         double data_prec, 
                                         double data_confidence, 
                                         size_type data_initsize) :
    MonteSampler(print_name, data_prec, data_confidence, data_initsize),
    m_property_name(_property_name),
    m_index(_index) {}
  
  
  /// \brief Sample data from a MonteCarlo calculation
  void VectorMonteSampler::sample(const MonteCarlo &mc) {
    data().push_back(mc.vector_property(m_property_name)(m_index));
  }
  
  /// \brief Clone this object
  std::unique_ptr<MonteSampler> VectorMonteSampler::clone() const {
    return std::unique_ptr<MonteSampler>(new VectorMonteSampler(*this));
  }
    
  
  
  // ---- CompMonteSampler Definitions ---------------------------------
  
  /// \brief Construct sampler that does not need to converge
  ///
  /// \param _index param composition index to sample, ex: 0 -> 'a', 1 -> 'b'
  /// \param _comp_converter CompositionConverter
  /// \param print_name Name to be printed, ex: "comp(a)"
  /// \param data_initsize For constructing MCData object
  ///
  CompMonteSampler::CompMonteSampler(size_type _index,
                                     const CompositionConverter& _comp_converter, 
                                     std::string print_name, 
                                     double data_confidence, 
                                     size_type data_initsize) :
    MonteSampler(print_name, data_confidence, data_initsize),
    m_index(_index),
    m_comp_converter(_comp_converter) {}
    
  /// \brief Construct sampler that must converge
  ///
  /// \param _index param composition index to sample, ex: 0 -> 'a', 1 -> 'b'
  /// \param _comp_converter CompositionConverter
  /// \param print_name Name to be printed, ex: "comp(a)"
  /// \param data_prec Required precision level
  /// \param data_confidence Required confidence level
  /// \param data_initsize For constructing MCData object
  ///
  CompMonteSampler::CompMonteSampler(size_type _index, 
                                     const CompositionConverter& _comp_converter, 
                                     std::string print_name, 
                                     double data_prec, 
                                     double data_confidence, 
                                     size_type data_initsize) :
    MonteSampler(print_name, data_prec, data_confidence, data_initsize),
    m_index(_index),
    m_comp_converter(_comp_converter) {}
  
  
  /// \brief Sample data from a MonteCarlo calculation
  void CompMonteSampler::sample(const MonteCarlo &mc) {
    
    auto comp_n = mc.vector_property("comp_n");
    auto comp = m_comp_converter.param_composition(comp_n);
    
    data().push_back(comp(m_index));
    
  }
  
  /// \brief Clone this object
  std::unique_ptr<MonteSampler> CompMonteSampler::clone() const {
    return std::unique_ptr<MonteSampler>(new CompMonteSampler(*this));
  }
  
  
  // ---- SiteFracMonteSampler Definitions ---------------------------------
  
  /// \brief Construct sampler that does not need to converge
  ///
  /// \param _index param composition index to sample, ex: 0 -> 'a', 1 -> 'b'
  /// \param _basis_size number of sites per primitive cell
  /// \param print_name Name to be printed, ex: "site_frac(Mg)"
  /// \param data_initsize For constructing MCData object
  ///
  SiteFracMonteSampler::SiteFracMonteSampler(size_type _index,
                                             size_type _basis_size, 
                                             std::string print_name, 
                                             double data_confidence, 
                                             size_type data_initsize) :
    MonteSampler(print_name, data_confidence, data_initsize),
    m_index(_index),
    m_basis_size(_basis_size) {}
    
  /// \brief Construct sampler that must converge
  ///
  /// \param _index param composition index to sample, ex: 0 -> 'a', 1 -> 'b'
  /// \param _basis_size number of sites per primitive cell
  /// \param print_name Name to be printed, ex: "site_frac(Mg)"
  /// \param data_prec Required precision level
  /// \param data_confidence Required confidence level
  /// \param data_initsize For constructing MCData object
  ///
  SiteFracMonteSampler::SiteFracMonteSampler(size_type _index, 
                                             size_type _basis_size, 
                                             std::string print_name, 
                                             double data_prec, 
                                             double data_confidence, 
                                             size_type data_initsize) :
    MonteSampler(print_name, data_prec, data_confidence, data_initsize),
    m_index(_index),
    m_basis_size(_basis_size) {}
  
  
  /// \brief Sample data from a MonteCarlo calculation
  void SiteFracMonteSampler::sample(const MonteCarlo &mc) {
    data().push_back(mc.vector_property("comp_n")(m_index)/m_basis_size);
  }
  
  /// \brief Clone this object
  std::unique_ptr<MonteSampler> SiteFracMonteSampler::clone() const {
    return std::unique_ptr<MonteSampler>(new SiteFracMonteSampler(*this));
  }
    
  
  // ---- AtomFracMonteSampler Definitions ---------------------------------
  
  /// \brief Construct sampler that does not need to converge
  ///
  /// \param _index param composition index to sample, ex: 0 -> 'a', 1 -> 'b'
  /// \param _vacancy_index index of vacancies in comp_n
  /// \param print_name Name to be printed, ex: "atom_frac(Zr)"
  /// \param data_initsize For constructing MCData object
  ///
  AtomFracMonteSampler::AtomFracMonteSampler(size_type _index,
                                             size_type _vacancy_index, 
                                             std::string print_name, 
                                             double data_confidence, 
                                             size_type data_initsize) :
    MonteSampler(print_name, data_confidence, data_initsize),
    m_index(_index),
    m_vacancy_index(_vacancy_index) {}
    
  /// \brief Construct sampler that must converge
  ///
  /// \param _index param composition index to sample, ex: 0 -> 'a', 1 -> 'b'
  /// \param _vacancy_index index of vacancies in comp_n
  /// \param print_name Name to be printed, ex: "atom_frac(Zr)"
  /// \param data_prec Required precision level
  /// \param data_confidence Required confidence level
  /// \param data_initsize For constructing MCData object
  ///
  AtomFracMonteSampler::AtomFracMonteSampler(size_type _index, 
                                             size_type _vacancy_index, 
                                             std::string print_name, 
                                             double data_prec, 
                                             double data_confidence, 
                                             size_type data_initsize) :
    MonteSampler(print_name, data_prec, data_confidence, data_initsize),
    m_index(_index),
    m_vacancy_index(_vacancy_index) {}
  
  
  /// \brief Sample data from a MonteCarlo calculation
  void AtomFracMonteSampler::sample(const MonteCarlo &mc) {
    
    double atom_sum = 0.0;
    for(size_type i=0; i<mc.vector_property("comp_n").size(); ++i) {
      if(i != m_vacancy_index) {
        atom_sum += mc.vector_property("comp_n")(i);
      }
    }
    
    data().push_back(mc.vector_property("comp_n")(m_index)/atom_sum);
  }
  
  /// \brief Clone this object
  std::unique_ptr<MonteSampler> AtomFracMonteSampler::clone() const {
    return std::unique_ptr<MonteSampler>(new AtomFracMonteSampler(*this));
  }
  

}

