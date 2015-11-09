#ifndef CASM_MonteSampler_HH
#define CASM_MonteSampler_HH

#include "casm/CASM_global_definitions.hh"
#include "casm/monte_carlo/MCData.hh"
#include "casm/clex/CompositionConverter.hh"

namespace CASM {
  
  class MonteCarlo;
  
  /// \brief An abstract base class for sampling and storing data observations
  ///
  /// - Derived classes "know" how to sample a particular property via implementation
  ///   of 'virtual void sample(const MonteCarlo &mc) = 0;'
  /// - Optionally may require and check for convergence to some level of precision 
  ///   given a particular confidence level
  ///
  class MonteSampler {
  
  public:
    
    typedef MCData::size_type size_type;
    
    /// \brief Construct sampler that does not need to converge
    MonteSampler(const std::string &print_name, 
                 double data_confidence, 
                 size_type data_initsize);

    /// \brief Construct sampler that must converge
    MonteSampler(const std::string &print_name, 
                 double data_prec, 
                 double data_confidence, 
                 size_type data_initsize);
    
    virtual ~MonteSampler() {}
    
    
    virtual void sample(const MonteCarlo &mc) = 0;
    
    /// \brief Clear all data observations
    void clear() {
      m_data.clear();
    }
    
    /// \brief Returns pair(true, equil_steps) if equilibration has occured to required precision
    ///
    /// - If must_converge() == false, always returns false
    ///
    std::pair<bool, size_type> is_equilibrated() const {
      
      if(!must_converge()) {
        return std::make_pair(false, m_data.size());
      }
      
      if(!m_convergence_uptodate) {
        m_equilibration = MCDataEquilibration(m_data.observations(), m_prec);
        m_equilibration_uptodate = true;
      }
      
      return std::make_pair(m_equilibration.is_equilibrated(), m_equilibration.equilibration_samples());
    }
    
    /// \brief Returns true if convergence criteria must be met for Monte Carlo calculation to be complete
    bool must_converge() const {
      return m_must_converge;
    }
    
    /// \brief Returns true if convergence criteria have been met for data sampled in range [equil_samples, end)
    ///
    /// - If must_converge() == false, always returns false
    bool is_converged(size_type equil_samples) const {
      
      if(!must_converge()) {
        return false;
      }
      
      // if not calculated, or calculated for a different range of data, re-calculate
      if(!m_convergence_uptodate || m_convergence_start_sample != equil_samples) {
        _check_convergence(equil_samples);
      }
      
      return m_convergence.is_converged(m_prec);
    }
    
    /// \brief Returns <X> for data sampled in range [equil_samples, end)
    double mean(size_type equil_samples) const {
      
      // if not calculated, or calculated for a different range of data, re-calculate
      if(!m_convergence_uptodate || m_convergence_start_sample != equil_samples) {
        _check_convergence(equil_samples);
      }
      
      return m_convergence.mean();
    }
    
    /// \brief Returns <X*X> for data sampled in range [equil_samples, end)
    double squared_norm(size_type equil_samples) const {
      
      // if not calculated, or calculated for a different range of data, re-calculate
      if(!m_convergence_uptodate || m_convergence_start_sample != equil_samples) {
        _check_convergence(equil_samples);
      }
      
      return m_convergence.squared_norm();
    }
    
    /// \brief Returns calculated precision on the mean for data sampled in range [equil_samples, end)
    double calculated_precision(size_type equil_samples) const {
      
      // if not calculated, or calculated for a different range of data, re-calculate
      if(!m_convergence_uptodate || m_convergence_start_sample != equil_samples) {
        _check_convergence(equil_samples);
      }
      
      return m_convergence.calculated_precision();
    }
    
    /// \brief const Access the raw data observation container
    const MCData &data() const {
      return m_data;
    }
    
    /// \brief Property name for printing
    std::string name() const {
      return m_name;
    }
    
    /// \brief Clone this object
    virtual std::unique_ptr<MonteSampler> clone() const = 0;
  
  
  protected:
    
    /// \brief Access the raw data observation container
    MCData &data() {
      m_convergence_uptodate = false;
      m_equilibration_uptodate = false;
      return m_data;
    }
    
  private:
    
    void _check_convergence(size_type equil_samples) const {
      
      m_convergence_start_sample = equil_samples;
      m_convergence = MCDataConvergence(m_data.observations().segment(equil_samples, m_data.observations().size() - equil_samples), m_conf);
      m_convergence_uptodate = true;
    }
    
    /// \brief If false, data observations are recorded but shouldn't be used to determine if Monte Carlo run has converged
    bool m_must_converge;
    
    /// \brief Requested precision, if must_converge
    double m_prec;
    
    /// \brief Requested confidence, if must_converge
    double m_conf;

    /// \brief An array of observations, one for every so many steps or passes
    MCData m_data;

    /// \brief Property name for printing
    std::string m_name;
    
    // enable storing equilibration and convergence info
    
    mutable bool m_equilibration_uptodate = false;
    mutable MCDataEquilibration m_equilibration;
    mutable bool m_convergence_uptodate = false;
    mutable size_type m_convergence_start_sample = 0;
    mutable MCDataConvergence m_convergence;
  };
  
  
  /// \brief Sampler for a scalar property
  ///
  class ScalarMonteSampler : public MonteSampler {
    
    public:
    
    /// \brief Construct sampler that does not need to converge
    ScalarMonteSampler(std::string _property_name, 
                       std::string print_name, 
                       double data_confidence, 
                       size_type data_initsize);
    
    /// \brief Construct sampler that must converge
    ScalarMonteSampler(std::string _property_name, 
                       std::string print_name, 
                       double data_prec, 
                       double data_confidence, 
                       size_type data_initsize);
    
    
    /// \brief Sample data from a MonteCarlo calculation
    void sample(const MonteCarlo &mc);
    
    /// \brief Clone this object
    std::unique_ptr<MonteSampler> clone() const;
    
    private:
    
    std::string m_property_name;
    
  };
  
  /// \brief Sampler for individual elements of a vector property
  ///
  class VectorMonteSampler : public MonteSampler {
    
    public:
    
    typedef Index size_type;

    /// \brief Construct sampler that does not need to converge
    VectorMonteSampler(std::string _property_name,
                       size_type _index, 
                       std::string print_name, 
                       double data_confidence, 
                       size_type data_initsize);
      
    /// \brief Construct sampler that must converge
    VectorMonteSampler(std::string _property_name,
                       size_type _index, 
                       std::string print_name, 
                       double data_prec, 
                       double data_confidence, 
                       size_type data_initsize);
    
    
    /// \brief Sample data from a MonteCarlo calculation
    void sample(const MonteCarlo &mc);
    
    /// \brief Clone this object
    std::unique_ptr<MonteSampler> clone() const;
    
    private:
    
    std::string m_property_name;
    
    size_type m_index;
    
  };
  
  /// \brief Sampler for parametric composition
  ///
  /// - Requires 'comp_n' property to exist in the sampled MonteCarlo object
  ///
  class CompMonteSampler : public MonteSampler {
    
    public:
    
    typedef Index size_type;

    /// \brief Construct sampler that does not need to converge
    CompMonteSampler(size_type _index,
                     const CompositionConverter& _comp_converter, 
                     std::string print_name, 
                     double data_confidence, 
                     size_type data_initsize);
      
    /// \brief Construct sampler that must converge
    CompMonteSampler(size_type _index, 
                     const CompositionConverter& _comp_converter, 
                     std::string print_name, 
                     double data_prec, 
                     double data_confidence, 
                     size_type data_initsize);
    
    
    /// \brief Sample data from a MonteCarlo calculation
    void sample(const MonteCarlo &mc);
    
    /// \brief Clone this object
    std::unique_ptr<MonteSampler> clone() const;
    
    private:
    
    size_type m_index;
    
    CompositionConverter m_comp_converter; 
                     
  };
  
  
  /// \brief Sampler for site fraction
  ///
  /// - Requires 'comp_n' property to exist in the sampled MonteCarlo object
  ///
  class SiteFracMonteSampler : public MonteSampler {
    
    public:
    
    typedef Index size_type;

    /// \brief Construct sampler that does not need to converge
    SiteFracMonteSampler(size_type _index,
                         size_type _basis_size, 
                         std::string print_name, 
                         double data_confidence, 
                         size_type data_initsize);
      
    /// \brief Construct sampler that must converge
    SiteFracMonteSampler(size_type _index, 
                         size_type _basis_size, 
                         std::string print_name, 
                         double data_prec, 
                         double data_confidence, 
                         size_type data_initsize);
    
    
    /// \brief Sample data from a MonteCarlo calculation
    void sample(const MonteCarlo &mc);
    
    /// \brief Clone this object
    std::unique_ptr<MonteSampler> clone() const;
    
    private:
    
    size_type m_index;
    
    size_type m_basis_size;
                     
  };
  
  
  /// \brief Sampler for atom fraction
  ///
  /// - Requires 'comp_n' property to exist in the sampled MonteCarlo object
  ///
  class AtomFracMonteSampler : public MonteSampler {
    
    public:
    
    typedef Index size_type;

    /// \brief Construct sampler that does not need to converge
    AtomFracMonteSampler(size_type _index,
                         size_type _vacancy_index, 
                         std::string print_name, 
                         double data_confidence, 
                         size_type data_initsize);
      
    /// \brief Construct sampler that must converge
    AtomFracMonteSampler(size_type _index, 
                         size_type _vacancy_index, 
                         std::string print_name, 
                         double data_prec, 
                         double data_confidence, 
                         size_type data_initsize);
    
    
    /// \brief Sample data from a MonteCarlo calculation
    void sample(const MonteCarlo &mc);
    
    /// \brief Clone this object
    std::unique_ptr<MonteSampler> clone() const;
    
    private:
    
    size_type m_index;
    
    size_type m_vacancy_index;
                     
  };

}

#endif


