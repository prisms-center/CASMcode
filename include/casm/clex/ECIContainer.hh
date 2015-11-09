#ifndef ECICONTAINER_HH
#define ECICONTAINER_HH

#include "casm/container/Array.hh"
#include "casm/clex/Correlation.hh"


namespace CASM {

  /**
   * Container class to package ScalarECI with their corresponding indexes. This way they can't get mixed up.
   * Holds only a ScalarECI (Array<double>) and a list of indexes (Array<Index>) that you would
   * need to know how to multiply with correlations.
   */

  class ECIContainer {
  public:
    //typedef Clexulator::size_type size_type;
    typedef unsigned int size_type;
    typedef Array<double> ScalarECI;
    ECIContainer() {};
    ECIContainer(const fs::path &eci_fit_file);

    const ScalarECI &eci_list() const {
      return m_eci_list;
    }
    const Array<size_type> &eci_index_list() const {
      return m_eci_index_list;
    }

  private:
    //Efective cluster interactions (typdef of Array<double>)
    ScalarECI m_eci_list;
    //These tell you which correlations m_eci correspond to
    Array<size_type> m_eci_index_list;

  };

  double operator*(const ECIContainer &_eci, const Correlation &_corr);
  
  /// \brief Evaluate property given an ECIContainer and pointer to beginning of range of correlation
  double operator*(const ECIContainer &_eci, double const* _corr_begin);
  
  namespace ECIContainer_impl {
    void populate_eci(const fs::path &filepath, ECIContainer::ScalarECI &mc_eci, Array<ECIContainer::size_type> &mc_eci_index);
  }
}
#endif
