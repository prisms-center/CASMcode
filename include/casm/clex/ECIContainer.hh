#ifndef CASM_ECICONTAINER_HH
#define CASM_ECICONTAINER_HH

#include "casm/clex/Correlation.hh"
#include "casm/clex/Clexulator.hh"

namespace CASM {

  /// \brief A sparse container of ECI values and their corresponding orbit indices.
  ///
  /// \ingroup ClexClex
  class ECIContainer {

  public:

    typedef Clexulator::size_type size_type;

    /// \brief Default constructor
    ECIContainer() {};

    /// \brief Construct from range of ECI values and corresponding orbit indices
    ///
    /// \param eci_begin,eci_end Iterators of range of ECI values
    /// \param index_begin Iterator to beginning of range (of size equal to
    ///        [eci_begin, eci_end)) containing the orbit indices of the corresponding
    ///        ECI values.
    ///
    template<typename SparseECIIterator, typename OrbitIndexIterator>
    ECIContainer(SparseECIIterator eci_begin, SparseECIIterator eci_end, OrbitIndexIterator index_begin) {
      auto eci_it = eci_begin;
      auto index_it = index_begin;
      for(; eci_it != eci_end; ++eci_it, ++index_it) {
        m_value.push_back(*eci_it);
        m_index.push_back(*index_it);
      }
    }

    /// \brief Number of eci specified (no guarentee they are all non-zero)
    size_type size() const {
      return m_value.size();
    }

    /// \brief const Access ECI values
    const std::vector<double> &value() const {
      return m_value;
    }

    /// \brief const Access orbit indices of ECI values
    const std::vector<size_type> &index() const {
      return m_index;
    }

  private:

    /// Efective cluster interaction values
    std::vector<double> m_value;

    /// Orbit index for each coefficient in m_value
    std::vector<size_type> m_index;

  };

  /// \brief Evaluate property given an ECIContainer and Correlation
  double operator*(const ECIContainer &_eci, const Correlation &_corr);

  /// \brief Evaluate property given an ECIContainer and pointer to beginning of range of correlation
  double operator*(const ECIContainer &_eci, double const *_corr_begin);

  /// \brief Read eci.out file from specified path (deprecated)
  ECIContainer read_eci_out(const fs::path &filepath);

  /// \brief Read eci.json file from specified path
  ECIContainer read_eci(const fs::path &filepath);

}
#endif
