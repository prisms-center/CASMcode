#include "casm/clexulator/BaseClexulator.hh"

#include "casm/clexulator/ClexParamPack.hh"

namespace CASM {
namespace clexulator {

BaseClexulator::BaseClexulator(size_type _nlist_size, size_type _corr_size,
                               size_type _n_point_corr)
    : m_nlist_size(_nlist_size),
      m_corr_size(_corr_size),
      m_n_point_corr(_n_point_corr),
      m_configdofvalues_ptr(nullptr),
      m_nlist_ptr(nullptr),
      m_occ_ptr(nullptr) {}

BaseClexulator::~BaseClexulator() {}

/// \brief Clone the Clexulator
std::unique_ptr<BaseClexulator> BaseClexulator::clone() const {
  return std::unique_ptr<BaseClexulator>(_clone());
}

/// \brief Obtain ClexParamKey for a particular parameter
ClexParamKey const &BaseClexulator::param_key(
    std::string const &_param_name) const {
  return param_pack().key(_param_name);
}

}  // namespace clexulator
}  // namespace CASM
