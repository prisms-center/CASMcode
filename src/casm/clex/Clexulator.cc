#include "casm/clex/Clexulator.hh"
#include "casm/clex/ClexParamPack.hh"
namespace CASM {

  /// \brief Obtain ClexParamKey for a particular parameter
  ClexParamKey const &Clexulator::param_key(std::string const &_param_name)const {
    return m_clex->param_key(_param_name);
  }

  /// \brief Alter evaluation of parameters specified by @param _param_key, using a custom double -> double function set
  void Clexulator::set_evaluation(ClexParamKey const _param_key, std::vector<std::function<double(ConfigDoF const &) > > const   &_basis_set) {
    m_clex->set_evaluation(_param_key, _basis_set);
  }

  /// \brief Alter evaluation of parameters specified by @param _param_key, using a custom int -> double function set
  void Clexulator::set_evaluation(ClexParamKey const _param_key, std::vector<std::function<double(std::vector<double> const &) > > const &_basis_set) {
    m_clex->set_evaluation(_param_key, _basis_set);
  }

  /// \brief Alter evaluation of parameters specified by @param _param_key, using the string  @param _eval_type,
  /// which can be at least either "READ" (i.e., read from ClexParamPack) or "DEFAULT" (i.e., the Clexulator's default implementation)
  void Clexulator::set_evaluation(ClexParamKey const _param_key, std::string _eval_type) {
    m_clex->set_evaluation(_param_key, _eval_type);
  }

  /// \brief Check evaluation mode of parameters specified by @param _param_key, which can be one of (at least)
  /// "READ" (i.e., read from ClexParamPack), "CUSTOM", or "DEFAULT" (i.e., the Clexulator's default implementation)
  std::string Clexulator::check_evaluation(ClexParamKey const _param_key) const {
    return m_clex->check_evaluation(_param_key);
  }

  namespace Clexulator_impl {
    /// \brief Alter evaluation of parameters specified by @param _param_key, using a custom double -> double function set
    void Base::set_evaluation(ClexParamKey const &_param_key, std::vector<std::function<double(ConfigDoF const &) > > const   &_basis_set) {}

    /// \brief Alter evaluation of parameters specified by @param _param_key, using a custom int -> double function set
    void Base::set_evaluation(ClexParamKey const &_param_key, std::vector<std::function<double(std::vector<double> const &) > > const &_basis_set) {}

    /// \brief Alter evaluation of parameters specified by @param _param_key, using the string  @param _eval_type,
    // which can be at least either "READ" (i.e., read from ClexParamPack) or "DEFAULT" (i.e., the Clexulator's default implementation)
    void Base::set_evaluation(ClexParamKey const &_param_key, std::string const &_eval_type) {}

    /// \brief Check evaluation mode of parameters specified by @param _param_key, which can be one of (at least)
    /// "READ" (i.e., read from ClexParamPack), "CUSTOM", or "DEFAULT" (i.e., the Clexulator's default implementation)
    std::string Base::check_evaluation(ClexParamKey const &_param_key) const {
      return "";
    }

    /// \brief Obtain ClexParamKey for a particular parameter
    ClexParamKey const  &Base::param_key(std::string const &_param_name)const {
      return param_pack().key(_param_name);
    }

  }

}
