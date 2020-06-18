#include "casm/clex/Clexulator.hh"
#include "casm/clex/ClexParamPack.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/app/LogRuntimeLibrary.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  /// \brief Construct a Clexulator
  ///
  /// \param name Class name for the Clexulator, typically 'X_Clexulator', with X
  ///             referring to the system of interest (i.e. 'NiAl_Clexulator')
  /// \param dirpath Directory containing the source code and compiled object file.
  /// \param nlist, A PrimNeighborList to be updated to include the neighborhood
  ///        of this Clexulator
  /// \param logging Print messages to inform users that compilation is occuring
  /// \param compile_options Compilation options
  /// \param so_options Shared library compilation options
  ///
  /// If 'name' is 'X_Clexulator', and 'dirpath' is '/path/to':
  /// - Looks for '/path/to/X_Clexulator.so' and tries to load it.
  /// - If not found, looks for 'path/to/X_Clexulator.cc' and tries to compile and load it.
  /// - If unsuccesful, will throw std::runtime_error.
  ///
  /// The Clexulator has shared ownership of the loaded library,
  /// so it is preferrable to duplicate the Clexulator using it's copy constructor rather
  /// than construct another using this constructor which will re-load the library.
  ///
  Clexulator::Clexulator(std::string name,
                         fs::path dirpath,
                         PrimNeighborList &nlist,
                         const Logging &logging,
                         std::string compile_options,
                         std::string so_options) {

    // Construct the RuntimeLibrary that will store the loaded clexulator library
    try {
      m_lib = log_make_shared_runtime_lib(
                (dirpath / name).string(),
                compile_options,
                so_options,
                "compile time depends on how many basis functions are included");
    }
    catch(std::exception &e) {
      logging.log() << "Clexulator construction failed: could not construct runtime library." << std::endl;
      throw;
    }

    // Get the Clexulator factory function
    std::function<Clexulator_impl::Base* (void)> factory;
    factory = m_lib->get_function<Clexulator_impl::Base* (void)>("make_" + name);

    // Use the factory to construct the clexulator and store it in m_clex
    m_clex.reset(factory());

    // Check nlist has the right weight_matrix
    if(nlist.weight_matrix() != m_clex->weight_matrix()) {
      std::cerr << "Error in Clexulator constructor: weight matrix of neighbor "
                "list does not match the weight matrix used to print the "
                "clexulator." << std::endl;
      std::cerr << "nlist weight matrix: \n" << nlist.weight_matrix() << std::endl;
      std::cerr << "clexulator weight matrix: \n" << m_clex->weight_matrix() << std::endl;
      throw std::runtime_error(
        "Error in Clexulator constructor: weight matrix of neighbor list does "
        "not match the weight matrix used to print the clexulator. Try 'casm bset -uf'.");
    }

    // Expand the given neighbor list as necessary
    nlist.expand(neighborhood().begin(), neighborhood().end());

  }


  /// \brief Copy constructor
  Clexulator::Clexulator(const Clexulator &B) :
    m_name(B.name()),
    m_lib(B.m_lib) {

    if(B.m_clex.get() != nullptr) {
      m_clex.reset(B.m_clex->clone().release());
    }
  }

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
