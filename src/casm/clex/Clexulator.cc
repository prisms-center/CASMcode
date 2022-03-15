#include "casm/clex/Clexulator.hh"

#include <boost/filesystem.hpp>

#include "casm/app/LogRuntimeLibrary.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clexulator/ClexParamPack.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {

/// \brief Construct a Clexulator
///
/// \param name Class name for the Clexulator, typically
///     'X_Clexulator_<basis_set_name>', with X referring to the system of
///     interest (i.e. 'NiAl_Clexulator'), or
///     'X_Clexulator_<basis_set_name>_<index>' for local clexulators, where
///     `<index>` is an equivalent local cluster expansion index.
/// \param dirpath Directory containing the source code and compiled object
/// file. \param nlist, A PrimNeighborList to be updated to include the
/// neighborhood
///        of this Clexulator
/// \param compile_options Compilation options
/// \param so_options Shared library compilation options
///
/// If 'name' is 'X_Clexulator_default', and 'dirpath' is '/path/to':
/// - Looks for '/path/to/X_Clexulator_default.so' and tries to load it.
/// - If not found, looks for 'path/to/X_Clexulator_default.cc' and tries to
//    compile and load it.
/// - If unsuccesful, will throw std::runtime_error.
///
/// The Clexulator has shared ownership of the loaded library,
/// so it is preferrable to duplicate the Clexulator using it's copy constructor
/// rather than construct another using this constructor which will re-load the
/// library.
///
Clexulator::Clexulator(std::string name, fs::path dirpath,
                       PrimNeighborList &nlist, std::string compile_options,
                       std::string so_options)
    : m_name(name) {
  // Construct the RuntimeLibrary that will store the loaded clexulator library
  try {
    m_lib = log_make_shared_runtime_lib(
        (dirpath / name).string(), compile_options, so_options,
        "compile time depends on how many basis functions are included");
  } catch (std::exception &e) {
    log() << "Clexulator construction failed: could not construct runtime "
             "library."
          << std::endl;
    throw;
  }

  // Get the Clexulator factory function
  std::function<clexulator::BaseClexulator *(void)> factory;
  factory =
      m_lib->get_function<clexulator::BaseClexulator *(void)>("make_" + name);

  // Use the factory to construct the clexulator and store it in m_clex
  m_clex.reset(factory());

  // Check nlist has the right weight_matrix
  if (nlist.weight_matrix() != m_clex->weight_matrix()) {
    std::cerr << "Error in Clexulator constructor: weight matrix of neighbor "
                 "list does not match the weight matrix used to print the "
                 "clexulator."
              << std::endl;
    std::cerr << "nlist weight matrix: \n"
              << nlist.weight_matrix() << std::endl;
    std::cerr << "clexulator weight matrix: \n"
              << m_clex->weight_matrix() << std::endl;
    throw std::runtime_error(
        "Error in Clexulator constructor: weight matrix of neighbor list does "
        "not match the weight matrix used to print the clexulator. Try 'casm "
        "bset -uf'.");
  }

  // Expand the given neighbor list as necessary
  nlist.expand(neighborhood().begin(), neighborhood().end());
}

/// \brief Copy constructor
Clexulator::Clexulator(const Clexulator &B) : m_name(B.name()), m_lib(B.m_lib) {
  if (B.m_clex != nullptr) {
    m_clex = B.m_clex->clone();
  }
}

/// \brief Move constructor
Clexulator::Clexulator(Clexulator &&B) { swap(*this, B); }

Clexulator::~Clexulator() {
  // ensure Clexulator is deleted before library
  m_clex.reset();
}

/// \brief Assignment operator
Clexulator &Clexulator::operator=(Clexulator B) {
  swap(*this, B);
  return *this;
}

/// \brief Obtain ClexParamKey for a particular parameter
clexulator::ClexParamKey const &Clexulator::param_key(
    std::string const &_param_name) const {
  return m_clex->param_key(_param_name);
}

/// \brief Clexulator factory function
Clexulator make_clexulator(std::string name, fs::path dirpath,
                           PrimNeighborList &nlist, std::string compile_options,
                           std::string so_options) {
  return Clexulator(name, dirpath, nlist, compile_options, so_options);
}

/// \brief Local Clexulator factory function
std::vector<Clexulator> make_local_clexulator(std::string name,
                                              fs::path dirpath,
                                              PrimNeighborList &nlist,
                                              std::string compile_options,
                                              std::string so_options) {
  std::vector<Clexulator> result;
  Index i = 0;
  fs::path equiv_dir = dirpath / fs::path(std::to_string(i));
  while (fs::exists(equiv_dir)) {
    std::string equiv_name = name + "_" + std::to_string(i);
    if (!fs::exists(equiv_dir / (equiv_name + ".cc"))) {
      break;
    }
    result.push_back(make_clexulator(equiv_name, equiv_dir, nlist,
                                     compile_options, so_options));
    ++i;
    equiv_dir = dirpath / fs::path(std::to_string(i));
  }
  return result;
}

}  // namespace CASM
