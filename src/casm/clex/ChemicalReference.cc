#include "casm/casm_io/container/stream_io.hh"
#include "casm/clex/ChemicalReference_impl.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace {
typedef std::vector<ChemicalReferenceState>::iterator RefStateIterator;
}
template ChemicalReference::ChemicalReference(const BasicStructure &,
                                              RefStateIterator,
                                              RefStateIterator, double);
template void ChemicalReference::set_config<RefStateIterator>(
    const std::string &, RefStateIterator, RefStateIterator, double);
template void ChemicalReference::set_supercell<RefStateIterator>(
    const std::string &, RefStateIterator, RefStateIterator, double);
template void ChemicalReference::set_global<RefStateIterator>(RefStateIterator,
                                                              RefStateIterator,
                                                              double);

/// \brief Construct using the results of n(config), and e(config)
///
/// The element of n(config) corresponding to vacancies is ignored, the
/// rest are mapped using 'species_num[molecule_name[i]] = n(i)'
ChemicalReferenceState::ChemicalReferenceState(
    const Configuration &config,
    std::function<Eigen::VectorXd(Configuration)> n,
    std::function<double(Configuration)> e) {
  auto names = xtal::struc_molecule_name(config.prim());
  auto vec = n(config);

  if (vec.size() != names.size()) {
    std::cerr << "Error constructing ChemicalReferenceState: Number of species "
                 "!= vec.size()\n";
    std::cerr << "  struc_molecule_name: " << names << "\n";
    std::cerr << "  n(config): " << vec << std::endl;
    throw std::runtime_error(
        "Error constructing ChemicalReferenceState: Number of species != "
        "vec.size()");
  }

  for (int i = 0; i < names.size(); ++i) {
    if (!xtal::is_vacancy(names[i])) {
      species_num[names[i]] = vec(i);
    }
  }
  energy_per_species = e(config);
}

// --- ChemicalReference implementations -----------

const std::string ChemicalReference::Name = "chem_ref";

const std::string ChemicalReference::Desc =
    "Returns a reference energy as interpolated via a composition-energy "
    "hyperplane.";

namespace {

Eigen::MatrixXd _species_frac_matrix(
    const PrimClex &primclex, const std::vector<std::string> &ref_config) {
  Eigen::MatrixXd _N(primclex.composition_axes().components().size(),
                     ref_config.size());
  for (int i = 0; i < ref_config.size(); ++i) {
    _N.col(i) = species_frac(*primclex.db<Configuration>().find(ref_config[i]));
  }
  return _N;
}

Eigen::MatrixXd _species_frac_space(const Eigen::MatrixXd &_N) {
  // The input space
  Eigen::MatrixXd N(_N.rows(), _N.cols() - 1);
  for (int i = 0; i < N.cols(); i++) {
    N.col(i) = _N.col(i + 1) - _N.col(0);
  }
  return N;
}

int _rank(const Eigen::MatrixXd &N, double lin_alg_tol) {
  auto Qr = N.transpose().fullPivHouseholderQr();
  Qr.setThreshold(lin_alg_tol);
  return Qr.rank();
}

// Eliminate ref config until the rank of the defined species_frac space
// is the same as the number of ref config.
void _prune_ref_config(const PrimClex &primclex,
                       std::vector<std::string> &ref_config,
                       double lin_alg_tol) {
  // Each column of _N is the species_frac of the corresponding ref config
  Eigen::MatrixXd _N = _species_frac_matrix(primclex, ref_config);

  // Contains vectors spanning the space defined by the ref config
  // N has _N.cols() - 1 columns, with N.col(i) being _N.col(i) - _N.col(0)
  Eigen::MatrixXd N = _species_frac_space(_N);

  while (_rank(N, lin_alg_tol) != N.cols()) {
    // remove ref_config whose summed distance in species_frac space to all
    // others is the minimum
    double min_dist = std::numeric_limits<double>::max();
    double min_dist_ref = -1;
    for (int i = 0; i < ref_config.size(); ++i) {
      double tot_dist = 0.0;
      for (int j = 0; j < ref_config.size(); ++j) {
        if (i != j) {
          tot_dist += (_N.col(i) - _N.col(j)).norm();
        }
      }
      if (tot_dist < min_dist) {
        min_dist_ref = i;
        min_dist = tot_dist;
      }
    }

    ref_config.erase(ref_config.begin() + min_dist_ref);
    _N = _species_frac_matrix(primclex, ref_config);
    N = _species_frac_space(_N);
  }
}

}  // namespace

/// \brief Constructor
///
/// \param _global_ref An Eigen::VectorXd giving the intercepts of the
///                    hyperplane used for the global reference
/// \param _supercell_ref A map of scelname to Eigen::VectorXd specializing the
///                       the reference value by Supercell
/// \param _config_ref A map of configname to Eigen::VectorXd specializing the
///                       the reference value by Configuration
///
/// A hyperplane reference, R, maps vector species_frac, x, to output
/// energy_per_species, y:
/// - \code y = R.dot(x) \endcode
///
/// The global reference, '_global_ref', is required, but may be specialized
/// to give a different R for a particular Supercell or Configuration via
/// optional '_supercell_ref' and '_config_ref'.
///
ChemicalReference::ChemicalReference(const BasicStructure &prim,
                                     const Eigen::VectorXd &_global_ref,
                                     SpecializedRef _supercell_ref,
                                     SpecializedRef _config_ref)
    : HyperPlaneReferenceBase(Name, Desc, _global_ref, ConfigIO::SpeciesFrac(),
                              _supercell_ref, _config_ref),
      m_prim(&prim) {}

/// \brief Clone
std::unique_ptr<ChemicalReference> ChemicalReference::clone() const {
  return notstd::make_unique<ChemicalReference>(*this->_clone());
}

/// \brief Get primitive BasicStructure
const BasicStructure &ChemicalReference::prim() const { return *m_prim; }

// --- Global reference ---

/// \brief const Access the global reference
///
const Eigen::VectorXd &ChemicalReference::global() const {
  return HyperPlaneReferenceBase::global();
}

/// \brief Set global hyperplane reference
///
/// - Erases associated RefStateVec
void ChemicalReference::set_global(const Eigen::VectorXd &ref) {
  _global() = ref;
  m_global_ref_vec.clear();
}

/// \brief const Access a map of configname to RefStateVec for Supercell
///        specialized references
///
/// - There may not be global reference states (maybe only the hyperplane is
///   known), in which case this is 'empty' / has 'size() == 0'
const ChemicalReference::RefStateVec &ChemicalReference::global_ref_states()
    const {
  return m_global_ref_vec;
}

// --- Supercell specialized references ---

/// \brief const Access a map of scelname to reference for Supercell
///        specialized references
///
const std::map<std::string, Eigen::VectorXd> &ChemicalReference::supercell()
    const {
  return HyperPlaneReferenceBase::supercell();
}

/// \brief Set hyperplane reference specialized for a Supercell
///
/// - Erases associated RefStateVec
void ChemicalReference::set_supercell(const std::string &scelname,
                                      const Eigen::VectorXd &ref) {
  _supercell()[scelname] = ref;
  m_supercell_ref_map.erase(scelname);
}

/// \brief Erase hyperplane reference specialized for a Supercell
ChemicalReference::size_type ChemicalReference::erase_supercell(
    const std::string &scelname) {
  auto result = _supercell().erase(scelname);
  m_supercell_ref_map.erase(scelname);
  return result;
}

/// \brief const Access a map of configname to RefStateVec for Supercell
///        specialized references
///
/// - A configuration with a specialized reference need not have an associated
///   RefStateVec
const ChemicalReference::RefStateMap &ChemicalReference::supercell_ref_states()
    const {
  return m_supercell_ref_map;
}

// --- Configuration specialized references ---

/// \brief const Access a map of configname to reference for Configuration
///        specialized references
///
const std::map<std::string, Eigen::VectorXd> &ChemicalReference::config()
    const {
  return HyperPlaneReferenceBase::config();
}

/// \brief Set hyperplane reference specialized for a Configuration
///
/// - Erases associated RefStateVec
void ChemicalReference::set_config(const std::string &configname,
                                   const Eigen::VectorXd &ref) {
  _config()[configname] = ref;
  m_config_ref_map.erase(configname);
}

/// \brief Erase hyperplane reference specialized for a Configuration
ChemicalReference::size_type ChemicalReference::erase_config(
    const std::string &configname) {
  auto result = _config().erase(configname);
  m_config_ref_map.erase(configname);
  return result;
}

/// \brief const Access a map of configname to RefStateVec for Configuration
///        specialized references
///
/// - A configuration with a specialized reference need not have an associated
///   RefStateVec
const ChemicalReference::RefStateMap &ChemicalReference::config_ref_states()
    const {
  return m_config_ref_map;
}

/// \brief Access the global reference
///
Eigen::VectorXd &ChemicalReference::_global() {
  return HyperPlaneReferenceBase::global();
}

/// \brief const Access a map of scelname to reference for Supercell
///        specialized references
///
std::map<std::string, Eigen::VectorXd> &ChemicalReference::_supercell() {
  return HyperPlaneReferenceBase::supercell();
}

/// \brief const Access a map of configname to reference for Configuration
///        specialized references
///
std::map<std::string, Eigen::VectorXd> &ChemicalReference::_config() {
  return HyperPlaneReferenceBase::config();
}

// --- For use in hyperplane() ------

/// \brief Clone
ChemicalReference *ChemicalReference::_clone() const {
  return new ChemicalReference(*this);
}

/// \brief Convert a set of ChemicalReferenceState to a hyperplane, including
/// checks
Eigen::VectorXd ChemicalReference::_calc_hyperplane(
    const BasicStructure &prim, const std::vector<std::string> &struc_mol_name,
    Eigen::MatrixXd _N, Eigen::VectorXd E, double tol) {
  // --- convert input compositions to atom_frac

  Index Va_index = find_index_if(struc_mol_name, [=](const std::string &str) {
    return xtal::is_vacancy(str);
  });
  bool has_Va = (Va_index != struc_mol_name.size());

  if (has_Va) {
    _N.row(Va_index) = Eigen::VectorXd::Zero(_N.cols());
  }
  for (int i = 0; i < _N.cols(); i++) {
    _N.col(i) /= _N.col(i).sum();
  }

  // --- check that the input space is full rank (excluding Va) --------------

  // The input space
  Eigen::MatrixXd N = _species_frac_space(_N);

  int r = _rank(N, tol);

  if (r != N.cols()) {
    std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
    std::cerr << "Input space (column vectors of atom_frac):\n"
              << N << std::endl;
    std::cerr << "Rows correspond to: " << jsonParser(struc_mol_name)
              << std::endl;
    std::cerr << "Input space rank: " << r << std::endl;
    throw std::runtime_error(
        "Error in ChemicalReference::hyperplane: Too many reference states "
        "specified");
  }

  //  --- check that the input space spans the prim space --------------------

  // get the prim composition space (column vectors are comp_n)
  Eigen::MatrixXd C = composition_space(allowed_molecule_names(prim), tol);

  Index prim_N_mol = struc_mol_name.size();

  // ensure that there is a solution X for:
  //             C = N.topRows(prim_N_mol) * X
  //   (prim_space = input_space involving prim species * X)
  Eigen::MatrixXd X = N.topRows(prim_N_mol).fullPivHouseholderQr().solve(C);

  double relative_error = (N.topRows(prim_N_mol) * X - C).norm() / C.norm();

  if (relative_error > tol) {
    std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
    std::cerr << "Input space does not span the composition space of your prim."
              << std::endl;

    std::cerr << "Input space (column vectors in atom_frac space):\n"
              << N.topRows(prim_N_mol) << std::endl;
    std::cerr << "End members:\n"
              << end_members(allowed_molecule_names(prim)) << std::endl;
    std::cerr << "Prim space (column vectors in atom_frac space):\n"
              << C << std::endl;
    std::cerr << "Rows correspond to: " << jsonParser(struc_mol_name)
              << std::endl;
    std::cerr << "X, prim_space = input_space*X: \n" << X << std::endl;
    std::cerr << "input_space*X: \n" << N.topRows(prim_N_mol) * X << std::endl;
    std::cerr << "relative_error: " << relative_error << std::endl;
    std::cerr << "tol: " << tol << std::endl;

    throw std::runtime_error(
        "Error in ChemicalReference::hyperplane: Input space does not span "
        "prim space");
  }

  // --- solve N.transpose() * P = E, for P, the hyperplane reference --------

  Eigen::VectorXd P = _N.transpose().fullPivHouseholderQr().solve(E);

  relative_error = (_N.transpose() * P - E).norm() / E.norm();

  if (relative_error > tol) {
    std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
    std::cerr << "Could not solve for hyperplane reference." << std::endl;

    std::cerr << "Input space (column vectors in atom_frac space), N:\n"
              << N << std::endl;
    std::cerr << "Rows correspond to: " << jsonParser(struc_mol_name)
              << std::endl;
    std::cerr << "Solve: _N.transpose()*P = E" << std::endl;
    std::cerr << "_N.transpose():\n" << _N.transpose() << std::endl;
    std::cerr << "Reference state energies, E:\n" << E << std::endl;
    std::cerr << "P:\n" << P.transpose() << std::endl;
    std::cerr << "relative_err: " << relative_error << std::endl;
    std::cerr << "tol: " << tol << std::endl;

    throw std::runtime_error(
        "Error in ChemicalReference::hyperplane: Could not solve for "
        "hyperplane reference");
  }

  if (has_Va && !almost_zero(P(Va_index))) {
    std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
    std::cerr << "Non-zero pure Va reference: " << P.transpose() << std::endl;
    std::cerr << "Elements correspond to: " << jsonParser(struc_mol_name)
              << std::endl;

    throw std::runtime_error(
        "Error in ChemicalReference::hyperplane: Input space does not span "
        "prim space");
  }

  return P.head(prim_N_mol);
}

/// \brief Automatically set ChemicalReference using calculated Configurations
///        with 'extreme' compositions
ChemicalReference auto_chemical_reference(const PrimClex &primclex,
                                          double lin_alg_tol) {
  auto closest_calculated_config = [&](const Eigen::VectorXd &target) {
    // return name of Configuration with param_comp closest to target_param_comp
    //   tie break goes to first Configuration with fewest atoms
    //
    //   must be Configurations for which the relaxed_energy has been calculated
    auto begin = primclex.db<Configuration>().begin();
    auto end = primclex.db<Configuration>().end();
    auto res = end;

    double close_dist = std::numeric_limits<double>::max();

    for (auto it = begin; it != end; ++it) {
      if (!it->calc_properties().has_scalar("relaxed_energy")) {
        continue;
      }
      double curr_dist = (target - comp(*it)).norm();
      if ((res == end) ||
          (almost_equal(curr_dist, close_dist, TOL) &&
           it->size() < res->size()) ||
          (curr_dist < close_dist)) {
        res = it;
        close_dist = curr_dist;
      }
    }

    if (res == end) {
      std::cerr << "Error in auto_chemical_reference: Could not find any "
                   "configurations\n";
      throw std::runtime_error(
          "Error in auto_chemical_reference: Could not find any "
          "configurations");
    }
    return res->name();
  };

  // get number of reference states needed minus one... we'll also look for the
  // 'origin'
  int Naxes = primclex.composition_axes().independent_compositions();
  Eigen::VectorXd target = Eigen::VectorXd::Zero(Naxes);

  // get 'origin' configuration
  std::vector<std::string> ref_config;
  ref_config.push_back(closest_calculated_config(target));

  for (int i = 0; i < Naxes; i++) {
    target(i) = 1.0;

    // get end member configurations
    std::string configname = closest_calculated_config(target);

    // ensure no repeats
    if (std::find(ref_config.begin(), ref_config.end(), configname) !=
        ref_config.end()) {
      std::cerr << "Error in auto_chemical_reference: Could not find enough "
                   "calculated configurations\n";
      throw std::runtime_error(
          "Error in auto_chemical_reference: Could not find enough calculated "
          "configurations");
    }

    ref_config.push_back(configname);

    target(i) = 0.0;
  }

  // make sure there are the right number of references in species_frac space
  // (may be 1 extra if pure Va configuration is possible)
  _prune_ref_config(primclex, ref_config, lin_alg_tol);

  // construct ChemicalReferenceState
  std::vector<ChemicalReferenceState> ref_states;
  for (auto it = ref_config.begin(); it != ref_config.end(); ++it) {
    ref_states.emplace_back(*primclex.db<Configuration>().find(*it),
                            ConfigIO::SpeciesFrac(),
                            ConfigIO::relaxed_energy_per_species());
  }

  return ChemicalReference(primclex.prim(), ref_states.begin(),
                           ref_states.end(), lin_alg_tol);
}

// --- ChemicalReferencePrinter implementations -----------

ChemicalReferencePrinter::ChemicalReferencePrinter(
    std::ostream &_stream, const ChemicalReference &_ref, int _indent,
    int _indent_incr)
    : stream(_stream),
      indent(_indent),
      indent_incr(_indent_incr),
      ref(_ref),
      struc_mol_name(xtal::struc_molecule_name(ref.prim())) {}

void ChemicalReferencePrinter::incr() { indent += indent_incr; }

void ChemicalReferencePrinter::decr() { indent -= indent_incr; }

// print regular string
void ChemicalReferencePrinter::print(const std::string &str) {
  stream << std::string(indent, ' ') << str << "\n";
}

// print plane as '<indent>mol_i(1): energy_per_species\n', for each molecule
void ChemicalReferencePrinter::print(const Eigen::VectorXd &plane) {
  for (int i = 0; i < plane.size(); ++i) {
    if (!xtal::is_vacancy(struc_mol_name[i])) {
      stream << std::string(indent, ' ') << struc_mol_name[i]
             << "(1): " << plane(i) << "\n";
    }
  }
}

// print plane as:
// \code
// <indent>mol_i(N_i)mol_j(N_j)...: energy_per_species //for each ref state
// ...
// \endcode
void ChemicalReferencePrinter::print(
    const std::vector<ChemicalReferenceState> &ref_state_vec) {
  for (auto it = ref_state_vec.begin(); it != ref_state_vec.end(); ++it) {
    const ChemicalReferenceState &ref_state = *it;

    stream << std::string(indent, ' ');
    for (auto it = ref_state.species_num.begin();
         it != ref_state.species_num.end(); ++it) {
      double num = it->second;
      if (almost_zero(num, 1e-14)) {
        continue;
      }

      stream << it->first << "(";
      // print integer if ~integer
      if (almost_equal(std::round(num), num, 1e-14)) {
        stream << std::lround(num);
      } else {
        stream << num;
      }
      stream << ")";
    }
    stream << ": " << ref_state.energy_per_species << "\n";
  }
}

// print supercell/config specific plane as:
// \code
// <indent>name:
// <indent>  mol_i: energy_per_species // for each molecule
// \endcode
void ChemicalReferencePrinter::print(
    const std::pair<std::string, Eigen::VectorXd> &_pair) {
  stream << std::string(indent, ' ') << _pair.first << ":\n";
  incr();
  print(_pair.second);
  decr();
}

// print supercell/config specific ref states as:
// \code
// <indent>name:
// <indent>  mol_i(N_i)mol_j(N_j)...: energy_per_species // for each ref state
// \endcode
void ChemicalReferencePrinter::print(
    const std::pair<std::string, std::vector<ChemicalReferenceState> > &_pair) {
  stream << std::string(indent, ' ') << _pair.first << ":\n";
  incr();
  print(_pair.second);
  decr();
}

void ChemicalReferencePrinter::print_global() {
  print("Global chemical reference:");
  incr();
  (ref.global_ref_states().empty()) ? print(ref.global())
                                    : print(ref.global_ref_states());
  decr();
}

void ChemicalReferencePrinter::print_supercell() {
  if (ref.supercell().size()) {
    print("Supercell specific chemical references:");
    for (auto it = ref.supercell().begin(); it != ref.supercell().end(); ++it) {
      print_supercell(it->first);
    }
  }
}

void ChemicalReferencePrinter::print_supercell(const std::string &name) {
  auto it = ref.supercell().find(name);
  auto res = ref.supercell_ref_states().find(name);
  auto ref_state_end = ref.supercell_ref_states().end();
  (res != ref_state_end) ? print(*res) : print(*it);
}

void ChemicalReferencePrinter::print_config() {
  if (ref.config().size()) {
    print("Configuration specific chemical references:");
    for (auto it = ref.config().begin(); it != ref.config().end(); ++it) {
      print_config(it->first);
    }
  }
}

void ChemicalReferencePrinter::print_config(const std::string &name) {
  auto it = ref.config().find(name);
  auto res = ref.config_ref_states().find(name);
  auto ref_state_end = ref.config_ref_states().end();
  (res != ref_state_end) ? print(*res) : print(*it);
}

void ChemicalReferencePrinter::print_all() {
  print_global();
  if (ref.supercell().size()) {
    print("");
    print_supercell();
  }
  if (ref.config().size()) {
    print("");
    print_config();
  }
}
}  // namespace CASM
