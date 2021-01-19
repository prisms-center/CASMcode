#ifndef CASM_ChemicalReference_impl
#define CASM_ChemicalReference_impl

#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {

/// \brief Construct global reference via range ChemicalReferenceState
///
template <typename RefStateIterator>
ChemicalReference::ChemicalReference(const BasicStructure &prim,
                                     RefStateIterator begin,
                                     RefStateIterator end, double tol)
    : HyperPlaneReferenceBase(Name, Desc, Eigen::VectorXd::Zero(0),
                              ConfigIO::AtomFrac()),
      m_prim(&prim) {
  set_global(begin, end, tol);
}

/// \brief Set global hyperplane reference
///
/// \param prim The Structure defining the composition space the reference
/// should span \param begin,end Iterators over a range of
/// ChemicalReferenceState \param tol Tolerance for checking that input spans
/// the prim composition
///            space and a solution for the hyperplane is found
///
/// sets global refrence to be the Eigen::VectorXd, R, that solves:
/// \code energy = R.dot(atom_frac) \endcode for each ChemicalReferenceState
///
///
/// - Inserts associated RefStateVec
template <typename RefStateIterator>
void ChemicalReference::set_global(RefStateIterator begin, RefStateIterator end,
                                   double tol) {
  Eigen::VectorXd ref = hyperplane(*m_prim, begin, end, tol);
  m_global_ref_vec = RefStateVec(begin, end);
  _global() = ref;
}

/// \brief Set hyperplane reference specialized for a Supercell
///
/// - Inserts associated RefStateVec
template <typename RefStateIterator>
void ChemicalReference::set_supercell(const std::string &scelname,
                                      RefStateIterator begin,
                                      RefStateIterator end, double tol) {
  Eigen::VectorXd ref = hyperplane(*m_prim, begin, end, tol);
  m_supercell_ref_map[scelname] = RefStateVec(begin, end);
  _supercell()[scelname] = ref;
}

/// \brief Set hyperplane reference specialized for a Configuration
///
/// - Inserts associated RefStateVec
template <typename RefStateIterator>
void ChemicalReference::set_config(const std::string &configname,
                                   RefStateIterator begin, RefStateIterator end,
                                   double tol) {
  Eigen::VectorXd ref = hyperplane(*m_prim, begin, end, tol);
  m_config_ref_map[configname] = RefStateVec(begin, end);
  _config()[configname] = ref;
}

/// \brief Convert a set of ChemicalReferenceState to a hyperplane, including
/// checks
///
/// \param prim The Structure defining the composition space the reference
/// should span \param begin,end Iterators over a range of
/// ChemicalReferenceState \param tol Tolerance for checking that input spans
/// the prim composition
///            space and a solution for the hyperplane is found
///
/// \returns Eigen::VectorXd, R, that solves: energy = R.dot(atom_frac) for
///          each ChemicalReferenceState
///
template <typename RefStateIterator>
Eigen::VectorXd ChemicalReference::hyperplane(const BasicStructure &prim,
                                              RefStateIterator begin,
                                              RefStateIterator end,
                                              double tol) {
  // store Molecule names in vector
  std::vector<std::string> struc_mol_name = struc_molecule_name(prim);

  // --- find any Molecule not in the prim, add to end of vector -------------

  // increase struc_mol_name to include all Molecule names in input
  // ensure no vacancies included
  for (auto it = begin; it != end; ++it) {
    for (auto mol_it = it->species_num.begin(); mol_it != it->species_num.end();
         ++mol_it) {
      if (xtal::is_vacancy(mol_it->first)) {
        throw std::runtime_error(
            "Error in ChemicalReference::hyperplane: Input should not include "
            "vacancies");
      }
      if (!contains(struc_mol_name, mol_it->first)) {
        struc_mol_name.push_back(mol_it->first);
      }
    }
  }

  // --- initialize vectors, matrices ----------------------------------------

  // reference 'relaxed_energy_per_species'
  Eigen::VectorXd E = Eigen::VectorXd::Zero(std::distance(begin, end));

  // column vector matrix of number of each Molecule in each reference state
  Eigen::MatrixXd N =
      Eigen::MatrixXd::Zero(struc_mol_name.size(), std::distance(begin, end));

  // --- get input values ---------------------------------------------------

  // populate E, N
  Index index = 0;
  for (auto it = begin; it != end; ++it, ++index) {
    E(index) = it->energy_per_species;
    for (auto mol_it = it->species_num.begin(); mol_it != it->species_num.end();
         ++mol_it) {
      N(find_index(struc_mol_name, mol_it->first), index) = mol_it->second;
    }
  }

  // use E, N to calculate hyperplane
  return _calc_hyperplane(prim, struc_mol_name, N, E, tol);
}
}  // namespace CASM

#endif
