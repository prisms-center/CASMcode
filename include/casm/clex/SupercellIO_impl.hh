#ifndef CASM_SupercellIO_impl
#define CASM_SupercellIO_impl

#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/SupercellIO.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

namespace ScelIO {

// --- template<typename Base> class SupercellCheckBase ---

template <typename Base>
SupercellCheckBase<Base>::SupercellCheckBase(std::string name, std::string desc)
    : Base(name, desc),
      m_refcell(nullptr),
      m_last_result(notstd::make_cloneable<result_type>()),
      m_last_scel(nullptr),
      m_last_unit(nullptr) {}

/// \brief Expects arguments of the form 'is_supercell_of(scelname)'
template <typename Base>
bool SupercellCheckBase<Base>::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
  if (args.empty()) {
    splt_vec.clear();
  }

  if (splt_vec.size() != 1) {
    std::stringstream ss;
    ss << this->name() << " expected 1 argument.  Received " << splt_vec.size()
       << ": '" << args << "'";
    throw std::runtime_error(ss.str());
  }

  m_refcell_name = args;
  return true;
}

/// \brief Set pointer to ref supercell
template <typename Base>
bool SupercellCheckBase<Base>::init(const Supercell &_tmplt) const {
  auto const &supercell_db = _tmplt.primclex().db<Supercell>();
  auto it = supercell_db.find(m_refcell_name);
  if (it == supercell_db.end()) {
    std::stringstream ss;
    ss << this->name() << " expected a supercell name.  Received '"
       << m_refcell_name << "': no supercell with this name was found";
    throw std::runtime_error(ss.str());
  }
  m_refcell = &*it;
  return true;
}

/// \brief col_header returns: {'name(refcell_name)'}
template <typename Base>
std::vector<std::string> SupercellCheckBase<Base>::col_header(
    const Supercell &_tmplt) const {
  return std::vector<std::string>{this->name() + "(" + m_refcell_name + ")"};
}

/// Call is_supercell using prim.factor_group() to try possible orientations
///
/// Returns (bool, SymOp op, Eigen::MatrixXi T) with:
/// - scel is supercell of unit?
/// - If true: scel.lattice() == apply(op, unit.lattice()) * T
///
template <typename Base>
const typename SupercellCheckBase<Base>::result_type &
SupercellCheckBase<Base>::_evaluate(const Supercell &scel,
                                    const Supercell &unit) const {
  if (&scel != m_last_scel || &unit != m_last_unit) {
    auto res = xtal::is_equivalent_superlattice(
        scel.lattice(), unit.lattice(), unit.prim().factor_group().begin(),
        unit.prim().factor_group().end(), unit.crystallography_tol());

    *m_last_result =
        std::make_tuple(res.first != unit.prim().factor_group().end(),
                        *res.first, iround(res.second));
  }
  return *m_last_result;
}

template <typename Base>
Supercell const &SupercellCheckBase<Base>::refcell() const {
  if (m_refcell == nullptr) {
    std::stringstream ss;
    ss << this->name() << " expected a supercell name.  Received '"
       << m_refcell_name << "': no supercell with this name was found";
    throw std::runtime_error(ss.str());
  }
  return *m_refcell;
}

}  // namespace ScelIO

}  // namespace CASM

#endif
