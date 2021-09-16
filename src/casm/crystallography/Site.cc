#include "casm/crystallography/Site.hh"

#include <exception>
#include <set>
#include <string>
#include <vector>

#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/SymTools.hh"

namespace {
using namespace CASM;
bool occupant_dof_are_equal(const std::vector<xtal::Molecule> &LHS,
                            const std::vector<xtal::Molecule> &RHS,
                            double tol) {
  if (RHS.size() != LHS.size()) return false;
  Index j;
  for (Index i = 0; i < LHS.size(); ++i) {
    for (j = 0; j < RHS.size(); ++j) {
      if (LHS[i].identical(RHS[j], tol)) {
        break;
      }
    }
    if (j == RHS.size()) {
      return false;
    }
  }
  return true;
}
}  // namespace

namespace CASM {
namespace xtal {

// TODO: Can we delete this?
Site::Site(const Lattice &init_home)
    : Coordinate(init_home), m_label(-1), m_type_ID(-1) {}

//****************************************************

Site::Site(const Coordinate &init_pos, const std::string &occ_name)
    : Site(init_pos, std::vector<Molecule>{Molecule::make_atom(occ_name)}) {}

/// \brief Construct site with initial position and the allowed Molecule
Site::Site(const Coordinate &init_pos, const std::vector<Molecule> &site_occ,
           const std::map<std::string, SiteDoFSet> &site_dof)
    : Coordinate(init_pos),
      m_label(-1),
      m_type_ID(-1),
      m_occupant_dof(site_occ),
      m_dof_map(site_dof) {}

Site::Site(const Coordinate &init_pos, const std::vector<Molecule> &site_occ)
    : Site(init_pos, site_occ, std::map<std::string, SiteDoFSet>()) {}

/// \brief Construct site with initial position, allowed molecules (occupants),
/// and local degrees of freedom
Site::Site(const Coordinate &init_pos, const std::vector<Molecule> &site_occ,
           const std::vector<SiteDoFSet> &dofset_vec)
    : Site(init_pos, site_occ, make_dofset_map(dofset_vec)) {}

Site::~Site() {}

//****************************************************

const std::vector<Molecule> &Site::occupant_dof() const {
  return m_occupant_dof;
}

//****************************************************

Index Site::dof_size() const { return m_dof_map.size(); }

//****************************************************

SiteDoFSet const &Site::dof(std::string const &_dof_type) const {
  auto it = m_dof_map.find(_dof_type);
  if (it != m_dof_map.end())
    return (it->second);
  else
    throw std::runtime_error(
        std::string("In Structure::dof(), this structure does not contain any "
                    "global DoF's of type ") +
        _dof_type);
}

//****************************************************
bool Site::has_dof(std::string const &_dof_type) const {
  return occupant_dof().size() > 1 ||
         m_dof_map.find(_dof_type) != m_dof_map.end();
}

//****************************************************

std::vector<std::string> Site::dof_types() const {
  std::vector<std::string> result;
  // if(occupant_dof().size() > 1)
  // result.push_back(occupant_dof().type_name());
  for (auto it = m_dof_map.begin(); it != m_dof_map.end(); ++it)
    result.push_back(it->first);
  return result;
}

//****************************************************

bool Site::time_reversal_active() const {
  for (auto const &_dof : m_dof_map)
    if (_dof.second.traits().time_reversal_active()) return true;

  for (auto const &mol : this->occupant_dof())
    if (mol.time_reversal_active()) return true;
  return false;
}

//****************************************************

Index Site::label() const { return m_label; };

//****************************************************

Site &Site::operator+=(const Coordinate &translation) {
  Coordinate::operator+=(translation);
  return *this;
}

//****************************************************
/**
 *
 */
//****************************************************

Site &Site::operator-=(const Coordinate &translation) {
  Coordinate::operator-=(translation);
  return *this;
}

//*******************************************************************************************
/**
 *
 */
//*******************************************************************************************

bool Site::compare(const Coordinate &test_coord) const {
  return (min_dist(test_coord) < lattice().tol());
}

//*******************************************************************************************
/**
 *
 */
//*******************************************************************************************

bool Site::compare(const Site &test_site) const {
  return (compare_type(test_site) && min_dist(test_site) < lattice().tol());
}

//*******************************************************************************************
/**
 *
 */
//*******************************************************************************************

bool Site::compare_type(const Site &test_site) const {
  return _type_ID() == test_site._type_ID();
}

//*******************************************************************************************

bool Site::operator==(const Site &test_site) const {
  return (compare_type(test_site) && Coordinate::operator==(test_site));
}

//*******************************************************************************************

bool Site::almost_equal(const Site &test_site) const {
  return (compare_type(test_site) && dist(test_site) < lattice().tol());
}

//****************************************************

bool Site::contains(const std::string &name) const {
  for (Index i = 0; i < occupant_dof().size(); i++)
    if (occupant_dof()[i].contains(name)) {
      return true;
    }
  return false;
}

//****************************************************

bool Site::contains(const std::string &name, int &index) const {
  for (Index i = 0; i < occupant_dof().size(); i++)
    if (occupant_dof()[i].contains(name)) {
      index = i;
      return true;
    }
  return false;
}

//****************************************************

void Site::set_allowed_occupants(std::vector<Molecule> const &_new_occ_domain) {
  m_occupant_dof = _new_occ_domain;
  m_type_ID = -1;
}

//****************************************************

void Site::set_dofs(std::map<std::string, SiteDoFSet> _dofs) {
  m_dof_map = std::move(_dofs);
  m_type_ID = -1;
}

//****************************************************

std::vector<std::string> Site::allowed_occupants() const {
  std::vector<std::string> occ_list;
  for (Index i = 0; i < occupant_dof().size(); i++) {
    occ_list.push_back(occupant_dof()[i].name());
  }
  return occ_list;
}

//****************************************************

void Site::set_label(Index new_ind) {
  if (new_ind == m_label) return;
  m_label = new_ind;

  m_type_ID = -1;
  return;
}

//****************************************************
//   read site, including all possible occupants
void Site::read(std::istream &stream, bool SD_is_on) {
  set_label(-1);

  char ch;

  Eigen::Vector3i SD_flag;

  Coordinate::read(stream, COORD_MODE::CHECK());
  if (SD_is_on) {
    for (int i = 0; i < 3; i++) {
      stream >> ch;
      if (ch == 'T') {
        SD_flag[i] = 1;
      } else if (ch == 'F') {
        SD_flag[i] = 0;
      }
    }
  }

  std::vector<Molecule> tocc;

  ch = stream.peek();
  while (ch != '\n' && ch != ':' && !stream.eof()) {
    while ((ch < 'a' || ch > 'z') && (ch < 'A' || ch > 'Z') && ch != '\n' &&
           ch != ':' && !stream.eof()) {
      stream.ignore();
      ch = stream.peek();
    }
    if (ch != '\n' && ch != ':' && !stream.eof()) {
      // Need to change this part for real molecules
      std::string mol_name;
      stream >> mol_name;
      tocc.push_back(Molecule::make_atom(mol_name));
    }
    ch = stream.peek();
  }

  if (ch == ':') {
    stream.ignore();
    stream.ignore();

    std::string mol_name;
    stream >> mol_name;

    if (tocc.size()) {
      m_occupant_dof = tocc;
      Index index = tocc.size();
      for (Index i = 0; i < tocc.size(); i++)
        if (tocc[i].name() == mol_name) {
          index = i;
          break;
        }
      /* if(index == tocc.size()) { */
      /*   std::stringstream ss; */
      /*   ss << "Error in Site::read(): Occupying molecule not listed in
       * possible occupants" << std::endl; */
      /*   ss << "  occupying molecule name: " << mol_name << "  index: " <<
       * index << std::endl; */
      /*   ss << "  possible occupants: "; */
      /*   for(Index i = 0; i < tocc.size(); i++) */
      /*     ss << tocc[i].name() << " "; */
      /*   ss << " " << std::endl; */
      /*   throw std::runtime_error(ss.str()); */
      /* } */
      /* else { */
      /*   m_occupant_dof->set_value(index); */
      /* } */

    } else {
      throw std::runtime_error(
          "Error in Site::read(): Trying to read Site info, but no valid input "
          "was received.");
    }
    m_type_ID = -1;
    return;
  }

  if (tocc.size()) {
    m_occupant_dof = tocc;
  } else {
    throw std::runtime_error(
        "Error in Site::read(): Trying to read Site info, but no valid input "
        "was received.");
  }
  stream.ignore(1000, '\n');  // ????

  m_type_ID = -1;
  return;
}

//****************************************************
// read site, using 'elem' as site occupant domain
void Site::read(std::istream &stream, std::string &elem, bool SD_is_on) {
  char ch;

  set_label(-1);

  Eigen::Vector3i SD_flag;

  Coordinate::read(stream, COORD_MODE::CHECK());
  if (SD_is_on) {
    for (int i = 0; i < 3; i++) {
      stream >> ch;
      if (ch == 'T') {
        SD_flag[i] = 1;
      } else if (ch == 'F') {
        SD_flag[i] = 0;
      }
    }
  }

  std::vector<Molecule> tocc;

  tocc.push_back(Molecule::make_atom(elem));

  if (tocc.size()) {
    m_occupant_dof = tocc;
  } else {
    throw std::runtime_error(
        "Error in Site::read(): Trying to read Site info, but no valid input "
        "was received.");
  }
  stream.ignore(1000, '\n');

  m_type_ID = -1;
  return;
}

//****************************************************
/**	Print coordinate of site with name of all possible
 *		occupying molecule
 */
//****************************************************

void Site::print(std::ostream &stream, Eigen::IOFormat format) const {
  Coordinate::print(stream, 0, format);
  stream << " ";
  for (const Molecule &m : this->occupant_dof()) {
    stream << m.name() << "  ";
  }
  stream << std::endl;

  return;
}

void Site::print_occupant_dof(const std::vector<Molecule> &allowed_occupants,
                              std::ostream &out_stream) {
  for (const Molecule &m : allowed_occupants) {
    out_stream << m.name() << ' ';
  }
  return;
}

//****************************************************

std::ostream &operator<<(std::ostream &stream, const Site &site) {
  site.print(stream);
  return stream;
}

//*******************************************************************************************

bool Site::_compare_type_no_ID(const Site &_other) const {
  // compare domain but not value
  if (!(label() == _other.label() &&
        ::occupant_dof_are_equal(occupant_dof(), _other.occupant_dof(), TOL)))
    return false;

  if (m_dof_map.size() != _other.m_dof_map.size()) return false;

  auto it1 = m_dof_map.begin(), it2 = _other.m_dof_map.begin();
  for (; it1 != m_dof_map.end(); ++it1, ++it2)
    if (!SiteDoFSetIsEquivalent_f(it1->second, CASM::TOL)(it2->second)) {
      return false;
    }

  return true;
}

//*******************************************************************************************

Index Site::_type_ID() const {
  if (!valid_index(m_type_ID)) {
    for (m_type_ID = 0; m_type_ID < _type_prototypes().size(); m_type_ID++) {
      if (_compare_type_no_ID(_type_prototypes()[m_type_ID])) return m_type_ID;
    }

    _type_prototypes().push_back(*this);
  }
  return m_type_ID;
}

//****************************************************

Site operator*(const SymOp &LHS, const Site &RHS) {
  return sym::copy_apply(LHS, RHS);
}

//****************************************************

Site operator+(const Site &LHS, const Coordinate &RHS) {
  return Site(LHS) += RHS;
}

//****************************************************

Site operator+(const Coordinate &LHS, const Site &RHS) {
  return Site(RHS) += LHS;
}

}  // namespace xtal

namespace sym {
xtal::Site &apply(const xtal::SymOp &op, xtal::Site &mutating_site) {
  xtal::Site transformed_site = copy_apply(op, mutating_site);
  std::swap(transformed_site, mutating_site);
  return mutating_site;
}

xtal::Site copy_apply(const xtal::SymOp &op, xtal::Site site) {
  xtal::Coordinate transformed_coord =
      copy_apply(op, static_cast<xtal::Coordinate &>(site));

  std::vector<xtal::Molecule> transformed_occupants = site.occupant_dof();
  for (xtal::Molecule &occ : transformed_occupants) {
    apply(op, occ);
  }

  std::map<std::string, xtal::SiteDoFSet> transformed_dof;
  for (const auto &name_dof_pr : site.dofs()) {
    transformed_dof.emplace(name_dof_pr.first,
                            sym::copy_apply(op, name_dof_pr.second));
  }

  return xtal::Site(transformed_coord, transformed_occupants, transformed_dof);
}
}  // namespace sym
}  // namespace CASM
