#include "casm/symmetry/SymGroup.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"

namespace CASM {

namespace Local {

//  count number of operations of each type
//  proper before improper for each order, in ascending order
//  [identity, inversion, 2-fold, mirror, 3-fold, improper 3-fold, etc.]
static std::vector<Index> _number_of_operation_types(SymGroup const &group) {
  std::vector<Index> result(20, 0);
  for (SymOp const &op : group) {
    SymInfo info(op, group.lattice());
    if (info.op_type == symmetry_type::identity_op) {
      result[0]++;
    } else if (info.op_type == symmetry_type::inversion_op) {
      result[1]++;
    } else {
      Index bit = 0;
      if (op.matrix().determinant() < 0) bit = 1;
      for (Index n = 2; n <= 20; ++n) {
        if ((round(std::abs(info.angle)) * n) % 360 == 0) {
          result[2 * (n - 1) + bit]++;
          break;
        }
      }
    }
  }
  return result;
}

// populates point group infor for centric point group (groups with inversion)
static std::map<std::string, std::string> _centric_point_group_info(
    std::vector<Index> const &op_types) {
  bool is_error = false;
  std::map<std::string, std::string> result;
  if (op_types[4] == 4) {  // 3-fold
    result["crystal_system"] = "Cubic";
    result["international_name"] = "m-3";
    result["name"] = "Th";
    result["latex_name"] = "T_h";
    result["space_group_range"] = "200-206";
  } else if (op_types[4] == 8) {  // 3-fold
    result["crystal_system"] = "Cubic";
    result["international_name"] = "m-3m";
    result["name"] = "Oh";
    result["latex_name"] = "O_h";
    result["space_group_range"] = "221-230";
  } else if (op_types[10] == 2) {  // 6-fold
    result["crystal_system"] = "Hexagonal";
    if (op_types[2] == 1) {  // 2-fold
      result["international_name"] = "6/m";
      result["name"] = "C6h";
      result["latex_name"] = "C_{6h}";
      result["space_group_range"] = "175-176";
    } else if (op_types[2] == 7) {  // 2-fold
      result["international_name"] = "6/mmm";
      result["name"] = "D6h";
      result["latex_name"] = "D_{6h}";
      result["space_group_range"] = "191-194";
    } else {
      is_error = true;
    }
  }                             // end hexagonal
  else if (op_types[4] == 2) {  // 3-fold
    result["crystal_system"] = "Rhombohedral";
    if (op_types[2] == 0) {  // 2-fold
      result["international_name"] = "-3";
      result["name"] = "S6";
      result["latex_name"] = "S_6";
      result["space_group_range"] = "147-148";
    } else if (op_types[2] == 3) {  // 2-fold
      result["international_name"] = "-3m";
      result["name"] = "D3d";
      result["latex_name"] = "D_{3d}";
      result["space_group_range"] = "162-167";
    } else {
      is_error = true;
    }
  }                             // end rhombohedral
  else if (op_types[6] == 2) {  // 4-fold
    result["crystal_system"] = "Tetragonal";
    if (op_types[2] == 1) {  // 2-fold
      result["international_name"] = "4/m";
      result["name"] = "C4h";
      result["latex_name"] = "C_{4h}";
      result["space_group_range"] = "83-88";
    } else if (op_types[2] == 5) {  // 2-fold
      result["international_name"] = "4/mmm";
      result["name"] = "D4h";
      result["latex_name"] = "D_{4h}";
      result["space_group_range"] = "123-142";
    } else {
      is_error = true;
    }
  }  // end tetragonal
  else if (op_types[2] == 3) {
    result["crystal_system"] = "Orthorhomic";
    result["international_name"] = "mmm";
    result["name"] = "D2h";
    result["latex_name"] = "D_{2h}";
    result["space_group_range"] = "47-84";
  }  // end orthorhombic
  else if (op_types[2] == 1) {
    result["crystal_system"] = "Monoclinic";
    result["international_name"] = "2/m";
    result["name"] = "C2h";
    result["latex_name"] = "C_{2h}";
    result["space_group_range"] = "10-15";
  }  // end monoclinic
  else {
    Index tot_ops = 0;
    for (Index m : op_types) {
      tot_ops += m;
    }

    result["crystal_system"] = "Triclinic";

    if (op_types[1] == 1 && tot_ops == 2) {
      result["international_name"] = "-1";
      result["name"] = "Ci";
      result["latex_name"] = "C_i";
      result["space_group_range"] = "2";
    } else {
      is_error = true;
    }
  }

  if (is_error || result["crystal_system"].empty()) {
    throw std::runtime_error(
        "Error finding centric point group type. Crystal system determined to "
        "be " +
        result["crystal_system"]);
  }
  return result;
}

std::map<std::string, std::string> _acentric_point_group_info(
    std::vector<Index> const &op_types) {
  std::map<std::string, std::string> result;
  bool is_error = false;
  if (op_types[4] == 8) {  // 3-fold
    result["crystal_system"] = "Cubic";
    if (op_types[6] == 0) {  // 4-fold
      result["international_name"] = "23";
      result["name"] = "T (Chiral)";
      result["latex_name"] = "T";
      result["space_group_range"] = "195-199";
    } else if (op_types[6] == 6) {  // 4-fold
      result["international_name"] = "432";
      result["name"] = "O (Chiral)";
      result["latex_name"] = "O";
      result["space_group_range"] = "207-214";
    } else if (op_types[7] == 6) {  // improper 4-fold
      result["international_name"] = "-43m";
      result["name"] = "Td";
      result["latex_name"] = "T_d";
      result["space_group_range"] = "215-220";
    } else {
      is_error = true;
    }
  }                             // end cubic;
  else if (op_types[5] == 2) {  // Improper 3-fold
    result["crystal_system"] = "Hexagonal";
    if (op_types[2] == 3) {  // 2-fold
      result["international_name"] = "-6m2";
      result["name"] = "D3h";
      result["latex_name"] = "D_{3h}";
      result["space_group_range"] = "187-190";
    }

  }                              // end hexagonal 1
  else if (op_types[10] == 2) {  // 6-fold
    result["crystal_system"] = "Hexagonal";
    if (op_types[2] == 7) {  // 2-fold
      result["international_name"] = "622";
      result["name"] = "D6 (Chiral)";
      result["latex_name"] = "D_{6h}";
      result["space_group_range"] = "177-182";
    } else if (op_types[3] == 6) {  // mirror
      result["international_name"] = "6mm";
      result["name"] = "C6v";
      result["latex_name"] = "C_{6v}";
      result["space_group_range"] = "183-186";
    } else if (op_types[3] == 0) {  // mirror
      result["international_name"] = "6";
      result["name"] = "C6 (Chiral)";
      result["latex_name"] = "C_6";
      result["space_group_range"] = "168-173";
    } else {
      is_error = true;
    }
  }                             // end hexagonal 2
  else if (op_types[4] == 2) {  // 3-fold
    result["crystal_system"] = "Rhombohedral";
    if (op_types[2] == 3) {  // 2-fold
      result["international_name"] = "32";
      result["name"] = "D3 (Chiral)";
      result["latex_name"] = "D_3";
      result["space_group_range"] = "149-155";
    } else if (op_types[3] == 3) {  // mirror
      result["international_name"] = "3m";
      result["name"] = "C3v";
      result["latex_name"] = "C_{3v}";
      result["space_group_range"] = "156-161";
    } else if ((op_types[2] + op_types[3]) == 0) {
      result["international_name"] = "3";
      result["name"] = "C3";
      result["latex_name"] = "C_{3}";
      result["space_group_range"] = "143-146";
    } else if (op_types[3] == 1) {  // mirror
      result["international_name"] = "-6";
      result["name"] = "C3h";
      result["latex_name"] = "C_{3h}";
      result["space_group_range"] = "174";
    } else {
      is_error = true;
    }
  }                             // end rhombohedral
  else if (op_types[6] == 2) {  // 4-fold
    result["crystal_system"] = "Tetragonal";
    if (op_types[2] == 5) {  // 2-fold
      result["international_name"] = "422";
      result["name"] = "D4 (Chiral)";
      result["latex_name"] = "D_4";
      result["space_group_range"] = "89-98";
    } else if (op_types[3] == 4) {  // mirror
      result["international_name"] = "4mm";
      result["name"] = "C4v";
      result["latex_name"] = "C_{4v}";
      result["space_group_range"] = "99-110";
    } else if (op_types[2] == 1) {  // two-fold
      result["international_name"] = "4";
      result["name"] = "C4 (Chiral)";
      result["latex_name"] = "C_4";
      result["space_group_range"] = "75-80";
    } else {
      is_error = true;
    }
  }                             // end tetragonal 1
  else if (op_types[7] == 2) {  // improper 4-fold
    result["crystal_system"] = "Tetragonal";
    if (op_types[2] == 3) {  // 2-fold
      result["international_name"] = "-42m";
      result["name"] = "D2d";
      result["latex_name"] = "D_{2d}";
      result["space_group_range"] = "111-122";
    } else if (op_types[2] == 1) {  // 2-fold
      result["international_name"] = "-4";
      result["name"] = "S4";
      result["latex_name"] = "S_4";
      result["space_group_range"] = "81-82";
    } else {
      is_error = true;
    }
  }                                             // end tetragonal 2
  else if ((op_types[2] + op_types[3]) == 3) {  // 2-fold and mirror
    result["crystal_system"] = "Orthorhomic";
    if (op_types[2] == 3) {  // 2-fold
      result["international_name"] = "222";
      result["name"] = "D2 (Chiral)";
      result["latex_name"] = "D_2";
      result["space_group_range"] = "16-24";
    } else if (op_types[3] == 2) {  // mirror
      result["international_name"] = "mm2";
      result["name"] = "C2v";
      result["latex_name"] = "C_{2v}";
      result["space_group_range"] = "25-46";
    } else {
      is_error = true;
    }
  }                                             // end orthorhombic
  else if ((op_types[2] + op_types[3]) == 1) {  // 2-fold and mirror
    result["crystal_system"] = "Monoclinic";
    if (op_types[2] == 1) {  // 2-fold
      result["international_name"] = "2";
      result["name"] = "C2 (Chiral)";
      result["latex_name"] = "C_2";
      result["space_group_range"] = "3-5";
    } else if (op_types[3] == 1) {  // mirror
      result["international_name"] = "m";
      result["name"] = "Cs";
      result["latex_name"] = "C_s";
      result["space_group_range"] = "6-9";
    } else {
      is_error = true;
    }
  }  // end Acentric monoclinic
  else {
    Index tot_ops = 0;
    for (Index m : op_types) {
      tot_ops += m;
    }

    result["crystal_system"] = "Triclinic";

    if (tot_ops == 1) {
      result["international_name"] = "1";
      result["name"] = "C1 (Chiral)";
      result["latex_name"] = "C_1";
      result["space_group_range"] = "1";
    } else {
      is_error = true;
    }
  }

  if (is_error || result["crystal_system"].empty()) {
    std::stringstream ss;
    ss << "op_types: {";
    for (Index op_type : op_types) {
      ss << op_type << ", ";
    }
    ss << "}";
    throw std::runtime_error(
        "Error finding acentric point group type. Crystal system determined to "
        "be " +
        result["crystal_system"] + " " + ss.str());
  }

  return result;
}

// Does not work for icosahedral groups
static std::map<std::string, std::string> _nonmagnetic_point_group_info(
    SymGroup const &g) {
  SymGroup pg = g.copy_no_trans(false);
  std::vector<Index> op_types = _number_of_operation_types(pg);
  std::map<std::string, std::string> result;
  // Calculate total number of rotation elements
  Index nm = op_types[2] + 1;
  for (Index k = 4; k < op_types.size(); k++) {
    nm += op_types[k];
  }
  if (op_types[1]) {
    result = _centric_point_group_info(op_types);
    result["centricity"] = "Centric";
  } else {
    result = _acentric_point_group_info(op_types);
    result["centricity"] = "Acentric";
  }
  return result;
}
}  // namespace Local

// INITIALIZE STATIC MEMBER MasterSymGroup::GROUP_COUNT
// THIS MUST OCCUR IN A .CC FILE; MAY CAUSE PROBLEMS IF WE
// CHANGE COMPILING/LINKING STRATEGY
Index MasterSymGroup::GROUP_COUNT(0);

//*******************************************************************************************

MasterSymGroup::MasterSymGroup(const MasterSymGroup &RHS)
    : SymGroup(RHS),
      m_group_index(RHS.m_group_index),
      m_coord_rep_ID(RHS.m_coord_rep_ID),
      m_reg_rep_ID(RHS.m_reg_rep_ID),
      m_identity_rep_IDs(RHS.m_identity_rep_IDs) {
  m_rep_array.reserve(RHS.m_rep_array.size());

  for (Index i = 0; i < RHS.m_rep_array.size(); i++) {
    _add_representation(RHS.m_rep_array[i]->copy());
  }

  for (Index i = 0; i < size(); i++) at(i).set_index(*this, i);
}

//*******************************************************************************************

MasterSymGroup::~MasterSymGroup() {
  clear();
  return;
}

//*******************************************************************************************
MasterSymGroup &MasterSymGroup::operator=(const MasterSymGroup &RHS) {
  SymGroup::operator=(RHS);
  m_coord_rep_ID = RHS.m_coord_rep_ID;
  m_reg_rep_ID = RHS.m_reg_rep_ID;
  m_rep_array.reserve(RHS.m_rep_array.size());
  for (Index i = 0; i < RHS.m_rep_array.size(); i++)
    _add_representation(RHS.m_rep_array[i]->copy());

  for (Index i = 0; i < size(); i++) at(i).set_index(*this, i);

  return *this;
}

//*******************************************************************************************

void MasterSymGroup::push_back(const SymOp &op) {
  SymGroup::push_back(op);

  back().set_index(*this, size() - 1);
  return;
}

//*******************************************************************************************

const SymGroup &MasterSymGroup::point_group() const {
  if (!m_point_group.size()) {
    m_point_group = copy_no_trans(false);
  }
  return m_point_group;
}

//*******************************************************************************************

void MasterSymGroup::clear() {
  SymGroup ::clear();
  m_point_group.clear();
  for (Index i = 0; i < m_rep_array.size(); i++) {
    if (m_rep_array[i]) delete m_rep_array[i];
  }
  m_rep_array.clear();

  m_reg_rep_ID = m_coord_rep_ID = SymGroupRepID();

  m_identity_rep_IDs.clear();

  // Yes, by the time you return GROUP_COUNT is greater than m_group_index, and
  // that's how we like it around here.
  this->m_group_index = ++MasterSymGroup::GROUP_COUNT;

  return;
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::coord_rep_ID() const {
  if (m_coord_rep_ID.empty()) _add_coord_rep();
  return m_coord_rep_ID;
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::reg_rep_ID() const {
  if (m_reg_rep_ID.empty()) _add_reg_rep();
  return m_reg_rep_ID;
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::identity_rep_ID(Index dim) const {
  if (m_identity_rep_IDs.size() < dim + 1) {
    auto tail = std::vector<SymGroupRepID>(dim + 1 - m_identity_rep_IDs.size());
    m_identity_rep_IDs.insert(m_identity_rep_IDs.end(), tail.begin(),
                              tail.end());
  }
  if (m_identity_rep_IDs[dim].empty()) {
    m_rep_array.push_back(new SymGroupRep(*this));
    for (Index i = 0; i < size(); i++) {
      m_rep_array.back()->set_rep(i, SymPermutation(Permutation(dim)));
    }
    m_identity_rep_IDs[dim] =
        SymGroupRepID(group_index(), m_rep_array.size() - 1);
  }
  return m_identity_rep_IDs[dim];
}

//*******************************************************************************************

SymGroupRep const &MasterSymGroup::coord_rep() const {
  if (m_coord_rep_ID.empty()) _add_coord_rep();
  return representation(m_coord_rep_ID);
}

//*******************************************************************************************

void MasterSymGroup::set_rep(SymGroupRepID rep_ID,
                             SymOpRepresentation const &_op_rep,
                             Index op_index) const {
  _representation_ptr(rep_ID)->set_rep(op_index, _op_rep);
}

//*******************************************************************************************

SymGroupRep const &MasterSymGroup::reg_rep() const {
  if (m_reg_rep_ID.empty()) _add_reg_rep();
  return representation(m_reg_rep_ID);
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::_add_coord_rep() const {
  SymGroupRep *coordrep(new SymGroupRep(*this));
  for (Index i = 0; i < size(); i++)
    coordrep->set_rep(i, SymMatrixXd(at(i).matrix()));

  m_coord_rep_ID = _add_representation(coordrep);
  return m_coord_rep_ID;
}

//*******************************************************************************************
SymGroupRepID MasterSymGroup::_add_reg_rep() const {
  SymGroupRep *regrep(new SymGroupRep(*this));
  Eigen::MatrixXd regrep_mat(size(), size());

  if (get_alt_multi_table().size() != size()) return SymGroupRepID();

  for (Index i = 0; i < size(); i++) {
    for (Index j = 0; j < size(); j++) {
      for (Index k = 0; k < size(); k++) {
        if (alt_multi_table[j][k] == i) {
          regrep_mat(j, k) = 1.0;
        } else {
          regrep_mat(j, k) = 0.0;
        }
      }
    }
    regrep->set_rep(i, SymMatrixXd(regrep_mat));
  }

  m_reg_rep_ID = _add_representation(regrep);

  return m_reg_rep_ID;
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::add_kronecker_rep(SymGroupRepID ID1,
                                                SymGroupRepID ID2) const {
  SymGroupRep const *rep1(_representation_ptr(ID1)),
      *rep2(_representation_ptr(ID2));
  if (!(rep1 && rep2)) return SymGroupRepID();

  SymGroupRep *new_rep(new SymGroupRep(*this));
  Eigen::MatrixXd tmat;
  for (Index i = 0; i < size(); i++) {
    // std::cout << "rep1 matrix:\n";
    // std::cout << *(rep1->MatrixXd(i)) << '\n';
    // std::cout << "rep2 matrix:\n";
    // std::cout << *(rep2->MatrixXd(i)) << '\n';

    kroneckerProduct(*(rep1->MatrixXd(i)), *(rep2->MatrixXd(i)), tmat);
    // std::cout << "Total matrix:\n" << tmat << '\n';
    new_rep->set_rep(i, SymMatrixXd(tmat));
  }
  return _add_representation(new_rep);
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::add_direct_sum_rep(
    const std::vector<SymGroupRepID> &rep_IDs) const {
  std::vector<SymGroupRep const *> treps;
  for (Index i = 0; i < rep_IDs.size(); i++) {
    treps.push_back(_representation_ptr(rep_IDs[i]));
    if (!treps.back()) return SymGroupRepID();
  }
  SymGroupRep *new_rep(new SymGroupRep(*this));

  int dim = 0;
  for (Index i = 0; i < treps.size(); i++) {
    if (treps[i]->size() != size()) return SymGroupRepID();
    if (!(treps[i]->MatrixXd(0))) return SymGroupRepID();

    dim += (treps[i]->MatrixXd(0))->cols();
  }

  Eigen::MatrixXd tmat(Eigen::MatrixXd::Zero(dim, dim));
  int corner = 0;
  for (Index i = 0; i < size(); i++) {
    corner = 0;
    for (Index j = 0; j < treps.size(); j++) {
      tmat.block(corner, corner, (treps[j]->MatrixXd(i))->cols(),
                 (treps[j]->MatrixXd(i))->cols()) = *(treps[j]->MatrixXd(i));
      corner += (treps[j]->MatrixXd(i))->cols();
    }
    new_rep->set_rep(i, SymMatrixXd(tmat));
  }
  return _add_representation(new_rep);
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::add_rotation_rep() const {
  // make sure coord_rep_ID exists?
  // Index coord_rep_index((*this).coord_rep_ID());
  // make a new symmetry representation
  SymGroupRep *new_rep(new SymGroupRep(*this));

  // we are going to use our knowledge of how rotation
  // matrices transform under symmetry (in the coord_rep)
  // in order to build
  // up a representation in the rotation basis.
  for (Index i = 0; i < size(); i++) {
    // build each 4d representaion for a = 1, b = 1, ...
    // rp = rotation parameter
    Eigen::Matrix3d coord_rep_mat = at(i).matrix();
    if (!almost_equal(coord_rep_mat.determinant(), 1.)) {
      // std::cout << "IN IF" << std::endl;
      coord_rep_mat *= coord_rep_mat.determinant();
      // continue;
    }

    Eigen::Quaterniond opq(coord_rep_mat);

    // std::cout << "Quaternion " << opq << std::endl;

    // this quaternion can be converted to a 4d rotation as follows
    // q = q0 q1 q2 q3
    // Right multiplication : r * q = u
    //
    // | q0 -q1 -q2 -q3 | | r0 |   | u0 |
    // | q1  q0  q3 -q2 | | r1 |   | u1 |
    // | q2 -q3  q0  q1 | | r2 | = | u2 |
    // | q3  q2 -q1  q0 | | q3 |   | u3 |
    //
    //  and Left multiplication : q * r = u
    //
    // | q0 -q1 -q2 -q3 | | r0 |   | u0 |
    // | q1  q0 -q3  q2 | | r1 |   | u1 |
    // | q2  q3  q0 -q1 | | r2 | = | u2 |
    // | q3 -q2  q1  q0 | | q3 |   | u3 |

    // building both for now
    Eigen::MatrixXd right_q(4, 4);
    Eigen::MatrixXd left_q(4, 4);

    right_q(0, 0) = opq.w();
    right_q(0, 1) = -opq.x();
    right_q(0, 2) = -opq.y();
    right_q(0, 3) = -opq.z();

    right_q(1, 0) = opq.x();
    right_q(1, 1) = opq.w();
    right_q(1, 2) = opq.z();
    right_q(1, 3) = -opq.y();

    right_q(2, 0) = opq.y();
    right_q(2, 1) = -opq.z();
    right_q(2, 2) = opq.w();
    right_q(2, 3) = opq.x();

    right_q(3, 0) = opq.z();
    right_q(3, 1) = opq.y();
    right_q(3, 2) = -opq.x();
    right_q(3, 3) = opq.w();

    left_q(0, 0) = opq.w();
    left_q(0, 1) = -opq.x();
    left_q(0, 2) = -opq.y();
    left_q(0, 3) = -opq.z();

    left_q(1, 0) = opq.x();
    left_q(1, 1) = opq.w();
    left_q(1, 2) = -opq.z();
    left_q(1, 3) = opq.y();

    left_q(2, 0) = opq.y();
    left_q(2, 1) = opq.z();
    left_q(2, 2) = opq.w();
    left_q(2, 3) = -opq.x();

    left_q(3, 0) = opq.z();
    left_q(3, 1) = -opq.y();
    left_q(3, 2) = opq.x();
    left_q(3, 3) = opq.w();

    // std::cout << right_q << std::endl;

    new_rep->set_rep(i, SymMatrixXd(left_q));
  }

  for (Index i = 0; i < new_rep->size(); i++) {
    // std::cout << "Rota Rep final mats " << i << std::endl <<
    // *(new_rep->MatrixXd(i)) << std::endl;
  }

  return _add_representation(new_rep);
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::add_transformed_rep(
    SymGroupRepID orig_ID, const Eigen::MatrixXd &trans_mat) const {
  SymGroupRep const *trep(_representation_ptr(orig_ID));
  if (!trep) return SymGroupRepID();
  return add_representation(coord_transformed_copy(*trep, trans_mat));
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::allocate_representation() const {
  SymGroupRepID new_ID(group_index(), m_rep_array.size());
  m_rep_array.push_back(new SymGroupRep(*this, new_ID));
  return new_ID;
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::add_representation(
    const SymGroupRep &new_rep) const {
  return _add_representation(new_rep.copy());
}

//*******************************************************************************************

SymGroupRepID MasterSymGroup::_add_representation(SymGroupRep *new_rep) const {
  SymGroupRepID new_ID(group_index(), m_rep_array.size());
  m_rep_array.push_back(new_rep);
  m_rep_array.back()->set_master_group(*this, new_ID);
  return new_ID;
}

//*******************************************************************************************
const SymGroupRep &MasterSymGroup::representation(SymGroupRepID _id) const {
  return *_representation_ptr(_id);
}
//*******************************************************************************************
SymGroupRep *MasterSymGroup::_representation_ptr(SymGroupRepID _id) const {
  if (_id.is_identity()) {
    // _id.rep_index() stores dimension of representation
    _id = identity_rep_ID(_id.rep_index());
  } else if (_id.group_index() != group_index()) {
    throw std::runtime_error(
        "Attempting to access representation from MasterGroup #" +
        std::to_string(group_index()) + " that resides in MasterGroup #" +
        std::to_string(_id.group_index()));
  }

  return m_rep_array[_id.rep_index()];
}

//*******************************************************************************************

void MasterSymGroup::sort() {
  SymGroup::sort();
  m_point_group.clear();
  bool broken_check(false);
  std::vector<Index> perm_array(size(), 0);
  for (Index i = 0; i < size(); i++) {
    perm_array[i] = at(i).index();
    if (at(i).index() != i) {
      at(i).set_index(*this, i);
      broken_check = true;
    }
  }
  if (broken_check && m_rep_array.size()) {
    // err_log() << "WARNING: Order of symmetry operations has been altered by
    // MasterSymGroup::sort_by_class(). Attempting to repair "
    //          << m_rep_array.size() << " symmetry representations.\n";
    for (Index i = 0; i < m_rep_array.size(); i++) {
      std::vector<SymOpRepresentation *> new_rep;
      std::transform(perm_array.begin(), perm_array.end(),
                     std::back_inserter(new_rep),
                     [&](Index &idx) { return (*m_rep_array[i])[idx]; });
      m_rep_array[i]->swap(new_rep);
      for (Index j = 0; j < m_rep_array[i]->size(); j++) {
        if (m_rep_array[i]->at(j)) {
          (m_rep_array[i]->at(j))
              ->set_identifiers(*this, m_rep_array[i]->symrep_ID(),
                                j);  // perm_array[j]);
        }
      }
    }
  }

  return;
}

//*******************************************************************************************
SymGroup::~SymGroup() {
  clear();
  return;
}

//*******************************************************************************************
SymGroup SymGroup::lattice_point_group(Lattice const &_lat) {
  xtal::SymOpVector lattice_point_group_operations =
      xtal::make_point_group(_lat);
  SymGroup point_group = adapter::Adapter<SymGroup, xtal::SymOpVector>()(
      lattice_point_group_operations, _lat);

  if (!point_group.is_group(_lat.tol())) {
    std::cerr << "*** WARNING *** \n"
              << "    SymGroup::lattice_point_group() has been called on an "
                 "ill-conditioned lattice \n"
              << "    (i.e., a well-defined point group could not be found "
                 "with the current tolerance of "
              << _lat.tol() << ").\n"
              << "    CASM will use the group closure of the symmetry "
                 "operations that were found.  Please consider using the \n"
              << "    CASM symmetrization tool on your input files.\n";
    std::cerr << "lat_column_mat:\n" << _lat.lat_column_mat() << "\n\n";

    point_group.enforce_group(_lat.tol());
  }
  // Sort point_group by trace/conjugacy class
  point_group.sort();

  return point_group;
}

//*******************************************************************************************
SymGroup::SymGroup(std::vector<SymOp> from_array, Lattice const *_lat_ptr,
                   PERIODICITY_TYPE init_type)
    : m_lat_ptr(_lat_ptr), m_group_periodicity(init_type), m_max_error(-1) {
  std::vector<SymOp>::swap(from_array);
}

//*******************************************************************************************

void SymGroup::set_lattice(const Lattice &lat) { m_lat_ptr = &lat; }

/*****************************************************************/

void SymGroup::push_back(const SymOp &new_op) {
  std::vector<SymOp>::push_back(new_op);
  if (back().map_error() > m_max_error) m_max_error = back().map_error();
  return;
}

//*******************************************************************************************

void SymGroup::clear_tables() {
  multi_table.clear();
  alt_multi_table.clear();
  conjugacy_classes.clear();
  class_names.clear();
  index2conjugacy_class.clear();
  irrep_IDs.clear();

  m_subgroups.clear();

  centralizer_table.clear();
  elem_order_table.clear();

  name.clear();
  latex_name.clear();
  comment.clear();

  return;
}

//*****************************************************************

/// Return the MasterSymGroup indices of the operations in this SymGroup
/// (equivalent to master_group_indices)
std::vector<Index> SymGroup::op_indices() const {
  return master_group_indices();
}

//*****************************************************************

/// Return the MasterSymGroup indices of the operations in this SymGroup
std::vector<Index> SymGroup::master_group_indices() const {
  std::vector<Index> ind_array(size());
  for (Index i = 0; i < size(); i++) {
    if (!valid_index(at(i).index()) && at(i).has_valid_master()) {
      ind_array[i] = at(i).master_group().find_periodic(at(i));
    } else {
      ind_array[i] = at(i).index();
    }
  }
  return ind_array;
}

//*****************************************************************

void SymGroup::clear() {
  std::vector<SymOp>::clear();

  clear_tables();

  m_max_error = -1;

  return;
}

//*****************************************************************
std::map<std::string, std::string> point_group_info(SymGroup const &g) {
  SymGroup pg = g.copy_no_trans(false);
  SymGroup nonmag;
  bool grey = false;
  for (SymOp const &op : pg) {
    if (!op.time_reversal()) {
      nonmag.push_back(op);
    } else if (!grey && op.matrix().isIdentity(TOL)) {
      grey = true;
    }
  }
  nonmag.set_lattice(g.lattice());
  auto nonmag_info = Local::_nonmagnetic_point_group_info(nonmag);

  // nonmagnetic group:
  if (nonmag.size() == pg.size()) {
    return nonmag_info;
  }

  // grey group:
  if (grey) {  //(2 * nonmag.size()) == pg.size()) {
    nonmag_info["space_group_range"] = "Magnetic group (not supported)";
    nonmag_info["international_name"] += "1'";
    nonmag_info["name"] += "'";
    nonmag_info["latex_name"].append("'");
    return nonmag_info;
  }

  auto tot_info = Local::_nonmagnetic_point_group_info(pg);

  // magnetic group:
  tot_info["international_name"] = "Magnetic group (not supported)";
  tot_info["space_group_range"] = "Magnetic group (not supported)";
  tot_info["name"] += ("(" + nonmag_info["name"] + ")");
  tot_info["latex_name"] += ("(" + nonmag_info["latex_name"] + ")");

  return tot_info;
}

//*******************************************************************************************
/* Loop through SymGroup and populate passed SymGroup
 * with all the operations but none of the shifts. Unle
 * explicitly requested, degenerate symmetry operations
 * will not be added.
 */
//*******************************************************************************************
SymGroup SymGroup::copy_no_trans(bool keep_repeated) const {
  SymGroup result;

  result.m_lat_ptr = m_lat_ptr;
  for (Index i = 0; i < size(); i++) {
    SymOp tsymop(at(i).no_trans());
    if (keep_repeated || !contains(result, tsymop)) {
      result.push_back(tsymop);
    }
  }
  return result;
}

//*******************************************************************************************
SymInfo SymGroup::info(Index i) const { return SymInfo(at(i), lattice()); }

//*******************************************************************************************

double SymGroup::max_error() { return m_max_error; }

//*******************************************************************************************

const Lattice &SymGroup::lattice() const {
  assert(m_lat_ptr &&
         "Attempting to access Lattice of a SymGroup, but it has not been "
         "initialized!");
  return *m_lat_ptr;
}

//*******************************************************************************************
/** This will name your conjugacy classes according to some
 *  crazy convention, which is similar, but not exactly
 *  identical to the Schoenflies notation.
 *
 *  We start by finding the principal axes. Often, there will
 *  only be one. However, in cases like fcc, there can be
 *  several. All rotations that are about the principal axis
 *  will be named "XCn" if they are proper and "XSn" if they
 *  are improper. Here, X is the size of their conjugacy
 *  class, and n is the foldedness of the rotation, i.e. 60
 *  degree rotations are 6-fold, while 120 degree rotations
 *  are 3-fold.
 */

void SymGroup::_generate_class_names() const {  // AAB
  if (!conjugacy_classes.size()) {
    _generate_conjugacy_classes();
  }

  class_names.resize(conjugacy_classes.size());

  std::vector<Eigen::Vector3d> highsym_axes;
  std::vector<int> mult;
  double angle = 360;
  double dprod;
  std::string symtype;
  Index to_name = conjugacy_classes.size();
  bool cubic = false;

  if ((name == "T") || (name == "Td") || (name == "Th") || (name == "O") ||
      (name == "Oh")) {
    cubic = true;
  }

  class_names[0] = "E";
  to_name--;

  std::vector<SymInfo> info;
  for (Index i = 0; i < size(); i++) info.push_back(SymInfo(at(i), lattice()));

  // If the SymGroup includes an inversion operation, name it "i"
  // We can just use back() here because inversion is always last
  if (info.back().op_type == symmetry_type::inversion_op) {
    class_names[class_names.size() - 1] = "i";
    to_name--;
  }

  if (to_name == 0) {
    return;
  }

  // C1h contains only the identity element and sigma h mirror plane, so we deal
  // with this specifically...
  if (name == "C1h") {
    class_names[1] = "h";
    to_name--;
  }

  if (to_name == 0) {
    return;
  }

  // This part will loop over all of the proper rotations and look for the axis
  // of highest rotational symmetry, which we are going to call the "principal
  // axis."
  highsym_axes.clear();
  mult.clear();

  for (Index i = 0; i < size(); i++) {
    if ((info[i].op_type == symmetry_type::rotation_op) ||
        (info[i].op_type == symmetry_type::screw_op)) {
      if (contains(highsym_axes,
                   info[i].axis.const_cart())) {  // Otherwise, check if the
                                                  // axis has been found;
        mult[find_index(highsym_axes, info[i].axis.const_cart())]++;
      } else {
        highsym_axes.push_back(info[i].axis.cart());
        mult.push_back(1);
      }
    }
  }

  Eigen::Vector3d hs_axis;
  bool hs_axis_set = false;
  double hangle = 360;

  for (Index i = 0; i < size(); i++) {
    if ((info[i].angle < hangle) && (info[i].angle > TOL)) {
      hangle = info[i].angle;
      hs_axis = info[i].axis.cart();
      hs_axis_set = true;
    }
  }

  int order = *std::max_element(mult.begin(), mult.end());

  for (int i = (int(mult.size()) - 1); i >= 0; i--) {
    if ((cubic == false) && (mult[i] != order)) {
      highsym_axes.erase(highsym_axes.begin() + i);
      mult.erase(mult.begin() + i);
    } else if (mult[i] < (order - 1)) {
      highsym_axes.erase(highsym_axes.begin() + i);
      mult.erase(mult.begin() + i);
    }
  }

  // This is kind of a hack... We might want to change this so that the
  // principal axis is actually the one with the highest fold rotation, and then
  // treat cases like cubic and orthorhombic separately....... but I am not sure
  // that's better...

  if (name == "D3d") {
    for (int i = int(mult.size()) - 1; i >= 0; i--) {
      if (!hs_axis_set) {
        throw std::runtime_error(
            "Error in _generate_class_names: using hs_axis unitialized");
      }
      if (!almost_zero(highsym_axes[i] - hs_axis)) {
        highsym_axes.erase(highsym_axes.begin() + i);
        mult.erase(mult.begin() + i);
      }
    }
  }

  int ind = 0;

  std::string cprime = "'";
  std::string sprime = "'";
  std::string vprime = "";

  for (Index i = 0; i < size(); i++) {  // Loop over all SymOps
    // ind = conjugacy_corr_table[i][0];
    ind = index2conjugacy_class[i];
    std::ostringstream s;
    if (!class_names[ind]
             .size()) {  // Check to see if this has already been named
      bool normal = false;
      for (Index j = 0; j < highsym_axes.size(); j++) {
        dprod = highsym_axes[j].dot(info[i].axis.const_cart());
        Eigen::Vector3d xprodvec =
            highsym_axes[j].cross(info[i].axis.const_cart());
        if (almost_zero(xprodvec.norm())) {
          normal = true;
          break;
        }
      }

      if (normal) {  // Check if the cross product with principal axis is zero
        if ((info[i].angle < 200) &&
            (info[i].angle >
             1)) {  // Only bother with angles that 360 is divisible by
          if ((info[i].op_type == symmetry_type::rotation_op) ||
              (info[i].op_type == symmetry_type::screw_op)) {
            angle = info[i].angle;
            symtype = "C";
            s << conjugacy_classes[ind].size() << symtype << int(360 / angle);
            class_names[ind] = s.str();
            to_name--;
            // std::cout << "  rotation\t==> " << s.str() << std::endl;
          } else if (info[i].op_type == symmetry_type::rotoinversion_op) {
            angle = info[i].angle;
            symtype = "S";
            s << conjugacy_classes[ind].size() << symtype << int(360 / angle);
            class_names[ind] = s.str();
            to_name--;
            // std::cout << "  rotoinversion\t==> " << s.str() << std::endl;
          }
        } else if ((info[i].op_type == symmetry_type::mirror_op) ||
                   (info[i].op_type == symmetry_type::glide_op)) {
          symtype = "h";
          s << conjugacy_classes[ind].size() << symtype;
          class_names[ind] = s.str();
          to_name--;
          // std::cout << "  mirror\t==> " << s.str() << std::endl;
        }
      } else {
        if ((info[i].angle < 200) && (info[i].angle > 1)) {
          if ((info[i].op_type == symmetry_type::rotation_op) ||
              (info[i].op_type == symmetry_type::screw_op)) {
            angle = info[i].angle;
            symtype = "C";
            s << conjugacy_classes[ind].size() << symtype << int(360 / angle)
              << cprime;
            class_names[ind] = s.str();
            to_name--;
            // std::cout << "  rotation\t==> " << s.str() << std::endl;
            cprime = "''";
          } else if (info[i].op_type == symmetry_type::rotoinversion_op) {
            angle = info[i].angle;
            symtype = "S";
            s << conjugacy_classes[ind].size() << symtype << int(360 / angle)
              << sprime;
            class_names[ind] = s.str();
            to_name--;
            // std::cout << "  rotoinversion\t==> " << s.str() << std::endl;
            sprime = "''";
          }
        } else if ((info[i].op_type == symmetry_type::mirror_op) ||
                   (info[i].op_type == symmetry_type::glide_op)) {
          if (almost_zero(dprod)) {
            symtype = "v";
            s << conjugacy_classes[ind].size() << symtype << vprime;
            class_names[ind] = s.str();
            to_name--;
            vprime = "'";
          } else {
            symtype = "d";
            s << conjugacy_classes[ind].size() << symtype;
            class_names[ind] = s.str();
            to_name--;
          }
          // std::cout << "  mirror\t==> " << s.str() << std::endl;
        }
      }
    }
  }

  for (Index i = 0; i < class_names.size(); i++) {
    std::string one = "1";
    if (class_names[i].find(one) != std::string::npos) {
      class_names[i].erase(class_names[i].find(one), 1);
    }
  }

  // This is also a hack... but again, not sure there is a better way...

  if (name == "D2h") {
    for (Index i = 0; i < conjugacy_classes.size(); i++) {
      if (info[conjugacy_classes[i][0]].op_type == symmetry_type::rotation_op) {
        if (almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -
                        Eigen::Vector3d::UnitX())) {
          class_names[i].append("(x)");
        } else if (almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -
                               Eigen::Vector3d::UnitY())) {
          class_names[i].append("(y)");
        } else if (almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -
                               Eigen::Vector3d::UnitY())) {
          class_names[i].append("(z)");
        }
      } else if (info[conjugacy_classes[i][0]].op_type ==
                 symmetry_type::mirror_op) {
        if (almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -
                        Eigen::Vector3d::UnitX())) {
          class_names[i].append("(yz)");
        } else if (almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -
                               Eigen::Vector3d::UnitY())) {
          class_names[i].append("(xz)");
        } else if (almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -
                               Eigen::Vector3d::UnitZ())) {
          class_names[i].append("(xy)");
        }
      }
    }
  }

  return;
}

//*******************************************************************************************

/** The centralizer of a subset S of a group G is the set of elements
 *  of G that commute with each element of S. The centralizer is a
 *  subgroup of G.
 *
 *  In our specific case, the subsets are the conjugacy classes, so
 *  we fill the centralizer_table with subgroups of G that commute
 *  with the corresponding conjugacy class.qq
 */

void SymGroup::_generate_centralizers() const {  // AAB
  if (!conjugacy_classes.size()) {
    _generate_conjugacy_classes();
  }

  bool all_commute = true;

  centralizer_table.resize(conjugacy_classes.size());

  for (Index i = 0; i < conjugacy_classes.size(); i++) {
    all_commute = true;
    for (Index j = 0; j < conjugacy_classes[i].size() && all_commute; j++) {
      for (Index k = 0; k < multi_table.size(); k++) {
        int ind = conjugacy_classes[i][j];
        if ((multi_table[ind][k] == multi_table[k][ind]) &&
            (!contains(centralizer_table[i], k))) {
          centralizer_table[i].push_back(k);
        } else {
          all_commute = false;
        }
      }
    }
  }
}

//*******************************************************************************************
void SymGroup::_generate_elem_order_table() const {
  if (!get_multi_table().size()) {
    return;
  }

  elem_order_table.resize(size());

  for (Index i = 0; i < size(); i++) {
    elem_order_table[i].push_back(i);
  }

  Index telem = 0;
  Index ind = 0;

  for (Index i = 0; i < size(); i++) {
    telem = i;
    ind = 0;
    do {
      telem = multi_table[i][telem];
      ind++;
    } while (telem != i);
    elem_order_table[i].push_back(ind);
  }

  return;
}

//*******************************************************************************************
std::vector<std::set<std::set<Index>>> SymGroup::_small_subgroups() const {
  std::vector<std::set<std::set<Index>>> result;

  int tempind = 0;
  Index i, j;
  if (get_alt_multi_table().size() != size()) {
    return result;
  }
  // identity is a small subgroup
  result.push_back({{0}});

  for (i = 1; i < multi_table.size(); i++) {
    std::set<Index> tgroup({0, i});
    j = i;  // ind_prod(i, i);
    while (j != 0) {
      j = ind_prod(i, j);
      tgroup.insert(j);
    }

    bool add = true;
    for (auto const &orbit : result) {
      if (orbit.begin()->size() == tgroup.size() && orbit.count(tgroup) > 0) {
        add = false;
        break;
      }
    }
    if (!add) continue;

    // use equiv_map to find the equivalent subgroups
    result.push_back({});
    // std::cout << "tgroup.size(): " << tgroup.size() << std::endl;
    // std::cout << "tgroup : ";
    // for(i : tgroup)
    // std::cout << i << "  ";
    // std::cout << std::endl;
    for (auto const &coset : left_cosets(tgroup.begin(), tgroup.end())) {
      std::set<Index> tequiv;
      for (Index op : tgroup) {
        tempind = ind_prod(op, ind_inverse(coset[0]));
        tequiv.insert(ind_prod(coset[0], tempind));
      }
      result.back().insert(std::move(tequiv));
    }
  }
  return result;
}

//*******************************************************************************************
/*  Start with m_subgroups = m_small_subgroups, then add new subgroups by
 *  finding the closure of a union of a large_group and a
 *  small_group. If the the new large_group is unique, add it as a
 *  large_group. Repeat for all (large_group, small_group) pairs,
 *  until no new m_subgroups are found. This is probably not the
 *  fastest algorithm, but it is complete
 */

void SymGroup::_generate_subgroups() const {
  auto small = _small_subgroups();
  m_subgroups = small;
  Index i, k, ii, jj, tempind;
  for (i = 0; i < m_subgroups.size(); i++) {
    // std::cout << "i is " << i << " and m_subgroups.size() is " <<
    // m_subgroups.size() << std::endl;
    for (auto const &orbit : small) {
      for (auto const &equiv : orbit) {
        std::set<Index> tgroup = *(m_subgroups[i].begin());
        Index init_size = tgroup.size();
        tgroup.insert(equiv.begin(), equiv.end());
        if (tgroup.size() == init_size) continue;

        // find group closure
        std::vector<Index> vgroup(tgroup.begin(), tgroup.end());
        for (ii = 0; ii < vgroup.size(); ii++) {
          for (jj = 0; jj < vgroup.size(); jj++) {
            Index prod = ind_prod(vgroup[ii], vgroup[jj]);
            if (tgroup.insert(prod).second) {
              vgroup.push_back(prod);
            }
          }
        }

        for (k = 0; k < m_subgroups.size(); k++) {
          if (m_subgroups[k].begin()->size() == tgroup.size() &&
              m_subgroups[k].count(tgroup) > 0)
            break;
        }
        // std::cout << " k is " << k << " and m_subgroups.size() is " <<
        // m_subgroups.size() << std::endl;
        if (k < m_subgroups.size()) continue;
        // add the new group

        // use equiv_map to find the equivalent subgroups
        m_subgroups.push_back({});

        for (auto const &coset : left_cosets(tgroup.begin(), tgroup.end())) {
          std::set<Index> tequiv;
          for (Index op : tgroup) {
            tempind = ind_prod(op, ind_inverse(coset[0]));
            tequiv.insert(ind_prod(coset[0], tempind));
          }
          m_subgroups.back().insert(std::move(tequiv));
        }
      }
    }
  }

  // Sort subgroups by number of elements and multiplicity
  for (Index i = 0; i < m_subgroups.size(); i++) {
    for (Index j = i + 1; j < m_subgroups.size(); j++) {
      if (m_subgroups[i].begin()->size() < m_subgroups[j].begin()->size())
        std::swap(m_subgroups[i], m_subgroups[j]);
      else if (m_subgroups[i].begin()->size() ==
                   m_subgroups[j].begin()->size() &&
               m_subgroups[i].size() < m_subgroups[j].size())
        std::swap(m_subgroups[i], m_subgroups[j]);
    }
  }
  return;
}

//*******************************************************************************************

std::vector<SymGroup> SymGroup::unique_subgroups() const {
  if (!m_subgroups.size()) _generate_subgroups();

  std::vector<std::string> sg_names, sg_names_limited;
  std::vector<bool> chosen_flag(m_subgroups.size(), false);
  for (Index i = 0; i < m_subgroups.size(); i++) {
    SymGroup sgroup;
    sgroup.m_lat_ptr = m_lat_ptr;
    for (Index op : *(m_subgroups[i].begin())) {
      sgroup.push_back(at(op));
    }

    sg_names.push_back(sgroup.get_name());
  }

  std::vector<std::vector<Index>> sg_tree(m_subgroups.size(),
                                          std::vector<Index>());
  for (Index i = 0; i < m_subgroups.size(); i++) {
    // std::cout << "Subgroup " << sg_names[i] << "-" << i << " is also a
    // subgroup of ";
    for (Index j = 0; j < m_subgroups.size(); j++) {
      for (auto const &equiv : m_subgroups[j]) {
        bool add = true;
        for (auto const &op : equiv) {
          if (m_subgroups[i].begin()->count(op) < 1) {
            add = false;
            break;
          }
        }
        if (add) {
          sg_tree[i].push_back(j);
          // std::cout << sg_names[j] << "-" << j << "-" << jj << "  ";
          break;
        }
      }
    }
    // std::cout << std::endl;
  }

  /* // commented out for now because datatype of m_subgroups has changed
  //attempt to maximize coincidence a
  int max_co, t_co;
  for(int i = int(m_subgroups.size()) - 1; i >= 0; i--) {
  if(chosen_flag[i]) continue;
  //chosen_flag[i]=true;
  for(int j = 0; j < i; j++) {
  if(sg_names[j] == sg_names[i]) {
  chosen_flag[j] = true;
  continue;
  }
  max_co = 0;
  for(Index jj = 0; jj < m_subgroups[j].size(); jj++) {
  t_co = m_subgroups[i][0].coincidence(m_subgroups[j][jj]);
  if(t_co > max_co) {
  max_co = t_co;
  m_subgroups[j].swap_elem(jj, 0);
  }
  }
  }
  }
  */
  std::vector<SymGroup> unique_sgroups;

  for (Index i = 0; i < m_subgroups.size(); i++) {
    if (chosen_flag[i]) continue;
    unique_sgroups.push_back(SymGroup());
    for (Index op : *(m_subgroups[i].begin())) {
      unique_sgroups.back().push_back(at(op));
    }
    unique_sgroups.back().sort();
  }

  return unique_sgroups;
}

//*******************************************************************************************
/** The number of irreducible representations is equal to
 *  the number of conjugacy classes.
 *
 *  The number of elements in each conjugacy class must be
 *  a divisor of the number of symmetry group operations.
 */
//*******************************************************************************************

void SymGroup::_generate_conjugacy_classes() const {
  if (get_multi_table().size() != size()) {
    return;
  }

  conjugacy_classes.clear();

  int k;

  for (Index i = 0; i < size(); i++) {
    bool dup_class(false);
    for (Index j = 0; j < conjugacy_classes.size(); j++) {
      if (contains(conjugacy_classes[j], i)) {
        dup_class = true;
        break;
      }
    }
    if (dup_class) continue;

    conjugacy_classes.push_back(std::vector<Index>());

    for (Index j = 0; j < size(); j++) {
      // std::cout << "for j=" << j << ", i=" << i << ": j-inverse= " <<
      // ind_inverse(j) << ", i*j-inverse= " << ind_prod(i, ind_inverse(j)); int
      // tk=alt_multi_table[j][i]; std::cout << "-- compare to amt[j][0]=" <<
      // alt_multi_table[j][0] << " and amt[j][i]=" << tk << " and result is
      // k=";
      k = ind_prod(j, ind_prod(i, ind_inverse(j)));
      // std::cout << k << " -- compare to explicit value " <<
      // multi_table[tk][j];

      if (!contains(conjugacy_classes.back(), k)) {
        // std::cout << " so " << k << " goes in class " <<
        // conjugacy_classes.size()-1;
        conjugacy_classes.back().push_back(k);
      }
      // std::cout << std::endl;
    }
    std::sort(conjugacy_classes.back().begin(), conjugacy_classes.back().end());
  }

  index2conjugacy_class.resize(size(), 0);
  for (Index i = 0; i < conjugacy_classes.size(); i++) {
    for (Index j = 0; j < conjugacy_classes[i].size(); j++) {
      index2conjugacy_class[conjugacy_classes[i][j]] = i;
    }
  }
  // std::cout << "index2conjugacy is " << index2conjugacy_class << '\n';
  return;
}

//*******************************************************************************************

Index SymGroup::ind_inverse(Index i) const {
  if (get_alt_multi_table().size() != size() || !valid_index(i) || i >= size())
    return -1;
  // std::cout << "Inside ind_inverse. 'i' is " << i << " and alt_multi_table
  // size is " << alt_multi_table.size() << std::endl;
  return alt_multi_table[i][0];
}

//*******************************************************************************************

Index SymGroup::ind_prod(Index i, Index j) const {
  if (get_multi_table().size() != size() || !valid_index(i) || i >= size() ||
      !valid_index(j) || j >= size()) {
    // err_log() << "WARNING: SymGroup::ind_prod() failed for " << i << ", " <<
    // j << std::endl;// and multi_table\n" << get_multi_table() << "\n\n";
    // assert(0);
    return -1;
  }
  // std::cout << "Inside ind_prod. 'i' is " << i << " and j is " << j << " and
  // multi_table size is " << multi_table.size() << std::endl;
  return multi_table[i][j];
}

//*******************************************************************************************

Index SymGroup::class_of_op(Index i) const {
  if (!get_conjugacy_classes().size() || i > size()) return -1;
  return index2conjugacy_class[i];
}

//*******************************************************************************************

void SymGroup::set_irrep_ID(Index i, SymGroupRepID ID) const {
  assert((valid_index(i) && i < irrep_IDs.size()) &&
         "Attempting to set ID for out-of-bounds irrep.");
  irrep_IDs[i] = ID;
  return;
}

//*******************************************************************************************

SymGroupRepID SymGroup::get_irrep_ID(Index i) const {
  if (!valid_index(i) || i >= irrep_IDs.size()) return SymGroupRepID();

  return irrep_IDs[i];
}

//*******************************************************************************************

SymGroupRepID SymGroup::coord_rep_ID() const {
  if (!size() || !at(0).has_valid_master()) {
    err_log() << "CRITICAL ERROR: In SymGroup::get_coord_rep_ID(), SymGroup is "
                 "improperly initialized.\n"
              << "                Exiting...\n";
    exit(1);
  }

  return at(0).master_group().coord_rep_ID();
}
//*******************************************************************************************

SymGroupRepID SymGroup::allocate_representation() const {
  if (!size() || !at(0).has_valid_master()) {
    err_log() << "CRITICAL ERROR: In SymGroup::allocate_representation(), "
                 "SymGroup is improperly initialized.\n"
              << "                Exiting...\n";
    exit(1);
  }

  return at(0).master_group().allocate_representation();
}

//*******************************************************************************************

SymGroupRep const &SymGroup::get_irrep(Index i) const {
  if (!size() || !valid_index(i) || i >= irrep_IDs.size())
    throw std::runtime_error(std::string("Cannot find irrep ") +
                             std::to_string(i) + " in the current SymGroup\n");

  return at(0).master_group().representation(irrep_IDs[i]);
}

//*******************************************************************************************
// The set of left cosets is identical to the equivalence_map formed by
// partitioning (*this) w.r.t. 'subgroup'
std::vector<std::vector<Index>> SymGroup::left_cosets(
    const std::vector<SymOp> &subgroup, double tol) const {
  std::vector<Index> sg_inds = find_all_periodic(subgroup, tol);
  return left_cosets(sg_inds.begin(), sg_inds.end());
}

//*******************************************************************************************

const std::vector<std::vector<Index>> &SymGroup::get_multi_table() const {
  if (multi_table.size() != size()) {
    // std::cout << "CALCULATING MULTI_TABLE for " << this <<  ": table size is
    // " << multi_table.size() << " and group size is " << size() << "!!\n";
    _generate_multi_table();
  }
  return multi_table;
}

//*******************************************************************************************

const std::vector<std::vector<Index>> &SymGroup::get_alt_multi_table() const {
  if (alt_multi_table.size() != size()) {
    // std::cout << "CALCULATING ALT_MULTI_TABLE " << this << ": table size is "
    // << alt_multi_table.size() << " and group size is " << size() << "!!\n";
    _generate_alt_multi_table();
  }
  return alt_multi_table;
}
//*******************************************************************************************

void SymGroup::invalidate_multi_tables() const {
  multi_table.resize(size(), std::vector<Index>(size(), -1));
  alt_multi_table.resize(size(), std::vector<Index>(size(), -1));
}

//*******************************************************************************************

const std::vector<std::vector<Index>> &SymGroup::get_conjugacy_classes() const {
  if (conjugacy_classes.size() != size()) _generate_conjugacy_classes();
  return conjugacy_classes;
}

//*******************************************************************************************

const std::string &SymGroup::get_name() const {
  if (!name.size()) {
    auto info = point_group_info(*this);
    name = info["name"];
    latex_name = info["latex_name"];
    if (!name.size()) {
      err_log()
          << "WARNING: In SymGroup::get_name(), unable to get symgroup type.\n";
      err_log() << "group size is " << size() << '\n';
      name = "unknown";
      latex_name = "unknown";
    }
  }

  return name;
}

//*******************************************************************************************

const std::string &SymGroup::get_latex_name() const {
  get_name();

  return latex_name;
}

//*******************************************************************************************

const std::vector<std::set<std::set<Index>>> &SymGroup::subgroups() const {
  if (!m_subgroups.size()) _generate_subgroups();
  return m_subgroups;
}

//*******************************************************************************************

bool SymGroup::_generate_multi_table() const {  // AAB
  Index i, j;
  multi_table.resize(size(), std::vector<Index>(size(), -1));

  for (i = 0; i < size(); i++) {
    for (j = 0; j < size(); j++) {
      multi_table[i][j] = find_periodic(at(i) * at(j));
      if (multi_table[i][j] >= size() ||
          find_index(multi_table[i], multi_table[i][j]) != j) {
        // this is a hack (sort of). If find_periodic doesn't work, we try
        // find_no trans, which *should* work. In other words, we are using
        // 'inuition' to determine that user doesn't really care about the
        // translational aspects. If our intuition is wrong, there will will
        // probably be an obvious failure later.
        multi_table[i][j] = find_no_trans(at(i) * at(j));
        if (multi_table[i][j] >= size() ||
            find_index(multi_table[i], multi_table[i][j]) != j) {
          // if(multi_table[i][j] >= size()) {
          // std::cout << "This SymGroup is not a group because the combination
          // of at least two of its elements is not contained in the set.\n";

          // ing a table of all 1's seems to make the most sense. This will
          // prevent weird recursion from happening.
          err_log() << "Failed to construc multiplication table!  Table in "
                       "progress:\n";
          for (Index m = 0; m < multi_table.size(); m++) {
            for (Index n = 0; n < multi_table[m].size(); n++) {
              err_log() << multi_table[m][n] << "   ";
            }
            err_log() << std::endl;
          }
          err_log() << std::endl;
          multi_table.resize(size(), std::vector<Index>(size(), -1));
          // multi_table.clear();
          return false;
        }
      }
    }
  }

  return true;
}

//*******************************************************************************************
void SymGroup::_generate_alt_multi_table() const {
  // by calling get_multi_table(), we ensure that multi_table is populated
  alt_multi_table.resize(get_multi_table().size());

  if (multi_table.size() && !valid_index(multi_table[0][0])) {
    alt_multi_table.resize(size(), std::vector<Index>(size(), -1));
    return;
  }
  for (Index i = 0; i < multi_table.size(); i++) {
    if (multi_table[i][i] != 0) {
      alt_multi_table[find_index(multi_table[i], 0)] = multi_table[i];
    } else {
      alt_multi_table[i] = multi_table[i];
    }
  }
}

//*******************************************************************************************
// Please keep in mind that this will only return the FIRST match that is found.
// This does not guarantee that it is the only match.

Index SymGroup::find_no_trans(const SymOp &test_op) const {
  for (Index i = 0; i < size(); i++) {
    if (almost_equal(at(i).matrix(), test_op.matrix())) {
      return i;
    }
  }

  return size();
}

//*******************************************************************************************

Index SymGroup::find_periodic(const SymOp &test_op, double tol) const {
  for (Index i = 0; i < size(); i++) {
    if (compare_periodic(at(i), test_op, lattice(), periodicity(), tol)) {
      return i;
    }
  }
  return size();
}

//*******************************************************************************************

std::vector<Index> SymGroup::find_all_periodic(
    const std::vector<SymOp> &subgroup, double tol) const {
  std::vector<Index> tarray;
  for (Index i = 0; i < subgroup.size(); i++) {
    tarray.push_back(find_periodic(subgroup[i], tol));
    if (tarray.back() == size()) {
      tarray.clear();
      break;
    }
  }
  return tarray;
}

//*******************************************************************************************
// Probably need to change this...
bool SymGroup::is_group(double tol) const {
  // This is important because it ensures that things like
  // within() and min_dist() work correctly for your group
  // regardless of what kind of group it is.

  for (Index i = 0; i < size(); i++) {
    for (Index j = 0; j < size(); j++) {
      if (!contains_periodic(at(i) * at(j), tol)) return false;
    }
    if (!contains_periodic(at(i).inverse(), tol)) return false;
  }
  return true;
}

//*******************************************************************************************

void SymGroup::enforce_group(double tol, Index max_size) {
  bool new_ops(true);

  while (new_ops && size() < max_size) {
    new_ops = false;
    for (Index i = 0; i < size() && size() < max_size; i++) {
      SymOp A_op(at(i).unregistered_copy());
      for (Index j = 0; j < size() && size() < max_size; j++) {
        SymOp B_op(at(j).unregistered_copy());

        SymOp tOp(A_op * B_op);

        if (!contains_periodic(tOp, tol)) {
          push_back(within_cell(tOp, lattice(), periodicity()));
          new_ops = true;
          //	  //std::cout << "Pushing back a SymOp due to multiplication
          // fail.\n";
        }
      }

      SymOp tOp(A_op.inverse());
      if (!contains_periodic(tOp, tol)) {
        push_back(within_cell(tOp, lattice(), periodicity()));
        new_ops = true;
        // std::cout << "Pushing back a SymOp due to inverse fail.\n";
      }
    }
  }
  if (size() >= max_size - 1) {
    err_log() << "In SymGroup::enforce_group() -- you have reached the maximum "
                 "allowed size you specified for your group (the default is "
                 "200). Unless you are generating a factor group in a large "
                 "supercell, you probably need to adjust your tolerances.\n";
    assert(0);
  }
  return;
}

//*******************************************************************************************

bool SymGroup::contains_periodic(const SymOp &test_op, double tol) const {
  return find_periodic(test_op, tol) != size();
}

//*******************************************************************************************

SymGroup &SymGroup::apply_sym(const SymOp &op) {
  for (Index i = 0; i < size(); i++) at(i).apply_sym(op);
  return *this;
}

//*******************************************************************************************

void SymGroup::print(std::ostream &out, COORD_TYPE mode) const {
  out << size() << " # " << xtal::COORD_MODE::NAME(mode)
      << " representation of group containing " << size() << " elements:\n\n";
  Eigen::Matrix3d c2f_mat(Eigen::Matrix3d::Identity());
  if (mode == FRAC) c2f_mat = lattice().inv_lat_column_mat();
  for (Index i = 0; i < size(); i++) {
    out << i + 1 << "  " << description(at(i), lattice(), mode) << std::endl;
    at(i).print(out, c2f_mat);
    out << std::endl;
  }
  return;
}

//*******************************************************************************************

void SymGroup::calc_space_group_in_cell(SymGroup &space_group_cell,
                                        const Lattice &_cell) const {
  if (!size()) return;

  Eigen::Vector3i max_trans(3, 3, 3);
  xtal::Coordinate trans(Eigen::Vector3d::Zero(), _cell, FRAC);
  space_group_cell.clear();

  std::vector<SymInfo> sg_info;
  for (Index i = 0; i < size(); i++) {
    EigenCounter<Eigen::Vector3i> lat_comb(-max_trans, max_trans,
                                           Eigen::Vector3i::Ones());
    do {
      trans.frac() = lat_comb().cast<double>();
      SymOp new_sym(SymOp::translation(trans.cart()) * at(i));
      SymInfo info(new_sym, lattice());
      trans = info.location;
      if (!trans.is_within()) {
        continue;
      }

      bool new_location = true;
      for (Index j = 0; j < space_group_cell.size(); j++) {
        if (almost_equal(new_sym.matrix(), space_group_cell[j].matrix()) &&
            almost_equal(info.location.const_cart(),
                         sg_info[j].location.const_cart())) {
          new_location = false;
          break;
        }
      }
      if (new_location) {
        space_group_cell.push_back(new_sym);
        sg_info.push_back(info);
      }
    } while (++lat_comb);
  }

  return;
}

//*******************************************************************************************

void SymGroup::calc_space_group_in_range(SymGroup &space_group,
                                         const Lattice &_cell,
                                         Eigen::Vector3i min_trans,
                                         Eigen::Vector3i max_trans) const {
  if (!size()) return;

  xtal::Coordinate trans(Eigen::Vector3d::Zero(), _cell, FRAC);

  for (Index i = 0; i < size(); i++) {
    EigenCounter<Eigen::Vector3i> lat_comb(min_trans, max_trans,
                                           Eigen::Vector3i::Ones());
    do {
      trans.frac() = lat_comb().cast<double>();

      SymOp new_sym(SymOp::translation(trans.cart()) * at(i));

      if (!contains(space_group, new_sym)) {
        space_group.push_back(new_sym);
      }

    } while (++lat_comb);
  }

  return;
}

//***************************************************
void SymGroup::print_locations(std::ostream &stream) const {
  // Assumes SymGroup is sorted with clumps of SymOps of common matrix type and
  // eigenvec sort();

  bool new_op = true;
  stream << "Locations for symmetry operations\n";
  SymInfo info(at(0), lattice());
  SymInfo next_info(info);
  Eigen::Matrix3d c2f_mat = lattice().inv_lat_column_mat();
  for (Index i = 0; i < size(); i++) {
    if (new_op) {
      at(i).print(stream, c2f_mat);

      stream << std::endl;
      at(i).print(stream, Eigen::Matrix3d::Identity());
      stream << std::endl;

      stream << "Location:" << std::endl;
      stream << "FRAC\t\t\t\t\tCART" << std::endl;
    }
    stream << info.location.const_frac();
    stream << "\t\t\t";
    stream << info.location.const_cart();
    stream << std::endl;

    if (i + 1 < size()) {
      next_info = SymInfo(at(i + 1), lattice());
      if (info.op_type == next_info.op_type &&
          almost_equal(info.axis.const_cart(), next_info.axis.const_cart())) {
        // Is this enough to know if it's a new symmetry or not?
        new_op = false;
      } else {
        new_op = true;
        stream << "------------------------------------------------------------"
                  "----------------------\n\n";
      }
    }
    info = next_info;
  }
  return;
}

//*******************************************************************************************

/// \brief Sort SymOp in the SymGroup
///
/// - If multiplication table can be generated:
///   - Generate conjugacy classes
///   - Sort SymOp in each conjugacy class
///   - Sort each conjugacy class by the first SymOp in the class
/// - Else:
///   - Sort all SymOp
///
/// SymOp are sorted by lexicographical comparison of: (-det, -trace, angle,
/// axis, tau)
/// - angle is positive
/// - axis[0] is positive
///
void SymGroup::sort() {
  // floating point comparison tolerance
  double tol = TOL;

  // COORD_TYPE print_mode = CART;

  // compare on vector of '-det', '-trace', 'angle', 'axis', 'tau'
  typedef Eigen::Matrix<double, 10, 1> key_type;
  auto make_key = [](const SymOp &op, const Lattice &lat) {
    key_type vec;
    int offset = 0;

    SymInfo info(op, lat);
    vec[offset] = double(op.time_reversal());
    offset++;

    vec[offset] = -op.matrix().determinant();
    offset++;

    vec[offset] = -op.matrix().trace();
    offset++;

    vec[offset] = info.angle;
    offset++;

    vec.segment<3>(offset) = info.axis.const_frac();
    offset += 3;

    vec.segment<3>(offset) = xtal::Coordinate(op.tau(), lat, CART).const_frac();
    offset += 3;

    return vec;
  };

  // define symop compare function
  auto op_compare = [tol](const key_type &A, const key_type &B) {
    return float_lexicographical_compare(A, B, tol);
  };

  typedef std::map<key_type, SymOp,
                   std::reference_wrapper<decltype(op_compare)>>
      map_type;

  // sort conjugacy class using the first symop in the sorted class
  auto cclass_compare = [tol](const map_type &A, const map_type &B) {
    return float_lexicographical_compare(A.begin()->first, B.begin()->first,
                                         tol);
  };

  // sort elements in each conjugracy class (or just put all elements in the
  // first map)
  std::set<map_type, std::reference_wrapper<decltype(cclass_compare)>> sorter(
      cclass_compare);

  // first put identity in position 0 in order to calculat multi_table correctly
  for (int i = 0; i < size(); ++i) {
    if (at(i).is_identity()) {
      std::swap(at(0), at(i));
      break;
    }
  }

  // if sorting by congujacy classes
  if (get_multi_table().size() == size()) {
    // get conjugacy classes
    get_conjugacy_classes().size();

    // insert elements into each conjugacy class to sort the class
    for (int i = 0; i < conjugacy_classes.size(); ++i) {
      map_type cclass(op_compare);
      for (int j = 0; j < conjugacy_classes[i].size(); ++j) {
        const SymOp &op = at(conjugacy_classes[i][j]);
        cclass.insert(std::make_pair(make_key(op, lattice()), op));
      }

      sorter.emplace(std::move(cclass));
    }
  } else {
    // else just sort element
    map_type all_op(op_compare);
    for (auto it = begin(); it != end(); ++it) {
      all_op.emplace(make_key(*it, lattice()), *it);
    }
    sorter.emplace(std::move(all_op));
  }

  // copy symop back into group
  int j = 0;
  for (auto const &cclass : sorter) {
    for (auto it = cclass.begin(); it != cclass.end(); ++it) {
      at(j) = it->second;
      ++j;
    }
  }

  clear_tables();
}

//*******************************************************************************************

bool SymGroup::is_irreducible() const {
  double tvalue = 0;

  for (Index i = 0; i < size(); i++) {
    tvalue += (at(i).matrix().trace()) * (at(i).matrix().trace());
  }

  if (Index(std::abs(tvalue)) == size()) {
    return true;
  } else {
    return false;
  }
}

//*******************************************************************************************
/// Translation operators for origin shift need to be defined
SymGroup &SymGroup::operator+=(
    const Eigen::Ref<const SymOp::vector_type> &shift) {
  for (Index ng = 0; ng < size(); ng++) at(ng) += shift;
  return (*this);
}

//*******************************************************************************************

SymGroup &SymGroup::operator-=(
    const Eigen::Ref<const SymOp::vector_type> &shift) {
  for (Index ng = 0; ng < size(); ng++) at(ng) -= shift;
  return (*this);
}

//*******************************************************************************************

MasterSymGroup make_master_sym_group(SymGroup const &_group,
                                     Lattice const &_lattice) {
  MasterSymGroup result;
  result.set_lattice(_lattice);
  for (auto const &op : _group) {
    result.push_back(op);
  }
  return result;
}

//**********************************************************

bool compare_periodic(const SymOp &a, const SymOp &b, const Lattice &lat,
                      PERIODICITY_TYPE periodicity, double _tol) {
  if (a.time_reversal() != b.time_reversal() ||
      !almost_equal(a.matrix(), b.matrix(), _tol))
    return false;
  // std::cout << "Operations:\n"
  //<< a.matrix() << "\n and \n"
  //<< b.matrix() << "\n are equal \n"
  //<< "a-tau " << a.tau().transpose() << std::endl
  //<< "b-tau " << b.tau().transpose() << std::endl;
  if (periodicity != PERIODIC) return almost_equal(a.tau(), b.tau(), _tol);

  return xtal::Coordinate(a.tau(), lat, CART)
             .min_dist(xtal::Coordinate(b.tau(), lat, CART)) < _tol;
}

//**********************************************************

SymOp within_cell(const SymOp &a, const Lattice &lat,
                  PERIODICITY_TYPE periodicity) {
  if (periodicity != PERIODIC) return a;

  xtal::Coordinate trans(a.tau(), lat, CART);
  trans.within();
  return SymOp(a.matrix(), trans.cart(), a.time_reversal(), a.map_error());
}

}  // namespace CASM
