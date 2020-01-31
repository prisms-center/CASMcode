#include "casm/symmetry/SymGroup.hh"
#include <boost/filesystem/fstream.hpp>
#include "casm/external/Eigen/CASM_AddOns"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/casm_io/Log.hh"

#include "casm/casm_io/container/json_io.hh"

namespace CASM {
  //INITIALIZE STATIC MEMBER MasterSymGroup::GROUP_COUNT
  //THIS MUST OCCUR IN A .CC FILE; MAY CAUSE PROBLEMS IF WE
  //CHANGE COMPILING/LINKING STRATEGY
  Index MasterSymGroup::GROUP_COUNT(0);

  //*******************************************************************************************

  MasterSymGroup::MasterSymGroup(const MasterSymGroup &RHS) :
    SymGroup(RHS),
    m_group_index(RHS.m_group_index),
    m_coord_rep_ID(RHS.m_coord_rep_ID),
    m_reg_rep_ID(RHS.m_reg_rep_ID),
    m_identity_rep_IDs(RHS.m_identity_rep_IDs) {

    m_rep_array.reserve(RHS.m_rep_array.size());

    for(Index i = 0; i < RHS.m_rep_array.size(); i++) {
      _add_representation(RHS.m_rep_array[i]->copy());
    }

    for(Index i = 0; i < size(); i++)
      at(i).set_index(*this, i);
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
    for(Index i = 0; i < RHS.m_rep_array.size(); i++)
      _add_representation(RHS.m_rep_array[i]->copy());

    for(Index i = 0; i < size(); i++)
      at(i).set_index(*this, i);

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
    if(!m_point_group.size()) {
      m_point_group = copy_no_trans(false);
    }
    return m_point_group;
  }

  //*******************************************************************************************

  void MasterSymGroup::clear() {
    SymGroup :: clear();
    m_point_group.clear();
    for(Index i = 0; i < m_rep_array.size(); i++) {
      if(m_rep_array[i])
        delete m_rep_array[i];
    }
    m_rep_array.clear();

    m_reg_rep_ID = m_coord_rep_ID = SymGroupRepID();

    m_identity_rep_IDs.clear();
    return;
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::coord_rep_ID() const {
    if(m_coord_rep_ID.empty())
      _add_coord_rep();
    return m_coord_rep_ID;
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::reg_rep_ID() const {
    if(m_reg_rep_ID.empty())
      _add_reg_rep();
    return m_reg_rep_ID;
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::identity_rep_ID(Index dim) const {
    if(m_identity_rep_IDs.size() < dim + 1) {
      auto tail = std::vector<SymGroupRepID>(dim + 1 - m_identity_rep_IDs.size());
      m_identity_rep_IDs.insert(m_identity_rep_IDs.end(), tail.begin(), tail.end());
    }
    if(m_identity_rep_IDs[dim].empty()) {
      m_rep_array.push_back(new SymGroupRep(*this));
      for(Index i = 0; i < size(); i++) {
        m_rep_array.back()->set_rep(i, SymPermutation(Permutation(dim)));
      }
      m_identity_rep_IDs[dim] = SymGroupRepID(group_index(), m_rep_array.size() - 1);
    }
    return m_identity_rep_IDs[dim];
  }

  //*******************************************************************************************

  SymGroupRep const &MasterSymGroup::coord_rep() const {
    if(m_coord_rep_ID.empty())
      _add_coord_rep();
    return representation(m_coord_rep_ID);
  }

  //*******************************************************************************************

  void MasterSymGroup::set_rep(SymGroupRepID rep_ID, SymOpRepresentation const &_op_rep, Index op_index) const {
    _representation_ptr(rep_ID)->set_rep(op_index, _op_rep);
  }

  //*******************************************************************************************

  SymGroupRep const &MasterSymGroup::reg_rep() const {
    if(m_reg_rep_ID.empty())
      _add_reg_rep();
    return representation(m_reg_rep_ID);
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::_add_coord_rep() const {
    SymGroupRep *coordrep(new SymGroupRep(*this));
    for(Index i = 0; i < size(); i++)
      coordrep->set_rep(i, SymMatrixXd(at(i).matrix()));

    m_coord_rep_ID = _add_representation(coordrep);
    return m_coord_rep_ID;
  }

  //*******************************************************************************************
  SymGroupRepID MasterSymGroup::_add_reg_rep() const {
    SymGroupRep *regrep(new SymGroupRep(*this));
    Eigen::MatrixXd regrep_mat(size(), size());

    if(get_alt_multi_table().size() != size())
      return SymGroupRepID();


    for(Index i = 0; i < size(); i++) {
      for(Index j = 0; j < size(); j++) {
        for(Index k = 0; k < size(); k++) {
          if(alt_multi_table[j][k] == i) {
            regrep_mat(j, k) = 1.0;
          }
          else {
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

  SymGroupRepID MasterSymGroup::add_kronecker_rep(SymGroupRepID ID1, SymGroupRepID ID2) const {
    SymGroupRep const *rep1(_representation_ptr(ID1)), *rep2(_representation_ptr(ID2));
    if(!(rep1 && rep2))
      return SymGroupRepID();

    SymGroupRep *new_rep(new SymGroupRep(*this));
    Eigen::MatrixXd tmat;
    for(Index i = 0; i < size(); i++) {
      //std::cout << "rep1 matrix:\n";
      //std::cout << *(rep1->MatrixXd(i)) << '\n';
      //std::cout << "rep2 matrix:\n";
      //std::cout << *(rep2->MatrixXd(i)) << '\n';

      kroneckerProduct(*(rep1->MatrixXd(i)), *(rep2->MatrixXd(i)), tmat);
      //std::cout << "Total matrix:\n" << tmat << '\n';
      new_rep->set_rep(i, SymMatrixXd(tmat));
    }
    return _add_representation(new_rep);
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::add_direct_sum_rep(const std::vector<SymGroupRepID> &rep_IDs) const {
    std::vector<SymGroupRep const *> treps;
    for(Index i = 0; i < rep_IDs.size(); i++) {
      treps.push_back(_representation_ptr(rep_IDs[i]));
      if(!treps.back())
        return SymGroupRepID();
    }
    SymGroupRep *new_rep(new SymGroupRep(*this));

    int dim = 0;
    for(Index i = 0; i < treps.size(); i++) {
      if(treps[i]->size() != size())
        return SymGroupRepID();
      if(!(treps[i]->MatrixXd(0)))
        return SymGroupRepID();

      dim += (treps[i]->MatrixXd(0))->cols();
    }

    Eigen::MatrixXd tmat(Eigen::MatrixXd::Zero(dim, dim));
    int corner = 0;
    for(Index i = 0; i < size(); i++) {
      corner = 0;
      for(Index j = 0; j < treps.size(); j++) {
        tmat.block(corner, corner, (treps[j]->MatrixXd(i))->cols(), (treps[j]->MatrixXd(i))->cols()) = *(treps[j]->MatrixXd(i));
        corner += (treps[j]->MatrixXd(i))->cols();
      }
      new_rep->set_rep(i, SymMatrixXd(tmat));
    }
    return _add_representation(new_rep);
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::add_rotation_rep() const {

    // make sure coord_rep_ID exists?
    //Index coord_rep_index((*this).coord_rep_ID());
    // make a new symmetry representation
    SymGroupRep *new_rep(new SymGroupRep(*this));

    // we are going to use our knowledge of how rotation
    // matrices transform under symmetry (in the coord_rep)
    // in order to build
    // up a representation in the rotation basis.
    for(Index i = 0; i < size(); i++) {

      // build each 4d representaion for a = 1, b = 1, ...
      // rp = rotation parameter
      Eigen::Matrix3d coord_rep_mat = at(i).matrix();
      if(!almost_equal(coord_rep_mat.determinant(), 1.)) {
        //std::cout << "IN IF" << std::endl;
        coord_rep_mat *= coord_rep_mat.determinant();
        //continue;
      }

      Eigen::Quaterniond opq(coord_rep_mat);

      //std::cout << "Quaternion " << opq << std::endl;

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

      right_q(0, 0) =  opq.w();
      right_q(0, 1) = -opq.x();
      right_q(0, 2) = -opq.y();
      right_q(0, 3) = -opq.z();

      right_q(1, 0) =  opq.x();
      right_q(1, 1) =  opq.w();
      right_q(1, 2) =  opq.z();
      right_q(1, 3) = -opq.y();

      right_q(2, 0) =  opq.y();
      right_q(2, 1) = -opq.z();
      right_q(2, 2) =  opq.w();
      right_q(2, 3) =  opq.x();

      right_q(3, 0) =  opq.z();
      right_q(3, 1) =  opq.y();
      right_q(3, 2) = -opq.x();
      right_q(3, 3) =  opq.w();

      left_q(0, 0) =  opq.w();
      left_q(0, 1) = -opq.x();
      left_q(0, 2) = -opq.y();
      left_q(0, 3) = -opq.z();

      left_q(1, 0) =  opq.x();
      left_q(1, 1) =  opq.w();
      left_q(1, 2) = -opq.z();
      left_q(1, 3) =  opq.y();

      left_q(2, 0) =  opq.y();
      left_q(2, 1) =  opq.z();
      left_q(2, 2) =  opq.w();
      left_q(2, 3) = -opq.x();

      left_q(3, 0) =  opq.z();
      left_q(3, 1) = -opq.y();
      left_q(3, 2) =  opq.x();
      left_q(3, 3) =  opq.w();

      //std::cout << right_q << std::endl;

      new_rep->set_rep(i, SymMatrixXd(left_q));
    }

    for(Index i = 0 ; i < new_rep->size() ; i++) {
      //std::cout << "Rota Rep final mats " << i << "\n" << *(new_rep->MatrixXd(i)) << std::endl;
    }

    return _add_representation(new_rep);
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::add_transformed_rep(SymGroupRepID orig_ID, const Eigen::MatrixXd &trans_mat) const {
    SymGroupRep const *trep(_representation_ptr(orig_ID));
    if(!trep)
      return SymGroupRepID();
    return add_representation(coord_transformed_copy(*trep, trans_mat));
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::allocate_representation() const {
    SymGroupRepID new_ID(group_index(), m_rep_array.size());
    m_rep_array.push_back(new SymGroupRep(*this, new_ID));
    return new_ID;
  }

  //*******************************************************************************************

  SymGroupRepID MasterSymGroup::add_representation(const SymGroupRep &new_rep) const {
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
    if(_id.is_identity()) {
      // _id.rep_index() stores dimension of representation
      _id = identity_rep_ID(_id.rep_index());
    }
    else if(_id.group_index() != group_index()) {
      throw std::runtime_error("Attempting to access representation from MasterGroup #" + std::to_string(group_index())
                               + " that resides in MasterGroup #" + std::to_string(_id.group_index()));
    }

    return m_rep_array[_id.rep_index()];
  }

  //*******************************************************************************************

  void MasterSymGroup::sort() {
    SymGroup::sort();
    m_point_group.clear();
    bool broken_check(false);
    std::vector<Index> perm_array(size(), 0);
    for(Index i = 0; i < size(); i++) {
      perm_array[i] = at(i).index();
      if(at(i).index() != i) {
        at(i).set_index(*this, i);
        broken_check = true;
      }
    }
    if(broken_check && m_rep_array.size()) {
      //default_err_log() << "WARNING: Order of symmetry operations has been altered by MasterSymGroup::sort_by_class(). Attempting to repair "
      //                << m_rep_array.size() << " symmetry representations.\n";
      for(Index i = 0; i < m_rep_array.size(); i++) {
        std::vector<SymOpRepresentation *> new_rep;
        std::transform(perm_array.begin(), perm_array.end(), std::back_inserter(new_rep), [&](Index & idx) {
          return (*m_rep_array[i])[idx];
        });
        m_rep_array[i]->swap(new_rep);
        for(Index j = 0; j < m_rep_array[i]->size(); j++) {
          if(m_rep_array[i]->at(j)) {
            (m_rep_array[i]->at(j))->set_identifiers(*this, m_rep_array[i]->symrep_ID(), j);//perm_array[j]);
          }
        }
      }
    }

    return;
  }

  //*******************************************************************************************
  SymGroup SymGroup::lattice_point_group(Lattice const &_lat) {

    xtal::SymOpVector lattice_point_group_operations = xtal::make_point_group(_lat);
    SymGroup point_group = adapter::Adapter<SymGroup, xtal::SymOpVector>()(lattice_point_group_operations, _lat);


    if(!point_group.is_group(_lat.tol())) {
      std::cerr << "*** WARNING *** \n"
                << "    SymGroup::lattice_point_group() has been called on an ill-conditioned lattice \n"
                << "    (i.e., a well-defined point group could not be found with the current tolerance of " << _lat.tol() << ").\n"
                << "    CASM will use the group closure of the symmetry operations that were found.  Please consider using the \n"
                << "    CASM symmetrization tool on your input files.\n";
      std::cout << "lat_column_mat:\n" << _lat.lat_column_mat() << "\n\n";

      point_group.enforce_group(_lat.tol());

    }
    //Sort point_group by trace/conjugacy class
    point_group.sort();

    return point_group;
  }

  //*******************************************************************************************
  SymGroup::SymGroup(std::vector<SymOp> from_array, Lattice const *_lat_ptr, PERIODICITY_TYPE init_type) :
    m_lat_ptr(_lat_ptr),
    m_group_periodicity(init_type),
    m_max_error(-1) {
    std::vector<SymOp>::swap(from_array);

  }

  //*******************************************************************************************

  void SymGroup::set_lattice(const Lattice &lat) {
    m_lat_ptr = &lat;
  }

  /*****************************************************************/

  void SymGroup::push_back(const SymOp &new_op) {
    std::vector<SymOp>::push_back(new_op);
    if(back().map_error() > m_max_error)
      m_max_error = back().map_error();
    return;
  }

  //*******************************************************************************************

  void SymGroup::clear_tables() {
    multi_table.clear();
    alt_multi_table.clear();
    conjugacy_classes.clear();
    class_names.clear();
    index2conjugacy_class.clear();
    m_character_table.clear();
    complex_irrep.clear();
    irrep_IDs.clear();
    irrep_names.clear();

    m_subgroups.clear();

    centralizer_table.clear();
    elem_order_table.clear();

    name.clear();
    latex_name.clear();
    comment.clear();

    return;
  }

  //*****************************************************************
  std::vector<Index> SymGroup::op_indices() const {
    std::vector<Index> ind_array(size());
    for(Index i = 0; i < size(); i++) {
      if(!valid_index(at(i).index()) && at(i).has_valid_master()) {
        ind_array[i] = at(i).master_group().find_periodic(at(i));
      }
      else {
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
  // Determination of point group type
  // By counting the number of unique type of rotation symmetry and
  // total number of rotation symmetry elements, the point group type
  // of a space group is determined.

  //
  // Donghee
  //*****************************************************************
  std::vector<Index> SymGroup::get_rotation_groups() const {
    std::vector<Index> result(12, 0);
    for(Index i = 0; i < size(); i++) {
      SymInfo info(at(i), lattice());
      if(info.op_type == symmetry_type::identity_op) {
        result[0]++;
      }
      else if(info.op_type == symmetry_type::inversion_op) {
        result[1]++;
      }
      else {
        Index bit = 0;
        if(at(i).matrix().determinant() < 0)
          bit = 1;
        for(Index n = 2; n <= 6; ++n) {
          if(almost_equal<double>(360. / n, info.angle)) {
            result[2 * (n - 1) + bit]++;
            break;
          }

        }
      }
    }
    return result;
  }
  //*****************************************************************
  std::map<std::string, std::string> SymGroup::point_group_info() const {
    std::vector<Index> rgroups = get_rotation_groups();
    std::map<std::string, std::string> result;
    // Calculate total number of rotation elements
    Index nm = rgroups[2] + 1;
    for(Index k = 4; k < rgroups.size(); k++) {
      nm += rgroups[k];
    }
    bool centric = true;
    result["centricity"] = "Centric";
    if(!rgroups[1]) {
      centric = false;
      result["centricity"] = "Acentric";
    }
    // naming point gorup
    if((rgroups[4] + rgroups[5]) == 8) {
      result["crystal_system"] = "Cubic";
      switch(nm) {
      case 12:
        if(centric) {
          result["international_name"] = "m-3";
          result["name"] = "Th";
          result["latex_name"] = "T_h";
          result["space_group_range"] = "200-206";
        }
        else {
          result["international_name"] = "23";
          result["name"] = "T (Chiral)";
          result["latex_name"] = "T";
          result["space_group_range"] = "195-199";
        }
        break;
      case 24:
        if(centric) {
          result["international_name"] = "m-3m";
          result["name"] = "Oh";
          result["latex_name"] = "O_h";
          result["space_group_range"] = "221-230";
        }
        else {
          if(rgroups[6] == 6) {
            result["international_name"] = "432";
            result["name"] = "O (Chiral)";
            result["latex_name"] = "O";
            result["space_group_range"] = "207-214";
          }
          else if(rgroups[7] == 6) {
            result["international_name"] = "-43m";
            result["name"] = "Td";
            result["latex_name"] = "T_d";
            result["space_group_range"] = "215-220";
          }
          else {
            default_err_log() << "\n Error Cubic case 24 Acentric \n ";
          }
        }
        break;
      default:
        default_err_log() << "\n Error Cubic \n";
      }
    } // for cubic;
    else if((rgroups[10] + rgroups[11]) == 2) {
      result["crystal_system"] = "Hexagonal";
      switch(nm) {
      case 6:
        if(centric) {
          result["international_name"] = "6/m";
          result["name"] = "C6h";
          result["latex_name"] = "C_{6h}";
          result["space_group_range"] = "175-176";
        }
        else {
          if(rgroups[10] == 2) {
            result["international_name"] = "6";
            result["name"] = "C6 (Chiral)";
            result["latex_name"] = "C_6";
            result["space_group_range"] = "168-173";
          }
          else if(rgroups[11] == 2) {
            result["international_name"] = "-6";
            result["name"] = "C3h";
            result["latex_name"] = "C_{3h}";
            result["space_group_range"] = "174";
          }
          else {
            default_err_log() << "\n Error Hexagonal case 6 Acentric \n ";
          }
        }
        break;
      case 12:
        if(centric) {
          result["international_name"] = "6/mmm";
          result["name"] = "D6h";
          result["latex_name"] = "D_{6h}";
          result["space_group_range"] = "191-194";
        }
        else {
          if(rgroups[10] == 2) {
            if(rgroups[3] == 7) {
              result["international_name"] = "622";
              result["name"] = "D6 (Chiral)";
              result["latex_name"] = "D_{6h}";
              result["space_group_range"] = "177-182";
            }
            else if(rgroups[4] == 6) {
              result["international_name"] = "6mm";
              result["name"] = "C6v";
              result["latex_name"] = "C_{6v}";
              result["space_group_range"] = "183-186";
            }
            else default_err_log() << "\n Error Hexagonal case 12 Ancentric #6 \n";
          }
          else if(rgroups[11] == 2) {
            result["international_name"] = "-6m2";
            result["name"] = "D3h";
            result["latex_name"] = "D_{3h}";
            result["space_group_range"] = "187-190";
          }
          else {
            default_err_log() << "\n Error Hexagonal case 12 Acentric \n ";
          }
        }
        break;
      default:
        default_err_log() << "\n Error Hexagonal \n";
      }
    } // for hexagonal
    else if((rgroups[4] + rgroups[5]) == 2) {
      result["crystal_system"] = "Trigonal";
      switch(nm) {
      case 3:
        if(centric) {
          result["international_name"] = "-3";
          result["name"] = "S6";
          result["latex_name"] = "S_6";
          result["space_group_range"] = "147-148";
        }
        else {
          result["international_name"] = "3";
          result["name"] = "C3 (Chiral)";
          result["latex_name"] = "C_3";
          result["space_group_range"] = "143-146";
        }
        break;
      case 6:
        if(centric) {
          result["international_name"] = "-3m";
          result["name"] = "D3d";
          result["latex_name"] = "D_{3d}";
          result["space_group_range"] = "162-167";
        }
        else {
          if(rgroups[2] == 3) {
            result["international_name"] = "32";
            result["name"] = "D3 (Chiral)";
            result["latex_name"] = "D_3";
            result["space_group_range"] = "149-155";
          }
          else if(rgroups[3] == 3) {
            result["international_name"] = "3m";
            result["name"] = "C3v";
            result["latex_name"] = "C_{3v}";
            result["space_group_range"] = "156-161";
          }
          else {
            default_err_log() << "\n Error Trigonal case 6 Acentric \n ";
          }
        }
        break;
      default:
        default_err_log() << "\n Error Trigonal \n";
      }
    } // for trigonal
    else if((rgroups[6] + rgroups[7]) == 2) {
      result["crystal_system"] = "Tetragonal";
      switch(nm) {
      case 4:
        if(centric) {
          result["international_name"] = "4/m";
          result["name"] = "C4h";
          result["latex_name"] = "C_{4h}";
          result["space_group_range"] = "83-88";
        }
        else {
          if(rgroups[6] == 2) {
            result["international_name"] = "4";
            result["name"] = "C4 (Chiral)";
            result["latex_name"] = "C_4";
            result["space_group_range"] = "75-80";
          }
          else if(rgroups[7] == 2) {
            result["international_name"] = "-4";
            result["name"] = "S4";
            result["latex_name"] = "S_4";
            result["space_group_range"] = "81-82";
          }
          else {
            default_err_log() << "\n Error Tetragonal case 4 Acentric \n ";
          }
        }

        break;
      case 8:
        if(centric) {
          result["international_name"] = "4/mmm";
          result["name"] = "D4h";
          result["latex_name"] = "D_{4h}";
          result["space_group_range"] = "123-142";
        }
        else {
          if(rgroups[6] == 2) {
            if(rgroups[3] == 5) {
              result["international_name"] = "422";
              result["name"] = "D4 (Chiral)";
              result["latex_name"] = "D_4";
              result["space_group_range"] = "89-98";
            }
            else if(rgroups[4] == 4) {
              result["international_name"] = "4mm";
              result["name"] = "C4v";
              result["latex_name"] = "C_{4v}";
              result["space_group_range"] = "99-110";
            }
            else default_err_log() << "\n Error Tetragonal case 8 Ancentric #4 \n";
          }
          else if(rgroups[7] == 2) {
            result["international_name"] = "-42m";
            result["name"] = "D2d";
            result["latex_name"] = "D_{2d}";
            result["space_group_range"] = "111-122";
          }
          else {
            default_err_log() << "\n Error Tetragonal case 8 Acentric \n ";
          }
        }
        break;
      default:
        default_err_log() << "\n Error Tetragonal \n";
      }
    } // for tetragonal
    else if((rgroups[2] + rgroups[3]) == 3 || ((rgroups[2] + rgroups[3]) == 6  && nm == 4)) {
      result["crystal_system"] = "Orthorhomic";
      if(centric) {
        result["international_name"] = "mmm";
        result["name"] = "D2h";
        result["latex_name"] = "D_{2h}";
        result["space_group_range"] = "47-84";
      }
      else {
        if(rgroups[3] == 3) {
          result["international_name"] = "222";
          result["name"] = "D2 (Chiral)";
          result["latex_name"] = "D_2";
          result["space_group_range"] = "16-24";
        }
        else if(rgroups[4] == 2) {
          result["international_name"] = "mm2";
          result["name"] = "C2v";
          result["latex_name"] = "C_{2v}";
          result["space_group_range"] = "25-46";
        }
        else {
          default_err_log() << "\n Error Orthorhombic Acentric \n ";
        }
      }
    } // for orthorhombic

    else if((rgroups[2] + rgroups[3]) == 1 || nm == 2) {
      result["crystal_system"] = "Monoclinic";
      if(centric) {
        result["international_name"] = "2/m";
        result["name"] = "C2h";
        result["latex_name"] = "C_{2h}";
        result["space_group_range"] = "10-15";
      }
      else {
        if(rgroups[3] == 1) {
          result["international_name"] = "2";
          result["name"] = "C2 (Chiral)";
          result["latex_name"] = "C_2";
          result["space_group_range"] = "3-5";
        }
        else if(rgroups[4] == 1) {
          result["international_name"] = "m";
          result["name"] = "Cs";
          result["latex_name"] = "C_s";
          result["space_group_range"] = "6-9";
        }
        else {
          default_err_log() << "\n Error Monoclinic Acentric \n ";
        }
      }
    } // for Acentric monoclinic
    else if(nm == 1) {
      result["crystal_system"] = "Triclinic";
      if(centric) {
        result["international_name"] = "-1";
        result["name"] = "Ci";
        result["latex_name"] = "C_i";
        result["space_group_range"] = "2";
      }
      else {
        result["international_name"] = "1";
        result["name"] = "C1 (Chiral)";
        result["latex_name"] = "C_1";
        result["space_group_range"] = "1";
      }
    }
    else {
      default_err_log() << "Error foind point group type \n";
    }
    return result;
  }

  //*******************************************************************************************

  void SymGroup::print_space_group_info(std::ostream &out) const {

    auto result = point_group_info();
    out << "  Crystal System : " << result["crystal_system"] << "\n";
    out << "  Space Group # : " << result["space_group_range"] << "\n";
    out << "  Point Group is " << result["centricity"] << " " << result["international_name"] << " ( " << result["name"] << " ) \n";
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
    for(Index i = 0; i < size(); i++) {
      SymOp tsymop(at(i).no_trans());
      if(keep_repeated || !contains(result, tsymop)) {
        result.push_back(tsymop);
      }
    }
    return result;
  }

  //*******************************************************************************************
  SymInfo SymGroup::info(Index i) const {
    return SymInfo(at(i), lattice());
  }

  //*******************************************************************************************

  double SymGroup::max_error() {

    return m_max_error;
  }

  //*******************************************************************************************

  const Lattice &SymGroup::lattice() const {
    assert(m_lat_ptr && "Attempting to access Lattice of a SymGroup, but it has not been initialized!");
    return *m_lat_ptr;
  }

  //***************************************************

  void SymGroup::_generate_irrep_names() const { //AAB
    irrep_names.resize(conjugacy_classes.size());

    std::vector<int> repeats;

    bool inversion = false;
    bool sigma_h = false;
    bool sigma_v = false;
    bool C2p = false;

    bool sym_wrt_paxis = true;
    bool sym_wrt_v = true;
    bool sym_wrt_C2p = true;
    bool sym_wrt_h = true;
    bool sym_wrt_inv = true;

    bool cubic = false;

    bool all_unique = false;

    _generate_class_names();

    if((name.find("O") != std::string::npos) || (name.find("T") != std::string::npos)) {
      cubic = true;
    }

    for(Index j = 0; j < conjugacy_classes.size(); j++) {
      if(class_names[j].find("i") != std::string::npos) {
        inversion = true;
      }
      if(class_names[j].find("h") != std::string::npos) {
        sigma_h = true;
      }
      if(class_names[j].find("v") != std::string::npos) {
        sigma_v = true;
      }
      if(class_names[j].find("C2'") != std::string::npos) {
        C2p = true;
      }
    }

    //std::cout << "STEP 1: Name according to representation dimension...\n";
    for(Index i = 0; i < irrep_names.size(); i++) {
      sym_wrt_paxis = true;
      if(m_character_table[i][0].real() == 1.0) {
        // Check for symmetry wrt principal axis...
        // This means finding unprimed proper rotations and checking their
        // dimension one characters' signs. If they are all +1, the IrRep
        // is symmetry wrt principal axis. If any of them are -1, then it is
        // antisymmetric wrt principal axis...
        // Symmetric => "A"
        // Antisymmetric => "B"


        // *** IF CUBIC, MUST ONLY CONSIDER THE C3 AXIS HERE!!
        // *** THERE SHOULD BE NO "B"

        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_paxis; j++) {
          if(class_names[j].find("C") != std::string::npos) {
            if(class_names[j].find("'") != std::string::npos) {
              continue;
            }
            else if(m_character_table[i][j].real() < 0) {
              sym_wrt_paxis = false;
            }
          }
        }

        if(sym_wrt_paxis == true) {
          irrep_names[i] = "A";
        }
        else if(!cubic) {
          irrep_names[i] = "B";
        }
        else {
          irrep_names[i] = "A";
        }

      } //Closes loop over 1-d representations

      else if(m_character_table[i][0].real() == 2.0) {
        irrep_names[i] = "E";
      }
      else if(m_character_table[i][0].real() == 3.0) {
        irrep_names[i] = "T";
      }
    }

    // This checks if each representation has a unique name...
    for(Index j = 0; j < irrep_names.size() && all_unique; j++) {
      for(Index k = 0; k < irrep_names.size() && all_unique; k++) {
        all_unique = true;
        if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0)) {
          //	    //std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }
    //std::cout << irrep_names << std::endl;
    //std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      repeats.clear();
      for(Index j = 0; j < irrep_names.size(); j++) {
        for(Index k = 0; k < irrep_names.size(); k++) {
          if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0) && (m_character_table[j][0].real() == m_character_table[k][0].real())) {
            int dim = int(m_character_table[k][0].real());
            if(!contains(repeats, dim)) {
              repeats.push_back(dim);
            }
          }
        }
      }
    }

    //    //std::cout << "Repeats in ..." << repeats << "\n";
    //std::cout << "STEP 2: Name according to sigma_v or C2p symmetry...\n";
    //Check for symmetry wrt vertical mirror plane or C2' axis...
    //I don't know how to do this for non-1D representations =(


    for(Index i = 0; i < irrep_names.size(); i++) {
      int C2ind = 0;
      if((sigma_v == true) && (!all_unique) && (contains(repeats, int(m_character_table[0][i].real()))) && (irrep_names[i] != "E") && (!cubic)) {
        sym_wrt_v = true;
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_v; j++) {
          if(class_names[j].find("v") != std::string::npos) {
            if(class_names[j].find("'") != std::string::npos) {
              continue;
            }
            else if(m_character_table[i][j].real() < 0) {
              sym_wrt_v = false;
            }
          }
        }
        if(sym_wrt_v == true) {
          irrep_names[i].append("1");
        }
        else {
          irrep_names[i].append("2");
        }
      }

      else if((sigma_v != true) && (C2p == true) && (!all_unique) && contains(repeats, int(m_character_table[0][i].real())) && (irrep_names[i] != "E") && (!cubic)) {
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_C2p; j++) {
          if(class_names[j].find("C2") != std::string::npos) {
            if((class_names[j].find("'") != std::string::npos) || (class_names[j].find("''") == std::string::npos)) {
              continue;
            }
            else if(m_character_table[i][j].real() < 0) {
              sym_wrt_C2p = false;
            }
          }
        }
        if(sym_wrt_C2p == true) {
          irrep_names[i].append("1");
        }
        else {
          irrep_names[i].append("2");
        }
      }

      else if(((sigma_v == true) || (C2p == true)) && (!all_unique) && (!cubic)) { //This only deals with 2d stuff...
        for(Index j = 0; j < conjugacy_classes.size(); j++) {
          if((class_names[j].find("C2") != std::string::npos) && (class_names[j].find("'") == std::string::npos) && !cubic) {
            C2ind = j;
          }
        }
        if(m_character_table[i][C2ind].real() > 0) {
          irrep_names[i].append("1");
        }
        else {
          irrep_names[i].append("2");
        }
      }

      else if(((sigma_v == true) || (C2p == true)) && (!all_unique) && (cubic)) {
        for(Index j = 0; j < conjugacy_classes.size(); j++) {
          if(((class_names[j].find("C2") != std::string::npos) && (class_names[j].find("'") != std::string::npos)) || (class_names[j].find("d") != std::string::npos)) {
            C2ind = j;
          }
        }

        if(m_character_table[i][C2ind].real() > 0) {
          if(irrep_names[i] == "A") {
            irrep_names[i].append("1");
          }
          else if(irrep_names[i] == "T") {
            irrep_names[i].append("2");
          }
        }
        else {
          if(irrep_names[i] == "A") {
            irrep_names[i].append("2");
          }
          else if(irrep_names[i] == "T") {
            irrep_names[i].append("1");
          }
        }

      }

      else if((name == "D2h") && (irrep_names[i] == "B")) {
        int c2z = 0, c2y = 0, c2x = 0;
        for(Index j = 0; j < conjugacy_classes.size(); j++) {
          if(class_names[j].find("(z)") != std::string::npos) {
            c2z = j;
          }
          else if(class_names[j].find("(y)") != std::string::npos) {
            c2y = j;
          }
          else if(class_names[j].find("(x)") != std::string::npos) {
            c2x = j;
          }
        }

        if(m_character_table[i][c2z].real() > 0) {
          irrep_names[i].append("1");
        }
        else if(m_character_table[i][c2y].real() > 0) {
          irrep_names[i].append("2");
        }
        else if(m_character_table[i][c2x].real() > 0) {
          irrep_names[i].append("3");
        }
      }
    }

    // This checks if each representation has a unique name...
    for(Index j = 0; j < irrep_names.size() && all_unique; j++) {
      for(Index k = 0; k < irrep_names.size() && all_unique; k++) {
        all_unique = true;
        if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0)) {
          //	    //std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }
    //    //std::cout << irrep_names << std::endl;
    //    //std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      repeats.clear();
      for(Index j = 0; j < irrep_names.size(); j++) {
        for(Index k = 0; k < irrep_names.size(); k++) {
          if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0) && (m_character_table[j][0].real() == m_character_table[k][0].real())) {
            int dim = int(m_character_table[k][0].real());
            if(!contains(repeats, dim)) {
              repeats.push_back(dim);
            }
          }
        }
      }
    }
    //    //std::cout << "Repeats in ..." << repeats << "\n";
    //std::cout << "STEP 3: Name according to inversion symmetry...\n";
    for(Index i = 0; i < irrep_names.size(); i++) {
      if((inversion == true) && (!all_unique) && (contains(repeats, int(m_character_table[0][i].real())))) {
        sym_wrt_inv = true;
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_inv; j++) {
          if((class_names[j].find("i") != std::string::npos) && (m_character_table[i][j].real() < 0)) {
            sym_wrt_inv = false;
          }
        }
        if(sym_wrt_inv == true) {
          irrep_names[i].append("g");
        }
        else {
          irrep_names[i].append("u");
        }
      }
    }

    for(Index j = 0; j < irrep_names.size() && all_unique; j++) {
      for(Index k = 0; k < irrep_names.size() && all_unique; k++) {
        all_unique = true;
        if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0)) {
          //	    //std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }
    //    //std::cout << irrep_names << std::endl;
    //std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      repeats.clear();
      for(Index j = 0; j < irrep_names.size(); j++) {
        for(Index k = 0; k < irrep_names.size(); k++) {
          if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0) && (m_character_table[j][0].real() == m_character_table[k][0].real())) {
            int dim = int(m_character_table[k][0].real());
            if(!contains(repeats, dim)) {
              repeats.push_back(dim);
            }
          }
        }
      }
    }

    //    //std::cout << "Repeats in ..." << repeats << "\n";

    //    //std::cout << "STEP 4: Name according to sigma_h symmetry...\n";
    for(Index i = 0; i < irrep_names.size(); i++) {
      if((sigma_h == true) && (!all_unique) && (contains(repeats, int(m_character_table[0][i].real())))) {
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_h; j++) {
          if((class_names[j].find("h") != std::string::npos) && (m_character_table[i][j].real() < 0)) {
            sym_wrt_h = false;
          }
        }

        if(sym_wrt_h == true) {
          irrep_names[i].append("'");
        }
        else {
          irrep_names[i].append("''");
        }
      }
    }//Close loop over representations

    //    //std::cout << "At the very end, check if everything is unique...\n";

    all_unique = true;
    for(Index j = 0; j < irrep_names.size() && all_unique; j++) {
      for(Index k = 0; k < irrep_names.size() && all_unique; k++) {
        all_unique = true;
        if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0)) {
          std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }

    //    std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      default_err_log() << "WARNING: Failed to name all irreps uniquely...  \n";
    }
    for(auto &irrep_name : irrep_names) {
      std::cout << irrep_name << "   ";
    }
    std::cout << std::endl;
    return;
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

  void SymGroup::_generate_class_names() const { //AAB
    if(!conjugacy_classes.size()) {
      _generate_conjugacy_classes();
    }

    class_names.resize(conjugacy_classes.size());

    std::vector<Eigen::Vector3d > highsym_axes;
    std::vector<int> mult;
    double angle = 360;
    double dprod;
    std::string symtype;
    Index to_name = conjugacy_classes.size();
    bool cubic = false;

    if((name == "T") || (name == "Td") || (name == "Th") || (name == "O") || (name == "Oh")) {
      cubic = true;
    }

    class_names[0] = "E";
    to_name--;

    std::vector<SymInfo> info;
    for(Index i = 0; i < size(); i++)
      info.push_back(SymInfo(at(i), lattice()));

    //If the SymGroup includes an inversion operation, name it "i"
    //We can just use back() here because inversion is always last
    if(info.back().op_type == symmetry_type::inversion_op) {
      class_names[class_names.size() - 1] = "i";
      to_name--;
    }

    if(to_name == 0) {
      return;
    }

    // C1h contains only the identity element and sigma h mirror plane, so we deal with this specifically...
    if(name == "C1h") {
      class_names[1] = "h";
      to_name--;
    }

    if(to_name == 0) {
      return;
    }

    // This part will loop over all of the proper rotations and look for the axis of highest rotational
    // symmetry, which we are going to call the "principal axis."
    highsym_axes.clear();
    mult.clear();

    for(Index i = 0; i < size(); i++) {
      if((info[i].op_type == symmetry_type::rotation_op) || (info[i].op_type == symmetry_type::screw_op)) {
        if(contains(highsym_axes, info[i].axis.const_cart())) { //Otherwise, check if the axis has been found;
          mult[ find_index(highsym_axes, info[i].axis.const_cart())]++;
        }
        else {
          highsym_axes.push_back(info[i].axis.cart());
          mult.push_back(1);
        }
      }
    }

    Eigen::Vector3d hs_axis;
    bool hs_axis_set = false;
    double hangle = 360;

    for(Index i = 0; i < size(); i++) {
      if((info[i].angle < hangle) && (info[i].angle > TOL)) {
        hangle = info[i].angle;
        hs_axis = info[i].axis.cart();
        hs_axis_set = true;
      }
    }


    int order = *std::max_element(mult.begin(), mult.end());

    for(int i = (int(mult.size()) - 1); i >= 0; i--) {
      if((cubic == false) && (mult[i] != order)) {
        highsym_axes.erase(highsym_axes.begin() + i);
        mult.erase(mult.begin() + i);
      }
      else if(mult[i] < (order - 1)) {
        highsym_axes.erase(highsym_axes.begin() + i);
        mult.erase(mult.begin() + i);
      }
    }

    //This is kind of a hack... We might want to change this so that the principal axis
    //is actually the one with the highest fold rotation, and then treat cases like
    //cubic and orthorhombic separately....... but I am not sure that's better...

    if(name == "D3d") {
      for(int i = int(mult.size()) - 1; i >= 0; i--) {
        if(!hs_axis_set) {
          throw std::runtime_error("Error in _generate_class_names: using hs_axis unitialized");
        }
        if(!almost_zero(highsym_axes[i] - hs_axis)) {
          highsym_axes.erase(highsym_axes.begin() + i);
          mult.erase(mult.begin() + i);
        }
      }
    }

    int ind = 0;

    std::string cprime = "'";
    std::string sprime = "'";
    std::string vprime = "";

    for(Index i = 0; i < size(); i++) { //Loop over all SymOps
      //ind = conjugacy_corr_table[i][0];
      ind = index2conjugacy_class[i];
      std::ostringstream s;
      if(!class_names[ind].size()) { //Check to see if this has already been named
        bool normal = false;
        for(Index j = 0; j < highsym_axes.size(); j++) {
          dprod = highsym_axes[j].dot(info[i].axis.const_cart());
          Eigen::Vector3d xprodvec = highsym_axes[j].cross(info[i].axis.const_cart());
          if(almost_zero(xprodvec.norm())) {
            normal = true;
            break;
          }
        }

        if(normal) { //Check if the cross product with principal axis is zero
          if((info[i].angle < 200) && (info[i].angle > 1)) { //Only bother with angles that 360 is divisible by
            if((info[i].op_type == symmetry_type::rotation_op) || (info[i].op_type == symmetry_type::screw_op)) {
              angle = info[i].angle;
              symtype = "C";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle);
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotation\t==> " << s.str() << std::endl;
            }
            else if(info[i].op_type == symmetry_type::rotoinversion_op) {
              angle = info[i].angle;
              symtype = "S";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle);
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotoinversion\t==> " << s.str() << std::endl;
            }
          }
          else if((info[i].op_type == symmetry_type::mirror_op) || (info[i].op_type == symmetry_type::glide_op)) {
            symtype = "h";
            s << conjugacy_classes[ind].size() << symtype;
            class_names[ind] = s.str();
            to_name--;
            //std::cout << "  mirror\t==> " << s.str() << std::endl;
          }
        }
        else {
          if((info[i].angle < 200) && (info[i].angle > 1)) {
            if((info[i].op_type == symmetry_type::rotation_op) || (info[i].op_type == symmetry_type::screw_op)) {
              angle = info[i].angle;
              symtype = "C";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle) << cprime;
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotation\t==> " << s.str() << std::endl;
              cprime = "''";
            }
            else if(info[i].op_type == symmetry_type::rotoinversion_op) {
              angle = info[i].angle;
              symtype = "S";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle) << sprime;
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotoinversion\t==> " << s.str() << std::endl;
              sprime = "''";
            }
          }
          else if((info[i].op_type == symmetry_type::mirror_op) || (info[i].op_type == symmetry_type::glide_op)) {
            if(almost_zero(dprod)) {
              symtype = "v";
              s << conjugacy_classes[ind].size() << symtype  << vprime;
              class_names[ind] = s.str();
              to_name--;
              vprime = "'";
            }
            else {
              symtype = "d";
              s << conjugacy_classes[ind].size() << symtype;
              class_names[ind] = s.str();
              to_name--;
            }
            //std::cout << "  mirror\t==> " << s.str() << std::endl;
          }
        }
      }
    }

    for(Index i = 0; i < class_names.size(); i++) {
      std::string one = "1";
      if(class_names[i].find(one) != std::string::npos) {
        class_names[i].erase(class_names[i].find(one), 1);
      }
    }

    //This is also a hack... but again, not sure there is a better way...

    if(name == "D2h") {
      for(Index i = 0; i < conjugacy_classes.size(); i++) {
        if(info[conjugacy_classes[i][0]].op_type == symmetry_type::rotation_op) {
          if(almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() - Eigen::Vector3d::UnitX())) {
            class_names[i].append("(x)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -  Eigen::Vector3d::UnitY())) {
            class_names[i].append("(y)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -  Eigen::Vector3d::UnitY())) {
            class_names[i].append("(z)");
          }
        }
        else if(info[conjugacy_classes[i][0]].op_type == symmetry_type::mirror_op) {
          if(almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -  Eigen::Vector3d::UnitX())) {
            class_names[i].append("(yz)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -  Eigen::Vector3d::UnitY())) {
            class_names[i].append("(xz)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis.const_cart() -  Eigen::Vector3d::UnitZ())) {
            class_names[i].append("(xy)");
          }
        }
      }
    }


    return;
  }


  //*******************************************************************************************
  void SymGroup::print_character_table(std::ostream &stream) {
    _generate_class_names();
    _generate_irrep_names();

    stream.precision(3);

    for(Index i = 0; i < m_character_table.size(); i++) {
      for(Index j = 0; j < m_character_table.size(); j++) {
        //This part cleans up the numbers in the table...
        if(almost_zero(std::abs(m_character_table[i][j].real()) - 0.5)) {
          if(m_character_table[i][j].real() > 0) {
            m_character_table[i][j] = std::complex<double>(0.5, m_character_table[i][j].imag());
          }
          else {
            m_character_table[i][j] = std::complex<double>(-0.5, m_character_table[i][j].imag());
          }
        }
        if(almost_zero(std::abs(m_character_table[i][j].imag()) - 0.5)) {
          if(m_character_table[i][j].imag() > 0) {
            m_character_table[i][j] = std::complex<double>(m_character_table[i][j].real(), 0.5);
          }
          else {
            m_character_table[i][j] = std::complex<double>(m_character_table[i][j].real(), -0.5);
          }
        }
        if(almost_zero(std::abs(m_character_table[i][j].real()) - 0.866)) {
          if(m_character_table[i][j].real() > 0) {
            m_character_table[i][j] = std::complex<double>(0.866, m_character_table[i][j].imag());
          }
          else {
            m_character_table[i][j] = std::complex<double>(-0.866, m_character_table[i][j].imag());
          }
        }
        if(almost_zero(std::abs(m_character_table[i][j].imag()) - 0.866)) {
          if(m_character_table[i][j].imag() > 0) {
            m_character_table[i][j] = std::complex<double>(m_character_table[i][j].real(), 0.866);
          }
          else {
            m_character_table[i][j] = std::complex<double>(m_character_table[i][j].real(), -0.866);
          }
        }
        if(almost_zero(m_character_table[i][j].real())) {
          m_character_table[i][j] = std::complex<double>(0, m_character_table[i][j].imag());
        }
        if(almost_zero(m_character_table[i][j].imag())) {
          m_character_table[i][j] = std::complex<double>(m_character_table[i][j].real(), 0);
        }
      }
    }

    // This will check if there are any non-integer values in the table.
    bool all_int = true;
    for(Index i = 0; i < m_character_table.size(); i++) {
      for(Index j = 0; j < m_character_table.size(); j++) {
        if(!almost_zero(m_character_table[i][j].real() - floor(m_character_table[i][j].real()))) {
          all_int = false;
        }
        if(!almost_zero(m_character_table[i][j].imag() - floor(m_character_table[i][j].imag()))) {
          all_int = false;
        }
      }
    }

    //This loop is for printing integer-only tables.
    if(all_int) {
      stream << " ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << "\n| Group: " << name << " (" << comment << ")";
      int space = (7 * (m_character_table.size() + 1)) - (name.size() + comment.size() + 11);
      for(int i = 0; i < space; i++) {
        stream << " ";
      }

      stream << "|\n ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << std::endl;

      stream << "|       |";
      for(Index i = 0; i < m_character_table.size(); i++) {
        double whitespace = (6 - double(class_names[i].size())) / 2.0;
        double leftspace = ceil(whitespace);
        double rightspace = floor(whitespace);
        for(int j = 0; j < int(leftspace); j++) {
          stream << " ";
        }
        stream << class_names[i];
        for(int j = 0; j < int(rightspace); j++) {
          stream << " ";
        }
        stream << "|";

      }
      stream << std::endl;

      stream << " ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << std::endl;

      for(Index i = 0; i < m_character_table.size(); i++) {
        stream << "|";
        double whitespace = (7 - double(irrep_names[i].size())) / 2.0;
        double leftspace = ceil(whitespace);
        double rightspace = floor(whitespace);
        for(int j = 0; j < int(leftspace); j++) {
          stream << " ";
        }
        stream << irrep_names[i];
        for(int j = 0; j < int(rightspace); j++) {
          stream << " ";
        }
        stream << "|";

        for(Index j = 0; j < m_character_table.size(); j++) {
          //This will make sure the positive and negative numbers align
          if(m_character_table[i][j].real() > 0) {
            stream << "   ";
          }
          else if(m_character_table[i][j].real() < 0) {
            stream << "  ";
          }
          else {
            if(m_character_table[i][j].imag() >= 0) {
              stream << "   ";
            }
            else {
              stream << "  ";
            }
          } //End alignment of negative numbers

          if((m_character_table[i][j].imag() > 0) && (m_character_table[i][j].real() > 0)) {
            if(m_character_table[i][j].imag() == 1) {
              stream << int(m_character_table[i][j].real()) <<  "i" << "  |";
            }
            else {
              stream << int(m_character_table[i][j].real()) << "+" << int(m_character_table[i][j].imag())  << "i" << "  |";
            }
          }
          else if((m_character_table[i][j].imag() < 0) && (m_character_table[i][j].real() > 0)) {
            if(m_character_table[i][j].imag() == -1) {
              stream << int(m_character_table[i][j].real()) << "-"  << "i" << "  |";
            }
            else {
              stream << int(m_character_table[i][j].real()) <<  int(m_character_table[i][j].imag())  << "i" << "  |";
            }
          }
          else if(m_character_table[i][j].imag() == 0) {
            stream << int(m_character_table[i][j].real()) << "  |";
          }
          else if(m_character_table[i][j].real() == 0) {
            if(m_character_table[i][j].imag() == 1) {
              stream  << "i" << "  |";
            }
            else if(m_character_table[i][j].imag() == -1) {
              stream <<  "-"  << "i" << "  |";
            }
            else {
              stream << int(m_character_table[i][j].imag()) << "i" << "  |";
            }
          }
        } //Closes loop over j
        stream << std::endl;
      } //Closes loop over i

      //Prints the line at the bottom of the table
      stream << " ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << std::endl;
    }//Closes the if all_int loop

    //If not everything is an integer, we use this loop instead...
    else {
      stream << " ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "----------------";
      }

      stream << "\n| Group: " << name << " (" << comment << ")";
      int space = (16 * (m_character_table.size() + 1)) - (name.size() + comment.size() + 12);
      for(int i = 0; i < space; i++) {
        stream << " ";
      }

      stream << "|\n ";

      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "----------------";
      }
      stream << std::endl;

      stream << "|               |";
      for(Index i = 0; i < m_character_table.size(); i++) {
        double whitespace = (15 - double(class_names[i].size())) / 2.0;
        double leftspace = ceil(whitespace);
        double rightspace = floor(whitespace);
        for(int j = 0; j < int(leftspace); j++) {
          stream << " ";
        }
        stream << class_names[i];
        for(int j = 0; j < int(rightspace); j++) {
          stream << " ";
        }
        stream << "|";
      }
      stream << std::endl;
      stream << " ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "----------------";
      }
      stream << std::endl;

      for(Index i = 0; i < m_character_table.size(); i++) {
        stream << "|";
        double whitespace = (15 - double(irrep_names[i].size())) / 2.0;
        double leftspace = ceil(whitespace);
        double rightspace = floor(whitespace);
        for(int j = 0; j < int(leftspace); j++) {
          stream << " ";
        }
        stream << irrep_names[i];
        for(int j = 0; j < int(rightspace); j++) {
          stream << " ";
        }
        stream << "|";

        for(Index j = 0; j < m_character_table.size(); j++) {
          //This will make sure the positive and negative numbers align
          if(m_character_table[i][j].real() > 0) {
            stream << "   ";
          }
          else if(m_character_table[i][j].real() < 0) {
            stream << "  ";
          }
          else {
            if(m_character_table[i][j].imag() >= 0) {
              stream << "   ";
            }
            else {
              stream << "  ";
            }
          } //End alignment of negative numbers

          if((m_character_table[i][j].imag() > 0) && (m_character_table[i][j].real() > 0)) { //Long
            if(m_character_table[i][j].imag() == 1) {
              stream << m_character_table[i][j].real() <<  "i" << "  |";
            }
            else {
              stream << m_character_table[i][j].real() << "+" << m_character_table[i][j].imag()  << "i" << "  |";
            }
          }
          else if((m_character_table[i][j].imag() < 0) && (m_character_table[i][j].real() > 0)) { //Long
            if(m_character_table[i][j].imag() == -1) {
              stream << m_character_table[i][j].real() << "-"  << "i" << "  |";
            }
            else {
              stream << m_character_table[i][j].real() << m_character_table[i][j].imag()  << "i" << "  |";
            }
          }
          else if(m_character_table[i][j].imag() == 0) { //Short
            stream << int(m_character_table[i][j].real()) << "           |";
          }
          else if(m_character_table[i][j].real() == 0) { //Short
            if(m_character_table[i][j].imag() == 1) {
              stream  << "i" << "           |";
            }
            else if(m_character_table[i][j].imag() == -1) {
              stream <<  "-"  << "i" << "           |";
            }
            else {
              stream << m_character_table[i][j].imag() << "i" << "          |";
            }
          }
          else if((m_character_table[i][j].imag() > 0) && (m_character_table[i][j].real() < 0)) { //Long
            stream << m_character_table[i][j].real() << "+" << m_character_table[i][j].imag()  << "i" << "  |";
          }
          else { //Long
            stream << m_character_table[i][j].real() << m_character_table[i][j].imag()  << "i" << "  |";
          }

        }//Closes the j loop

        //stream << "\n|-";
        //for(int x=0; x<m_character_table.size(); x++){
        //  stream << "--------------|";
        //}

        stream << "\n";
      } //Closes the i loop
      stream << " ";
      for(Index i = 0; i < m_character_table.size() + 1; i++) {
        stream << "----------------";
      }
      stream << "\n";
    }//Closes else loop (non-integers)

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

  void SymGroup::_generate_centralizers() const { //AAB
    if(!conjugacy_classes.size()) {
      _generate_conjugacy_classes();
    }

    bool all_commute = true;

    centralizer_table.resize(conjugacy_classes.size());

    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      all_commute = true;
      for(Index j = 0; j < conjugacy_classes[i].size() && all_commute; j++) {
        for(Index k = 0; k < multi_table.size(); k++) {
          int ind = conjugacy_classes[i][j];
          if((multi_table[ind][k] == multi_table[k][ind]) && (!contains(centralizer_table[i], k))) {
            centralizer_table[i].push_back(k);
          }
          else {
            all_commute = false;
          }
        }
      }
    }
  }

  //*******************************************************************************************

  void SymGroup::_generate_character_table()const { //AAB
    //std::cout << "Calculating character table of SymGroup " << this << '\n';
    if(!size()) {
      default_err_log() << "WARNING: This is an empty group. It has no character.\n";
      name = "EMPTY";
      return;
    }

    int sigma_h_ind = 0;
    Index d1 = 0, d2 = 0, d3 = 0;
    Index order1 = 0, order2 = 0;
    bool mirror(false);
    bool inversion(false);
    SymGroup subgroup;
    subgroup.m_lat_ptr = m_lat_ptr;

    std::vector<SymInfo> info;
    for(Index i = 0; i < size(); i++)
      info.push_back(SymInfo(at(i), lattice()));

    if(!size() || get_multi_table().size() != size() || !valid_index(get_multi_table()[0][0])) {
      default_err_log() << "WARNING: This is not a group!!!\n";
      name = "NG";
      return;
    }

    _generate_conjugacy_classes();
    //_generate_conjugacy_corr_table();
    _generate_centralizers();

    Index h = multi_table.size(); //This is the dimensionality of the group.
    Index nc = conjugacy_classes.size(); //This is the number of conjugacy classes, which is also the number of irreducible representations.

    m_character_table.resize(nc, std::vector<std::complex<double> >(nc, -7));
    complex_irrep.resize(nc, false);
    irrep_IDs.resize(nc);


    /** We need to figure out the dimensionality of each irreducible representation.
     *  We know that the total number of irreducible representations is the same as
     *  the total number of conjugacy classes, which we are calling nc. Also, we know
     *  that the sum of the squares of the dimensionalities of all of theirreducible
     *  representations must be equal to the size of the group, h.
     */

    for(Index i = 0; i <= nc; i++) {
      for(Index j = 0; j <= nc; j++) {
        for(Index k = 0; k <= nc; k++) {
          if(((1 * i + 4 * j + 9 * k) == h) && ((i + j + k) == nc) && (d1 >= d2) && (d1 >= d3)) {
            d1 = i;
            d2 = j;
            d3 = k;
          }
        }
      }
    }

    // Set characters related to identity
    for(Index i = 0; i < m_character_table.size(); i++) {
      m_character_table[0][i] = 1;
      if(i < d1) {
        m_character_table[i][0] = 1;
      }
      else if(i >= d1 && i < (d1 + d2)) {
        m_character_table[i][0] = 2;
      }
      else {
        m_character_table[i][0] = 3;
      }
    }

    // count over all possible 1d representations
    Counter<std::vector<int> > count(std::vector<int>(nc, -1), std::vector<int>(nc, 1), std::vector<int>(nc, 2));
    int sum = 0;
    order1 = 1;
    std::complex<double> ortho = 0;
    int telem = 0, elem1 = 0, elem2 = 0;
    int try1 = 0, try2 = 0;

    int try1_ind = 0, try2_ind = 0;
    std::vector<bool> check(nc, false);


    do {
      sum = 0;
      ortho = 0;
      bool satisfy_multi_table = true;

      if(count[0] == 1) {

        for(Index i = 0; i < nc; i++) {
          sum += count[i] * conjugacy_classes[i].size();
          check[i] = true;
        }

        if(sum == 0) {
          for(Index i = 0; i < nc && check[i]; i++) {
            for(Index j = 0; j < nc && check[i]; j++) {
              telem = (count[i] * count[j]);

              for(Index x = 0; x < conjugacy_classes[i].size() && check[i]; x++) {
                for(Index y = 0; y < conjugacy_classes[j].size() && check[i]; y++) {
                  elem1 = conjugacy_classes[i][x];
                  elem2 = conjugacy_classes[j][y];

                  try1 = multi_table[elem1][elem2];
                  try2 = multi_table[elem2][elem1];

                  //try1_ind = conjugacy_corr_table[0][try1];
                  //try2_ind = conjugacy_corr_table[0][try2];
                  try1_ind = index2conjugacy_class[try1];
                  try2_ind = index2conjugacy_class[try2];

                  if((telem != count[try1_ind]) || (telem != count[try2_ind])) {
                    check[i] = false;
                    break;
                  }
                }
              }
            }
          }

          std::vector<std::complex<double> > ortharray;

          for(Index i = 0; i < nc; i++) {
            ortho = 0.0;
            ortharray.resize(0);
            for(Index j = 0; j < nc; j++) {
              if((m_character_table[i][1].real() != -7) && (m_character_table[i][1].imag() != 0)) {
                std::complex<double> temp = std::complex<double>(double(count[j] * conjugacy_classes[j].size()), 0.0);
                ortho += (temp * m_character_table[i][j]);
              }
              ortharray.push_back(ortho);
            }
          }

          ortho = 0;

          for(Index i = 0; i < ortharray.size(); i++) {
            if(!almost_zero(ortharray[i])) {
              ortho = std::complex<double>(1, 1);
            }
          }

          Index wtf = 0;
          for(Index i = 0; i < nc; i++) {
            wtf += check[i];
          }

          if(wtf == nc) {
            satisfy_multi_table = true;
          }
          else {
            satisfy_multi_table = false;
          }

          if((ortho == std::complex<double>(0.0, 0.0)) && (satisfy_multi_table == true)) {
            for(Index i = 1; i < nc; i++) {
              if(order1 < nc) {
                m_character_table[order1][i] = count[i];
              }
            }
            order1++;
          }
        }
      }
    }
    while(++count);

    if(order1 == nc) {
      /** If you're in this block, you have finished your character table using only real
       *  +/- 1 characters because all of your representations are dimension 1. Life is
       *  good. Also, your group is Abelian.
       *
       *  The group you are dealing with is one of the following:
       *  C1, C2, C2v, C1h (S1), C2h, S2 (Ci), D2, D2h.
       */

      /** This group is one of the following...
       *  If nc = 1: C1
       *  If nc = 2: C2, C1h, S2.
       *  If nc = 4: C2v, C2h, D2
       *  If nc = 8: D2h
       */

      if(nc == 1) {
        name = "C1";
        latex_name = "C_1";
        comment = "#1";
        return;
      }
      else if(nc == 8) {
        name = "D2h";
        latex_name = "D_{2h}";
        comment = "#47-74";
        return;
      }

      bool inv = false;
      bool mir = false;

      for(Index i = 0; i < size(); i++) {
        if(info[i].op_type == symmetry_type::inversion_op) {
          inv = true;
        }
        else if(info[i].op_type == symmetry_type::mirror_op || info[i].op_type == symmetry_type::glide_op) {
          mir = true;
        }
      }

      if(nc == 2) {
        if(inv == true) {
          name = "S2"; //This can also be called Ci
          latex_name = "S_2";
          comment = "#2";
        }
        else if(mir == true) {
          name = "C1h"; //This can also be called S1 or Cs
          latex_name = "C_{1h}";
          comment = "#6-9";
        }
        else {
          name = "C2";
          latex_name = "C_2";
          comment = "#3-5";
        }
        return;
      }

      else if(nc == 4) {
        if(inv == true) {
          name = "C2h";
          latex_name = "C_{2h}";
          comment = "#10-15";
        }
        else if(mir == true) {
          name = "C2v";
          latex_name = "C_{2v}";
          comment = "#25-46";
        }
        else {
          name = "D2";
          latex_name = "D_2";
          comment = "#16-24";
        }
        return;
      }

      //You're done. Probably don't need to check column ortho, but why not?
      return;
    }

    if(!centralizer_table.size()) {
      _generate_centralizers();
    }

    std::vector<std::complex<double> > centralizers(nc, 0);

    // //std::cout << "----------------------------------------------------------------------------\n";
    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      centralizers[i] = centralizer_table[i].size();
      //   //std::cout << centralizer_table[i].size() << "\t";
    }
    //std::cout << std::endl;

    //std::cout << "----------------------------------------------------------------------------\n";
    //std::cout << "----------------------------------------------------------------------------\n";

    //std::cout << order1 << std::endl;

    if(order1 != nc) {
      /** If you're in this loop, then you're not done yet because you either have
       *  complex dimension 1 characters or higher dimension real characters still
       *  to look for. In this case, we need to check our group for inversion and
       *  mirror symmetries to determine if it can be broken up into subgroups in
       *  a convenient way that will help us construct the rest of the character
       *  table.
       */

      //We need to find sigma_h in our group so we can use it to hit stuff

      for(Index i = 0; i < size(); i++) {
        if(info[i].op_type == symmetry_type::inversion_op) {
          inversion = true;
        }
        else if((info[i].op_type == symmetry_type::mirror_op) || (info[i].op_type == symmetry_type::glide_op)) {
          //if(conjugacy_classes[conjugacy_corr_table[0][i]].size() == 1) {
          if(conjugacy_classes[index2conjugacy_class[i]].size() == 1) {
            sigma_h_ind = i;
            mirror = true;

          }
        }
      }
    }

    /** This part will handle groups that do not have special symmetries (sigma_h mirror
     *  planes or inversion). We then need to know if the group is Abelian. We can check
     *  this by looking at the number of conjugacy classes and comparing to the number of
     *  dimension 1 representations we are supposed to find. So if d1 = nc, the group
     *  must be Abelian. If this is not the case, the group must have higher dimension
     *  representations, which will be dealt with separately.
     */

    if((inversion == false) && (mirror == false) && (d1 == conjugacy_classes.size())) {

      /** Then the group is cyclic & Abelian. This means that we can proceed by determining
       *  the smallest rotation angle in the group, which will give us omega because
       *  omega is just 2*pi*i/3 for 120 degree rotations, 2*pi*i/4 for 90 degree rotations,
       *  or 2*pi*i/6 for 60 degree rotations. We then "guess" powers of omega for our
       *  one of the complex characters if our representation and determine the rest of the
       *  characters in that representation by enforcing the cycle (which is equivalent to
       *  enforcing the multiplication table).
       *
       *  The possibilities here are:
       *  If nc = 3: C3
       *  If nc = 4: C4
       *  If nc = 6: C6
       */

      double angle = 360;

      int generator = 0;
      bool generator_found = false;

      for(Index i = 0; i < size(); i++) {
        if((info[i].angle < angle) && (info[i].angle > TOL)) {
          angle = info[i].angle;
          generator = i;
          generator_found = true;
        }
      }
      if(!generator_found) {
        throw std::runtime_error("Error in _generate_character_table: generator not found");
      }

      /** We need to keep in mind that the angle returned by get_rotation_angle() is
       *  real and is given in degrees. We need to convert it to complex radian form,
       *  which we will do using Euler's formula: e^(ix) = cos(x) + i*sin(x).
       */

      angle = angle * (M_PI / 180);

      std::complex<double> iangle(cos(angle), sin(angle));

      std::vector<std::complex<double> > tchar;
      std::vector<std::complex<double> > tcharconj;
      tchar.resize(nc);
      tcharconj.resize(nc);
      tchar[0] = std::complex<double>(1, 0);
      tchar[generator] = iangle;

      Index cycle = 0;
      Index ind2 = generator;

      while(order1 < nc) {
        std::complex<double> tangle = iangle;
        do {
          tangle *= iangle;
          tchar[multi_table[ind2][generator]] = tangle;
          ind2 = multi_table[ind2][generator];
          cycle++;
        }
        while(cycle < nc);

        /** For every complex representation in the character table, there is
         *  also a complex conjugate representation, so once we have found a
         *  full set of complex characters, we can conjugate them to get another
         *  row in the table.
         */

        for(Index i = 0; i < tchar.size(); i++) {
          tcharconj[i] = conj(tchar[i]);
        }

        for(Index i = 0; i < tchar.size(); i++) {
          if(std::abs(tchar[i].real()) < TOL) {
            tchar[i] = std::complex<double>(0, tchar[i].imag());
            tcharconj[i] = std::complex<double>(0, tcharconj[i].imag());
          }
          if(std::abs(tchar[i].imag()) < TOL) {
            tchar[i] = std::complex<double>(tchar[i].real(), 0);
            tcharconj[i] = std::complex<double>(tcharconj[i].real(), 0);
          }
        }
        m_character_table[order1] = tchar;
        complex_irrep[order1] = true;
        order1++;
        m_character_table[order1] = tcharconj;
        complex_irrep[order1] = true;
        order1++;
        iangle *= iangle;
      }

      // Should check for orthogonality here...

      if(nc == 3) {
        name = "C3";
        latex_name = "C_3";
        comment = "#143-146";
      }
      else if(nc == 4) {
        bool inv = false;
        for(Index g = 0; g < size() && !inv; g++) {
          if(info[g].op_type == symmetry_type::rotoinversion_op) {
            inv = true;
          }
        }

        if(!inv) {
          name = "C4";
          latex_name = "C_4";
          comment = "#75-80";
        }
        else {
          name = "S4";
          latex_name = "S_4";
          comment = "81-82";
        }
      }
      else if(nc == 6) {
        name = "C6";
        latex_name = "C_6";
        comment = "#168-173";
      }

      return;
    }

    else if((mirror == true) || (inversion == true) /*&& (size()>3)*/) {
      /** If we are in this loop, then our group has mirror or inversion symmetry,
       *  which means that we can determine its character table by breaking it up
       *  into subgroups and computing their character tables first. It turns out
       *  that all groups containing i or sigma_h symmetry can be broken up into
       *  the product of some other group and C1h or S2 group, for which the
       *  character table is really simple.
       *
       *  We can find this subgroup by looking at the sign representation, or
       *  equivalently, the determinants of our SymOp matrices. If we collect
       *  all SymOps that have positive determinant matrices, we will have found
       *  the large subgroup we are looking for.
       */

      for(Index i = 0; i < size(); i++) {
        if(at(i).matrix().determinant() > 0) {
          subgroup.push_back(at(i));
        }
      }

      /** Yes, this part is going to be recursive and therefore scary. But don't
       *  worry, we thought this through. Any subgroup that is created by only
       *  including group elements with positive determinant matrices cannot
       *  have either inversion or sigma_h symmetry. Thus, we can never end
       *  up in this loop more than once...
       */

      subgroup._generate_character_table();
      subgroup._generate_class_names();
      //subgroup.print_character_table(std::cout);

      /** This next part is a little tricky. We're effectively going to paste
       *  our sub-character-table into our big character table four times. We
       *  have to be careful because we need to make sure that the columns
       *  or conjugacy classes that we are copying correspond correctly to
       *  the columns we are copying them into. Luckily, we know which SymOps
       *  the columns correspond do, we know a little bit about the order
       *  that they appear in, and we can check the weird ones by just hiting
       *  the normal ones with inversion or sigma_h. We're going to use the
       *  inversion symmetry whenever it is possible. There is at least one
       *  case where we are forced to use mirror symmetry, but this will
       *  become obvious.
       */

      //std::cout << "Subgroup name is : " << subgroup.name << std::endl;

      for(Index i = 0; i < subgroup.size(); i++) { //Loop over subgroup SymOps
        for(Index j = 0; j < size(); j++) { //Loop over big group SymOps
          //This loop is for the normal SymOps
          if(almost_equal(subgroup[i].matrix(), at(j).matrix(), TOL)) {
            //Index ind_i = subgroup.conjugacy_corr_table[0][i];
            //Index ind_j = conjugacy_corr_table[0][j];
            Index ind_i = subgroup.index2conjugacy_class[i];
            Index ind_j = index2conjugacy_class[j];

            for(Index k = 0; k < subgroup.m_character_table.size(); k++) { //rows
              m_character_table[k][ind_j] = subgroup.m_character_table[k][ind_i];
              m_character_table[subgroup.m_character_table.size() + k][ind_j] = subgroup.m_character_table[k][ind_i];
            }
          }
          //This loop is for the inverted SymOps
          if(almost_equal(subgroup[i].matrix()*back().matrix(), at(j).matrix(), TOL) && inversion) {
            //Index ind_i = subgroup.conjugacy_corr_table[0][i];
            //Index ind_j = conjugacy_corr_table[0][j];
            Index ind_i = subgroup.index2conjugacy_class[i];
            Index ind_j = index2conjugacy_class[j];
            for(Index k = 0; k < subgroup.m_character_table.size(); k++) {
              m_character_table[k][ind_j] = subgroup.m_character_table[k][ind_i];
              m_character_table[subgroup.m_character_table.size() + k][ind_j] = -subgroup.m_character_table[k][ind_i];
            }

            if(nc == 12) {
              if(h == 12) {
                name = "C6h";
                latex_name = "C_{6h}";
                comment = "#175-176";
              }
              else if(h == 24) {
                name = "D6h";
                latex_name = "D_{6h}";
                comment = "#191-194";
              }
            }
            else if(nc == 10) {
              if(h == 16) {
                name = "D4h";
                latex_name = "D_{4h}";
                comment = "#123-142";
              }
              else if(h == 48) {
                name = "Oh";
                latex_name = "O_h";
                comment = "#221-230";
              }
            }

            else if(nc == 8) {
              if(subgroup.name == "C4") {
                name = "C4h";
                latex_name = "C_{4h}";
                comment = "#83-88";
              }
              else if(subgroup.name == "D2") {
                name = "D2h";
                latex_name = "D_{2h}";
                comment = "#47-74";
              }
              else if(subgroup.name == "T") {
                name = "Th";
                latex_name = "T_h";
                comment = "#200-206";
              }
            }
            else if(nc == 6) {
              if(subgroup.name == "D3") {
                name = "D3d";
                latex_name = "D_{3d}";
                comment = "#162-167";
              }
              else if(subgroup.name == "C3") {
                name = "S6";
                latex_name = "S_6";
                comment = "#147-148";
              }
            }
          }

          //This loop is for mirrored SymOps when there is no inversion
          //This only happens in groups C3h and D3h

          else if(almost_equal(subgroup[i].matrix()*at(sigma_h_ind).matrix(), at(j).matrix(), TOL) && mirror && !inversion) {
            //Index ind_i = subgroup.conjugacy_corr_table[0][i];
            //Index ind_j = conjugacy_corr_table[0][j];
            Index ind_i = subgroup.index2conjugacy_class[i];
            Index ind_j = index2conjugacy_class[j];
            for(Index k = 0; k < subgroup.m_character_table.size(); k++) {
              m_character_table[k][ind_j] = subgroup.m_character_table[k][ind_i];
              m_character_table[subgroup.m_character_table.size() + k][ind_j] = -subgroup.m_character_table[k][ind_i];
            }

            if(h == 12) {
              name = "D3h";
              latex_name = "D_{3h}";
              comment = "#187-190";
            }
            else if(h == 6) {
              name = "C3h";
              latex_name = "C_{3h}";
              comment = "#174";
            }
          }
        }
      }

      //You're done! Check orthogonality and centralizers just in case.
      return;
    }

    else if((d3 == 0) && (d2 > 0)) {
      std::vector<std::complex<double> > centrcheck;
      centrcheck.resize(nc);

      for(Index i = 0; i < d1; i++) {
        for(Index j = 0; j < m_character_table.size(); j++) {
          centrcheck[j] += (m_character_table[i][j] * m_character_table[i][j]);
        }
      }

      if(d2 == 1) {
        for(Index i = 1; i < nc; i++) {
          if(abs(centrcheck[i] - centralizers[i]) < TOL) {
            m_character_table[nc - 1][i] = std::complex<double>(0, 0);
          }
          else {
            //This is a bit of a hack, but it turns out to always be true...
            m_character_table[nc - 1][i] = -sqrt(abs(centrcheck[i] - centralizers[i]));
          }
        }

        bool mir = false;
        int mircount = 0;

        for(Index i = 0; i < size(); i++) {
          if(info[i].op_type == symmetry_type::mirror_op || info[i].op_type == symmetry_type::glide_op) {
            mir = true;
            mircount++;
          }
        }

        if(nc == 3) {
          if(mir == true) {
            name = "C3v";
            latex_name = "C_{3v}";
            comment = "#156-161";
          }
          else {
            name = "D3";
            latex_name = "D_3";
            comment = "#149-155";
          }
        }
        else if(nc == 5) {
          if(mircount == 0) {
            name = "D4";
            latex_name = "D_4";
            comment = "#89-98";
          }
          else if(mircount == 2) {
            name = "D2d";
            latex_name = "D_{2d}";
            comment = "#111-122";
          }
          else if(mircount == 4) {
            name = "C4v";
            latex_name = "C_{4v}";
            comment = "#99-110";
          }
        }

        //You're done! Check orthogonality just in case.
        return;
      }

      else if((d2 == 2) && !subgroup.size()) {
        std::vector<std::complex<double> > remainder;
        int nr = 0;
        remainder.resize(nc);
        for(Index i = 1; i < nc; i++) {
          if(abs(centrcheck[i] - centralizers[i]) < TOL) {
            for(Index j = 0; j < d2; j++) {
              m_character_table[(nc - j) - 1][i] = std::complex<double>(0, 0);
              remainder[i] = std::complex<double>(0, 0);
            }
          }
          else {
            remainder[i] = sqrt(abs(centralizers[i] - centrcheck[i]) / 2);
            nr++;
          }
        }

        Counter<std::vector<int> > signcount(std::vector<int>(nr, -1), std::vector<int>(nr, 1), std::vector<int> (nr, 2));
        std::vector<std::complex<double> > trow;
        order2 = (nc - d2);
        std::vector<std::vector<std::complex<double> >> tset;

        do {
          trow.resize(nc);
          int j = 0;
          for(Index i = 0; i < nc; i++) {
            if(m_character_table[nc - 1][i].real() != -7) {
              trow[i] = m_character_table[nc - 1][i];
            }
            if(m_character_table[order2][i].real() == -7) {
              trow[i] = signcount[j] * remainder[i].real();
              j++;
            }
          }

          int sum = 0;

          for(Index i = 0; i < nc; i++) {
            sum += trow[i].real() * conjugacy_classes[i].size();
          }
          double orthcheck;
          std::vector<double> ortharray;
          if(sum == 0) {
            ortharray.resize(0);
            for(Index i = 0; i < d1; i++) {
              orthcheck = 0.0;
              for(Index j = 0; j < nc; j++) {
                double temp = (trow[j].real()) * (conjugacy_classes[j].size());
                orthcheck += (temp * m_character_table[i][j].real());
              }
              ortharray.push_back(orthcheck);
            }

            ortho = 0;

            for(Index i = 0; i < ortharray.size(); i++) {
              if(!almost_zero(ortharray[i])) {
                ortho = std::complex<double>(1, 1);
              }
            }
          }

          if((sum == 0) && (ortho.real() == 0) && (ortho.imag() == 0)) {
            tset.push_back(trow);
          }

        }
        while(++signcount);

        if(tset.size() == d2) {
          for(Index i = 0; i < d2; i++) {
            m_character_table[order2 + i] = tset[i];
          }
        }

        bool mir = false;
        for(Index i = 0; i < size(); i++) {
          if(info[i].op_type == symmetry_type::mirror_op || info[i].op_type == symmetry_type::glide_op) {
            mir = true;
          }
        }

        if(mir == true) {
          name = "C6v";
          latex_name = "C_{6v}";
          comment = "#183-186";
        }
        else {
          name = "D6";
          latex_name = "D_6";
          comment = "#177-182";
        }

        //You're done. Check orthogonality of columns.
        return;
      }
    }

    else if(d3 != 0) {
      /** At this point, we gave up on elegance entirely and decided to
       *  tabulate the three remaining character tables that would be
       *  stupidly impractical to attempt to calculate. So here they are.
       *  Point groups: Td, T, O
       */

      if(nc == 4) {
        name = "T";
        latex_name = "T";
        comment = "#195-199";
        complex_irrep[1] = true;
        complex_irrep[2] = true;
        std::complex<double> omega = std::complex<double>(-0.5, 0.866);
        for(Index j = 1; j < nc; j++) {
          if(conjugacy_classes[j].size() == 3) {
            m_character_table[1][j] = std::complex<double>(1, 0);
            m_character_table[2][j] = std::complex<double>(1, 0);
            m_character_table[3][j] = std::complex<double>(-1, 0);
          }
          else {
            m_character_table[1][j] = omega;
            m_character_table[2][j] = omega * omega;
            m_character_table[3][j] = std::complex<double>(0, 0);
            omega *= omega;
          }
        }
      }

      bool sigma_d = false;

      for(Index i = 0; i < size(); i++) {
        if((info[i].op_type == symmetry_type::mirror_op) || (info[i].op_type == symmetry_type::glide_op)) {
          sigma_d = true;
        }
      }

      if(sigma_d == true) {
        name = "Td";
        latex_name = "T_d";
        comment = "#215-220";

        for(Index j = 1; j < nc; j++) { //Loop over columns
          if(conjugacy_classes[j].size() == 3) {
            m_character_table[1][j] = std::complex<double>(1, 0);
            m_character_table[2][j] = std::complex<double>(2, 0);
            m_character_table[3][j] = std::complex<double>(-1, 0);
            m_character_table[4][j] = std::complex<double>(-1, 0);
          }
          else if(conjugacy_classes[j].size() == 8) {
            m_character_table[1][j] = std::complex<double>(1, 0);
            m_character_table[2][j] = std::complex<double>(-1, 0);
            m_character_table[3][j] = std::complex<double>(0, 0);
            m_character_table[4][j] = std::complex<double>(0, 0);
          }
          else if((conjugacy_classes[j].size() == 6) && (info[conjugacy_classes[j][0]].angle > 10)) {
            m_character_table[1][j] = std::complex<double>(-1, 0);
            m_character_table[2][j] = std::complex<double>(0, 0);
            m_character_table[3][j] = std::complex<double>(1, 0);
            m_character_table[4][j] = std::complex<double>(-1, 0);
          }
          else {
            m_character_table[1][j] = std::complex<double>(-1, 0);
            m_character_table[2][j] = std::complex<double>(0, 0);
            m_character_table[3][j] = std::complex<double>(-1, 0);
            m_character_table[4][j] = std::complex<double>(1, 0);
          }
        }
      }

      else if((sigma_d == false) && (nc == 5)) {
        name = "O";
        latex_name = "O";
        comment = "#207-214";
        for(Index j = 1; j < nc; j++) { //Loop over columns
          if(conjugacy_classes[j].size() == 3) {
            m_character_table[1][j] = std::complex<double>(1, 0);
            m_character_table[2][j] = std::complex<double>(2, 0);
            m_character_table[3][j] = std::complex<double>(-1, 0);
            m_character_table[4][j] = std::complex<double>(-1, 0);
          }
          else if(conjugacy_classes[j].size() == 8) {
            m_character_table[1][j] = std::complex<double>(1, 0);
            m_character_table[2][j] = std::complex<double>(-1, 0);
            m_character_table[3][j] = std::complex<double>(0, 0);
            m_character_table[4][j] = std::complex<double>(0, 0);
          }
          else if((conjugacy_classes[j].size() == 6) && (almost_equal(info[conjugacy_classes[j][0]].angle, 180.0))) {
            m_character_table[1][j] = std::complex<double>(-1, 0);
            m_character_table[2][j] = std::complex<double>(0, 0);
            m_character_table[3][j] = std::complex<double>(-1, 0);
            m_character_table[4][j] = std::complex<double>(1, 0);
          }
          else {
            m_character_table[1][j] = std::complex<double>(-1, 0);
            m_character_table[2][j] = std::complex<double>(0, 0);
            m_character_table[3][j] = std::complex<double>(1, 0);
            m_character_table[4][j] = std::complex<double>(-1, 0);
          }
        }
      }
    }

    return;
  }

  //*******************************************************************************************
  void SymGroup::_generate_elem_order_table() const {

    if(!get_multi_table().size()) {
      return;
    }

    elem_order_table.resize(size());

    for(Index i = 0; i < size(); i++) {
      elem_order_table[i].push_back(i);
    }

    Index telem = 0;
    Index ind = 0;

    for(Index i = 0; i < size(); i++) {
      telem = i;
      ind = 0;
      do {
        telem = multi_table[i][telem];
        ind++;
      }
      while(telem != i);
      elem_order_table[i].push_back(ind);
    }

    return;
  }



  //*******************************************************************************************
  std::vector< std::set<std::set<Index> > > SymGroup::_small_subgroups() const {
    std::vector< std::set<std::set<Index> > > result;

    int tempind = 0;
    Index i, j;
    if(get_alt_multi_table().size() != size()) {
      return result;
    }
    // identity is a small subgroup
    result.push_back({{0}});

    for(i = 1; i < multi_table.size(); i++) {
      std::set<Index> tgroup({0, i});
      j = i;//ind_prod(i, i);
      while(j != 0) {
        j = ind_prod(i, j);
        tgroup.insert(j);
      }

      bool add = true;
      for(auto const &orbit : result) {
        if(orbit.begin()->size() == tgroup.size()
           && orbit.count(tgroup) > 0) {
          add = false;
          break;
        }
      }
      if(!add)
        continue;

      // use equiv_map to find the equivalent subgroups
      result.push_back({});
      //std::cout << "tgroup.size(): " << tgroup.size() << "\n";
      //std::cout << "tgroup : ";
      //for(i : tgroup)
      //std::cout << i << "  ";
      //std::cout << "\n";
      for(auto const &coset : left_cosets(tgroup.begin(), tgroup.end())) {
        std::set<Index> tequiv;
        for(Index op : tgroup) {
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
    for(i = 0; i < m_subgroups.size(); i++) {
      //std::cout << "i is " << i << " and m_subgroups.size() is " << m_subgroups.size() << std::endl;
      for(auto const &orbit : small) {
        for(auto const &equiv : orbit) {
          std::set<Index> tgroup = *(m_subgroups[i].begin());
          Index init_size = tgroup.size();
          tgroup.insert(equiv.begin(), equiv.end());
          if(tgroup.size() == init_size)
            continue;


          // find group closure
          std::vector<Index> vgroup(tgroup.begin(), tgroup.end());
          for(ii = 0; ii < vgroup.size(); ii++) {
            for(jj = 0; jj < vgroup.size(); jj++) {
              Index prod = ind_prod(vgroup[ii], vgroup[jj]);
              if(tgroup.insert(prod).second) {
                vgroup.push_back(prod);
              }
            }
          }

          for(k = 0; k < m_subgroups.size(); k++) {
            if(m_subgroups[k].begin()->size() == tgroup.size()
               && m_subgroups[k].count(tgroup) > 0)
              break;
          }
          //std::cout << " k is " << k << " and m_subgroups.size() is " << m_subgroups.size() << std::endl;
          if(k < m_subgroups.size())
            continue;
          // add the new group

          // use equiv_map to find the equivalent subgroups
          m_subgroups.push_back({});

          for(auto const &coset : left_cosets(tgroup.begin(), tgroup.end())) {
            std::set<Index> tequiv;
            for(Index op : tgroup) {
              tempind = ind_prod(op, ind_inverse(coset[0]));
              tequiv.insert(ind_prod(coset[0], tempind));
            }
            m_subgroups.back().insert(std::move(tequiv));
          }
        }
      }
    }

    // Sort subgroups by number of elements and multiplicity
    for(Index i = 0; i < m_subgroups.size(); i++) {
      for(Index j = i + 1; j < m_subgroups.size(); j++) {
        if(m_subgroups[i].begin()->size() < m_subgroups[j].begin()->size())
          std::swap(m_subgroups[i], m_subgroups[j]);
        else if(m_subgroups[i].begin()->size() == m_subgroups[j].begin()->size()
                && m_subgroups[i].size() < m_subgroups[j].size())
          std::swap(m_subgroups[i], m_subgroups[j]);
      }
    }
    return;
  }

  //*******************************************************************************************


  std::vector<SymGroup> SymGroup::unique_subgroups() const {
    if(!m_subgroups.size()) _generate_subgroups();

    std::vector<std::string> sg_names, sg_names_limited;
    std::vector<bool> chosen_flag(m_subgroups.size(), false);
    for(Index i = 0; i < m_subgroups.size(); i++) {
      SymGroup sgroup;
      sgroup.m_lat_ptr = m_lat_ptr;
      for(Index op : * (m_subgroups[i].begin())) {
        sgroup.push_back(at(op));
      }
      sgroup._generate_character_table();
      sg_names.push_back(sgroup.name);
      //std::cout << sgroup.name << "-" << i << " has equivalencies:\n";
      for(Index j = 0; j < m_subgroups[i].size(); j++) {
        //std::cout << "  " << m_subgroups[i][j] << "\n";
      }
      //std::cout << "\n";
    }

    std::vector< std::vector< Index > > sg_tree(m_subgroups.size(), std::vector<Index>());
    for(Index i = 0; i < m_subgroups.size(); i++) {
      //std::cout << "Subgroup " << sg_names[i] << "-" << i << " is also a subgroup of ";
      for(Index j = 0; j < m_subgroups.size(); j++) {
        for(auto const &equiv : m_subgroups[j]) {
          bool add = true;
          for(auto const &op : equiv) {
            if(m_subgroups[i].begin()->count(op) < 1) {
              add = false;
              break;
            }
          }
          if(add) {
            sg_tree[i].push_back(j);
            //std::cout << sg_names[j] << "-" << j << "-" << jj << "  ";
            break;
          }
        }
      }
      //std::cout << "\n";
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

    for(Index i = 0; i < m_subgroups.size(); i++) {
      if(chosen_flag[i]) continue;
      unique_sgroups.push_back(SymGroup());
      for(Index op :  * (m_subgroups[i].begin())) {
        unique_sgroups.back().push_back(at(op));
      }
      unique_sgroups.back().sort();
      unique_sgroups.back()._generate_character_table();

      //std::cout << "Added group " << unique_sgroups.back().name  << " having bit-string " << m_subgroups[i][0] << std::endl;
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

    if(get_multi_table().size() != size()) {
      return;
    }


    conjugacy_classes.clear();

    int k;

    for(Index i = 0; i < size(); i++) {
      bool dup_class(false);
      for(Index j = 0; j < conjugacy_classes.size(); j++) {
        if(contains(conjugacy_classes[j], i)) {
          dup_class = true;
          break;
        }
      }
      if(dup_class) continue;

      conjugacy_classes.push_back(std::vector<Index>());

      for(Index j = 0; j < size(); j++) {
        //std::cout << "for j=" << j << ", i=" << i << ": j-inverse= " << ind_inverse(j) << ", i*j-inverse= " << ind_prod(i, ind_inverse(j));
        //int tk=alt_multi_table[j][i];
        //std::cout << "-- compare to amt[j][0]=" << alt_multi_table[j][0] << " and amt[j][i]=" << tk << " and result is k=";
        k = ind_prod(j, ind_prod(i, ind_inverse(j)));
        //std::cout << k << " -- compare to explicit value " << multi_table[tk][j];

        if(!contains(conjugacy_classes.back(), k)) {
          //std::cout << " so " << k << " goes in class " << conjugacy_classes.size()-1;
          conjugacy_classes.back().push_back(k);
        }
        //std::cout << "\n";
      }
      std::sort(conjugacy_classes.back().begin(), conjugacy_classes.back().end());
    }

    index2conjugacy_class.resize(size(), 0);
    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      for(Index j = 0; j < conjugacy_classes[i].size(); j++) {
        index2conjugacy_class[conjugacy_classes[i][j]] = i;
      }
    }
    //std::cout << "index2conjugacy is " << index2conjugacy_class << '\n';
    return;
  }

  //*******************************************************************************************

  Index SymGroup::ind_inverse(Index i) const {
    if(get_alt_multi_table().size() != size() || !valid_index(i) || i >= size())
      return -1;
    //std::cout << "Inside ind_inverse. 'i' is " << i << " and alt_multi_table size is " << alt_multi_table.size() << "\n";
    return alt_multi_table[i][0];
  }

  //*******************************************************************************************

  Index SymGroup::ind_prod(Index i, Index j) const {
    if(get_multi_table().size() != size()
       || !valid_index(i) || i >= size()
       || !valid_index(j) || j >= size()) {


      //default_err_log() << "WARNING: SymGroup::ind_prod() failed for " << i << ", " << j << "\n";// and multi_table\n" << get_multi_table() << "\n\n";
      //assert(0);
      return -1;
    }
    //std::cout << "Inside ind_prod. 'i' is " << i << " and j is " << j << " and multi_table size is " << multi_table.size() << "\n";
    return multi_table[i][j];
  }

  //*******************************************************************************************

  Index SymGroup::class_of_op(Index i) const {
    if(!get_conjugacy_classes().size() || i > size())
      return -1;
    return index2conjugacy_class[i];
  }

  //*******************************************************************************************

  void SymGroup::set_irrep_ID(Index i, SymGroupRepID ID) const {
    assert((valid_index(i) && i < irrep_IDs.size()) && "Attempting to set ID for out-of-bounds irrep.");
    irrep_IDs[i] = ID;
    return;

  }

  //*******************************************************************************************

  SymGroupRepID SymGroup::get_irrep_ID(Index i) const {
    if(!valid_index(i) || i >= irrep_IDs.size())
      return SymGroupRepID();

    return irrep_IDs[i];

  }

  //*******************************************************************************************

  SymGroupRepID SymGroup::coord_rep_ID() const {
    if(!size() || !at(0).has_valid_master()) {
      default_err_log() << "CRITICAL ERROR: In SymGroup::get_coord_rep_ID(), SymGroup is improperly initialized.\n"
                        << "                Exiting...\n";
      exit(1);
    }

    return at(0).master_group().coord_rep_ID();
  }
  //*******************************************************************************************

  SymGroupRepID SymGroup::allocate_representation() const {
    if(!size() || !at(0).has_valid_master()) {
      default_err_log() << "CRITICAL ERROR: In SymGroup::allocate_representation(), SymGroup is improperly initialized.\n"
                        << "                Exiting...\n";
      exit(1);
    }

    return at(0).master_group().allocate_representation();
  }

  //*******************************************************************************************

  SymGroupRep const &SymGroup::get_irrep(Index i) const {
    if(!size() || !valid_index(i) || i >= irrep_IDs.size())
      throw std::runtime_error(std::string("Cannot find irrep ") + std::to_string(i) + " in the current SymGroup\n");

    return at(0).master_group().representation(irrep_IDs[i]);

  }

  //*******************************************************************************************
  // The set of left cosets is identical to the equivalence_map formed by partitioning (*this) w.r.t. 'subgroup'
  std::vector<std::vector<Index> > SymGroup::left_cosets(const std::vector<SymOp> &subgroup, double tol) const {
    std::vector<Index> sg_inds = find_all_periodic(subgroup, tol);
    return left_cosets(sg_inds.begin(), sg_inds.end());
  }



  //*******************************************************************************************

  const std::vector<std::vector<Index>> &SymGroup::get_multi_table() const {
    if(multi_table.size() != size()) {
      //std::cout << "CALCULATING MULTI_TABLE for " << this <<  ": table size is " << multi_table.size() << " and group size is " << size() << "!!\n";
      _generate_multi_table();
    }
    return multi_table;
  }

  //*******************************************************************************************

  const std::vector<std::vector<Index>> &SymGroup::get_alt_multi_table() const {
    if(alt_multi_table.size() != size()) {
      //std::cout << "CALCULATING ALT_MULTI_TABLE " << this << ": table size is " << alt_multi_table.size() << " and group size is " << size() << "!!\n";
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
    if(conjugacy_classes.size() != size())
      _generate_conjugacy_classes();
    return conjugacy_classes;
  }
  //*******************************************************************************************

  const std::vector<bool > &SymGroup::get_complex_irrep_list() const {
    if(!m_character_table.size())
      _generate_character_table();
    return complex_irrep;
  }

  //*******************************************************************************************

  const std::vector<std::vector<std::complex<double> >> &SymGroup::character_table() const {
    if(!m_character_table.size())
      _generate_character_table();
    return m_character_table;
  }

  //*******************************************************************************************

  const std::string &SymGroup::get_name() const {
    if(!name.size()) {
      _generate_character_table();
      if(!name.size()) {

        default_err_log() << "WARNING: In SymGroup::get_name(), unable to get symgroup type.\n";
        default_err_log() << "group size is " << size() << '\n';
        name = "unknown";
      }
    }

    return name;
  }

  //*******************************************************************************************

  const std::string &SymGroup::get_latex_name() const {
    character_table();
    return latex_name;
  }

  //*******************************************************************************************

  const std::vector<std::set<std::set<Index> > > &SymGroup::subgroups() const {
    if(!m_subgroups.size())
      _generate_subgroups();
    return m_subgroups;
  }

  //*******************************************************************************************

  bool SymGroup::_generate_multi_table() const { //AAB
    Index i, j;
    multi_table.resize(size(), std::vector<Index>(size(), -1));

    for(i = 0; i < size(); i++) {
      for(j = 0; j < size(); j++) {
        multi_table[i][j] = find_periodic(at(i) * at(j));
        if(multi_table[i][j] >= size() || find_index(multi_table[i], multi_table[i][j]) != j) {
          // this is a hack (sort of). If find_periodic doesn't work, we try find_no trans, which *should* work.
          // In other words, we are using 'inuition' to determine that user doesn't really care about the translational aspects.
          // If our intuition is wrong, there will will probably be an obvious failure later.
          multi_table[i][j] = find_no_trans(at(i) * at(j));
          if(multi_table[i][j] >= size() || find_index(multi_table[i], multi_table[i][j]) != j) {

            //if(multi_table[i][j] >= size()) {
            //std::cout << "This SymGroup is not a group because the combination of at least two of its elements is not contained in the set.\n";

            //ing a table of all 1's seems to make the most sense. This will prevent weird recursion from happening.
            default_err_log() << "Failed to construc multiplication table!  Table in progress:\n";
            for(Index m = 0; m < multi_table.size(); m++) {
              for(Index n = 0; n < multi_table[m].size(); n++) {
                default_err_log() << multi_table[m][n] << "   ";
              }
              default_err_log() << "\n";
            }
            default_err_log() << "\n";
            multi_table.resize(size(), std::vector<Index>(size(), -1));
            //multi_table.clear();
            return false;
          }
        }
      }
    }

    return true;

  }

  //*******************************************************************************************
  void SymGroup::_generate_alt_multi_table() const {
    //by calling get_multi_table(), we ensure that multi_table is populated
    alt_multi_table.resize(get_multi_table().size());

    if(multi_table.size() && !valid_index(multi_table[0][0])) {
      alt_multi_table.resize(size(), std::vector<Index>(size(), -1));
      return;
    }
    for(Index i = 0; i < multi_table.size(); i++) {
      if(multi_table[i][i] != 0) {
        alt_multi_table[find_index(multi_table[i], 0)] = multi_table[i];
      }
      else {
        alt_multi_table[i] = multi_table[i];
      }
    }

  }

  //*******************************************************************************************
  // Please keep in mind that this will only return the FIRST match that is found.
  // This does not guarantee that it is the only match.

  Index SymGroup::find_no_trans(const SymOp &test_op) const {
    for(Index i = 0; i < size(); i++) {
      if(almost_equal(at(i).matrix(), test_op.matrix())) {
        return i;
      }
    }

    return size();
  }

  //*******************************************************************************************

  Index SymGroup::find_periodic(const SymOp &test_op, double tol) const {
    for(Index i = 0; i < size(); i++) {
      if(compare_periodic(at(i), test_op, lattice(), periodicity(), tol)) {
        return i;
      }
    }
    return size();
  }

  //*******************************************************************************************

  std::vector<Index> SymGroup::find_all_periodic(const std::vector<SymOp> &subgroup, double tol) const {
    std::vector<Index> tarray;
    for(Index i = 0; i < subgroup.size(); i++) {
      tarray.push_back(find_periodic(subgroup[i], tol));
      if(tarray.back() == size()) {
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

    for(Index i = 0; i < size(); i++) {
      for(Index j = 0; j < size(); j++) {
        if(!contains_periodic(at(i)*at(j), tol))
          return false;
      }
      if(!contains_periodic(at(i).inverse(), tol))
        return false;
    }
    return true;
  }

  //*******************************************************************************************

  void SymGroup::enforce_group(double tol, Index max_size) {
    bool new_ops(true);

    while(new_ops && size() < max_size) {
      new_ops = false;
      for(Index i = 0; i < size() && size() < max_size; i++) {
        SymOp A_op(at(i).unregistered_copy());
        for(Index j = 0; j < size() && size() < max_size; j++) {
          SymOp B_op(at(j).unregistered_copy());

          SymOp tOp(A_op * B_op);

          if(!contains_periodic(tOp, tol)) {
            push_back(within_cell(tOp, lattice(), periodicity()));
            new_ops = true;
            //	  //std::cout << "Pushing back a SymOp due to multiplication fail.\n";
          }
        }

        SymOp tOp(A_op.inverse());
        if(!contains_periodic(tOp, tol)) {
          push_back(within_cell(tOp, lattice(), periodicity()));
          new_ops = true;
          //std::cout << "Pushing back a SymOp due to inverse fail.\n";
        }
      }
    }
    if(size() >= max_size - 1) {
      default_err_log() << "In SymGroup::enforce_group() -- you have reached the maximum allowed size you specified for your group (the default is 200). Unless you are generating a factor group in a large supercell, you probably need to adjust your tolerances.\n";
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
    for(Index i = 0; i < size(); i++)
      at(i).apply_sym(op);
    return *this;
  }

  //*******************************************************************************************

  void SymGroup::write(std::string filename, COORD_TYPE mode) const {
    std::ofstream outfile;
    outfile.open(filename.c_str());
    print(outfile, mode);
    outfile.close();
    return;
  }

  //*******************************************************************************************

  void SymGroup::print(std::ostream &out, COORD_TYPE mode) const {
    out << size() << " # " << xtal::COORD_MODE::NAME(mode) << " representation of group containing " << size() << " elements:\n\n";
    Eigen::Matrix3d c2f_mat(Eigen::Matrix3d::Identity());
    if(mode == FRAC)
      c2f_mat = lattice().inv_lat_column_mat();
    for(Index i = 0; i < size(); i++) {
      out << i << "  " << description(at(i), lattice(), mode) << "\n";
      at(i).print(out, c2f_mat);
      out << std::endl;
    }
    return;
  }


  //*******************************************************************************************

  void SymGroup::calc_space_group_in_cell(SymGroup &space_group_cell, const Lattice &_cell) const {
    if(!size()) return;

    Eigen::Vector3i max_trans(3, 3, 3);
    Coordinate trans(Eigen::Vector3d::Zero(), _cell, FRAC);
    space_group_cell.clear();

    std::vector<SymInfo> sg_info;
    for(Index i = 0; i < size(); i++) {
      EigenCounter<Eigen::Vector3i> lat_comb(-max_trans, max_trans, Eigen::Vector3i::Ones());
      do {
        trans.frac() = lat_comb().cast<double>();
        SymOp new_sym(SymOp::translation(trans.cart())*at(i));
        SymInfo info(new_sym, lattice());
        trans = info.location;
        if(!trans.is_within()) {
          continue;
        }


        bool new_location = true;
        for(Index j = 0; j < space_group_cell.size(); j++) {

          if(almost_equal(new_sym.matrix(), space_group_cell[j].matrix()) && almost_equal(info.location.const_cart(), sg_info[j].location.const_cart())) {
            new_location = false;
            break;

          }
        }
        if(new_location) {
          space_group_cell.push_back(new_sym);
          sg_info.push_back(info);
        }
      }
      while(++lat_comb);
    }

    return;
  }

  //*******************************************************************************************

  void SymGroup::calc_space_group_in_range(SymGroup &space_group,
                                           const Lattice &_cell,
                                           Eigen::Vector3i min_trans,
                                           Eigen::Vector3i max_trans) const {
    if(!size()) return;


    Coordinate trans(Eigen::Vector3d::Zero(), _cell, FRAC);

    for(Index i = 0; i < size(); i++) {
      EigenCounter<Eigen::Vector3i> lat_comb(min_trans, max_trans, Eigen::Vector3i::Ones());
      do {
        trans.frac() = lat_comb().cast<double>();

        SymOp new_sym(SymOp::translation(trans.cart())*at(i));

        if(!contains(space_group, new_sym)) {
          space_group.push_back(new_sym);

        }

      }
      while(++lat_comb);
    }

    return;
  }

  //***************************************************
  void SymGroup::print_locations(std::ostream &stream) const {
    //Assumes SymGroup is sorted with clumps of SymOps of common matrix type and eigenvec
    //sort();

    bool new_op = true;
    stream << "Locations for symmetry operations\n";
    SymInfo info(at(0), lattice());
    SymInfo next_info(info);
    Eigen::Matrix3d c2f_mat = lattice().inv_lat_column_mat();
    for(Index i = 0; i < size(); i++) {
      if(new_op) {
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

      if(i + 1 < size()) {
        next_info = SymInfo(at(i + 1), lattice());
        if(info.op_type == next_info.op_type && almost_equal(info.axis.const_cart(), next_info.axis.const_cart())) {
          //Is this enough to know if it's a new symmetry or not?
          new_op = false;
        }
        else {
          new_op = true;
          stream << "----------------------------------------------------------------------------------\n\n";
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
  /// SymOp are sorted by lexicographical comparison of: (-det, -trace, angle, axis, tau)
  /// - angle is positive
  /// - axis[0] is positive
  ///
  void SymGroup::sort() {

    // floating point comparison tolerance
    double tol = TOL;

    //COORD_TYPE print_mode = CART;

    // compare on vector of '-det', '-trace', 'angle', 'axis', 'tau'
    typedef Eigen::Matrix<double, 10, 1> key_type;
    auto make_key = [](const SymOp & op, const Lattice & lat) {

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

      vec.segment<3>(offset) =  Coordinate(op.tau(), lat, CART).const_frac();
      offset += 3;

      return vec;
    };

    // define symop compare function
    auto op_compare = [tol](const key_type & A, const key_type & B) {
      return float_lexicographical_compare(A, B, tol);
    };

    typedef std::map<key_type, SymOp, std::reference_wrapper<decltype(op_compare)> > map_type;

    // sort conjugacy class using the first symop in the sorted class
    auto cclass_compare = [tol](const map_type & A, const map_type & B) {
      return float_lexicographical_compare(A.begin()->first, B.begin()->first, tol);
    };

    // sort elements in each conjugracy class (or just put all elements in the first map)
    std::set<map_type, std::reference_wrapper<decltype(cclass_compare)> > sorter(cclass_compare);

    // first put identity in position 0 in order to calculat multi_table correctly
    for(int i = 0; i < size(); ++i) {
      if(at(i).is_identity()) {
        std::swap(at(0), at(i));
        break;
      }
    }

    // if sorting by congujacy classes
    if(get_multi_table().size() == size()) {

      // get conjugacy classes
      get_conjugacy_classes().size();

      // insert elements into each conjugacy class to sort the class
      for(int i = 0; i < conjugacy_classes.size(); ++i) {
        map_type cclass(op_compare);
        for(int j = 0; j < conjugacy_classes[i].size(); ++j) {
          const SymOp &op = at(conjugacy_classes[i][j]);
          cclass.insert(std::make_pair(make_key(op, lattice()), op));
        }

        sorter.emplace(std::move(cclass));
      }
    }
    else {
      // else just sort element
      map_type all_op(op_compare);
      for(auto it = begin(); it != end(); ++it) {
        all_op.emplace(make_key(*it, lattice()), *it);
      }
      sorter.emplace(std::move(all_op));
    }

    // copy symop back into group
    int j = 0;
    for(auto const &cclass : sorter) {
      for(auto it = cclass.begin(); it != cclass.end(); ++it) {
        at(j) = it->second;
        ++j;
      }
    }

    clear_tables();
  }

  //*******************************************************************************************

  bool SymGroup::is_irreducible() const {
    double tvalue = 0;

    for(Index i = 0; i < size(); i++) {
      tvalue += (at(i).matrix().trace()) * (at(i).matrix().trace());
    }

    if(Index(std::abs(tvalue)) == size()) {
      return true;
    }
    else {
      return false;
    }
  }

  //*******************************************************************************************
  /// Translation operators for origin shift need to be defined
  SymGroup &SymGroup::operator+=(const Eigen::Ref<const SymOp::vector_type> &shift) {
    for(Index ng = 0; ng < size(); ng++)
      at(ng) += shift;
    return (*this);
  }

  //*******************************************************************************************

  SymGroup &SymGroup::operator-=(const Eigen::Ref<const SymOp::vector_type> &shift) {
    for(Index ng = 0; ng < size(); ng++)
      at(ng) -= shift;
    return (*this);
  }

  //*******************************************************************************************

  std::vector<Index> SymGroup::get_irrep_decomposition() const {
    std::vector<Index> tdec;
    std::vector<double> repchar;

    tdec.resize(conjugacy_classes.size());
    repchar.resize(conjugacy_classes.size());

    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      repchar[i] += at(conjugacy_classes[i][0]).matrix().trace();
    }

    //    //std::cout << "The reducible characters are:\n";

    //    //std::cout << repchar << "\n\n";

    for(Index i = 0; i < m_character_table.size(); i++) { // Loop over irreducible representations
      std::complex<double> temp = std::complex<double>(0, 0);
      for(Index j = 0; j < m_character_table.size(); j++) { // Loop over conjugacy classes
        temp += (m_character_table[i][j] * std::complex<double>(conjugacy_classes[j].size(), 0) * std::complex<double>(repchar[j], 0));
      }
      //      //std::cout << temp << "\t";

      tdec[i] = round(temp.real() / size());
    }

    //    //std::cout << "\n\nThe irreducible projection is:\n";
    //std::cout << tdec << "\n\n";

    return tdec;
  }

  //*******************************************************************************************

  jsonParser &SymGroup::to_json(jsonParser &json) const {
    json.put_obj();

    json["symop"].put<std::vector<SymOp> >(*this);

    for(int i = 0; i < json["symop"].size(); ++i) {
      SymInfo info(at(i), lattice());
      auto &j = json["symop"][i];
      add_sym_info(info, j);
    }

    // PERIODICITY_TYPE group_periodicity;
    json["group_periodicity"] = periodicity();

    // mutable std::vector<std::vector<int> > multi_table;
    json["multi_table"] = multi_table;

    // mutable std::vector<std::vector<int> > alt_multi_table;
    json["alt_multi_table"] = alt_multi_table;

    // mutable std::vector<std::vector<int> > conjugacy_classes;
    json["conjugacy_classes"] = conjugacy_classes;

    // mutable std::vector<std::string> class_names;
    json["class_names"] = class_names;

    // mutable std::vector<int> index2conjugacy_class;
    json["index2conjugacy_class"] = index2conjugacy_class;

    // mutable std::vector<std::vector<std::complex<double> > > m_character_table;
    json["character_table"] = m_character_table;
    for(int i = 0; i < json["character_table"].size(); i++) {
      json["character_table"][i].set_force_row();
      for(int j = 0; j < json["character_table"][i].size(); j++) {
        json["character_table"][i][j].set_force_row();
        json["character_table"][i][j]["real"].set_remove_trailing_zeros();
        json["character_table"][i][j]["imag"].set_remove_trailing_zeros();
      }
    }

    // mutable std::vector<int> irrep_IDs;
    json["irrep_IDs"] = irrep_IDs;

    // mutable std::vector<bool> complex_irrep;
    json["complex_irrep"] = complex_irrep;

    // mutable std::vector<std::string> irrep_names;
    json["irrep_names"] = irrep_names;

    //json["unique_subgroups"] = unique_subgroups();
    json["unique_subgroups"] = std::string("TODO: fix SymGroup::unique_subgroups()");

    // mutable std::vector<std::vector<std::vector<int> > > m_subgroups;
    json["m_subgroups"] = m_subgroups;

    // mutable std::vector<std::vector<int> > centralizer_table;
    json["centralizer_table"] = centralizer_table;

    // mutable std::vector<std::vector<int> > elem_order_table;
    json["elem_order_table"] = elem_order_table;

    // mutable double max_error;
    json["max_error"] = m_max_error;

    for(auto const &info : point_group_info())
      json[info.first] = info.second;

    return json;
  }

  //*******************************************************************************************

  void SymGroup::from_json(const jsonParser &json) {
    try {

      // class SymGroup : public std::vector<SymOp>
      clear();

      for(int i = 0; i < json["symop"].size(); i++) {
        push_back(json["symop"][i].get<SymOp>());
      }

      // PERIODICITY_TYPE group_periodicity;
      CASM::from_json(m_group_periodicity, json["group_periodicity"]);

      // mutable std::vector<std::vector<int> > multi_table;
      CASM::from_json(multi_table, json["multi_table"]);

      // mutable std::vector<std::vector<int> > alt_multi_table;
      CASM::from_json(alt_multi_table, json["alt_multi_table"]);

      // mutable std::vector<std::vector<int> > conjugacy_classes;
      CASM::from_json(conjugacy_classes, json["conjugacy_classes"]);

      // mutable std::vector<std::string> class_names;
      CASM::from_json(class_names, json["class_names"]);

      // mutable std::vector<int> index2conjugacy_class;
      CASM::from_json(index2conjugacy_class, json["index2conjugacy_class"]);

      // mutable std::vector<std::vector<std::complex<double> > > m_character_table;
      CASM::from_json(m_character_table, json["character_table"]);

      // mutable std::vector<int> irrep_IDs;
      CASM::from_json(irrep_IDs, json["irrep_IDs"]);

      // mutable std::vector<bool> complex_irrep;
      CASM::from_json(complex_irrep, json["complex_irrep"]);

      // mutable std::vector<std::string> irrep_names;
      CASM::from_json(irrep_names, json["irrep_names"]);

      // mutable std::vector<std::vector<std::vector<int> > > m_subgroups;
      CASM::from_json(m_subgroups, json["m_subgroups"]);

      // mutable std::vector<std::vector<int> > centralizer_table;
      CASM::from_json(centralizer_table, json["centralizer_table"]);

      // mutable std::vector<std::vector<int> > elem_order_table;
      CASM::from_json(elem_order_table, json["elem_order_table"]);

      // mutable double max_error;
      CASM::from_json(m_max_error, json["max_error"]);

    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //*******************************************************************************************

  jsonParser &to_json(const SymGroup &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  void from_json(SymGroup &sym, const jsonParser &json) {
    try {
      sym.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }
  //*******************************************************************************************

  jsonParser &MasterSymGroup::to_json(jsonParser &json) const {

    // class MasterSymGroup : public SymGroup
    SymGroup::to_json(json);

    // mutable std::vector<SymGroupRep *> m_rep_array;
    json["m_rep_array"].put_array();
    for(Index i = 0; i < m_rep_array.size(); i++) {
      json["rep_array"].push_back(m_rep_array[i]);
    }

    // mutable int coord_rep_ID, reg_rep_ID;
    json["coord_rep_ID"] = m_coord_rep_ID;
    json["reg_rep_ID"] = m_reg_rep_ID;

    // mutable SymGroup point_group_internal;
    json["point_group"] = m_point_group;

    return json;
  }

  //*******************************************************************************************

  // Note: as a hack this expects at(0) to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the SymOp
  void MasterSymGroup::from_json(const jsonParser &json) {
    try {
      clear();

      // class MasterSymGroup : public SymGroup
      SymGroup::from_json(json);

      // mutable std::vector<SymGroupRep *> m_rep_array;
      // destruct exisiting
      for(Index i = 0; i < m_rep_array.size(); i++) {
        delete m_rep_array[i];
      }
      m_rep_array.resize(json["rep_array"].size());
      for(int i = 0; i < json["rep_array"].size(); i++) {
        m_rep_array[i] = new SymGroupRep(*this);
        m_rep_array[i]->from_json(json["rep_array"][i]);
      }

      // mutable int coord_rep_ID, reg_rep_ID;
      CASM::from_json(m_coord_rep_ID, json["coord_rep_ID"]);
      CASM::from_json(m_reg_rep_ID, json["reg_rep_ID"]);

      // mutable SymGroup m_point_group;
      m_point_group.clear();
      m_point_group.from_json(json["point_group"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //*******************************************************************************************

  jsonParser &to_json(const MasterSymGroup &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  void from_json(MasterSymGroup &sym, const jsonParser &json) {
    try {
      sym.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //**********************************************************

  bool compare_periodic(const SymOp &a,
                        const SymOp &b,
                        const Lattice &lat,
                        PERIODICITY_TYPE periodicity,
                        double _tol) {


    if(a.time_reversal() != b.time_reversal() || !almost_equal(a.matrix(), b.matrix(), _tol))
      return false;
    //std::cout << "Operations:\n"
    //<< a.matrix() << "\n and \n"
    //<< b.matrix() << "\n are equal \n"
    //<< "a-tau " << a.tau().transpose() << "\n"
    //<< "b-tau " << b.tau().transpose() << "\n";
    if(periodicity != PERIODIC)
      return almost_equal(a.tau(), b.tau(), _tol);

    return Coordinate(a.tau(), lat, CART).min_dist(Coordinate(b.tau(), lat, CART)) < _tol;

  }

  //**********************************************************

  SymOp within_cell(const SymOp &a,
                    const Lattice &lat,
                    PERIODICITY_TYPE periodicity) {
    if(periodicity != PERIODIC)
      return a;

    Coordinate trans(a.tau(), lat, CART);
    trans.within();
    return SymOp(a.matrix(), trans.cart(), a.time_reversal(), a.map_error());
  }


}
