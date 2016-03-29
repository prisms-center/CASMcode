#include "casm/symmetry/SymGroup.hh"

#include "casm/external/Eigen/CASM_AddOns"
#include "casm/misc/CASM_math.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"

#include "casm/casm_io/json_io/container.hh"
#include "casm/casm_io/json_io/global.hh"

namespace CASM {

  MasterSymGroup::MasterSymGroup(const MasterSymGroup &RHS) :
    SymGroup(RHS), coord_rep_ID(RHS.coord_rep_ID),
    reg_rep_ID(RHS.reg_rep_ID) {
    rep_array.reserve(RHS.rep_array.size());
    for(Index i = 0; i < RHS.rep_array.size(); i++) {
      add_representation(RHS.rep_array[i]->copy());
      rep_array.back()->set_master_group(*this);
    }

    for(Index i = 0; i < size(); i++)
      at(i).set_index(*this, i);
  }

  //***************************************************

  MasterSymGroup::~MasterSymGroup() {
    for(Index i = 0; i < rep_array.size(); i++)
      delete rep_array[i];
    return;
  }

  //***************************************************
  MasterSymGroup &MasterSymGroup::operator=(const MasterSymGroup &RHS) {
    SymGroup::operator=(RHS);
    coord_rep_ID = RHS.coord_rep_ID;
    reg_rep_ID = RHS.reg_rep_ID;
    rep_array.reserve(RHS.rep_array.size());
    for(Index i = 0; i < RHS.rep_array.size(); i++)
      add_representation(RHS.rep_array[i]->copy());

    for(Index i = 0; i < size(); i++)
      at(i).set_index(*this, i);

    return *this;
  }

  //***************************************************

  void MasterSymGroup::push_back(const SymOp &op) {
    SymGroup::push_back(op);

    back().set_index(*this, size() - 1);
    if(!size())
      std::cout << "PUSHBACK TO MASTERGROUP " << this << "\n";
    //std::cout << "Added SymOp and home_group is " << back().master_group() << " for back and " << at(0).master_group() << " for front." << std::endl;
    return;
  }

  //***************************************************

  const SymGroup &MasterSymGroup::point_group() const {
    if(!point_group_internal.size()) {
      copy_no_trans(point_group_internal, false);
    }
    return point_group_internal;
  }

  //***************************************************

  void MasterSymGroup::clear() {
    SymGroup :: clear();
    point_group_internal.clear();
    for(Index i = 0; i < rep_array.size(); i++) {
      delete rep_array[i];
    }
    rep_array.clear();

    reg_rep_ID = coord_rep_ID = -1;

    return;
  }

  //***************************************************

  Index MasterSymGroup::get_coord_rep_ID() const {
    if(coord_rep_ID == Index(-1)) add_coord_rep();
    return coord_rep_ID;
  }

  //***************************************************

  Index MasterSymGroup::get_reg_rep_ID() const {
    if(reg_rep_ID == Index(-1)) add_reg_rep();
    return reg_rep_ID;
  }

  //***************************************************

  SymGroupRep const *MasterSymGroup::get_coord_rep() const {
    if(coord_rep_ID == Index(-1)) add_coord_rep();
    return representation(coord_rep_ID);
  }

  //***************************************************

  SymGroupRep const *MasterSymGroup::get_reg_rep() const {
    if(reg_rep_ID == Index(-1)) add_reg_rep();
    return representation(reg_rep_ID);
  }

  //***************************************************

  Index MasterSymGroup::add_coord_rep() const {
    SymGroupRep *coordrep(new SymGroupRep(*this));
    Eigen::MatrixXd coordrep_mat(3, 3);
    for(Index i = 0; i < size(); i++) {
      Eigen::Matrix3d tmat = at(i).matrix();
      for(int j = 0; j < 3; j++) {
        for(int k = 0; k < 3; k++) {
          coordrep_mat(j, k) = tmat(j, k);
        }
      }
      coordrep->push_back(SymMatrixXd(coordrep_mat));
    }
    coord_rep_ID = add_representation(coordrep);
    return coord_rep_ID;
  }

  //***************************************************
  Index MasterSymGroup::add_reg_rep() const {
    SymGroupRep *regrep(new SymGroupRep(*this));
    Eigen::MatrixXd regrep_mat(size(), size());

    if(get_alt_multi_table().size() != size())
      return -1;


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
      regrep->push_back(SymMatrixXd(regrep_mat));
    }

    reg_rep_ID = add_representation(regrep);

    return reg_rep_ID;
  }

  //***************************************************

  Index MasterSymGroup::add_kronecker_rep(Index ID1, Index ID2) const {
    SymGroupRep const *rep1(representation(ID1)), *rep2(representation(ID2));
    if(!(rep1 && rep2)) return -1;
    SymGroupRep *new_rep(new SymGroupRep(*this));
    Eigen::MatrixXd tmat;
    for(Index i = 0; i < size(); i++) {
      //std::cout << "rep1 matrix:\n";
      //std::cout << *(rep1->get_MatrixXd(i)) << '\n';
      //std::cout << "rep2 matrix:\n";
      //std::cout << *(rep2->get_MatrixXd(i)) << '\n';

      kroneckerProduct(*(rep1->get_MatrixXd(i)), *(rep2->get_MatrixXd(i)), tmat);
      //std::cout << "Total matrix:\n" << tmat << '\n';
      new_rep->push_back(SymMatrixXd(tmat));
    }
    return add_representation(new_rep);
  }

  //***************************************************

  Index MasterSymGroup::add_direct_sum_rep(const Array<Index> &rep_IDs) const {
    Array<SymGroupRep const *> treps;
    for(Index i = 0; i < rep_IDs.size(); i++) {
      treps.push_back(representation(rep_IDs[i]));
      if(!treps.back())
        return -1;
    }
    SymGroupRep *new_rep(new SymGroupRep(*this));

    int dim = 0;
    for(Index i = 0; i < treps.size(); i++) {
      if(treps[i]->size() != size())
        return -1;
      if(!(treps[i]->get_MatrixXd(0)))
        return -1;

      dim += (treps[i]->get_MatrixXd(0))->cols();
    }

    Eigen::MatrixXd tmat(Eigen::MatrixXd::Zero(dim, dim));
    int corner = 0;
    for(Index i = 0; i < size(); i++) {
      corner = 0;
      for(Index j = 0; j < treps.size(); j++) {
        tmat.block(corner, corner, (treps[j]->get_MatrixXd(i))->cols(), (treps[j]->get_MatrixXd(i))->cols()) = *(treps[j]->get_MatrixXd(i));
        corner += (treps[j]->get_MatrixXd(i))->cols();
      }
      new_rep->push_back(SymMatrixXd(tmat));
    }
    return add_representation(new_rep);
  }

  //***************************************************
  Index MasterSymGroup::make_empty_representation() const {
    SymGroupRep *t_rep(new SymGroupRep(*this));
    t_rep->resize(size(), nullptr);
    rep_array.push_back(t_rep);

    return rep_array.back()->get_ID();
  }
  //***************************************************

  Index MasterSymGroup::add_representation(SymGroupRep *new_ptr) const {
    rep_array.push_back(new_ptr);
    return rep_array.back()->get_ID();
  }

  //***************************************************

  Index MasterSymGroup::add_representation(const SymGroupRep &new_rep) const {
    rep_array.push_back(new_rep.copy());
    return rep_array.back()->get_ID();
  }

  //***************************************************
  SymGroupRep const *MasterSymGroup::representation(Index i) const {
    for(Index j = 0; j < rep_array.size(); j++) {
      if(rep_array[j]->get_ID() == i)
        return rep_array[j];
    }

    return nullptr;
  }
  //***************************************************

  void MasterSymGroup::sort() {
    SymGroup::sort();
    point_group_internal.clear();
    bool broken_check(false);
    Array<Index> perm_array(size(), 0);
    for(Index i = 0; i < size(); i++) {
      perm_array[i] = at(i).index();
      if(at(i).index() != i) {
        at(i).set_index(*this, i);
        broken_check = true;
      }
    }
    if(broken_check && rep_array.size()) {
      std::cerr << "WARNING: Order of symmetry operations has been altered by MasterSymGroup::sort(). Attempting to repair "
                << rep_array.size() << " symmetry representations.\n";
      for(Index i = 0; i < rep_array.size(); i++) {
        rep_array[i]->permute(perm_array);
        for(Index j = 0; j < rep_array[i]->size(); j++) {
          if(rep_array[i]->at(j)) {
            (rep_array[i]->at(j))->set_identifiers(*this, rep_array[i]->get_ID(), perm_array[j]);
          }
        }
      }
    }
    return;
  }

  //***************************************************

  void MasterSymGroup::sort_by_class() {
    SymGroup::sort_by_class();
    point_group_internal.clear();
    bool broken_check(false);
    Array<Index> perm_array(size(), 0);
    for(Index i = 0; i < size(); i++) {
      perm_array[i] = at(i).index();
      if(at(i).index() != i) {
        at(i).set_index(*this, i);
        broken_check = true;
      }
    }
    if(broken_check && rep_array.size()) {
      std::cerr << "WARNING: Order of symmetry operations has been altered by MasterSymGroup::sort_by_class(). Attempting to repair "
                << rep_array.size() << " symmetry representations.\n";
      for(Index i = 0; i < rep_array.size(); i++) {
        rep_array[i]->permute(perm_array);
        for(Index j = 0; j < rep_array[i]->size(); j++) {
          if(rep_array[i]->at(j)) {
            (rep_array[i]->at(j))->set_identifiers(*this, rep_array[i]->get_ID(), perm_array[j]);
          }
        }
      }
    }

    return;
  }


  /*****************************************************************/

  void SymGroup::push_back(const SymOp &new_op) {
    Array<SymOp>::push_back(new_op);
    if(!size())
      std::cout << "PUSHING BACK INTO SYMGROUP " << this << "\n";
    if(back().map_error() > max_error)
      max_error = back().map_error();

    return;
  }

  //***************************************************

  void SymGroup::clear_tables() {
    multi_table.clear();
    alt_multi_table.clear();
    conjugacy_classes.clear();
    class_names.clear();
    index2conjugacy_class.clear();
    character_table.clear();
    complex_irrep.clear();
    irrep_IDs.clear();
    irrep_names.clear();

    unique_subgroups.clear();
    small_groups.clear();
    large_groups.clear();

    centralizer_table.clear();
    elem_order_table.clear();

    name.clear();
    latex_name.clear();
    comment.clear();

    return;
  }

  /*****************************************************************/

  void SymGroup::clear() {
    Array<SymOp>::clear();

    clear_tables();

    max_error = -1;

    return;
  }

  /*****************************************************************/
  // Determination of point group type
  // By counting the number of unique type of rotation symmetry and
  // total number of rotation symmetry elements, the point group type
  // of a space group is determined.

  //
  // Donghee
  /*****************************************************************/
  void SymGroup::get_rotation_groups() const {
    Array<SymOp> E, TF, THF, FF, SF;
    Array<SymOp> I, ITF, ITHF, IFF, ISF;
    for(Index i = 0; i < size(); i++) {
      auto info = at(i).info();
      if(info.op_type == SymOp::identity_op) {
        E.push_back(at(i));
      }
      if(info.op_type == SymOp::inversion_op) {
        I.push_back(at(i));
      }
      if(info.op_type == SymOp::rotation_op) {
        if(std::abs(180 - info.angle) < TOL) {
          TF.push_back(at(i));
        }
        else if(std::abs(120 - info.angle) < TOL) {
          THF.push_back(at(i));
        }
        else if(std::abs(90 - info.angle) < TOL) {
          FF.push_back(at(i));
        }
        else if(std::abs(60 - info.angle) < TOL) {
          SF.push_back(at(i));
        }
        else continue;
      }
      if(info.op_type == SymOp::rotoinversion_op || info.op_type == SymOp::mirror_op) {
        if(info.op_type == SymOp::mirror_op) {
          ITF.push_back(at(i));
        }
        else if(std::abs(120 - info.angle) < TOL) {
          ITHF.push_back(at(i));
        }
        else if(std::abs(90 - info.angle) < TOL) {
          IFF.push_back(at(i));
        }
        else if(std::abs(60 - info.angle) < TOL) {
          ISF.push_back(at(i));
        }
        else continue;
      }
    }
    rotation_groups.push_back(E);
    rotation_groups.push_back(I);
    rotation_groups.push_back(TF);
    rotation_groups.push_back(ITF);
    rotation_groups.push_back(THF);
    rotation_groups.push_back(ITHF);
    rotation_groups.push_back(FF);
    rotation_groups.push_back(IFF);
    rotation_groups.push_back(SF);
    rotation_groups.push_back(ISF);

    return;
  }
  /*****************************************************************/
  void SymGroup::get_point_group_type() const {
    if(!rotation_groups.size()) {
      get_rotation_groups();
    }
    // Calculate total number of rotation elements
    Index nm = rotation_groups[2].size() + 1;
    for(Index k = 4; k < rotation_groups.size(); k++) {
      nm += rotation_groups[k].size();
    }
    centric = true;
    if(!rotation_groups[1].size()) {
      centric = false;
    }
    // naming point gorup
    group_name.resize(2);
    group_number.resize(2);
    if((rotation_groups[4].size() + rotation_groups[5].size()) == 8) {
      crystal_system = "Cubic";
      switch(nm) {
      case 12:
        if(centric) {
          group_name[0] = "m-3";
          group_name[1] = "Th";
          group_number[0] = 200;
          group_number[1] = 206;
        }
        else {
          group_name[0] = "23";
          group_name[1] = "T (Chiral)";
          group_number[0] = 195;
          group_number[1] = 199;
        }
        break;
      case 24:
        if(centric) {
          group_name[0] = "m-3m";
          group_name[1] = "Oh";
          group_number[0] = 221;
          group_number[1] = 230;
        }
        else {
          if(rotation_groups[6].size() == 6) {
            group_name[0] = "432";
            group_name[1] = "O (Chiral)";
            group_number[0] = 207;
            group_number[1] = 214;
          }
          else if(rotation_groups[7].size() == 6) {
            group_name[0] = "-43m";
            group_name[1] = "Td";
            group_number[0] = 215;
            group_number[1] = 220;
          }
          else {
            std::cout << "\n Error Cubic case 24 Acentric \n ";
          }
        }
        break;
      default:
        std::cout << "\n Error Cubic \n";
      }
    } // for cubic;
    else if((rotation_groups[8].size() + rotation_groups[9].size()) == 2) {
      crystal_system = "Hexagonal";
      switch(nm) {
      case 6:
        if(centric) {
          group_name[0] = "6/m";
          group_name[1] = "C6h";
          group_number[0] = 175;
          group_number[1] = 176;
        }
        else {
          if(rotation_groups[8].size() == 2) {
            group_name[0] = "6";
            group_name[1] = "C6 (Chiral)";
            group_number[0] = 168;
            group_number[1] = 173;
          }
          else if(rotation_groups[9].size() == 2) {
            group_name[0] = "-6";
            group_name[1] = "C3h";
            group_number[0] = 174;
            group_number[1] = 174;
          }
          else {
            std::cout << "\n Error Hexagonal case 6 Acentric \n ";
          }
        }
        break;
      case 12:
        if(centric) {
          group_name[0] = "6/mmm";
          group_name[1] = "D6h";
          group_number[0] = 191;
          group_number[1] = 194;

        }
        else {
          if(rotation_groups[8].size() == 2) {
            if(rotation_groups[3].size() == 7) {
              group_name[0] = "622";
              group_name[1] = "D6 (Chiral)";
              group_number[0] = 177;
              group_number[1] = 182;
            }
            else if(rotation_groups[4].size() == 6) {
              group_name[0] = "6mm";
              group_name[1] = "C6v";
              group_number[0] = 183;
              group_number[1] = 186;
            }
            else std::cout << "\n Error Hexagonal case 12 Ancentric #6 \n";
          }
          else if(rotation_groups[9].size() == 2) {
            group_name[0] = "-6m2";
            group_name[1] = "D3h";
            group_number[0] = 187;
            group_number[1] = 190;
          }
          else {
            std::cout << "\n Error Hexagonal case 12 Acentric \n ";
          }
        }
        break;
      default:
        std::cout << "\n Error Hexagonal \n";
      }
    } // for hexagonal
    else if((rotation_groups[4].size() + rotation_groups[5].size()) == 2) {
      crystal_system = "Trigonal";
      switch(nm) {
      case 3:
        if(centric) {
          group_name[0] = "-3";
          group_name[1] = "S6";
          group_number[0] = 147;
          group_number[1] = 148;
        }
        else {
          group_name[0] = "3";
          group_name[1] = "C3 (Chiral)";
          group_number[0] = 143;
          group_number[1] = 146;
        }
        break;
      case 6:
        if(centric) {
          group_name[0] = "-3m";
          group_name[1] = "D3d";
          group_number[0] = 162;
          group_number[1] = 167;
        }
        else {
          if(rotation_groups[2].size() == 3) {
            group_name[0] = "32";
            group_name[1] = "D3 (Chiral)";
            group_number[0] = 149;
            group_number[1] = 155;
          }
          else if(rotation_groups[3].size() == 3) {
            group_name[0] = "3m";
            group_name[1] = "C3v";
            group_number[0] = 156;
            group_number[1] = 161;
          }
          else {
            std::cout << "\n Error Trigonal case 6 Acentric \n ";
          }
        }
        break;
      default:
        std::cout << "\n Error Trigonal \n";
      }
    } // for trigonal
    else if((rotation_groups[6].size() + rotation_groups[7].size()) == 2) {
      crystal_system = "Tetragonal";
      switch(nm) {
      case 4:
        if(centric) {
          group_name[0] = "4/m";
          group_name[1] = "C4h";
          group_number[0] = 83;
          group_number[1] = 88;
        }
        else {
          if(rotation_groups[6].size() == 2) {
            group_name[0] = "4";
            group_name[1] = "C4 (Chiral)";
            group_number[0] = 75;
            group_number[1] = 80;
          }
          else if(rotation_groups[7].size() == 2) {
            group_name[0] = "-4";
            group_name[1] = "S4";
            group_number[0] = 81;
            group_number[1] = 82;
          }
          else {
            std::cout << "\n Error Tetragonal case 4 Acentric \n ";
          }
        }

        break;
      case 8:
        if(centric) {
          group_name[0] = "4/mmm";
          group_name[1] = "D4h";
          group_number[0] = 123;
          group_number[1] = 142;
        }
        else {
          if(rotation_groups[6].size() == 2) {
            if(rotation_groups[3].size() == 5) {
              group_name[0] = "422";
              group_name[1] = "D4 (Chiral)";
              group_number[0] = 89;
              group_number[1] = 98;
            }
            else if(rotation_groups[4].size() == 4) {
              group_name[0] = "4mm";
              group_name[1] = "C4v";
              group_number[0] = 99;
              group_number[1] = 110;
            }
            else std::cout << "\n Error Tetragonal case 8 Ancentric #4 \n";
          }
          else if(rotation_groups[7].size() == 2) {
            group_name[0] = "-42m";
            group_name[1] = "D2d";
            group_number[0] = 111;
            group_number[1] = 122;
          }
          else {
            std::cout << "\n Error Tetragonal case 8 Acentric \n ";
          }
        }
        break;
      default:
        std::cout << "\n Error Tetragonal \n";
      }
    } // for tetragonal
    else if((rotation_groups[2].size() + rotation_groups[3].size()) == 3 || ((rotation_groups[2].size() + rotation_groups[3].size()) == 6  && nm == 4)) {
      crystal_system = "Orthorhomic";
      if(centric) {
        group_name[0] = "mmm";
        group_name[1] = "D2h";
        group_number[0] = 47;
        group_number[1] = 74;
      }
      else {
        if(rotation_groups[3].size() == 3) {
          group_name[0] = "222";
          group_name[1] = "D2 (Chiral)";
          group_number[0] = 16;
          group_number[1] = 24;
        }
        else if(rotation_groups[4].size() == 2) {
          group_name[0] = "mm2";
          group_name[1] = "C2v";
          group_number[0] = 25;
          group_number[1] = 46;
        }
        else {
          std::cout << "\n Error Orthorhombic Acentric \n ";
        }
      }
    } // for orthorhombic

    else if((rotation_groups[2].size() + rotation_groups[3].size()) == 1 || nm == 2) {
      crystal_system = "Monoclinic";
      if(centric) {
        group_name[0] = "2/m";
        group_name[1] = "C2h";
        group_number[0] = 10;
        group_number[1] = 15;
      }
      else {
        if(rotation_groups[3].size() == 1) {
          group_name[0] = "2";
          group_name[1] = "C2 (Chiral)";
          group_number[0] = 3;
          group_number[1] = 5;
        }
        else if(rotation_groups[4].size() == 1) {
          group_name[0] = "m";
          group_name[1] = "Cs";
          group_number[0] = 6;
          group_number[1] = 9;
        }
        else {
          std::cout << "\n Error Monoclinic Acentric \n ";
        }
      }
    } // for Acentric monoclinic
    else if(nm == 1) {
      crystal_system = "Triclinic";
      if(centric) {
        group_name[0] = "-1";
        group_name[1] = "Ci";
        group_number[0] = 2;
        group_number[1] = 2;
      }
      else {
        group_name[0] = "1";
        group_name[1] = "C1 (Chiral)";
        group_number[0] = 1;
        group_number[1] = 1;
      }
    }
    else {
      std::cout << "Error foind point group type \n";
    }
    return;
  }

  //***************************************************

  void SymGroup::print_space_group_info(std::ostream &out) const {
    if(!rotation_groups.size())
      get_rotation_groups();

    if(!crystal_system.size())
      get_point_group_type();

    std::cout << "  Crystal System : " << crystal_system << "\n";
    std::cout << "  Point Group # : " << group_number[0] << " - " << group_number[1] << "\n";
    std::cout << "  Point Group is " ;
    if(!centric) {
      std::cout << " Acentric ";
    }
    else std::cout << " Centric ";
    std::cout << " " << group_name[0] << " ( " << group_name[1] << " ) \n";
  }

  //***************************************************
  /* Loop through SymGroup and populate passed SymGroup
   * with all the operations but none of the shifts. Unle
   * explicitly requested, degenerate symmetry operations
   * will not be added.
   */
  //***************************************************
  void SymGroup::copy_no_trans(SymGroup &shiftless, bool keep_repeated) const {
    if(shiftless.size() != 0) {
      std::cerr << "WARNING in SymGroup::copy_no_trans" << std::endl;
      std::cerr << "The provided SymGroup wasn't empty and it's about to be erased. Say goodbye." << std::endl;
      shiftless.clear();
    }

    for(Index i = 0; i < size(); i++) {
      SymOp tsymop(at(i).no_trans());
      if(keep_repeated) {
        shiftless.push_back(tsymop);
      }
      else if(!shiftless.contains(tsymop)) {
        shiftless.push_back(tsymop);
      }
    }
    return;
  }

  //***************************************************

  double SymGroup::get_max_error() {

    return max_error;
  }

  //***************************************************

  void SymGroup::generate_irrep_names() const { //AAB
    irrep_names.resize(conjugacy_classes.size());

    Array<int> repeats;

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

    generate_class_names();

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

    // std::cout << "The group is: " << name << "\n\n";
    //std::cout << "\nFound these characteristics:\n";
    //std::cout << "inversion:\t" << inversion << "\n";
    //std::cout << "sigma_h:\t" << sigma_h << "\n";
    //std::cout << "sigma_v:\t" << sigma_v << "\n";
    //std::cout << "C2p:\t" << C2p << "\n";
    //std::cout << "Cubic:\t" << cubic << "\n\n";

    //std::cout << "STEP 1: Name according to representation dimension...\n";
    for(Index i = 0; i < irrep_names.size(); i++) {
      sym_wrt_paxis = true;
      if(character_table[i][0].real() == 1.0) {
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
            else if(character_table[i][j].real() < 0) {
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

      else if(character_table[i][0].real() == 2.0) {
        irrep_names[i] = "E";
      }
      else if(character_table[i][0].real() == 3.0) {
        irrep_names[i] = "T";
      }
    }

    // This checks if each representation has a unique name...
    for(Index j = 0; j < irrep_names.size() && all_unique; j++) {
      for(Index k = 0; k < irrep_names.size() && all_unique; k++) {
        all_unique = true;
        if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0)) {
          //	    std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }
    std::cout << irrep_names << std::endl;
    //std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      repeats.clear();
      for(Index j = 0; j < irrep_names.size(); j++) {
        for(Index k = 0; k < irrep_names.size(); k++) {
          if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0) && (character_table[j][0].real() == character_table[k][0].real())) {
            int dim = int(character_table[k][0].real());
            if(!(repeats.contains(dim))) {
              repeats.push_back(dim);
            }
          }
        }
      }
    }

    //    std::cout << "Repeats in ..." << repeats << "\n";
    //std::cout << "STEP 2: Name according to sigma_v or C2p symmetry...\n";
    //Check for symmetry wrt vertical mirror plane or C2' axis...
    //I don't know how to do this for non-1D representations =(


    for(Index i = 0; i < irrep_names.size(); i++) {
      int C2ind = 0;
      if((sigma_v == true) && (!all_unique) && (repeats.contains(int(character_table[0][i].real()))) && (irrep_names[i] != "E") && (!cubic)) {
        sym_wrt_v = true;
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_v; j++) {
          if(class_names[j].find("v") != std::string::npos) {
            if(class_names[j].find("'") != std::string::npos) {
              continue;
            }
            else if(character_table[i][j].real() < 0) {
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

      else if((sigma_v != true) && (C2p == true) && (!all_unique) && (repeats.contains(int(character_table[0][i].real()))) && (irrep_names[i] != "E") && (!cubic)) {
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_C2p; j++) {
          if(class_names[j].find("C2") != std::string::npos) {
            if((class_names[j].find("'") != std::string::npos) || (class_names[j].find("''") == std::string::npos)) {
              continue;
            }
            else if(character_table[i][j].real() < 0) {
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
        if(character_table[i][C2ind].real() > 0) {
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

        if(character_table[i][C2ind].real() > 0) {
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

        if(character_table[i][c2z].real() > 0) {
          irrep_names[i].append("1");
        }
        else if(character_table[i][c2y].real() > 0) {
          irrep_names[i].append("2");
        }
        else if(character_table[i][c2x].real() > 0) {
          irrep_names[i].append("3");
        }
      }
    }

    // This checks if each representation has a unique name...
    for(Index j = 0; j < irrep_names.size() && all_unique; j++) {
      for(Index k = 0; k < irrep_names.size() && all_unique; k++) {
        all_unique = true;
        if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0)) {
          //	    std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }
    //    std::cout << irrep_names << std::endl;
    //    std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      repeats.clear();
      for(Index j = 0; j < irrep_names.size(); j++) {
        for(Index k = 0; k < irrep_names.size(); k++) {
          if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0) && (character_table[j][0].real() == character_table[k][0].real())) {
            int dim = int(character_table[k][0].real());
            if(!(repeats.contains(dim))) {
              repeats.push_back(dim);
            }
          }
        }
      }
    }
    //    std::cout << "Repeats in ..." << repeats << "\n";
    //std::cout << "STEP 3: Name according to inversion symmetry...\n";
    for(Index i = 0; i < irrep_names.size(); i++) {
      if((inversion == true) && (!all_unique) && (repeats.contains(int(character_table[0][i].real())))) {
        sym_wrt_inv = true;
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_inv; j++) {
          if((class_names[j].find("i") != std::string::npos) && (character_table[i][j].real() < 0)) {
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
          //	    std::cout << irrep_names[j] << " = " << irrep_names[k] << std::endl;
          all_unique = false;
        }
      }
    }
    //    std::cout << irrep_names << std::endl;
    //std::cout << "Check if all of the current names are unique... " << all_unique << "\n";

    if(!all_unique) {
      repeats.clear();
      for(Index j = 0; j < irrep_names.size(); j++) {
        for(Index k = 0; k < irrep_names.size(); k++) {
          if((j != k) && (irrep_names[j].compare(irrep_names[k]) == 0) && (character_table[j][0].real() == character_table[k][0].real())) {
            int dim = int(character_table[k][0].real());
            if(!(repeats.contains(dim))) {
              repeats.push_back(dim);
            }
          }
        }
      }
    }

    //    std::cout << "Repeats in ..." << repeats << "\n";

    //    std::cout << "STEP 4: Name according to sigma_h symmetry...\n";
    for(Index i = 0; i < irrep_names.size(); i++) {
      if((sigma_h == true) && (!all_unique) && (repeats.contains(int(character_table[0][i].real())))) {
        for(Index j = 0; j < conjugacy_classes.size() && sym_wrt_h; j++) {
          if((class_names[j].find("h") != std::string::npos) && (character_table[i][j].real() < 0)) {
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

    //    std::cout << "At the very end, check if everything is unique...\n";

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
      std::cout << "Failed to name them all uniquely... =( \n";
    }

    std::cout << irrep_names << std::endl;

    return;
  }



  //***************************************************
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

  void SymGroup::generate_class_names() const { //AAB
    if(!conjugacy_classes.size()) {
      calc_conjugacy_classes();
    }

    class_names.resize(conjugacy_classes.size());

    Array<Eigen::Vector3d > highsym_axes;
    Array<int> mult;
    double angle = 360;
    Eigen::Vector3d xprodvec;
    double dprod;
    std::string symtype;
    Index to_name = conjugacy_classes.size();
    bool cubic = false;

    if((name == "T") || (name == "Td") || (name == "Th") || (name == "O") || (name == "Oh")) {
      cubic = true;
    }

    class_names[0] = "E";
    to_name--;

    std::vector<SymOp::SymInfo> info;
    for(Index i = 0; i < size(); i++)
      info.push_back(at(i).info());

    //If the SymGroup includes an inversion operation, name it "i"
    //We can just use back() here because inversion is always last
    if(info.back().op_type == SymOp::inversion_op) {
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
      if((info[i].op_type == SymOp::rotation_op) || (info[i].op_type == SymOp::screw_op)) {
        if(highsym_axes.contains(info[i].axis)) { //Otherwise, check if the axis has been found;
          mult[highsym_axes.find(info[i].axis)]++;
        }
        else {
          highsym_axes.push_back(info[i].axis);
          mult.push_back(1);
        }
      }
    }

    Eigen::Vector3d hs_axis;
    double hangle = 360;

    for(Index i = 0; i < size(); i++) {
      if((info[i].angle < hangle) && (info[i].angle > TOL)) {
        hangle = info[i].angle;
        hs_axis = info[i].axis;
      }
    }


    int order = mult.max();

    for(int i = (int(mult.size()) - 1); i >= 0; i--) {
      if((cubic == false) && (mult[i] != order)) {
        highsym_axes.remove(i);
        mult.remove(i);
      }
      else if(mult[i] < (order - 1)) {
        highsym_axes.remove(i);
        mult.remove(i);
      }
    }

    //This is kind of a hack... We might want to change this so that the principal axis
    //is actually the one with the highest fold rotation, and then treat cases like
    //cubic and orthorhombic separately....... but I am not sure that's better...

    if(name == "D3d") {
      for(int i = int(mult.size()) - 1; i >= 0; i--) {
        if(!almost_zero(highsym_axes[i] - hs_axis)) {
          highsym_axes.remove(i);
          mult.remove(i);
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
        for(Index j = 0; j < highsym_axes.size() && !normal; j++) {
          dprod = highsym_axes[j].dot(info[i].axis);
          xprodvec = highsym_axes[j].cross(info[i].axis);
          if(almost_zero(xprodvec.norm())) {
            normal = true;
          }
        }

        if(almost_zero(xprodvec.norm())) { //Check if the cross product with principal axis is zero
          if((info[i].angle < 200) && (info[i].angle > 1)) { //Only bother with angles that 360 is divisible by
            if((info[i].op_type == SymOp::rotation_op) || (info[i].op_type == SymOp::screw_op)) {
              angle = info[i].angle;
              symtype = "C";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle);
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotation\t==> " << s.str() << std::endl;
            }
            else if(info[i].op_type == SymOp::rotoinversion_op) {
              angle = info[i].angle;
              symtype = "S";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle);
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotoinversion\t==> " << s.str() << std::endl;
            }
          }
          else if((info[i].op_type == SymOp::mirror_op) || (info[i].op_type == SymOp::glide_op)) {
            symtype = "h";
            s << conjugacy_classes[ind].size() << symtype;
            class_names[ind] = s.str();
            to_name--;
            //std::cout << "  mirror\t==> " << s.str() << std::endl;
          }
        }
        else {
          if((info[i].angle < 200) && (info[i].angle > 1)) {
            if((info[i].op_type == SymOp::rotation_op) || (info[i].op_type == SymOp::screw_op)) {
              angle = info[i].angle;
              symtype = "C";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle) << cprime;
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotation\t==> " << s.str() << std::endl;
              cprime = "''";
            }
            else if(info[i].op_type == SymOp::rotoinversion_op) {
              angle = info[i].angle;
              symtype = "S";
              s << conjugacy_classes[ind].size() << symtype << int(360 / angle) << sprime;
              class_names[ind] = s.str();
              to_name--;
              //std::cout << "  rotoinversion\t==> " << s.str() << std::endl;
              sprime = "''";
            }
          }
          else if((info[i].op_type == SymOp::mirror_op) || (info[i].op_type == SymOp::glide_op)) {
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
        if(info[conjugacy_classes[i][0]].op_type == SymOp::rotation_op) {
          if(almost_zero(info[conjugacy_classes[i][0]].axis - Eigen::Vector3d::UnitX())) {
            class_names[i].append("(x)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis -  Eigen::Vector3d::UnitY())) {
            class_names[i].append("(y)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis -  Eigen::Vector3d::UnitY())) {
            class_names[i].append("(z)");
          }
        }
        else if(info[conjugacy_classes[i][0]].op_type == SymOp::mirror_op) {
          if(almost_zero(info[conjugacy_classes[i][0]].axis -  Eigen::Vector3d::UnitX())) {
            class_names[i].append("(yz)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis -  Eigen::Vector3d::UnitY())) {
            class_names[i].append("(xz)");
          }
          else if(almost_zero(info[conjugacy_classes[i][0]].axis -  Eigen::Vector3d::UnitZ())) {
            class_names[i].append("(xy)");
          }
        }
      }
    }


    return;
  }


  //***************************************************
  void SymGroup::print_character_table(std::ostream &stream) {
    generate_class_names();
    generate_irrep_names();

    stream.precision(3);

    for(Index i = 0; i < character_table.size(); i++) {
      for(Index j = 0; j < character_table.size(); j++) {
        //This part cleans up the numbers in the table...
        if(almost_zero(std::abs(character_table[i][j].real()) - 0.5)) {
          if(character_table[i][j].real() > 0) {
            character_table[i][j] = std::complex<double>(0.5, character_table[i][j].imag());
          }
          else {
            character_table[i][j] = std::complex<double>(-0.5, character_table[i][j].imag());
          }
        }
        if(almost_zero(std::abs(character_table[i][j].imag()) - 0.5)) {
          if(character_table[i][j].imag() > 0) {
            character_table[i][j] = std::complex<double>(character_table[i][j].real(), 0.5);
          }
          else {
            character_table[i][j] = std::complex<double>(character_table[i][j].real(), -0.5);
          }
        }
        if(almost_zero(std::abs(character_table[i][j].real()) - 0.866)) {
          if(character_table[i][j].real() > 0) {
            character_table[i][j] = std::complex<double>(0.866, character_table[i][j].imag());
          }
          else {
            character_table[i][j] = std::complex<double>(-0.866, character_table[i][j].imag());
          }
        }
        if(almost_zero(std::abs(character_table[i][j].imag()) - 0.866)) {
          if(character_table[i][j].imag() > 0) {
            character_table[i][j] = std::complex<double>(character_table[i][j].real(), 0.866);
          }
          else {
            character_table[i][j] = std::complex<double>(character_table[i][j].real(), -0.866);
          }
        }
        if(almost_zero(character_table[i][j].real())) {
          character_table[i][j] = std::complex<double>(0, character_table[i][j].imag());
        }
        if(almost_zero(character_table[i][j].imag())) {
          character_table[i][j] = std::complex<double>(character_table[i][j].real(), 0);
        }
      }
    }

    // This will check if there are any non-integer values in the table.
    bool all_int = true;
    for(Index i = 0; i < character_table.size(); i++) {
      for(Index j = 0; j < character_table.size(); j++) {
        if(!almost_zero(character_table[i][j].real() - floor(character_table[i][j].real()))) {
          all_int = false;
        }
        if(!almost_zero(character_table[i][j].imag() - floor(character_table[i][j].imag()))) {
          all_int = false;
        }
      }
    }

    //This loop is for printing integer-only tables.
    if(all_int) {
      stream << " ";
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << "\n| Group: " << name << " (" << comment << ")";
      int space = (7 * (character_table.size() + 1)) - (name.size() + comment.size() + 11);
      for(int i = 0; i < space; i++) {
        stream << " ";
      }

      stream << "|\n ";
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << std::endl;

      stream << "|       |";
      for(Index i = 0; i < character_table.size(); i++) {
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
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << std::endl;

      for(Index i = 0; i < character_table.size(); i++) {
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

        for(Index j = 0; j < character_table.size(); j++) {
          //This will make sure the positive and negative numbers align
          if(character_table[i][j].real() > 0) {
            stream << "   ";
          }
          else if(character_table[i][j].real() < 0) {
            stream << "  ";
          }
          else {
            if(character_table[i][j].imag() >= 0) {
              stream << "   ";
            }
            else {
              stream << "  ";
            }
          } //End alignment of negative numbers

          if((character_table[i][j].imag() > 0) && (character_table[i][j].real() > 0)) {
            if(character_table[i][j].imag() == 1) {
              stream << int(character_table[i][j].real()) <<  "i" << "  |";
            }
            else {
              stream << int(character_table[i][j].real()) << "+" << int(character_table[i][j].imag())  << "i" << "  |";
            }
          }
          else if((character_table[i][j].imag() < 0) && (character_table[i][j].real() > 0)) {
            if(character_table[i][j].imag() == -1) {
              stream << int(character_table[i][j].real()) << "-"  << "i" << "  |";
            }
            else {
              stream << int(character_table[i][j].real()) <<  int(character_table[i][j].imag())  << "i" << "  |";
            }
          }
          else if(character_table[i][j].imag() == 0) {
            stream << int(character_table[i][j].real()) << "  |";
          }
          else if(character_table[i][j].real() == 0) {
            if(character_table[i][j].imag() == 1) {
              stream  << "i" << "  |";
            }
            else if(character_table[i][j].imag() == -1) {
              stream <<  "-"  << "i" << "  |";
            }
            else {
              stream << int(character_table[i][j].imag()) << "i" << "  |";
            }
          }
        } //Closes loop over j
        stream << std::endl;
      } //Closes loop over i

      //Prints the line at the bottom of the table
      stream << " ";
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "-------";
      }
      stream << std::endl;
    }//Closes the if all_int loop

    //If not everything is an integer, we use this loop instead...
    else {
      stream << " ";
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "----------------";
      }

      stream << "\n| Group: " << name << " (" << comment << ")";
      int space = (16 * (character_table.size() + 1)) - (name.size() + comment.size() + 12);
      for(int i = 0; i < space; i++) {
        stream << " ";
      }

      stream << "|\n ";

      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "----------------";
      }
      stream << std::endl;

      stream << "|               |";
      for(Index i = 0; i < character_table.size(); i++) {
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
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "----------------";
      }
      stream << std::endl;

      for(Index i = 0; i < character_table.size(); i++) {
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

        for(Index j = 0; j < character_table.size(); j++) {
          //This will make sure the positive and negative numbers align
          if(character_table[i][j].real() > 0) {
            stream << "   ";
          }
          else if(character_table[i][j].real() < 0) {
            stream << "  ";
          }
          else {
            if(character_table[i][j].imag() >= 0) {
              stream << "   ";
            }
            else {
              stream << "  ";
            }
          } //End alignment of negative numbers

          if((character_table[i][j].imag() > 0) && (character_table[i][j].real() > 0)) { //Long
            if(character_table[i][j].imag() == 1) {
              stream << character_table[i][j].real() <<  "i" << "  |";
            }
            else {
              stream << character_table[i][j].real() << "+" << character_table[i][j].imag()  << "i" << "  |";
            }
          }
          else if((character_table[i][j].imag() < 0) && (character_table[i][j].real() > 0)) { //Long
            if(character_table[i][j].imag() == -1) {
              stream << character_table[i][j].real() << "-"  << "i" << "  |";
            }
            else {
              stream << character_table[i][j].real() << character_table[i][j].imag()  << "i" << "  |";
            }
          }
          else if(character_table[i][j].imag() == 0) { //Short
            stream << int(character_table[i][j].real()) << "           |";
          }
          else if(character_table[i][j].real() == 0) { //Short
            if(character_table[i][j].imag() == 1) {
              stream  << "i" << "           |";
            }
            else if(character_table[i][j].imag() == -1) {
              stream <<  "-"  << "i" << "           |";
            }
            else {
              stream << character_table[i][j].imag() << "i" << "          |";
            }
          }
          else if((character_table[i][j].imag() > 0) && (character_table[i][j].real() < 0)) { //Long
            stream << character_table[i][j].real() << "+" << character_table[i][j].imag()  << "i" << "  |";
          }
          else { //Long
            stream << character_table[i][j].real() << character_table[i][j].imag()  << "i" << "  |";
          }

        }//Closes the j loop

        //stream << "\n|-";
        //for(int x=0; x<character_table.size(); x++){
        //  stream << "--------------|";
        //}

        stream << "\n";
      } //Closes the i loop
      stream << " ";
      for(Index i = 0; i < character_table.size() + 1; i++) {
        stream << "----------------";
      }
      stream << "\n";
    }//Closes else loop (non-integers)

    return;
  }

  //***************************************************

  /** The centralizer of a subset S of a group G is the set of elements
   *  of G that commute with each element of S. The centralizer is a
   *  subgroup of G.
   *
   *  In our specific case, the subsets are the conjugacy classes, so
   *  we fill the centralizer_table with subgroups of G that commute
   *  with the corresponding conjugacy class.qq
   */

  void SymGroup::calc_centralizers() const { //AAB
    if(!conjugacy_classes.size()) {
      calc_conjugacy_classes();
    }

    bool all_commute = true;

    centralizer_table.resize(conjugacy_classes.size());

    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      all_commute = true;
      for(Index j = 0; j < conjugacy_classes[i].size() && all_commute; j++) {
        for(Index k = 0; k < multi_table.size(); k++) {
          int ind = conjugacy_classes[i][j];
          if((multi_table[ind][k] == multi_table[k][ind]) && (!centralizer_table[i].contains(k))) {
            centralizer_table[i].push_back(k);
          }
          else {
            all_commute = false;
          }
        }
      }
    }
  }

  //***************************************************

  void SymGroup::calc_character_table()const { //AAB
    //std::cout << "Calculating character table of SymGroup " << this << '\n';
    if(!size()) {
      std::cout << "This is an empty group. It has no character.\n";
      name = "EMPTY";
      return;
    }

    int sigma_h_ind = 0;
    Index d1 = 0, d2 = 0, d3 = 0;
    Index order1 = 0, order2 = 0;
    bool mirror(false);
    bool inversion(false);
    SymGroup subgroup;

    std::vector<SymOp::SymInfo> info;
    for(Index i = 0; i < size(); i++)
      info.push_back(at(i).info());

    if(!size() || get_multi_table().size() != size() || !valid_index(get_multi_table()[0][0])) {
      std::cout << "This is not a group!!!\n";
      name = "NG";
      return;
    }

    calc_conjugacy_classes();
    //calc_conjugacy_corr_table();
    calc_centralizers();

    Index h = multi_table.size(); //This is the dimensionality of the group.
    Index nc = conjugacy_classes.size(); //This is the number of conjugacy classes, which is also the number of irreducible representations.

    character_table.resize(nc, Array<std::complex<double> >(nc, -7));
    complex_irrep.resize(nc, false);
    irrep_IDs.resize(nc, -1);


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

    //std::cout << "d1 = " << d1 << "\t" << "d2 = " << d2 << "\t" << "d3 = " << d3 << "\n\n";

    // Set characters related to identity
    for(Index i = 0; i < character_table.size(); i++) {
      character_table[0][i] = 1;
      if(i < d1) {
        character_table[i][0] = 1;
      }
      else if(i >= d1 && i < (d1 + d2)) {
        character_table[i][0] = 2;
      }
      else {
        character_table[i][0] = 3;
      }
    }

    // count over all possible 1d representations
    Counter<Array<int> > count(Array<int>(nc, -1), Array<int>(nc, 1), Array<int>(nc, 2));
    int sum = 0;
    order1 = 1;
    std::complex<double> ortho = 0;
    int telem = 0, elem1 = 0, elem2 = 0;
    int try1 = 0, try2 = 0;

    int try1_ind = 0, try2_ind = 0;
    Array<bool> check(nc, false);


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

          Array<std::complex<double> > ortharray;

          for(Index i = 0; i < nc; i++) {
            ortho = 0.0;
            ortharray.resize(0);
            for(Index j = 0; j < nc; j++) {
              if((character_table[i][1].real() != -7) && (character_table[i][1].imag() != 0)) {
                std::complex<double> temp = std::complex<double>(double(count[j] * conjugacy_classes[j].size()), 0.0);
                ortho += (temp * character_table[i][j]);
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
                character_table[order1][i] = count[i];
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
        if(info[i].op_type == SymOp::inversion_op) {
          inv = true;
        }
        else if(info[i].op_type == SymOp::mirror_op || info[i].op_type == SymOp::glide_op) {
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
      calc_centralizers();
    }

    Array<std::complex<double> > centralizers(nc, 0);

    // std::cout << "----------------------------------------------------------------------------\n";
    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      centralizers[i] = centralizer_table[i].size();
      //   std::cout << centralizer_table[i].size() << "\t";
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
        if(info[i].op_type == SymOp::inversion_op) {
          inversion = true;
        }
        else if((info[i].op_type == SymOp::mirror_op) || (info[i].op_type == SymOp::glide_op)) {
          //if(conjugacy_classes[conjugacy_corr_table[0][i]].size() == 1) {
          if(conjugacy_classes[index2conjugacy_class[i]].size() == 1) {
            sigma_h_ind = i;
            mirror = true;

          }
        }
      }
    }
    /*
    if((inversion == true) && ((mirror == true) && (sigma_h_ind !=0))){
      std::cout << "This group has INVERSION AND MIRROR symmetry.\n";
    }
    else if(inversion == true){
      std::cout << "This group has INVERSION SYMMETRY.\n";
    }
    else if((mirror == true) && (sigma_h_ind != 0)){
      std::cout << "This group has MIRROR symmetry.\n";
    }
    else{
      std::cout << "This group does NOT have any special symmetries.\n";
    }
    */

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
      int generator;

      for(Index i = 0; i < size(); i++) {
        if((info[i].angle < angle) && (info[i].angle > TOL)) {
          angle = info[i].angle;
          generator = i;
        }
      }

      /** We need to keep in mind that the angle returned by get_rotation_angle() is
       *  real and is given in degrees. We need to convert it to complex radian form,
       *  which we will do using Euler's formula: e^(ix) = cos(x) + i*sin(x).
       */

      angle = angle * (M_PI / 180);

      std::complex<double> iangle(cos(angle), sin(angle));

      Array<std::complex<double> > tchar;
      Array<std::complex<double> > tcharconj;
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
        character_table[order1] = tchar;
        complex_irrep[order1] = true;
        order1++;
        character_table[order1] = tcharconj;
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
          if(info[g].op_type == SymOp::rotoinversion_op) {
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

      subgroup.calc_character_table();
      subgroup.generate_class_names();
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

            for(Index k = 0; k < subgroup.character_table.size(); k++) { //rows
              character_table[k][ind_j] = subgroup.character_table[k][ind_i];
              character_table[subgroup.character_table.size() + k][ind_j] = subgroup.character_table[k][ind_i];
            }
          }
          //This loop is for the inverted SymOps
          if(almost_equal(subgroup[i].matrix()*back().matrix(), at(j).matrix(), TOL) && inversion) {
            //Index ind_i = subgroup.conjugacy_corr_table[0][i];
            //Index ind_j = conjugacy_corr_table[0][j];
            Index ind_i = subgroup.index2conjugacy_class[i];
            Index ind_j = index2conjugacy_class[j];
            for(Index k = 0; k < subgroup.character_table.size(); k++) {
              character_table[k][ind_j] = subgroup.character_table[k][ind_i];
              character_table[subgroup.character_table.size() + k][ind_j] = -subgroup.character_table[k][ind_i];
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
            for(Index k = 0; k < subgroup.character_table.size(); k++) {
              character_table[k][ind_j] = subgroup.character_table[k][ind_i];
              character_table[subgroup.character_table.size() + k][ind_j] = -subgroup.character_table[k][ind_i];
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
      Array<std::complex<double> > centrcheck;
      centrcheck.resize(nc);

      for(Index i = 0; i < d1; i++) {
        for(Index j = 0; j < character_table.size(); j++) {
          centrcheck[j] += (character_table[i][j] * character_table[i][j]);
        }
      }

      /*
      std::cout << "----------------------------------------------------------------------------\n";
      for(Index i=0; i<nc; i++){
      std::cout << centrcheck[i] << "\t";
      }
      std::cout << std::endl;

      std::cout << "----------------------------------------------------------------------------\n";
      */

      if(d2 == 1) {
        for(Index i = 1; i < nc; i++) {
          if(abs(centrcheck[i] - centralizers[i]) < TOL) {
            character_table[nc - 1][i] = std::complex<double>(0, 0);
          }
          else {
            //This is a bit of a hack, but it turns out to always be true...
            character_table[nc - 1][i] = -sqrt(abs(centrcheck[i] - centralizers[i]));
          }
        }

        bool mir = false;
        int mircount = 0;

        for(Index i = 0; i < size(); i++) {
          if(info[i].op_type == SymOp::mirror_op || info[i].op_type == SymOp::glide_op) {
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
        Array<std::complex<double> > remainder;
        int nr = 0;
        remainder.resize(nc);
        for(Index i = 1; i < nc; i++) {
          if(abs(centrcheck[i] - centralizers[i]) < TOL) {
            for(Index j = 0; j < d2; j++) {
              character_table[(nc - j) - 1][i] = std::complex<double>(0, 0);
              remainder[i] = std::complex<double>(0, 0);
            }
          }
          else {
            remainder[i] = sqrt(abs(centralizers[i] - centrcheck[i]) / 2);
            nr++;
          }
        }

        Counter<Array<int> > signcount(Array<int>(nr, -1), Array<int>(nr, 1), Array<int> (nr, 2));
        Array<std::complex<double> > trow;
        order2 = (nc - d2);
        Array<Array<std::complex<double> > > tset;

        do {
          trow.resize(nc);
          int j = 0;
          for(Index i = 0; i < nc; i++) {
            if(character_table[nc - 1][i].real() != -7) {
              trow[i] = character_table[nc - 1][i];
            }
            if(character_table[order2][i].real() == -7) {
              trow[i] = signcount[j] * remainder[i].real();
              j++;
            }
          }

          int sum = 0;

          for(Index i = 0; i < nc; i++) {
            sum += trow[i].real() * conjugacy_classes[i].size();
          }
          double orthcheck;
          Array<double> ortharray;
          if(sum == 0) {
            ortharray.resize(0);
            for(Index i = 0; i < d1; i++) {
              orthcheck = 0.0;
              for(Index j = 0; j < nc; j++) {
                double temp = (trow[j].real()) * (conjugacy_classes[j].size());
                orthcheck += (temp * character_table[i][j].real());
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
            character_table[order2 + i] = tset[i];
          }
        }

        bool mir = false;
        for(Index i = 0; i < size(); i++) {
          if(info[i].op_type == SymOp::mirror_op || info[i].op_type == SymOp::glide_op) {
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
            character_table[1][j] = std::complex<double>(1, 0);
            character_table[2][j] = std::complex<double>(1, 0);
            character_table[3][j] = std::complex<double>(-1, 0);
          }
          else {
            character_table[1][j] = omega;
            character_table[2][j] = omega * omega;
            character_table[3][j] = std::complex<double>(0, 0);
            omega *= omega;
          }
        }
      }

      bool sigma_d = false;

      for(Index i = 0; i < size(); i++) {
        if((info[i].op_type == SymOp::mirror_op) || (info[i].op_type == SymOp::glide_op)) {
          sigma_d = true;
        }
      }

      if(sigma_d == true) {
        name = "Td";
        latex_name = "T_d";
        comment = "#215-220";

        for(Index j = 1; j < nc; j++) { //Loop over columns
          if(conjugacy_classes[j].size() == 3) {
            character_table[1][j] = std::complex<double>(1, 0);
            character_table[2][j] = std::complex<double>(2, 0);
            character_table[3][j] = std::complex<double>(-1, 0);
            character_table[4][j] = std::complex<double>(-1, 0);
          }
          else if(conjugacy_classes[j].size() == 8) {
            character_table[1][j] = std::complex<double>(1, 0);
            character_table[2][j] = std::complex<double>(-1, 0);
            character_table[3][j] = std::complex<double>(0, 0);
            character_table[4][j] = std::complex<double>(0, 0);
          }
          else if((conjugacy_classes[j].size() == 6) && (info[conjugacy_classes[j][0]].angle > 10)) {
            character_table[1][j] = std::complex<double>(-1, 0);
            character_table[2][j] = std::complex<double>(0, 0);
            character_table[3][j] = std::complex<double>(1, 0);
            character_table[4][j] = std::complex<double>(-1, 0);
          }
          else {
            character_table[1][j] = std::complex<double>(-1, 0);
            character_table[2][j] = std::complex<double>(0, 0);
            character_table[3][j] = std::complex<double>(-1, 0);
            character_table[4][j] = std::complex<double>(1, 0);
          }
        }
      }

      else if((sigma_d == false) && (nc == 5)) {
        name = "O";
        latex_name = "O";
        comment = "#207-214";
        for(Index j = 1; j < nc; j++) { //Loop over columns
          if(conjugacy_classes[j].size() == 3) {
            character_table[1][j] = std::complex<double>(1, 0);
            character_table[2][j] = std::complex<double>(2, 0);
            character_table[3][j] = std::complex<double>(-1, 0);
            character_table[4][j] = std::complex<double>(-1, 0);
          }
          else if(conjugacy_classes[j].size() == 8) {
            character_table[1][j] = std::complex<double>(1, 0);
            character_table[2][j] = std::complex<double>(-1, 0);
            character_table[3][j] = std::complex<double>(0, 0);
            character_table[4][j] = std::complex<double>(0, 0);
          }
          else if((conjugacy_classes[j].size() == 6) && (almost_equal(info[conjugacy_classes[j][0]].angle, 180.0))) {
            character_table[1][j] = std::complex<double>(-1, 0);
            character_table[2][j] = std::complex<double>(0, 0);
            character_table[3][j] = std::complex<double>(-1, 0);
            character_table[4][j] = std::complex<double>(1, 0);
          }
          else {
            character_table[1][j] = std::complex<double>(-1, 0);
            character_table[2][j] = std::complex<double>(0, 0);
            character_table[3][j] = std::complex<double>(1, 0);
            character_table[4][j] = std::complex<double>(-1, 0);
          }
        }
      }
    }

    return;
  }

  //***************************************************
  void SymGroup::calc_elem_order_table() const {

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



  //***************************************************
  void SymGroup::calc_small_subgroups() const { //AAB

    int tempind = 0;
    Index i, j, k, l;
    if(get_alt_multi_table().size() != size()) {
      return;
    }
    small_groups.clear();

    for(i = 0; i < multi_table.size(); i++) {
      Array<Index> tgroup;
      tgroup.push_back(0);
      j = i;
      while(j != 0) {
        if(!tgroup.contains(j)) {
          tgroup.push_back(j);
        }
        if(!tgroup.contains(alt_multi_table[j][0])) {
          tgroup.push_back(alt_multi_table[j][0]);
        }
        j = multi_table[i][j];
      }
      //if(tgroup.size() == 1) continue;
      for(j = 0; j < small_groups.size(); j++) {
        if(small_groups[j][0].size() == tgroup.size()
           && which_unique_combination(tgroup, small_groups[j]) != small_groups[j].size()) {
          break;
        }
      }
      //std::cout << " j is " << j << " and small_groups.size() is " << small_groups.size() << std::endl;
      if(j == small_groups.size()) { // add the new group
        //std::cout << "tgroup " << tgroup << " is unique, so we add it " << std::endl;
        Array<Index> tconj(tgroup);
        small_groups.push_back(Array< Array<Index> >());
        for(k = 0; k < alt_multi_table.size(); k++) {
          for(l = 0; l < tgroup.size(); l++) {
            tempind = alt_multi_table[k][tgroup[l]];
            tconj[l] = multi_table[tempind][k];
          }
          if(which_unique_combination(tconj, small_groups.back()) == small_groups.back().size()) {
            small_groups.back().push_back(tconj);
          }

        }
      }
    }

    /*
    for(i = 0; i < small_groups.size(); i++) {
      std::cout << "Equivalent small_groups of class " << i + 1 << ":\n";
      for(j = 0; j < small_groups[i].size(); j++) {
        std::cout << small_groups[i][j]
                  << "  has inverse: ";
        for(k = 0; k < small_groups[i][j].size(); k++)
          std::cout << alt_multi_table[small_groups[i][j][k]][0] << "  ";
        std::cout << "\n";
      }
      std::cout << "\n";
    }
    */
  }

  //***************************************************
  /*  Start with large_groups = small_groups, then add new subgroups by
   *  finding the closure of a union of a large_group and a
   *  small_group. If the the new large_group is unique, add it as a
   *  large_group. Repeat for all (large_group, small_group) pairs,
   *  until no new large_groups are found. This is probably not the
   *  fastest algorithm, but it is complete
   */

  void SymGroup::calc_large_subgroups() const {
    if(!small_groups.size()) {
      calc_small_subgroups();
    }
    Index i, j, k, l, m, jj, iii, jjj;

    large_groups = small_groups;


    for(i = 0; i < large_groups.size(); i++) {
      //std::cout << "i is " << i << " and large_groups.size() is " << large_groups.size() << std::endl;
      for(j = 0; j < small_groups.size(); j++) {
        //for(ii=0; ii<large_groups[i].size(); ii++){
        for(jj = 0; jj < small_groups[j].size(); jj++) {
          Array<Index> tgroup1(large_groups[i][0]);

          // append unique elements from small_groups[j][jj]
          for(jjj = 0; jjj < small_groups[j][jj].size(); jjj++) {
            if(!tgroup1.contains(small_groups[j][jj][jjj])) {
              tgroup1.push_back(small_groups[j][jj][jjj]);
            }
          }

          // find group closure
          for(iii = 0; iii < tgroup1.size(); iii++) {
            for(jjj = 0; jjj < tgroup1.size(); jjj++) {
              if(!tgroup1.contains(multi_table[tgroup1[iii]][tgroup1[jjj]])) {
                tgroup1.push_back(multi_table[tgroup1[iii]][tgroup1[jjj]]);
              }
            }
          }

          for(k = 0; k < large_groups.size(); k++) {
            if(large_groups[k][0].size() == tgroup1.size()
               && which_unique_combination(tgroup1, large_groups[k]) != large_groups[k].size()) {
              break;
            }
          }
          //std::cout << " k is " << k << " and large_groups.size() is " << large_groups.size() << std::endl;
          if(k < large_groups.size())
            continue;
          // add the new group
          //std::cout << "tgroup1 " << tgroup1 << " is unique, so we add it " << std::endl;
          //std::cout << "tgroup2 " << tgroup2 << " yields all_in = " << tgroup1.all_in(tgroup2) << std::endl;
          Array<Index> tconj(tgroup1);
          large_groups.push_back(Array< Array<Index> >());
          for(l = 0; l < alt_multi_table.size(); l++) {
            for(m = 0; m < tgroup1.size(); m++) {
              iii = alt_multi_table[l][tgroup1[m]];
              tconj[m] = multi_table[iii][l];
            }
            if(which_unique_combination(tconj, large_groups.back()) == large_groups.back().size()) {
              tconj.sort();
              large_groups.back().push_back(tconj);
            }

          }

        }
      }
    }

    // Sort subgroups by number of elements and multiplicity
    for(Index i = 0; i < large_groups.size(); i++) {
      for(Index j = i + 1; j < large_groups.size(); j++) {
        if(large_groups[i][0].size() < large_groups[j][0].size())
          large_groups.swap_elem(i, j);
        if(large_groups[i][0].size() == large_groups[j][0].size()
           && large_groups[i].size() < large_groups[j].size())
          large_groups.swap_elem(i, j);
      }
    }
    return;
  }

  //**********************************************


  void SymGroup::make_unique_subgroups() {
    if(!large_groups.size()) calc_large_subgroups();

    Array<std::string> sg_names, sg_names_limited;
    Array<bool> chosen_flag(large_groups.size(), false);
    for(Index i = 0; i < large_groups.size(); i++) {
      SymGroup sgroup;
      for(Index j = 0; j < large_groups[i][0].size(); j++) {
        sgroup.push_back(at(large_groups[i][0][j]));
      }
      sgroup.calc_character_table();
      sg_names.push_back(sgroup.name);
      std::cout << sgroup.name << "-" << i << " has equivalencies:\n";
      for(Index j = 0; j < large_groups[i].size(); j++) {
        std::cout << "  " << large_groups[i][j] << "\n";
      }
      std::cout << "\n";
    }

    Array< Array< Index > > sg_tree(large_groups.size(), Array<Index>());
    for(Index i = 0; i < large_groups.size(); i++) {
      std::cout << "Subgroup " << sg_names[i] << "-" << i << " is also a subgroup of ";
      for(Index j = 0; j < large_groups.size(); j++) {
        for(Index jj = 0; jj < large_groups[j].size(); jj++) {
          if(large_groups[i][0].all_in(large_groups[j][jj])) {
            sg_tree[i].push_back(j);
            std::cout << sg_names[j] << "-" << j << "-" << jj << "  ";
            break;
          }
        }
      }
      std::cout << "\n";
    }

    //attempt to maximize coincidence a
    int max_co, t_co;
    for(int i = int(large_groups.size()) - 1; i >= 0; i--) {
      if(chosen_flag[i]) continue;
      //chosen_flag[i]=true;
      for(int j = 0; j < i; j++) {
        if(sg_names[j] == sg_names[i]) {
          chosen_flag[j] = true;
          continue;
        }
        max_co = 0;
        for(Index jj = 0; jj < large_groups[j].size(); jj++) {
          t_co = large_groups[i][0].coincidence(large_groups[j][jj]);
          if(t_co > max_co) {
            max_co = t_co;
            large_groups[j].swap_elem(jj, 0);
          }
        }
      }
    }

    unique_subgroups.clear();
    for(Index i = 0; i < large_groups.size(); i++) {
      if(chosen_flag[i]) continue;
      unique_subgroups.push_back(SymGroup());
      for(Index ii = 0; ii < large_groups[i][0].size(); ii++) {
        unique_subgroups.back().push_back(at(large_groups[i][0][ii]));
      }
      unique_subgroups.back().sort();
      unique_subgroups.back().calc_character_table();

      std::cout << "Added group " << unique_subgroups.back().name  << " having bit-string " << large_groups[i][0] << std::endl;
    }

    return;
  }

  //***************************************************
  /** The number of irreducible representations is equal to
   *  the number of conjugacy classes.
   *
   *  The number of elements in each conjugacy class must be
   *  a divisor of the number of symmetry group operations.
   ********************************************************/

  void SymGroup::calc_conjugacy_classes() const {

    if(get_multi_table().size() != size()) {
      return;
    }


    conjugacy_classes.clear();

    int k;

    for(Index i = 0; i < size(); i++) {
      bool dup_class(false);
      for(Index j = 0; j < conjugacy_classes.size(); j++) {
        if(conjugacy_classes[j].contains(i)) {
          dup_class = true;
          break;
        }
      }
      if(dup_class) continue;

      conjugacy_classes.push_back(Array<Index>());

      for(Index j = 0; j < size(); j++) {
        //std::cout << "for j=" << j << ", i=" << i << ": j-inverse= " << ind_inverse(j) << ", i*j-inverse= " << ind_prod(i, ind_inverse(j));
        //int tk=alt_multi_table[j][i];
        //std::cout << "-- compare to amt[j][0]=" << alt_multi_table[j][0] << " and amt[j][i]=" << tk << " and result is k=";
        k = ind_prod(j, ind_prod(i, ind_inverse(j)));
        //std::cout << k << " -- compare to explicit value " << multi_table[tk][j];

        if(!conjugacy_classes.back().contains(k)) {
          //std::cout << " so " << k << " goes in class " << conjugacy_classes.size()-1;
          conjugacy_classes.back().push_back(k);
        }
        //std::cout << "\n";
      }
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


  //***************************************************

  Index SymGroup::ind_inverse(Index i) const {
    if(get_alt_multi_table().size() != size() || !valid_index(i) || i >= size())
      return -1;
    //std::cout << "Inside ind_inverse. 'i' is " << i << " and alt_multi_table size is " << alt_multi_table.size() << "\n";
    return alt_multi_table[i][0];
  }

  //***************************************************

  Index SymGroup::ind_prod(Index i, Index j) const {
    if(get_multi_table().size() != size()
       || !valid_index(i) || i >= size()
       || !valid_index(j) || j >= size()) {

      return -1;

    }
    //std::cout << "Inside ind_prod. 'i' is " << i << " and j is " << j << " and multi_table size is " << multi_table.size() << "\n";
    return multi_table[i][j];
  }

  //***************************************************

  Index SymGroup::class_of_op(Index i) const {
    if(!get_conjugacy_classes().size() || i > size())
      return -1;
    return index2conjugacy_class[i];
  }

  //***************************************************

  void SymGroup::set_irrep_ID(Index i, Index ID) const {
    assert((valid_index(i) && i < irrep_IDs.size()) && "Attempting to set ID for out-of-bounds irrep.");
    std::cout << "\n setting irrep ID of SymGroup " << this << "\n";
    std::cout << "IDs before: " << irrep_IDs << '\n';

    irrep_IDs[i] = ID;
    std::cout << "IDs after: " << irrep_IDs << '\n';
    return;

  }

  //***************************************************

  Index SymGroup::get_irrep_ID(Index i) const {
    if(!valid_index(i) || i >= irrep_IDs.size())
      return -1;

    return irrep_IDs[i];

  }

  //***************************************************

  Index SymGroup::make_empty_representation() const {
    if(!size() || !at(0).has_valid_master()) {
      std::cerr << "CRITICAL ERROR: In SymGroup::make_empty_representation(), SymGroup is improperly initialized.\n"
                << "                Exiting...\n";
      exit(1);
    }

    return at(0).master_group().make_empty_representation();
  }

  //***************************************************

  SymGroupRep const *SymGroup::get_irrep(Index i) const {
    if(!size() || !valid_index(i) || i >= irrep_IDs.size())
      return nullptr;

    return at(0).master_group().representation(irrep_IDs[i]);

  }
  //***************************************************

  const Array<Array<Index> > &SymGroup::get_multi_table() const {
    if(multi_table.size() != size()) {
      //std::cout << "CALCULATING MULTI_TABLE for " << this <<  ": table size is " << multi_table.size() << " and group size is " << size() << "!!\n";
      calc_multi_table();
    }
    return multi_table;
  }

  //***************************************************

  const Array<Array<Index> > &SymGroup::get_alt_multi_table() const {
    if(alt_multi_table.size() != size()) {
      //std::cout << "CALCULATING ALT_MULTI_TABLE " << this << ": table size is " << alt_multi_table.size() << " and group size is " << size() << "!!\n";
      calc_alt_multi_table();
    }
    return alt_multi_table;
  }
  //***************************************************

  void SymGroup::invalidate_multi_tables() const {
    multi_table.resize(size(), Array<Index>(size(), -1));
    alt_multi_table.resize(size(), Array<Index>(size(), -1));

  }

  //***************************************************

  const Array<Array<Index> > &SymGroup::get_conjugacy_classes() const {
    if(conjugacy_classes.size() != size())
      calc_conjugacy_classes();
    return conjugacy_classes;
  }
  //***************************************************

  const Array<bool > &SymGroup::get_complex_irrep_list() const {
    if(!character_table.size())
      calc_character_table();
    return complex_irrep;
  }

  //***************************************************

  const Array<Array<std::complex<double> > > &SymGroup::get_character_table() const {
    if(!character_table.size())
      calc_character_table();
    return character_table;
  }

  //***************************************************

  const std::string &SymGroup::get_name() const {
    if(!name.size()) {
      calc_character_table();
      if(!name.size()) {

        std::cerr << "WARNING: In SymGroup::get_name(), unable to get symgroup type.\n";
        std::cout << "group size is " << size() << '\n';
        name = "unknown";
      }
    }

    return name;
  }

  //***************************************************

  const std::string &SymGroup::get_latex_name() const {
    get_character_table();
    return latex_name;
  }

  //***************************************************

  const Array< Array<Array<Index> > > &SymGroup::get_small_subgroups() const {
    if(!small_groups.size())
      calc_small_subgroups();
    return small_groups;
  }

  //***************************************************

  const Array< Array<Array<Index> > > &SymGroup::get_large_subgroups() const {
    if(!large_groups.size())
      calc_large_subgroups();
    return large_groups;
  }

  //***************************************************

  bool SymGroup::calc_multi_table() const { //AAB
    Index i, j;

    multi_table.resize(size(), Array<Index>(size(), -1));

    for(i = 0; i < size(); i++) {
      for(j = 0; j < size(); j++) {
        multi_table[i][j] = find_periodic(at(i) * at(j));
        if(multi_table[i][j] >= size() || multi_table[i].find(multi_table[i][j]) != j) {

          //if(multi_table[i][j] >= size()) {
          //std::cout << "This SymGroup is not a group because the combination of at least two of its elements is not contained in the set.\n";

          //Returning a table of all 1's seems to make the most sense. This will prevent weird recursion from happening.
          multi_table.resize(size(), Array<Index>(size(), -1));
          //multi_table.clear();
          return false;
        }
      }
    }

    return true;

  }

  //***************************************************
  void SymGroup::calc_alt_multi_table() const {
    //by calling get_multi_table(), we ensure that multi_table is populated
    alt_multi_table.resize(get_multi_table().size());

    /*std::cout << "INSIDE CALC_ALT_MULTI_TABLE ---- Multiplication table is:\n";
    for(Index i=0; i<multi_table.size(); i++){
      std::cout << multi_table[i] << '\n';
      }*/

    if(multi_table.size() && !valid_index(multi_table[0][0])) {
      alt_multi_table.resize(size(), Array<Index>(size(), -1));
      return;
    }
    for(Index i = 0; i < multi_table.size(); i++) {
      if(multi_table[i][i] != 0) {
        alt_multi_table[multi_table[i].find(0)] = multi_table[i];
      }
      else {
        alt_multi_table[i] = multi_table[i];
      }
    }

    //     std::cout << "INSIDE CALC_ALT_MULTI_TABLE ---- Alt Multiplication table is:\n";
    //     for(Index i=0; i<alt_multi_table.size(); i++){
    //       std::cout << alt_multi_table[i] << '\n';
    //     }


  }

  //***************************************************

  /*Index SymGroup::find(const SymOp &test_op) const {
    for(Index i = 0; i < size(); i++) {
      if(at(i) == test_op) {
        return i;
      }
    }
    return size();
  }
  */
  //***************************************************
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

  //***************************************************

  Index SymGroup::find_periodic(const SymOp &test_op, double tol) const {
    for(Index i = 0; i < size(); i++) {
      if(compare_periodic(at(i), test_op, _lattice(), periodicity(), tol)) {
        return i;
      }
    }
    return size();
  }

  //***************************************************
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

  //***************************************************

  void SymGroup::enforce_group(double tol, Index max_size) {
    bool new_ops(true);

    while(new_ops && size() < max_size) {
      new_ops = false;
      for(Index i = 0; i < size() && size() < max_size; i++) {
        for(Index j = 0; j < size() && size() < max_size; j++) {
          SymOp tOp(at(i)*at(j));
          if(!contains_periodic(tOp, tol)) {
            push_back(within_cell(tOp, _lattice(), periodicity()));
            new_ops = true;
            //	  std::cout << "Pushing back a SymOp due to multiplication fail.\n";
          }
        }

        SymOp tOp(at(i).inverse());
        if(!contains_periodic(tOp, tol)) {
          push_back(within_cell(tOp, _lattice(), periodicity()));
          new_ops = true;
          //std::cout << "Pushing back a SymOp due to inverse fail.\n";
        }
      }
    }
    if(size() == max_size) {
      std::cerr << "In SymGroup::enforce_group() -- you have reached the maximum allowed size you specified for your group (the default is 1000). Unless you are generating a factor group in a large supercell, you probably need to adjust your tolerances.\n";
    }
    return;
  }

  //***************************************************

  /*bool SymGroup::contains(const SymOp &test_op) const {
    return Array<SymOp> :: contains(test_op);
    }*/

  //***************************************************

  bool SymGroup::contains_periodic(const SymOp &test_op, double tol) const {
    return find_periodic(test_op, tol) != size();
  }

  //***************************************************

  SymGroup &SymGroup::apply_sym(const SymOp &op) {
    for(Index i = 0; i < size(); i++)
      at(i).apply_sym(op);
    return *this;
  }

  //***************************************************

  void SymGroup::write(std::string filename, COORD_TYPE mode, const Eigen::Ref<const Eigen::Matrix3d> &c2f_mat) const {
    std::ofstream outfile;
    outfile.open(filename.c_str());
    print(outfile, mode, c2f_mat);
    outfile.close();
    return;
  }

  //***************************************************

  void SymGroup::print(std::ostream &out, COORD_TYPE mode, const Eigen::Ref<const Eigen::Matrix3d> &c2f_mat) const {
    out << size() << " # " << COORD_MODE::NAME(mode) << " representation of group containing " << size() << " elements:\n\n";
    for(Index i = 0; i < size(); i++) {
      out << i << "  ";
      at(i).print(out, c2f_mat);
      out << std::endl;
    }
    return;
  }


  //***************************************************

  void SymGroup::calc_space_group_in_cell(SymGroup &space_group_cell, const Lattice &_cell) const {
    if(!size()) return;

    Eigen::Vector3i max_trans(3, 3, 3);
    Coordinate trans(Eigen::Vector3d::Zero(), _cell, FRAC);
    space_group_cell.clear();

    std::vector<SymOp::SymInfo> sg_info;
    for(Index i = 0; i < size(); i++) {
      EigenCounter<Eigen::Vector3i> lat_comb(-max_trans, max_trans, Eigen::Vector3i::Ones());
      do {
        trans.frac() = lat_comb().cast<double>();
        SymOp new_sym(SymOp::translation(trans.cart())*at(i));
        auto info = new_sym.info();
        trans.cart() = info.location;
        if(!trans.is_within()) {
          continue;
        }


        bool new_location = true;
        for(Index j = 0; j < space_group_cell.size(); j++) {

          if(almost_equal(new_sym.matrix(), space_group_cell[j].matrix()) && almost_equal(info.location, sg_info[j].location)) {
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

  //***************************************************

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

        if(!space_group.contains(new_sym)) {
          space_group.push_back(new_sym);

        }

      }
      while(++lat_comb);
    }

    return;
  }

  //***************************************************
  void SymGroup::print_locations(std::ostream &stream, const Eigen::Ref<const Eigen::Matrix3d> &c2f_mat) const {
    //Assumes SymGroup is sorted with clumps of SymOps of common matrix type and eigenvec
    //sort();

    bool new_op = true;
    stream << "Locations for symmetry operations\n";
    SymOp::SymInfo info = at(0).info(), next_info;

    for(Index i = 0; i < size(); i++) {
      if(new_op) {
        at(i).print(stream, c2f_mat);
        stream << std::endl;
        at(i).print(stream, Eigen::Matrix3d::Identity());
        stream << std::endl;

        stream << "Location:" << std::endl;
        stream << "FRAC\t\t\t\t\tCART" << std::endl;
      }
      stream << c2f_mat *info.location;
      stream << "\t\t\t";
      stream << info.location;
      stream << std::endl;

      if(i + 1 < size()) {
        next_info = at(i + 1).info();
        if(info.op_type == next_info.op_type && almost_equal(info.axis, next_info.axis)) {
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

  //***************************************************

  //Sorts SymOp array by trace, with largest trace first.
  //Identity will be first, inversion is last, if it is in the group.
  //Rotations are sorted by angle magnitude, where angle is between -Pi and +Pi
  //Mirror operations appear together
  void SymGroup::sort() {
    if(!size())
      return;
    SymOp t_op(at(0));
    for(Index i = 0; i < size(); i++) {
      for(Index j = i + 1; j < size(); j++) {

        // added by Donghee - always non-translation symmetry showed up first;
        if(almost_zero(at(j).tau())) {
          t_op = at(j);
          at(j) = at(i);
          at(i) = t_op;
        }
        if((at(j).matrix().determinant() - at(i).matrix().determinant()) > TOL) {
          t_op = at(j);
          at(j) = at(i);
          at(i) = t_op;
          continue;
        }

        if(std::abs(at(j).matrix().determinant() - at(i).matrix().determinant()) < TOL
           && (at(j).matrix().trace() - at(i).matrix().trace()) > TOL) {
          t_op = at(j);
          at(j) = at(i);
          at(i) = t_op;
        }
      }
    }
    clear_tables();
    return;
  }

  //***************************************************

  void SymGroup::sort_by_class() {
    if(!size()) {
      return;
    }

    sort();
    if(get_multi_table().size() != size()) {
      return;
    }

    calc_conjugacy_classes();
    //calc_conjugacy_corr_table();

    SymOp t_op(at(0));
    Array<int> tarray;
    int tind;

    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      for(Index j = 0; j < conjugacy_classes[i].size(); j++) {
        tarray.push_back(conjugacy_classes[i][j]);
      }
    }

    for(Index i = 0; i < tarray.size(); i++) {
      for(Index j = i + 1; j < tarray.size(); j++) {
        if(tarray[j] < tarray[i]) {
          t_op = at(j);
          at(j) = at(i);
          at(i) = t_op;
          tind = tarray[j];
          tarray[j] = tarray[i];
          tarray[i] = tind;
        }
      }
    }

    clear_tables();

    return;
  }

  //***************************************************

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

  //***************************************************
  /// Translation operators for origin shift need to be defined
  SymGroup &SymGroup::operator+=(const Eigen::Ref<const SymOp::vector_type> &shift) {
    for(Index ng = 0; ng < size(); ng++)
      at(ng) += shift;
    return (*this);
  }

  //***************************************************

  SymGroup &SymGroup::operator-=(const Eigen::Ref<const SymOp::vector_type> &shift) {
    for(Index ng = 0; ng < size(); ng++)
      at(ng) -= shift;
    return (*this);
  }

  //***************************************************

  ReturnArray<Index> SymGroup::get_irrep_decomposition() const {
    Array<Index> tdec;
    Array<double> repchar;

    tdec.resize(conjugacy_classes.size());
    repchar.resize(conjugacy_classes.size());

    for(Index i = 0; i < conjugacy_classes.size(); i++) {
      repchar[i] += at(conjugacy_classes[i][0]).matrix().trace();
    }

    //    std::cout << "The reducible characters are:\n";

    //    std::cout << repchar << "\n\n";

    for(Index i = 0; i < character_table.size(); i++) { // Loop over irreducible representations
      std::complex<double> temp = std::complex<double>(0, 0);
      for(Index j = 0; j < character_table.size(); j++) { // Loop over conjugacy classes
        temp += (character_table[i][j] * std::complex<double>(conjugacy_classes[j].size(), 0) * std::complex<double>(repchar[j], 0));
      }
      //      std::cout << temp << "\t";

      tdec[i] = round(temp.real() / size());
    }

    //    std::cout << "\n\nThe irreducible projection is:\n";
    //std::cout << tdec << "\n\n";

    return tdec;
  }

  //***************************************************

  jsonParser &SymGroup::to_json(jsonParser &json) const {
    json.put_obj();

    json["symop"].put<Array<SymOp> >(*this);

    // PERIODICITY_TYPE group_periodicity;
    json["group_periodicity"] = periodicity();

    // mutable Array<Array<int> > multi_table;
    json["multi_table"] = multi_table;

    // mutable Array<Array<int> > alt_multi_table;
    json["alt_multi_table"] = alt_multi_table;

    // mutable Array<Array<int> > conjugacy_classes;
    json["conjugacy_classes"] = conjugacy_classes;

    // mutable Array<std::string> class_names;
    json["class_names"] = class_names;

    // mutable Array<int> index2conjugacy_class;
    json["index2conjugacy_class"] = index2conjugacy_class;

    // mutable Array<Array<std::complex<double> > > character_table;
    json["character_table"] = character_table;
    for(int i = 0; i < json["character_table"].size(); i++) {
      json["character_table"][i].set_force_row();
      for(int j = 0; j < json["character_table"][i].size(); j++) {
        json["character_table"][i][j].set_force_row();
        json["character_table"][i][j]["real"].set_remove_trailing_zeros();
        json["character_table"][i][j]["imag"].set_remove_trailing_zeros();
      }
    }

    // mutable Array<int> irrep_IDs;
    json["irrep_IDs"] = irrep_IDs;

    // mutable Array<bool> complex_irrep;
    json["complex_irrep"] = complex_irrep;

    // mutable Array<std::string> irrep_names;
    json["irrep_names"] = irrep_names;

    // Array<SymGroup> unique_subgroups;
    json["unique_subgroups"] = unique_subgroups;

    // mutable Array<Array<Array<int> > > small_groups;
    json["small_groups"] = small_groups;

    // mutable Array<Array<Array<int> > > large_groups;
    json["large_groups"] = large_groups;

    // mutable Array<Array<int> > centralizer_table;
    json["centralizer_table"] = centralizer_table;

    // mutable Array<Array<int> > elem_order_table;
    json["elem_order_table"] = elem_order_table;

    // mutable std::string name;
    json["name"] = name;

    // mutable std::string latex_name;
    json["latex_name"] = latex_name;

    //mutable std::string comment;
    json["comment"] = comment;

    // mutable double max_error;
    json["max_error"] = max_error;

    // mutable Array<Array<SymOp> > rotation_groups;
    json["rotation_groups"] = rotation_groups;

    // mutable std::string crystal_system;
    json["crystal_system"] = crystal_system;

    // mutable bool centric;
    json["centric"] = centric;

    // mutable Array<int> group_number; // space group number (min and max)
    json["group_number"] = group_number;

    // mutable Array<std::string> group_name; // 0: International 1: Schonflies
    json["group_name"] = group_name;

    return json;
  }

  //***************************************************

  void SymGroup::from_json(const jsonParser &json) {
    try {

      // class SymGroup : public Array<SymOp>
      std::cout << "In SymGroup::from_json" << std::endl;

      clear();
      std::cout << "Reading in the symops" << std::endl;
      for(int i = 0; i < json["symop"].size(); i++) {
        push_back(json["symop"][i].get<SymOp>());
      }

      std::cout << "Reading in group_periodicity" << std::endl;
      // PERIODICITY_TYPE group_periodicity;
      CASM::from_json(m_group_periodicity, json["group_periodicity"]);

      std::cout << "Reading in multi_table" << std::endl;
      // mutable Array<Array<int> > multi_table;
      CASM::from_json(multi_table, json["multi_table"]);

      std::cout << "Reading in alt_multi_table" << std::endl;
      // mutable Array<Array<int> > alt_multi_table;
      CASM::from_json(alt_multi_table, json["alt_multi_table"]);

      // mutable Array<Array<int> > conjugacy_classes;
      CASM::from_json(conjugacy_classes, json["conjugacy_classes"]);

      // mutable Array<std::string> class_names;
      CASM::from_json(class_names, json["class_names"]);

      // mutable Array<int> index2conjugacy_class;
      CASM::from_json(index2conjugacy_class, json["index2conjugacy_class"]);

      // mutable Array<Array<std::complex<double> > > character_table;
      CASM::from_json(character_table, json["character_table"]);

      // mutable Array<int> irrep_IDs;
      CASM::from_json(irrep_IDs, json["irrep_IDs"]);

      // mutable Array<bool> complex_irrep;
      CASM::from_json(complex_irrep, json["complex_irrep"]);

      // mutable Array<std::string> irrep_names;
      CASM::from_json(irrep_names, json["irrep_names"]);

      // Array<SymGroup> unique_subgroups;
      unique_subgroups.clear();
      for(int i = 0; i < json["unique_subgroups"].size(); i++) {

        SymGroup group(periodicity());

        CASM::from_json(group, json["unique_subgroups"][i]);
        unique_subgroups.push_back(group);
      }

      // mutable Array<Array<Array<int> > > small_groups;
      CASM::from_json(small_groups, json["small_groups"]);

      // mutable Array<Array<Array<int> > > large_groups;
      CASM::from_json(large_groups, json["large_groups"]);

      // mutable Array<Array<int> > centralizer_table;
      CASM::from_json(centralizer_table, json["centralizer_table"]);

      // mutable Array<Array<int> > elem_order_table;
      CASM::from_json(elem_order_table, json["elem_order_table"]);

      // mutable std::string name;
      CASM::from_json(name, json["name"]);

      // mutable std::string latex_name;
      CASM::from_json(latex_name, json["latex_name"]);

      //mutable std::string comment;
      CASM::from_json(comment, json["comment"]);

      // mutable double max_error;
      CASM::from_json(max_error, json["max_error"]);

      // mutable Array<Array<SymOp> > rotation_groups;
      //CASM::from_json( rotation_groups, json["rotation_groups"]);
      rotation_groups.clear();
      Array<SymOp> group;
      for(int i = 0; i < json["rotation_groups"].size(); i++) {
        group.clear();
        for(int j = 0; i < json["rotation_groups"][i].size(); j++) {
          group.push_back(json["rotation_groups"][i][j].get<SymOp>());
        }
        rotation_groups.push_back(group);
      }

      // mutable std::string crystal_system;
      CASM::from_json(crystal_system, json["crystal_system"]);

      // mutable bool centric;
      CASM::from_json(centric, json["centric"]);

      // mutable Array<int> group_number; // space group number (min and max)
      CASM::from_json(group_number, json["group_number"]);

      // mutable Array<std::string> group_name; // 0: International 1: Schonflies
      CASM::from_json(group_name, json["group_name"]);

    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //**********************************************************

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

  //**********************************************************

  jsonParser &MasterSymGroup::to_json(jsonParser &json) const {

    // class MasterSymGroup : public SymGroup
    SymGroup::to_json(json);

    // mutable Array<SymGroupRep *> rep_array;
    json["rep_array"].put_array();
    for(Index i = 0; i < rep_array.size(); i++) {
      json["rep_array"].push_back(rep_array[i]);
    }

    // mutable int coord_rep_ID, reg_rep_ID;
    json["coord_rep_ID"] = coord_rep_ID;
    json["reg_rep_ID"] = reg_rep_ID;

    // mutable SymGroup point_group_internal;
    json["point_group"] = point_group_internal;

    return json;
  }

  //***************************************************

  // Note: as a hack this expects at(0) to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the SymOp
  void MasterSymGroup::from_json(const jsonParser &json) {
    try {
      clear();

      // class MasterSymGroup : public SymGroup
      SymGroup::from_json(json);

      // mutable Array<SymGroupRep *> rep_array;
      // destruct exisiting
      for(Index i = 0; i < rep_array.size(); i++) {
        delete rep_array[i];
      }
      rep_array.resize(json["rep_array"].size());
      for(int i = 0; i < json["rep_array"].size(); i++) {
        rep_array[i] = new SymGroupRep(*this);
        rep_array[i]->from_json(json["rep_array"][i]);
      }

      // mutable int coord_rep_ID, reg_rep_ID;
      CASM::from_json(coord_rep_ID, json["coord_rep_ID"]);
      CASM::from_json(reg_rep_ID, json["reg_rep_ID"]);

      // mutable SymGroup point_group_internal;
      point_group_internal.clear();
      point_group_internal.from_json(json["point_group"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //**********************************************************

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
    if(!almost_equal(a.matrix(), b.matrix()))
      return false;

    if(periodicity != PERIODIC)
      return almost_equal(a.tau(), b.tau());

    return almost_zero(Coordinate(a.tau(), lat, CART).min_dist(Coordinate(b.tau(), lat, CART)) < _tol);

  }

  SymOp within_cell(const SymOp &a,
                    const Lattice &lat,
                    PERIODICITY_TYPE periodicity) {
    if(periodicity != PERIODIC)
      return a;

    Coordinate trans(a.tau(), lat, CART);
    trans.within();
    return SymOp(a.matrix(), trans.cart(), a.map_error());
  }


};


