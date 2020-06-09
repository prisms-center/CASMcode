#include "casm/symmetry/SymOp.hh"

#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"

namespace CASM {

  const double &SymOp::map_error() const {
    return m_map_error;
  }

  //*******************************************************************************************

  void SymOp::set_map_error(const double &value) {
    m_map_error = value;
    return;
  }


  //*******************************************************************************************

  void SymOp::set_index(const MasterSymGroup &new_group, Index new_index) {
    if((valid_index(new_index) && new_index < new_group.size())
       && (this == &(new_group[new_index]) ||
           almost_equal(matrix(), new_group[new_index].matrix()))) {
      m_master_group = &new_group;
      m_op_index = new_index;
      _set_integral_tau();
    }
    else {
      m_master_group = &new_group;
      //m_master_group = NULL;
      m_op_index = -1;
      _set_integral_tau();
    }
  }

  //*******************************************************************************************

  SymOp SymOp::operator*(const SymOp &RHS) const {
    SymOp t_op(matrix() * RHS.matrix(),
               tau() + matrix() * RHS.tau(),
               time_reversal() != RHS.time_reversal(),  //time reversal set by XOR operation
               sqrt(map_error()*map_error() + RHS.map_error()*RHS.map_error()),
               -1,
               nullptr);

    if(m_master_group && (m_master_group == RHS.m_master_group)) {
      t_op.set_index(master_group(), master_group().ind_prod(index(), RHS.index()));
    }
    else if(RHS.m_master_group && !m_master_group && is_identity()) {
      t_op.set_index(RHS.master_group(), RHS.index());
    }
    else if(m_master_group && !RHS.m_master_group && RHS.is_identity()) {
      t_op.set_index(master_group(), index());
    }
    //The following blocks caused problems at some point (mainly for non-primitive structures)
    //else if(is_identity() && RHS.is_identity()) {
    //t_op.m_op_index = 0;
    //}
    //else{
    //std::cout << "This symmetry is " << symmetry << " with head " << m_master_group << " and RHS symmetry is " << RHS.symmetry << " with head " << RHS.m_master_group << "\n";
    //}

    return t_op;

  }

  //*******************************************************************************************

  SymOp &SymOp::operator+=(const Eigen::Ref<const SymOp::vector_type> &RHS) {
    m_tau += RHS - matrix() * RHS;
    return (*this);
  }

  //*******************************************************************************************

  SymOp &SymOp::operator-=(const Eigen::Ref<const SymOp::vector_type> &RHS) {
    m_tau -= RHS - matrix() * RHS;
    return (*this);
  }

  //*******************************************************************************************
  // SymOp matrix is unitary, so inverse is equivalent to transpose.
  // To do inverse of translation, you must perform
  // inverse matrix operaton on translation and subtract
  SymOp SymOp::inverse() const {
    SymOp t_op(matrix().transpose(),
               -(matrix().transpose() * tau()),
               time_reversal(),
               map_error(),
               -1,
               nullptr);
    if(m_master_group) {
      t_op.set_index(master_group(), master_group().ind_inverse(index()));
    }
    else if(is_identity()) {
      t_op.m_op_index = 0;
      t_op._set_integral_tau();
    }

    return t_op;
  }

  //*******************************************************************************************

  SymOp SymOp::no_trans() const {
    return SymOp(matrix(), vector_type::Zero(), time_reversal(), map_error(), index(), m_master_group);
  }

  //*******************************************************************************************

  bool SymOp::operator==(const SymOp &RHS) const {
    return
      almost_equal(matrix(), RHS.matrix()) &&
      almost_equal(tau(), RHS.tau()) &&
      time_reversal() == RHS.time_reversal();
  };

  //*******************************************************************************************

  SymOp &SymOp::apply_sym(const SymOp &op) {
    (*this) = op * (*this) * (op.inverse());
    return *this;
  }

  //*******************************************************************************************

  void SymOp::print_short(std::ostream &stream, const Eigen::Ref<const SymOp::matrix_type> &c2f_mat) const {
    print(stream, c2f_mat);
  }

  //*******************************************************************************************

  void SymOp::print(std::ostream &stream, const Eigen::Ref<const SymOp::matrix_type> &c2f_mat) const {

    int tprec = stream.precision();
    std::ios::fmtflags tflags = stream.flags();

    stream.precision(3);

    stream.flags(std::ios::left);
    stream << std::setw(53) << "Symmetry Operation Matrix" << "Shift \n"; //SOM has 25 chars, width of those 3 columns are 14 each, so 42 width.  Shift width is 22, so spacing of 9-11 extra characters, so add 5 more to get in the middle

    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream.precision(9);
    matrix_type tmat(c2f_mat * matrix()*c2f_mat.inverse());
    vector_type ttau(c2f_mat * tau());
    for(int i = 0; i < 3; i++) {
      //Print each row of the symmetry matrix separately
      for(int j = 0; j < 3; j++) {
        stream << std::setw(14) << tmat(i, j);
      }

      stream << std::setw(22) << ttau(i) << "\n";
    }

    stream << "Time Reversal: " << (time_reversal() ? -1 : 1) << "\n";
    stream.precision(tprec);
    stream.flags(tflags);
    return;
  }

  //*******************************************************************************************

  const SymOp::matrix_type &get_matrix(const SymOp &op) {
    return op.matrix();
  }

  const SymOp::vector_type &get_translation(const SymOp &op) {
    return op.tau();
  }

  bool get_time_reversal(const SymOp &op) {
    return op.time_reversal();
  }

  //*******************************************************************************************

  jsonParser &SymOp::to_json(jsonParser &json) const {
    json.put_obj();

    // Members not included:
    //
    // From SymOpRepresentation:
    //   MasterSymGroup const master_group();
    //
    // From SymOp:
    //   Lattice const *home;
    //   Array<SymOpRepresentation *> representation;

    json["SymOpRep_type"] = "SymOp";

    json["op_index"] = m_op_index;
    json["rep_ID"] = m_rep_ID;

    // mutable SymOp::matrix_type symmetry_mat[2];
    json["symmetry_mat"] = matrix();

    // mutable Coordinate tau_vec;
    to_json_array(tau(), json["tau"]);

    // double map_error;
    json["map_error"] = map_error();

    return json;
  }

  //*******************************************************************************************

  void SymOp::from_json(const jsonParser &json) {
    try {
      //std::cout<<"Inside of SymOp::from_json"<<std::endl;
      //std::cout<<"Reading in op_index"<<std::endl;
      CASM::from_json(m_op_index, json["op_index"]);
      //std::cout<<"Reading in rep_id"<<std::endl;
      CASM::from_json(m_rep_ID, json["rep_ID"]);

      // mutable SymOp::matrix_type symmetry_mat[2];
      //std::cout<<"Reading in symmetry_mat"<<std::endl;
      CASM::from_json(m_mat, json["symmetry_mat"]);

      // mutable Coordinate tau_vec;
      //std::cout<<"Reading in tau_vec"<<std::endl;
      if(json.contains("tau"))
        CASM::from_json(m_tau, json["tau"]);
      _set_integral_tau();
      //std::cout<<"Reading in map_error"<<std::endl;
      // double map_error;
      CASM::from_json(m_map_error, json["map_error"]);
      //std::cout<<"Done Reading in the SymOp"<<std::endl;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


  void SymOp::_set_integral_tau() {
    if(has_valid_master() && valid_index(index())) {
      m_integral_tau = tau() - master_group()[index()].tau();
      m_valid_integral_tau = true;
    }
    else if(is_identity()) {
      m_integral_tau = tau();
      m_valid_integral_tau = true;
    }
    else {
      m_valid_integral_tau = false;
    }
  }




  //*******************************************************************************************

  jsonParser &to_json(const SymOp &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  //*******************************************************************************************
  void from_json(SymOp &sym, const jsonParser &json) {
    try {
      sym.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  /// \brief Print formatted SymOp matrix and tau
  ///
  /// \code
  ///                       Matrix |       Tau
  ///      0.5  -0.86603         0 |     1.617
  ///  0.86603       0.5         0 |   0.93357
  ///        0         0         1 |    2.5843
  /// \endcode
  ///
  /// - Includes Log indentation level
  ///
  void print_matrix_tau_col(Log &log, const SymOp &op, Index prec) {
    Index w1 = 3;
    w1 = print_matrix_width(log, op.matrix(), w1);
    Eigen::IOFormat fmt1(prec, w1);

    Index w2 = 3;
    w2 = print_matrix_width(log, op.tau(), w2);
    Eigen::IOFormat fmt2(prec, w2);

    std::string header = std::string(2 + w1 * 3 - 6, ' ') + "Matrix | " + std::string(w2 - 3, ' ') + "Tau";
    log << log.indent_str() << header << std::endl;
    log << log.indent_str() << op.matrix().row(0).format(fmt1) << " | " << op.tau().row(0).format(fmt2) << std::endl;
    log << log.indent_str() << op.matrix().row(1).format(fmt1) << " | " << op.tau().row(1).format(fmt2) << std::endl;
    log << log.indent_str() << op.matrix().row(2).format(fmt1) << " | " << op.tau().row(2).format(fmt2) << std::endl;
  }

}
