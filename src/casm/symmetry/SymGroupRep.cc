#include "casm/symmetry/SymGroupRep.hh"

#include <numeric>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/container/Permutation.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/VectorSymCompare.hh"

namespace CASM {

SymGroupRep::SymGroupRep(const SymGroupRep &RHS)
    : std::vector<SymOpRepresentation *>(RHS.size(), NULL) {
  (*this) = RHS;
}

//*******************************************************************************************

SymGroupRep::~SymGroupRep() { clear(); }

//*******************************************************************************************

SymGroupRep &SymGroupRep::operator=(const SymGroupRep &RHS) {
  m_master_group = RHS.m_master_group;
  clear();
  if (RHS.size() > 0 && has_valid_master() &&
      master_group().size() != RHS.size()) {
    throw std::runtime_error(
        "Invalid assignment of SymGroupRep.  Sizes are incompatible.\n");
  }
  resize(RHS.size(), NULL);
  for (Index i = 0; i < RHS.size(); i++) {
    if (RHS[i]) set_rep(i, *RHS[i]);
  }
  return *this;
}

//*******************************************************************************************

void SymGroupRep::set_master_group(const MasterSymGroup &master,
                                   const SymGroupRepID &_rep_ID) {
  m_master_group = &master;
  if (_rep_ID.empty() || &(master.representation(_rep_ID)) != this) {
    throw std::runtime_error(std::string(
        "SymGroupRep::set_master_group() attempted to assign SymGroupRepID "
        "that does not match the current SymGroupRep!\n"));
  }
  m_rep_ID = _rep_ID;
  if (size() == 0)
    std::vector<SymOpRepresentation *>::resize(master.size());
  else if (size() == master.size()) {
    for (Index i = 0; i < size(); i++) {
      if (at(i)) at(i)->set_identifiers(master, m_rep_ID, i);
    }
  } else {
    throw std::runtime_error(
        std::string("SymGroupRep::set_master_group() passed new master whose "
                    "size is incompatible with the size\n") +
        "of the current representation.\n");
  }
}

//*******************************************************************************************

void SymGroupRep::set_rep(const SymOp &base_op,
                          const SymOpRepresentation &new_rep) {
  if (!has_valid_master()) {
    err_log() << "CRITICAL ERROR: In SymGroupRep::set_rep(), you are trying to "
                 "assign the representation of a SymOp whose factor_group is "
                 "not specified!\n"
              << "                Exiting...\n";
    assert(0);
    exit(1);
  }
  if (valid_index(base_op.index())) return set_rep(base_op.index(), new_rep);

  // else:
  return set_rep(master_group().find_periodic(base_op), new_rep);
}

//*******************************************************************************************

void SymGroupRep::set_rep(const SymOpRepresentation &base_op,
                          const SymOpRepresentation &new_rep) {
  if (!has_valid_master()) {
    err_log() << "CRITICAL ERROR: In SymGroupRep::set_rep(), you are trying to "
                 "assign the representation of a SymOpRepresentation whose "
                 "factor_group is not specified!\n"
              << "                Exiting...\n";
    assert(0);
    exit(1);
  }
  if (valid_index(base_op.index()))
    set_rep(base_op.index(), new_rep);
  else {
    err_log() << "CRITICAL ERROR: In SymGroupRep::set_rep(), you are trying to "
                 "assign the representation of a SymOpRepresentation whose "
                 "index is not specified!\n"
              << "                Exiting...\n";
    assert(0);
    exit(1);
  }
}

//*******************************************************************************************

void SymGroupRep::set_rep(Index op_index, const SymOpRepresentation &new_rep) {
  assert(valid_index(op_index) && op_index < size() &&
         "In SymGroupRep::set_rep(), reference representation is improperly "
         "initialized.");
  if (at(op_index)) {
    err_log() << "CRITICAL ERROR: In SymGroupRep::set_rep(), representation "
                 "already exists for operation "
              << op_index << ".\n"
              << "                Exiting...\n";
    assert(0);
    exit(1);
  }

  SymOpRepresentation *tcopy = new_rep.copy();
  if (has_valid_master())
    tcopy->set_identifiers(master_group(), m_rep_ID, op_index);

  // avoid doing this elsewhere in CASM
  const_cast<SymOpRepresentation *&>(at(op_index)) = tcopy;
}

//*******************************************************************************************

void SymGroupRep::clear() {
  for (Index i = 0; i < size(); i++) delete at(i);
  std::vector<SymOpRepresentation *>::clear();
}

//*******************************************************************************************

Eigen::MatrixXd const *SymGroupRep::MatrixXd(Index i) const {
  return at(i)->MatrixXd();
}

//*******************************************************************************************

Eigen::MatrixXd const *SymGroupRep::MatrixXd(
    const SymOpRepresentation &op) const {
  return at(op.index())->MatrixXd();
}

//*******************************************************************************************

Permutation const *SymGroupRep::permutation(Index i) const {
  return at(i)->permutation();
}

//*******************************************************************************************

Permutation const *SymGroupRep::permutation(
    const SymOpRepresentation &op) const {
  return at(op.index())->permutation();
}

//*******************************************************************************************

//*******************************************************************************************
// Calculates new SymGroupRep that is the results of performing coordinate
// transformation specified by trans_mat The ROWS of trans_mat are the new basis
// vectors in terms of the old such that new_symrep_matrix = trans_mat *
// old_symrep_matrix * trans_mat.transpose();
SymGroupRep coord_transformed_copy(SymGroupRep const &_rep,
                                   const Eigen::MatrixXd &trans_mat) {
  SymGroupRep new_rep(_rep.master_group());
  if (!_rep.size()) return new_rep;
  if (_rep[0] && !(_rep.MatrixXd(0))) {
    err_log() << "CRITICAL ERROR: Trying to perform matrix transformation on a "
                 "non-matrix SymRep. Exiting...\n";
    assert(0);
    exit(1);
  }
  Eigen::MatrixXd rightmat;
  rightmat =
      trans_mat.transpose()
          .jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(Eigen::MatrixXd::Identity(trans_mat.cols(), trans_mat.cols()))
          .transpose();

  for (Index i = 0; i < _rep.size(); i++) {
    if (!_rep[i]) continue;

    new_rep.set_rep(i,
                    SymMatrixXd(trans_mat * (*(_rep.MatrixXd(i))) * rightmat));
  }
  return new_rep;
}

}  // namespace CASM
