#include "casm/symmetry/SymGroupRep.hh"

#include "casm/external/Eigen/CASM_AddOns"
#include "casm/misc/CASM_math.hh"

#include "casm/container/Permutation.hh"
#include "casm/symmetry/SymMatrixXd.hh"

namespace CASM {
  //INITIALIZE STATIC MEMBER SymGroupRep::REP_COUNT
  //THIS MUST OCCUR IN A .CC FILE; MAY CAUSE PROBLEMS IF WE
  //CHANGE COMPILING/LINKING STRATEGY
  Index SymGroupRep::REP_COUNT(0);

  //***************************************************

  SymGroupRep::SymGroupRep(const SymGroupRep &RHS) :
    Array<SymOpRepresentation * > (0), m_rep_ID(REP_COUNT++) {
    (*this) = RHS;
  }

  //***************************************************

  SymGroupRep::~SymGroupRep() {
    for(Index i = 0; i < size(); i++) delete at(i);
  };

  //***************************************************

  SymGroupRep &SymGroupRep::operator=(const SymGroupRep &RHS) {
    m_home_group = RHS.m_home_group;
    clear();
    for(Index i = 0; i < RHS.size(); i++)
      push_back(*RHS[i]);
    return *this;
  }

  //***************************************************

  void SymGroupRep::push_back(const SymOpRepresentation &new_op) {
    Array<SymOpRepresentation *>::push_back(new_op.copy());
    back()->set_identifiers(*m_home_group, m_rep_ID, size() - 1);
  }

  //***************************************************

  void SymGroupRep::set_master_group(const MasterSymGroup &master) {
    m_home_group = &master;
    for(Index i = 0; i < size(); i++)
      at(i)->set_identifiers(master, m_rep_ID, i);
  }

  //***************************************************

  void SymGroupRep::set_rep(Index i, const SymOpRepresentation &new_rep) const {
    assert(valid_index(i) && i < size() && "In SymGroupRep::set_rep(), index out of bounds.");
    if(at(i)) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::set_rep(), representation already exists for operation " << i << ".\n"
                << "                Exiting...\n";
      exit(1);
    }

    // avoid doing this elsewhere in CASM
    const_cast<SymOpRepresentation *&>(at(i)) = new_rep.copy();

  }

  //***************************************************

  Index SymGroupRep::add_self_to_home() {
    if(!m_home_group) return -1;

    const SymGroupRep *t_rep(m_home_group->representation(m_rep_ID));
    if(t_rep == this)
      return m_rep_ID;

    if(t_rep == nullptr)
      return m_home_group->add_representation(this);

    //else - update rep_ID and add this
    m_rep_ID = REP_COUNT++;
    return m_home_group->add_representation(this);

  }

  //***************************************************

  void SymGroupRep::clear() {
    for(Index i = 0; i < size(); i++) delete at(i);
    Array<SymOpRepresentation *>::clear();
  }

  //John G 050513

  //***************************************************

  void SymGroupRep::print_permutation(std::ostream &stream) const {
    for(Index i = 0; i < this->size(); i++) {
      if(at(i)->get_permutation())
        stream << *(at(i)->get_permutation()) << '\n';
      else
        stream << "nullptr\n";
    }
    return;
  }

  //\John G

  //***************************************************

  void SymGroupRep::print_MatrixXd(std::ostream &stream) const {
    for(Index i = 0; i < size(); i++) {
      stream << "SymRep Matrix " << i + 1 << " of " << size() << ":\n";
      if(at(i)->get_MatrixXd())
        stream << *(at(i)->get_MatrixXd()) << '\n';
      else
        stream << "nullptr\n";
      stream << "\n";
    }
    return;
  }

  //***************************************************

  void SymGroupRep::print_MatrixXd(std::ostream &stream, const SymGroup &subgroup) const {
    for(Index i = 0; i < subgroup.size(); i++) {
      stream << "SymRep Matrix " << i + 1 << " of " << subgroup.size() << ":\n";
      if(at(subgroup[i].index())->get_MatrixXd())
        stream << *(at(subgroup[i].index())->get_MatrixXd()) << '\n';
      else
        stream << "nullptr\n";
      stream << "\n";
    }
    return;
  }

  //***************************************************
  Eigen::MatrixXd SymGroupRep::block_shape_matrix() const {
    if(!size() || !get_MatrixXd(0))
      return Eigen::MatrixXd();

    Eigen::MatrixXd block_shape(get_MatrixXd(0)->cwiseProduct(*get_MatrixXd(0)));

    for(Index i = 1; i < size(); i++) {
      block_shape += get_MatrixXd(i)->cwiseProduct(*get_MatrixXd(i));

    }
    return block_shape;
  }
  //***************************************************
  Eigen::MatrixXd SymGroupRep::block_shape_matrix(const SymGroup &subgroup) const {
    if(!size() || !subgroup.size() || !get_MatrixXd(subgroup[0].index()))
      return Eigen::MatrixXd();

    Eigen::MatrixXd block_shape(get_MatrixXd(subgroup[0].index())->cwiseProduct(*get_MatrixXd(subgroup[0].index())));

    for(Index i = 1; i < subgroup.size(); i++) {
      block_shape += get_MatrixXd(subgroup[i].index())->cwiseProduct(*get_MatrixXd(subgroup[i].index()));
    }

    return block_shape;
  }

  //***************************************************
  Index SymGroupRep::num_blocks() const {
    Eigen::MatrixXd bmat(block_shape_matrix());

    if(bmat.cols() == 0)
      return 0;

    Index Nb = 1;

    //count zeros on first superdiagonal
    for(EigenIndex i = 0; i + 1 < bmat.rows(); i++) {
      if(almost_zero(bmat(i, i + 1)))
        Nb++;
    }

    return Nb;

  }
  //***************************************************

  Index SymGroupRep::num_blocks(const SymGroup &subgroup) const {
    Eigen::MatrixXd bmat(block_shape_matrix(subgroup));

    if(bmat.cols() == 0)
      return 0;

    Index Nb = 1;

    //count zeros on first superdiagonal
    for(EigenIndex i = 0; i + 1 < bmat.rows(); i++) {
      if(almost_zero(bmat(i, i + 1)))
        Nb++;
    }

    return Nb;
  }
  //***************************************************

  Eigen::MatrixXd const *SymGroupRep::get_MatrixXd(Index i) const {
    return at(i)->get_MatrixXd();
  }

  //***************************************************

  Eigen::MatrixXd const *SymGroupRep::get_MatrixXd(const SymOpRepresentation &op) const {
    return at(op.index())->get_MatrixXd();
  }

  //***************************************************

  Permutation const *SymGroupRep::get_permutation(Index i) const {
    return at(i)->get_permutation();
  }

  //***************************************************

  Permutation const *SymGroupRep::get_permutation(const SymOpRepresentation &op) const {
    return at(op.index())->get_permutation();
  }



  //***************************************************

  jsonParser &SymGroupRep::to_json(jsonParser &json) const {
    json.put_obj();

    // Member not included in json:
    //
    //   Pointer to the m_home_group that generated this SymGroupRep
    //   MasterSymGroup const *m_home_group;

    // class SymGroupRep : public Array<SymOpRepresentation *>
    json["symop_representations"].put_array();
    for(Index i = 0; i < size(); i++) {
      at(i)->to_json(json["symop_representations"]);
    }

    // static int REP_COUNT;
    json["REP_COUNT"] = REP_COUNT;

    // int m_rep_ID;
    json["m_rep_ID"] = m_rep_ID;

    return json;
  }

  //***************************************************

  // If 'm_home_group' is not nullptr, should be initialized accordingly
  // Lattice may be necessary for constructing SymOps
  void SymGroupRep::from_json(const jsonParser &json, const Lattice &lat) {
    try {

      // Member not included in json:
      //
      //   Pointer to the m_home_group that generated this SymGroupRep
      //   MasterSymGroup const *m_home_group;

      // class SymGroupRep : public Array<SymOpRepresentation *>

      for(Index i = 0; i < size(); i++) {
        delete at(i);
      }
      clear();
      std::cout << "Resizing SymGroupRep to " << json["symop_representations"].size() << std::endl;
      resize(json["symop_representations"].size());
      std::cout << "Reading in the symmetry operations" << std::endl;
      for(int i = 0; i < json["symop_representations"].size(); i++) {
        /// This allocates a new object to 'at(i)'.
        ///   It might need a Lattice
        std::cout << "Working on:" << i << std::endl;
        CASM::from_json(at(i), json["symop_representations"][i], lat);
      }

      std::cout << "Reading in REP_COUNT" << std::endl;
      // static int REP_COUNT;
      CASM::from_json(REP_COUNT, json["REP_COUNT"]);

      std::cout << "Reading in m_rep_id" << std::endl;
      // int m_rep_ID;
      CASM::from_json(m_rep_ID, json["m_rep_ID"]);

      std::cout << "Done reading in the permute_group" << std::endl;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //**********************************************************

  jsonParser &to_json(const SymGroupRep &rep, jsonParser &json) {
    return rep.to_json(json);
  }

  // If 'm_home_group' is not nullptr, should be initialized accordingly
  // Lattice may be necessary for constructing SymOps
  void from_json(SymGroupRep &rep, const jsonParser &json, const Lattice &lat) {
    try {
      rep.from_json(json, lat);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


};

