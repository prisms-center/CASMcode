#include "casm/basis_set/DoFSet.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

  DoFSet::DoFSet(BasicTraits const &_traits) :
    m_traits(_traits),
    m_info(SymGroupRepID(), Eigen::MatrixXd::Zero(m_traits.dim(), 0)) {

  }


  //********************************************************************
  void DoFSet::allocate_symrep(SymGroup const &_group) const {
    if(!m_info.symrep_ID().empty())
      throw std::runtime_error("In DoFSet::allocate_symrep(), representation has already been allocated for this symrep.");
    //std::cout << "Allocating symrep...\n"
    m_info.set_symrep_ID(_group.allocate_representation());

  }
  //********************************************************************

  bool DoFSet::identical(DoFSet const &rhs)const {
    if(!std::equal(begin(), end(), rhs.begin(), compare_no_value)) {
      return false;
    }

    return almost_equal(basis(), rhs.basis());
  }

  //********************************************************************

  void DoFSet::transform_basis(Eigen::Ref<const Eigen::MatrixXd> const &trans_mat) {
    m_info.set_basis(trans_mat * basis());
    m_info.set_symrep_ID(SymGroupRepID());
  }

  //********************************************************************

  bool DoFSet::update_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs) {

    Index ID_ind;
    bool is_updated(false);
    for(Index i = 0; i < m_components.size(); i++) {
      // IMPORTANT: Do before_IDs.find(), NOT m_components().find() (if such a thing existed)
      ID_ind = find_index(before_IDs, m_components[i].ID());
      // Only set ID if DoF doesn't have an ID lock
      if(ID_ind < after_IDs.size() && !m_components[i].is_locked()) {
        m_components[i].set_ID(after_IDs[ID_ind]);
        // The new ID only changes the formula if the corresponding coeff is nonzero
        //if(!almost_zero(m_coeffs[i]))
        is_updated = true;
      }
    }
    return is_updated;
  }

  //********************************************************************

  void DoFSet::from_json(jsonParser const &json) {
    /*
      if(json.contains("coordinate_space"))
      m_info.basis=json["coordinate_space"].get<Eigen::MatrixXd>();
    */
    //std::vector<ContinuousDoF> tdof;

    m_components.clear();
    m_excluded_occs.clear();
    m_info.set_symrep_ID(SymGroupRepID());

    bool is_error = false;

    json.get_if(m_excluded_occs, "excluded_occupants");
    auto it = json.find("axes");
    if(it != json.end()) {
      if(it->is_array()) {
        std::vector<std::string> anames;
        json.get_if(anames, "axis_names");
        if(anames.size() != it->size()) {
          throw std::runtime_error("Parsing DoF " + type_name() + ", field \"axes\" must have the same number of elements as field \"axes\"");
        }

        for(std::string const &aname : anames)
          m_components.push_back(ContinuousDoF(traits(),
                                               aname,
                                               -1, // ID
                                               -std::numeric_limits<double>::infinity(),
                                               std::numeric_limits<double>::infinity()));

        Eigen::MatrixXd tbasis = it->get<Eigen::MatrixXd>().transpose();
        if(tbasis.rows() != traits().dim()) {
          throw std::runtime_error("Parsing DoF " + type_name() + ", number of columns in field \"axes\" must be " + std::to_string(traits().dim()));
        }
        m_info.set_basis(tbasis);
      }
      else
        is_error = true;

      if(is_error) {
        std::stringstream ss;
        ss << json;
        throw std::runtime_error("Parsing malformed JSON DoF object " + ss.str() + " each element of object \"basis\" must be an array containing " + std::to_string(traits().dim() + 1) + " elements.\n The first element of each sub-array must be a string, and the remaining elements must be numbers. ");
      }

    }
    else {
      m_info.set_basis(Eigen::MatrixXd::Identity(traits().dim(), traits().dim()));
      for(std::string var_name : traits().standard_var_names()) {
        //std::cout << "Adding var_name " << var_name << "\n";
        m_components.push_back(ContinuousDoF(traits(),
                                             var_name,
                                             -1, // ID
                                             -std::numeric_limits<double>::infinity(),
                                             std::numeric_limits<double>::infinity()));
        if(traits().global())
          m_components.back().lock_ID();

      }
    }
    //traits().from_json(*this, json);

  }

  //********************************************************************

  jsonParser &DoFSet::to_json(jsonParser &json) const {
    /*
    if(!m_info.basis.isIdentity(1e-5))
      json["coordinate_space"]=coordinate_space();
    */
    if(!m_excluded_occs.empty())
      json["excluded_occupants"] = m_excluded_occs;

    //DoF::traits(type_name()).to_json(*this, json);
    return json;
  }

  //********************************************************************

  /// \brief Apply SymOp to a DoFSet

  DoFSet &apply(const xtal::SymOp &op, DoFSet &_dof) {
    _dof.transform_basis(_dof.traits().symop_to_matrix(op.matrix(), op.tau(), op.time_reversal()));
    return _dof;
  }

  //********************************************************************

  /// \brief Copy and apply SymOp to a DoFSet
  DoFSet copy_apply(const xtal::SymOp &op, const DoFSet &_dof) {
    DoFSet result(_dof);
    return apply(op, result);
  }

  //********************************************************************

  DoFSet jsonConstructor<DoFSet>::from_json(const jsonParser &json, DoFSet::BasicTraits const &_type) {
    DoFSet value(_type);
    value.from_json(json);
    return value;
  }

  //********************************************************************

  DoFSet DoFSet::make_default(DoFSet::BasicTraits const &_type) {
    DoFSet result(_type);
    result.m_info.set_basis(Eigen::MatrixXd::Identity(_type.dim(), _type.dim()));
    for(std::string var_name : _type.standard_var_names()) {
      //std::cout << "Adding var_name " << var_name << "\n";
      result.m_components.push_back(ContinuousDoF(_type,
                                                  var_name,
                                                  -1, // ID
                                                  -std::numeric_limits<double>::infinity(),
                                                  std::numeric_limits<double>::infinity()));
      if(_type.global())
        result.m_components.back().lock_ID();

    }
    return result;
  }

}
