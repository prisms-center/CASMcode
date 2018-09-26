#include "casm/basis_set/DoFSet.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace DoFType {

    /// \brief implements json parsing of a specialized DoFSet.
    void BasicTraits::from_json(DoFSet &_dof, jsonParser const &json) const {
      return;
    }

    /// \brief implements json parsing of a specialized DoFSet.
    void BasicTraits::to_json(DoFSet const &_dof, jsonParser &json) const {
      return;
    }
  }

  //********************************************************************

  DoFType::Traits const &DoFSet::traits() const {
    return DoFType::traits(type_name());
  }

  //********************************************************************
  void DoFSet::allocate_symrep(SymGroup const &_group) const {
    if(!m_symrep_ID.empty())
      throw std::runtime_error("In DoFSet::allocate_symrep(), representation has already been allocated for this symrep.");

    m_symrep_ID = _group.allocate_representation();

  }
  //********************************************************************

  bool DoFSet::identical(DoFSet const &rhs)const {
    if(!std::equal(begin(), end(), rhs.begin(), compare_no_value)) {
      return false;
    }

    return almost_equal(m_basis, rhs.m_basis);
  }

  //********************************************************************

  void DoFSet::transform_basis(Eigen::Ref<const Eigen::MatrixXd> const &trans_mat) {
    m_basis = trans_mat * m_basis;
    m_symrep_ID = SymGroupRepID();
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
      m_basis=json["coordinate_space"].get<Eigen::MatrixXd>();
    */
    std::vector<ContinuousDoF> tdof;
    bool is_error = false;

    auto const &traits_ref = DoFType::traits(type_name());

    json.get_if(m_excluded_occs, "excluded_occupants");
    auto it = json.find("basis");
    if(it != json.end()) {
      if(it->is_array()) {
        m_basis.resize(traits_ref.dim(), it->size());
        int col = 0;
        for(auto const &el : *it) {
          if(is_error || el.size() != (traits_ref.dim() + 1) || !el[0].is_string()) {
            is_error = true;
            break;
          }

          tdof.push_back(ContinuousDoF(traits_ref,
                                       el[0].get<std::string>(), //var_name
                                       -1, // ID
                                       -std::numeric_limits<double>::infinity(),
                                       std::numeric_limits<double>::infinity()));
          for(Index row = 0; row < traits_ref.dim() + 1; row++) {
            if(!el[row].is_number()) {
              is_error = true;
              break;
            }
            m_basis(row, col) = el[row].get<double>();

          }
        }
      }
      else
        is_error = true;

      if(is_error) {
        std::stringstream ss;
        ss << json;
        throw std::runtime_error("Parsing malformed JSON DoF object " + ss.str() + " each element of object \"basis\" must be an array containing " + std::to_string(traits_ref.dim() + 1) + " elements.\n The first element of each sub-array must be a string, and the remaining elements must be numbers. ");
      }

    }
    else {
      m_basis.setIdentity(traits_ref.dim(), traits_ref.dim());
      for(std::string var_name : traits_ref.standard_var_names()) {
        tdof.push_back(ContinuousDoF(traits_ref,
                                     var_name,
                                     -1, // ID
                                     -std::numeric_limits<double>::infinity(),
                                     std::numeric_limits<double>::infinity()));

      }
    }
    traits_ref.from_json(*this, json);

  }

  //********************************************************************

  jsonParser &DoFSet::to_json(jsonParser &json) const {
    /*
    if(!m_basis.isIdentity(1e-5))
      json["coordinate_space"]=coordinate_space();
    */
    if(!m_excluded_occs.empty())
      json["excluded_occupants"] = m_excluded_occs;

    DoF::traits(type_name()).to_json(*this, json);
    return json;
  }

  //********************************************************************

  /// \brief Apply SymOp to a DoFSet
  DoFSet &apply(const SymOp &op, DoFSet &_dof) {
    _dof.transform_basis(DoFType::traits(_dof.type_name()).symop_to_matrix(op));
    return _dof;
  }

  //********************************************************************

  /// \brief Copy and apply SymOp to a DoFSet
  DoFSet copy_apply(const SymOp &op, const DoFSet &_dof) {
    DoFSet result(_dof);
    return apply(op, result);
  }

  //********************************************************************

  DoFSet jsonConstructor<DoFSet>::from_json(const jsonParser &json, DoF::BasicTraits const &_type) {
    DoFSet value(_type);
    value.from_json(json);
    return value;
  }
}
