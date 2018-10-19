#ifndef ConfigDoF_HH
#define ConfigDoF_HH

#include <vector>
#include "casm/CASM_global_definitions.hh"
#include "casm/CASM_global_Eigen.hh"
#include "casm/clex/ConfigDoFValues.hh"
namespace CASM {

  class PermuteIterator;
  class jsonParser;

  /// \brief A container class for the different degrees of freedom a Configuration
  /// might have
  ///
  /// Contains an std::vector<int> that tells you the current occupant of each
  /// site, an Eigen::MatrixXd that tells you the displacements at each site, and
  /// a LatticeStrain that tells you the strain of the Configuration. Everything
  /// is public.
  ///
  /// \ingroup Configuration
  ///
  class ConfigDoF {

  public:
    using GlobalValueType = GlobalContinuousConfigDoFValues;
    using LocalValueType = LocalContinuousConfigDoFValues;

    // Can treat as a Eigen::VectorXd
    //typedef displacement_matrix_t::ColXpr displacement_t;
    //typedef displacement_matrix_t::ConstColXpr const_displacement_t;

    /// Initialize with number of sites, and dimensionality of global and local DoFs
    /// GlobalContainerType is an iterable container of value_type std::pair<DoFKey,ContinuousDoFInfo>
    /// LocalContainerType is an iterable container of value_type std::pair<DoFKey,std::vector<ContinuousDoFInfo>  >
    template<typename GlobalContainerType, typename LocalContainerType>
    ConfigDoF(Index _N_sublat, Index _N_vol, GlobalContainerType const &global_dof_dims, LocalContainerType const &local_dof_dims, double _tol);

    ///Initialize with explicit occupation
    //ConfigDoF(const std::vector<int> &_occupation, double _tol = TOL);

    Index size() const {
      return m_occupation.size();
    }

    double tol() const {
      return m_tol;
    }

    void clear();


    // -- Occupation ------------------

    int &occ(Index i) {
      return m_occupation[i];
    };

    const int &occ(Index i) const {
      return m_occupation[i];
    };

    /// set_occupation ensures that ConfigDoF::size() is compatible with _occupation.size()
    /// or if ConfigDoF::size()==0, sets ConfigDoF::size() to _occupation.size()
    void set_occupation(const std::vector<int> &_occupation);

    const std::vector<int> &occupation() const {
      return m_occupation;
    };

    bool has_occupation() const {
      return size() != 0 && occupation().size() == size();
    }

    //void resize(Index _N);

    std::map<DoFKey, GlobalValueType> const &global_dofs() const {
      return m_global_dofs;
    }

    GlobalValueType const &global_dof(DoFKey const &_key) const {
      return global_dofs().at(_key);
    }

    GlobalValueType &global_dof(DoFKey const &_key) {
      auto it = m_global_dofs.find(_key);
      if(it == m_global_dofs.end())
        throw std::runtime_error("Attempting to access unitialized ConfigDoF value for '" + _key + "'");
      return it->second;
    }

    bool has_global_dof(DoFKey const &_key) const {
      return global_dofs().count(_key);
    }

    //void reset_global_dof(DoFKey const &_key);

    void set_global_dof(DoFKey const &_key, Eigen::Ref<const Eigen::VectorXd> const &_val);

    std::map<DoFKey, LocalValueType> const &local_dofs() const {
      return m_local_dofs;
    }

    LocalValueType const &local_dof(DoFKey const &_key) const {
      return local_dofs().at(_key);
    }

    LocalValueType &local_dof(DoFKey const &_key) {
      auto it = m_local_dofs.find(_key);
      if(it == m_local_dofs.end())
        throw std::runtime_error("Attempting to access unitialized ConfigDoF value for '" + _key + "'");
      return it->second;
    }



    bool has_local_dof(DoFKey const &_key) const {
      return local_dofs().count(_key);
    }

    //void reset_local_dof(DoFKey const &_key);

    void set_local_dof(DoFKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val);

    ConfigDoF &apply_sym(PermuteIterator const &it);

    /// \brief Calculate transformed ConfigDoF from PermuteIterator,
    /// using only the effect of symmetry on the value at each site
    /// sites are not permuted as a result
    ConfigDoF &apply_sym_no_permute(SymOp const &_op);

    void swap(ConfigDoF &RHS);

    //**** I/O ****
    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);//, Index NB);


  private:

    // ***DON'T FORGET: If you add something here, also update ConfigDoF::swap!!


    /// \brief Number of sites in the Configuration
    Index m_N;

    /// With one value for each site in the Configuration, this std::vector
    /// describes which occupant is at each of the 'N' sites of the configuration
    ///
    /// 'occupation' is a list of the indices describing the occupants in each crystal site.
    ///   prim().basis()[ sublat(i) ].site_occupant[ occupation[i]] -> Molecule on site i
    ///   This means that for the background structure, 'occupation' is all 0
    ///
    /// Configuration sites are arranged by basis, and then prim:
    ///   occupation: [basis0                |basis1               |basis2          |...] up to prim.basis().size()
    ///       basis0: [prim0|prim1|prim2|...] up to supercell.volume()
    ///
    ///
    std::vector<int> m_occupation;

    std::map<std::string, GlobalContinuousConfigDoFValues> m_global_dofs;

    std::map<std::string, LocalContinuousConfigDoFValues> m_local_dofs;

    /// Tolerance used for transformation to canonical form -- used also for comparisons, since
    /// Since comparisons are only meaningful to within the tolerance used for finding the canonical form
    /// (This is relevant only for displacement and deformation degrees of freedom
    mutable double m_tol;

  };

  template<>
  struct jsonConstructor<ConfigDoF> {

    static ConfigDoF from_json(const jsonParser &json, Index NB, std::map<DoFKey, DoFSetInfo> const &global_info, std::map<DoFKey, std::vector<DoFSetInfo> > const &local_info, double _tol);
    //from_json(const jsonParser &json, Index NB);
  };

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json);

  void from_json(ConfigDoF &value, const jsonParser &json);//, Index NB);

  void swap(ConfigDoF &A, ConfigDoF &B);

  inline
  void reset_properties(ConfigDoF &_dof) {
    return;
  }

  /// Initialize with number of sites, and dimensionality of global and local DoFs
  /// GlobalContainerType is an iterable container of value_type std::pair<DoFKey,ContinuousDoFInfo>
  /// LocalContainerType is an iterable container of value_type std::pair<DoFKey,std::vector<ContinuousDoFInfo>  >
  template<typename GlobalContainerType, typename LocalContainerType>
  ConfigDoF::ConfigDoF(Index _N_sublat, Index _N_vol, GlobalContainerType const &global_dof_info, LocalContainerType const &local_dof_info, double _tol) :
    m_occupation(_N_sublat * _N_vol, 0),
    m_tol(_tol) {
    for(auto const &dof : global_dof_info) {
      DoF::BasicTraits const &ttraits = DoF::traits(dof.first);

      if(!ttraits.global())
        throw std::runtime_error("Attempting to initialize ConfigDoF global value using local DoF " + dof.first);
      m_global_dofs[dof.first] = GlobalContinuousConfigDoFValues(ttraits, _N_sublat, _N_vol, Eigen::VectorXd::Zero(dof.second.dim()), dof.second);
    }
    for(auto const &dof : local_dof_info) {
      DoF::BasicTraits const &ttraits = DoF::traits(dof.first);
      if(_N_sublat == 0)
        continue;
      if(ttraits.global())
        throw std::runtime_error("Attempting to initialize ConfigDoF local value using global DoF " + dof.first);
      if(_N_sublat != dof.second.size()) {
        throw std::runtime_error("Attempting to initialize ConfigDoF local value '"
                                 + dof.first
                                 + "' with improperly initialized parameter 'local_dof_info'.");
      }
      Index dim = 0;
      for(auto const &info : dof.second)
        dim = max(dim, info.dim());

      m_local_dofs[dof.first] = LocalContinuousConfigDoFValues(ttraits, _N_sublat, _N_vol, Eigen::MatrixXd::Zero(dim, _N_sublat * _N_vol), dof.second);
    }
  }
}

#endif
