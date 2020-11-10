#ifndef ConfigDoF_HH
#define ConfigDoF_HH

#include <vector>
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/clex/ConfigDoFValues.hh"
#include "casm/container/ContainerTraits.hh"
namespace CASM {

  class PermuteIterator;
  class jsonParser;
  class SymGroupRepID;
  class SymOp;

  /// A container class for the different degrees of freedom values a Configuration might have
  ///
  /// \ingroup Configuration
  ///
  /// ConfigDoF has:
  /// - The values of the discrete site DoF within the supercell ("occupation", Eigen::VectorXi)
  ///
  ///   Example: Occupation values, accessed via `Eigen::VectorXi const &ConfigDoF::occupation()`,
  ///   with integer value corresponding to which Molecule in the the a Site::occupant_dof vector
  ///   is occupying a particular site:
  ///
  ///       [<- sublattice 0 "occ" values -> | <- sublattice 1 "occ" values -> | ... ]
  ///
  /// - The values of the continuous site DoF within the supercell ("local_dofs",
  ///   std::map<DoFKey, LocalContinuousConfigDoFValues>)
  ///
  ///   Example: Displacement values, with prim DoF basis equal to the standard basis (dx, dy, dz),
  ///   accessed via `Eigen::MatrixXd const &ConfigDoF::local_dofs("disp").values()`:
  ///
  ///       [<- sublattice 0 dx values -> | <- sublattice 1 dx values -> | ... ]
  ///       [<- sublattice 0 dy values -> | <- sublattice 1 dy values -> | ... ]
  ///       [<- sublattice 0 dz values -> | <- sublattice 1 dz values -> | ... ]
  ///
  ///   Example: Displacement values, with non-standard prim DoF basis:
  ///
  ///       "basis" : [ {
  ///           "coordinate": [c0x, c0y, c0z],
  ///           "occupants": [...],
  ///           "dofs": {
  ///             "disp" : {
  ///               "axis_names" : ["dxy", "dz"],
  ///               "axes" : [[1.0, 1.0, 0.0],
  ///                       [0.0, 0.0, 1.0]]}}
  ///         },
  ///         {
  ///           "coordinate": [c1x, c1y, c1z],
  ///           "occupants": [...],
  ///           "dofs": {
  ///             "disp" : {
  ///               "axis_names" : ["d\bar{x}y", "dz"],
  ///               "axes" : [[-1.0, 1.0, 0.0],
  ///                       [0.0, 0.0, 1.0]]}}
  //          },
  ///         ...
  ///       }
  ///
  ///       [<- sublattice 0 dxy values -> | <- sublattice 1 d\bar{x}y values ->| ... ]
  ///       [<- sublattice 0 dz values  -> | <- sublattice 1 dz values ->       | ... ]
  ///       [<- 0.0 values ->              | <- 0.0 values ->                   | ... ]
  ///
  /// - The values of the continuous global DoF ("global_dofs",
  ///   std::map<DoFKey, GlobalDoFContainerType>).
  ///
  ///   Example: GLstrain values, with prim DoF basis equal to the standard basis, accessed via
  ///   `Eigen::VectorXd const &ConfigDoF::global_dofs("GLstrain").values()`:
  ///
  ///       [e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy]
  ///
  /// Note: Continuous DoF values stored in memory in ConfigDoF are coordinates in the
  /// prim DoF basis (the "axes" given in "prim.json" which set xtal::SiteDoFSet::basis() or
  /// xtal::DoFSet::basis()), but when saved to file (i.e. `.casm/config/config_list.json`) they are
  /// saved as coordinates in the standard DoF basis (with axes meaning as described by
  /// AnisoValTraits::standard_var_names()).
  ///
  class ConfigDoF {

  public:
    using GlobalDoFContainerType = GlobalContinuousConfigDoFValues;
    using LocalDoFContainerType = LocalContinuousConfigDoFValues;
    using OccDoFContainerType = LocalDiscreteConfigDoFValues;

    using OccValueType = OccDoFContainerType::ValueType;

    // Can treat as a Eigen::VectorXd
    //typedef displacement_matrix_t::ColXpr displacement_t;
    //typedef displacement_matrix_t::ConstColXpr const_displacement_t;

    /// Initialize with number of sites, and dimensionality of global and local DoFs
    /// GlobalInfoContainerType is an iterable container having value_type std::pair<DoFKey,ContinuousDoFInfo>
    /// LocalInfoContainerType is an iterable container having value_type std::pair<DoFKey,std::vector<ContinuousDoFInfo>  >
    /// OccInfoContainerType is an iterable container having value_type SymGroupRepID
    template<typename GlobalInfoContainerType, typename LocalInfoContainerType>
    ConfigDoF(Index _N_sublat,
              Index _N_vol,
              GlobalInfoContainerType const &global_dof_info,
              LocalInfoContainerType const &local_dof_info,
              std::vector<SymGroupRepID> const &occ_symrep_IDs,
              double _tol);

    ///\brief Number of sites in the ConfigDoF
    Index size() const {
      return occupation().size();
    }

    ///\brief Integer volume of ConfigDoF
    Index n_vol() const {
      return m_occupation.n_vol();
    }

    /// \brief Number of sublattices in ConfigDoF
    Index n_sublat() const {
      return m_occupation.n_sublat();
    }

    /// \brief tolerance for comparison of continuous DoF values
    double tol() const {
      return m_tol;
    }

    void clear();


    // -- Occupation ------------------

    int &occ(Index i) {
      return m_occupation.values()[i];
    }

    const int &occ(Index i) const {
      return m_occupation.values()[i];
    }

    /// set_occupation ensures that ConfigDoF::size() is compatible with _occupation.size()
    /// or if ConfigDoF::size()==0, sets ConfigDoF::size() to _occupation.size()
    template <typename OtherOccContainerType, typename std::enable_if<std::is_integral<typename ContainerTraits<OtherOccContainerType>::value_type>::value>::type * = nullptr>
    void set_occupation(const OtherOccContainerType &_occupation) {
      if(occupation().size() != _occupation.size())
        throw std::runtime_error("Size mismatch in ConfigDoF::set_occupation()");
      for(Index i = 0; i < occupation().size(); ++i)
        occ(i) = _occupation[i];
    }

    OccValueType const &occupation() const {
      return m_occupation.values();
    }

    bool has_occupation() const {
      return size() != 0 && occupation().size() == size();
    }

    std::map<DoFKey, GlobalDoFContainerType> const &global_dofs() const {
      return m_global_dofs;
    }

    GlobalDoFContainerType const &global_dof(DoFKey const &_key) const {
      auto it = m_global_dofs.find(_key);
      if(it == m_global_dofs.end())
        throw std::runtime_error("Attempting to access uninitialized ConfigDoF value for '" + _key + "'");
      return it->second;
    }

    GlobalDoFContainerType &global_dof(DoFKey const &_key) {
      auto it = m_global_dofs.find(_key);
      if(it == m_global_dofs.end())
        throw std::runtime_error("Attempting to access uninitialized ConfigDoF value for '" + _key + "'");
      return it->second;
    }

    bool has_global_dof(DoFKey const &_key) const {
      return global_dofs().count(_key);
    }

    void set_global_dof(DoFKey const &_key, Eigen::Ref<const Eigen::VectorXd> const &_val);

    std::map<DoFKey, LocalDoFContainerType> const &local_dofs() const {
      return m_local_dofs;
    }

    LocalDoFContainerType const &local_dof(DoFKey const &_key) const {
      auto it = m_local_dofs.find(_key);
      if(it == m_local_dofs.end())
        throw std::runtime_error("Attempting to access uninitialized ConfigDoF value for '" + _key + "'");
      return it->second;
    }

    LocalDoFContainerType &local_dof(DoFKey const &_key) {
      auto it = m_local_dofs.find(_key);
      if(it == m_local_dofs.end())
        throw std::runtime_error("Attempting to access uninitialized ConfigDoF value for '" + _key + "'");
      return it->second;
    }

    bool has_local_dof(DoFKey const &_key) const {
      return local_dofs().count(_key);
    }

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
    OccDoFContainerType m_occupation;

    std::map<std::string, GlobalDoFContainerType> m_global_dofs;

    std::map<std::string, LocalDoFContainerType> m_local_dofs;

    /// Tolerance used for transformation to canonical form -- used also for continuous DoF
    /// comparisons, since comparisons are only meaningful to within the tolerance used for finding
    /// the canonical form.
    mutable double m_tol;

  };

  void swap(ConfigDoF &A, ConfigDoF &B);

  inline
  void reset_properties(ConfigDoF &_dof) {
    return;
  }

  /// Initialize with number of sites, and dimensionality of global and local DoFs
  /// GlobalInfoContainerType is an iterable container of value_type std::pair<DoFKey,ContinuousDoFInfo>
  /// LocalInfoContainerType is an iterable container of value_type std::pair<DoFKey,std::vector<ContinuousDoFInfo>  >
  template<typename GlobalInfoContainerType, typename LocalInfoContainerType>
  ConfigDoF::ConfigDoF(Index _N_sublat,
                       Index _N_vol,
                       GlobalInfoContainerType const &global_dof_info,
                       LocalInfoContainerType const &local_dof_info,
                       std::vector<SymGroupRepID> const &occ_symrep_IDs,
                       double _tol) :
    m_occupation(DoF::BasicTraits("occ"), _N_sublat, _N_vol, OccValueType::Zero(_N_sublat * _N_vol), occ_symrep_IDs),
    m_tol(_tol) {
    for(auto const &dof : global_dof_info) {
      DoF::BasicTraits ttraits(dof.first);

      if(!ttraits.global())
        throw std::runtime_error("Attempting to initialize ConfigDoF global value using local DoF " + dof.first);
      m_global_dofs[dof.first] = GlobalContinuousConfigDoFValues(ttraits, _N_sublat, _N_vol, Eigen::VectorXd::Zero(dof.second.dim()), dof.second);
    }
    for(auto const &dof : local_dof_info) {
      DoF::BasicTraits ttraits(dof.first);
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
