#ifndef CASM_clexulator_BaseClexulator
#define CASM_clexulator_BaseClexulator
#include <algorithm>
#include <cstddef>
#include <sstream>

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace clexulator {

class ClexParamPack;
class ClexParamKey;

/// \brief Abstract base class for cluster expansion correlation calculations
class BaseClexulator {
 public:
  typedef unsigned int size_type;

  BaseClexulator(size_type _nlist_size, size_type _corr_size,
                 size_type _n_point_corr);

  virtual ~BaseClexulator();

  /// \brief Neighbor list size
  size_type nlist_size() const { return m_nlist_size; }

  /// \brief Number of correlations
  size_type corr_size() const { return m_corr_size; }

  /// \brief Number of distinct point correlations per unit cell
  ///
  /// Note:
  /// - This is not the number of point functions
  /// - [0, n_point_corr()) is the valid range for the `neighbor_ind`
  /// argument to for the point correlation calculating member functions
  /// (calc_point_corr, calc_restricted_point_corr, calc_delta_point_corr, and
  /// calc_restricted_delta_point_corr)
  /// - For periodic clex, this is number of prim basis sites with site degrees
  /// of freedom included in the cluster functions.
  /// - For local clex, this is the number of sites in the local neighborhood
  /// with site degrees of freedom included in the cluster functions  .
  size_type n_point_corr() const { return m_n_point_corr; }

  /// \brief Clone the Clexulator
  std::unique_ptr<BaseClexulator> clone() const;

  /// \brief Obtain const reference to abstract ClexParamPack object
  virtual ClexParamPack const &param_pack() const = 0;

  /// \brief Obtain reference to abstract ClexParamPack object
  virtual ClexParamPack &param_pack() = 0;

  /// \brief Obtain ClexParamKey for a particular parameter
  ClexParamKey const &param_key(std::string const &_param_name) const;

  /// \brief The UnitCellCoord involved in calculating the basis functions,
  /// relative origin UnitCell
  const std::set<xtal::UnitCell> &neighborhood() const {
    return m_neighborhood;
  }

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  const std::set<xtal::UnitCell> &neighborhood(
      size_type linear_orbit_index) const {
    return m_orbit_neighborhood[linear_orbit_index];
  }

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  const std::set<xtal::UnitCellCoord> &site_neighborhood(
      size_type linear_orbit_index) const {
    return m_orbit_site_neighborhood[linear_orbit_index];
  }

  /// \brief The weight matrix used for ordering the neighbor list
  const PrimNeighborList::Matrix3Type &weight_matrix() const {
    return m_weight_matrix;
  }

  /// \brief The sublattice indices included in the neighbor list
  std::set<int> const &sublat_indices() const { return m_sublat_indices; }

  /// \brief The total number of sublattices in the prim
  size_type n_sublattices() const { return m_n_sublattices; }

  // Note: The following all require setting the DoFValues and NeighborList
  // prior to use. This can be done with:
  //
  // myclexulatorbase.set_configdofvalues(myconfigdofvalues);
  // auto nlist_begin =
  //     supercell_neighbor_list.sites(linear_unitcell_index).data();
  // myclexulatorbase.set_nlist(nlist_begin);

  /// \brief Set internal pointers to correct DoF values
  ///
  /// \param _configdofvalues DoF values to be used as input to the basis
  ///      functions
  /// \param _force
  ///
  /// Notes:
  /// - In the vast majority of cases this is handled by the `calc_X` method
  void set_configdofvalues(ConfigDoFValues const &_configdofvalues,
                           bool _force = false) const {
    if (m_configdofvalues_ptr != &_configdofvalues || _force) {
      m_configdofvalues_ptr = &_configdofvalues;
      m_occ_ptr = _configdofvalues.occupation.data();
      for (auto const &dof : m_local_dof_registry) {
        auto it = _configdofvalues.local_dof_values.find(dof.first);
        if (it == _configdofvalues.local_dof_values.end()) {
          std::stringstream msg;
          msg << "Clexulator error: ConfigDoFValues missing required local DoF "
                 "type '"
              << dof.first << "'";
          throw std::runtime_error(msg.str());
        }
        m_local_dof_ptrs[dof.second] = &it->second;
      }

      for (auto const &dof : m_global_dof_registry) {
        auto it = _configdofvalues.global_dof_values.find(dof.first);
        if (it == _configdofvalues.global_dof_values.end()) {
          std::stringstream msg;
          msg << "Clexulator error: ConfigDoFValues missing required global "
                 "DoF type '"
              << dof.first << "'";
          throw std::runtime_error(msg.str());
        }
        m_global_dof_ptrs[dof.second] = &it->second;
      }
    }
  }

  /// \brief Set pointer to neighbor list
  ///
  /// Notes:
  /// - A SuperNeighborList instance contains one array of neighbor indices for
  ///   each unit cell in the supercell.
  /// - Correlations are evaluated for a particular unit cell, scanning over
  ///   all unit cells in the supercell if necessary.
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // UnitCellCoord of site in
  /// Configuration int l_index = my_supercell.find(bijk); // Linear index of
  /// site in Configuration
  ///
  /// auto nlist_begin =
  ///     supercell_neighbor_list.sites(linear_unitcell_index).data();
  /// myclexulatorbase.set_nlist(nlist_begin);
  /// \endcode
  ///
  void set_nlist(const long int *_nlist_begin) const {
    m_nlist_ptr = _nlist_begin;
    return;
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  void calc_global_corr_contribution(double *_corr_begin) const {
    _calc_global_corr_contribution(_corr_begin);
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  ///
  /// Notes:
  /// - This overload writes values to the ParamPack, but not to a correlations
  ///   vector
  ///
  void calc_global_corr_contribution() const {
    _calc_global_corr_contribution();
  }

  /// \brief Calculate contribution to select global correlations from one unit
  /// cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  /// \param _corr_ind_begin,_corr_ind_end Pointers to range indicating which
  /// correlations should be calculated
  ///
  void calc_restricted_global_corr_contribution(
      double *_corr_begin, size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const {
    _calc_restricted_global_corr_contribution(_corr_begin, _corr_ind_begin,
                                              _corr_ind_end);
  }

  /// \brief Calculate point correlations about basis site 'neighbor_ind'
  ///
  /// \param neighbor_ind Basis site index about which to calculate correlations
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  void calc_point_corr(int neighbor_ind, double *_corr_begin) const {
    _calc_point_corr(neighbor_ind, _corr_begin);
  }

  /// \brief Calculate select point correlations about basis site 'neighbor_ind'
  ///
  /// \param neighbor_ind Neighbor within specified neighborhood about which to
  ///     calculate correlations
  /// \param _corr_begin Pointer to beginning of data structure where
  ///     correlations are written
  /// \param _corr_ind_begin,_corr_ind_end Pointers to range indicating which
  ///     correlations should be calculated
  ///
  void calc_restricted_point_corr(int neighbor_ind, double *_corr_begin,
                                  size_type const *_corr_ind_begin,
                                  size_type const *_corr_ind_end) const {
    _calc_restricted_point_corr(neighbor_ind, _corr_begin, _corr_ind_begin,
                                _corr_ind_end);
  }

  /// \brief Calculate the change in point correlations due to changing an
  /// occupant
  ///
  /// \param neighbor_ind Basis site index about which to calculate correlations
  /// \param occ_i,occ_f Initial and final occupant variable
  /// \param _corr_begin Pointer to beginning of data structure where difference
  /// in correlations are written
  ///
  void calc_delta_point_corr(int neighbor_ind, int occ_i, int occ_f,
                             double *_corr_begin) const {
    _calc_delta_point_corr(neighbor_ind, occ_i, occ_f, _corr_begin);
  }

  /// \brief Calculate the change in select point correlations due to changing
  /// an occupant
  ///
  /// \param neighbor_ind Basis site index about which to calculate correlations
  /// \param occ_i,occ_f Initial and final occupant variable
  /// \param _corr_begin Pointer to beginning of data structure where difference
  /// in correlations are written \param _corr_ind_begin,_corr_ind_end Pointers
  /// to range indicating which correlations should be calculated
  ///
  void calc_restricted_delta_point_corr(int neighbor_ind, int occ_i, int occ_f,
                                        double *_corr_begin,
                                        size_type const *_corr_ind_begin,
                                        size_type const *_corr_ind_end) const {
    _calc_restricted_delta_point_corr(neighbor_ind, occ_i, occ_f, _corr_begin,
                                      _corr_ind_begin, _corr_ind_end);
  }

 private:
  /// \brief Clone the Clexulator
  virtual BaseClexulator *_clone() const = 0;

  /// \brief The neighbor list size
  size_type m_nlist_size;

  /// \brief The number of correlations
  size_type m_corr_size;

  /// \brief Valid range for `neighbor_ind` argument to calc_point_corr
  size_type m_n_point_corr;

  std::map<std::string, Index> m_local_dof_registry;

  std::map<std::string, Index> m_global_dof_registry;

 protected:
  virtual void _calc_global_corr_contribution() const = 0;

  virtual void _calc_global_corr_contribution(double *_corr_begin) const = 0;

  virtual void _calc_restricted_global_corr_contribution(
      size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const = 0;

  virtual void _calc_restricted_global_corr_contribution(
      double *_corr_begin, size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const = 0;

  virtual void _calc_point_corr(int neighbor_ind) const = 0;

  virtual void _calc_point_corr(int neighbor_ind,
                                double *_corr_begin) const = 0;

  virtual void _calc_restricted_point_corr(
      int neighbor_ind, size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const = 0;

  virtual void _calc_restricted_point_corr(
      int neighbor_ind, double *_corr_begin, size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const = 0;

  virtual void _calc_delta_point_corr(int neighbor_ind, int occ_i,
                                      int occ_f) const = 0;

  virtual void _calc_delta_point_corr(int neighbor_ind, int occ_i, int occ_f,
                                      double *_corr_begin) const = 0;

  virtual void _calc_restricted_delta_point_corr(
      int neighbor_ind, int occ_i, int occ_f, size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const = 0;

  virtual void _calc_restricted_delta_point_corr(
      int neighbor_ind, int occ_i, int occ_f, double *_corr_begin,
      size_type const *_corr_ind_begin,
      size_type const *_corr_ind_end) const = 0;

  void _register_local_dof(std::string const &_type_name, Index _ind) {
    Index new_size = std::max(Index(_ind), Index(m_local_dof_ptrs.size())) + 1;
    m_local_dof_ptrs.resize(new_size, nullptr);
    m_global_dof_ptrs.resize(new_size, nullptr);
    m_local_dof_registry[_type_name] = _ind;
  }

  void _register_global_dof(std::string const &_type_name, Index _ind) {
    Index new_size = std::max(Index(_ind), Index(m_local_dof_ptrs.size())) + 1;
    m_local_dof_ptrs.resize(new_size, nullptr);
    m_global_dof_ptrs.resize(new_size, nullptr);
    m_global_dof_registry[_type_name] = _ind;
  }

  /// \brief access reference to internally pointed ConfigDoF
  ConfigDoFValues const &_configdofvalues() const {
    return *m_configdofvalues_ptr;
  }

  /// \brief access reference to internally pointed ConfigDoF
  Index const &_l(Index nlist_ind) const { return *(m_nlist_ptr + nlist_ind); }

  /// \brief access reference to internally pointed occupation list
  int const &_occ(Index nlist_ind) const {
    return *(m_occ_ptr + *(m_nlist_ptr + nlist_ind));
  }

  /// \brief The UnitCell involved in calculating the basis functions,
  /// relative origin UnitCell
  std::set<xtal::UnitCell> m_neighborhood;

  /// \brief The UnitCell involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  std::vector<std::set<xtal::UnitCell> > m_orbit_neighborhood;

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  std::vector<std::set<xtal::UnitCellCoord> > m_orbit_site_neighborhood;

  /// \brief The weight matrix used for ordering the neighbor list
  PrimNeighborList::Matrix3Type m_weight_matrix;

  /// \brief The sublattice indices included in the neighbor list
  std::set<int> m_sublat_indices;

  /// \brief The total number of sublattices in the prim
  size_type m_n_sublattices;

  mutable std::vector<Eigen::MatrixXd const *> m_local_dof_ptrs;

  mutable std::vector<Eigen::VectorXd const *> m_global_dof_ptrs;

 private:
  /// \brief Pointer to ConfigDoFValues for which evaluation is occuring
  mutable ConfigDoFValues const *m_configdofvalues_ptr;

  /// \brief Pointer to neighbor list
  mutable long int const *m_nlist_ptr;

  /// \brief Pointer to occupation list
  mutable int const *m_occ_ptr;
};

}  // namespace clexulator
}  // namespace CASM

#endif
