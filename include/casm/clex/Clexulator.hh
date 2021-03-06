#ifndef CLEXULATOR_HH
#define CLEXULATOR_HH
#include <boost/filesystem/path.hpp>
#include <cstddef>

#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
/** \defgroup Clexulator
 *  \ingroup ClexClex
 *  \brief Clexulator is a class for evaluating crystal basis functions and
 * cluster expansions
 *
 *  @{
 */

class ClexParamPack;
class ClexParamKey;

namespace Clexulator_impl {

/// \brief Abstract base class for cluster expansion correlation calculations
class Base {
 public:
  typedef unsigned int size_type;

  Base(size_type _nlist_size, size_type _corr_size, size_type _n_point_corr);

  virtual ~Base();

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
  std::unique_ptr<Base> clone() const;

  /// \brief Obtain const reference to abstract ClexParamPack object
  virtual ClexParamPack const &param_pack() const = 0;

  /// \brief Obtain reference to abstract ClexParamPack object
  virtual ClexParamPack &param_pack() = 0;

  /// \brief Obtain ClexParamKey for a particular parameter
  ClexParamKey const &param_key(std::string const &_param_name) const;

  /// \brief Alter evaluation of parameters specified by @param _param_key,
  /// using a custom double -> double function set
  void set_evaluation(
      ClexParamKey const &_param_key,
      std::vector<std::function<double(ConfigDoF const &)> > const &_basis_set);

  /// \brief Alter evaluation of parameters specified by @param _param_key,
  /// using a custom int -> double function set
  void set_evaluation(
      ClexParamKey const &_param_key,
      std::vector<std::function<double(std::vector<double> const &)> > const
          &_basis_set);

  /// \brief Alter evaluation of parameters specified by @param _param_key,
  /// using the string  @param _eval_type,
  // which can be at least either "READ" (i.e., read from ClexParamPack) or
  // "DEFAULT" (i.e., the Clexulator's default implementation)
  void set_evaluation(ClexParamKey const &_param_key,
                      std::string const &_eval_type);

  /// \brief Check evaluation mode of parameters specified by @param _param_key,
  /// which can be one of (at least) "READ" (i.e., read from ClexParamPack),
  /// "CUSTOM", or "DEFAULT" (i.e., the Clexulator's default implementation)
  std::string check_evaluation(ClexParamKey const &_param_key) const;

  /// \brief The UnitCellCoord involved in calculating the basis functions,
  /// relative origin UnitCell
  const std::set<UnitCell> &neighborhood() const { return m_neighborhood; }

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  const std::set<UnitCell> &neighborhood(size_type linear_orbit_index) const {
    return m_orbit_neighborhood[linear_orbit_index];
  }

  /// \brief The weight matrix used for ordering the neighbor list
  const PrimNeighborList::Matrix3Type &weight_matrix() const {
    return m_weight_matrix;
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get
  /// contribution from int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration
  /// myclexulator.calc_global_corr_contribution(my_configdof,
  ///                                            my_supercell.get_nlist(l_index).begin(),
  ///                                            correlation_array.begin());
  /// \endcode
  ///
  void calc_global_corr_contribution(ConfigDoF const &_input_configdof,
                                     long int const *_nlist_begin,
                                     long int const *_nlist_end,
                                     double *_corr_begin,
                                     double *_corr_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_global_corr_contribution(_corr_begin);
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  ///
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get
  /// contribution from int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration
  /// myclexulator.calc_global_corr_contribution(my_configdof,
  ///                                            my_supercell.get_nlist(l_index).begin(),
  ///                                            my_supercell.get_nlist(l_index).end())
  /// \endcode
  ///
  void calc_global_corr_contribution(ConfigDoF const &_input_configdof,
                                     long int const *_nlist_begin,
                                     long int const *_nlist_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_global_corr_contribution();
  }

  /// \brief Calculate contribution to select global correlations from one unit
  /// cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written \param _corr_ind_begin,_corr_ind_end Pointers to
  /// range indicating which correlations should be calculated
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get
  /// contribution from int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration std::vector<int> _corr_ind = {0, 2, 4, 6}; // Get
  /// contribution to correlations 0, 2, 4, and 6
  /// myclexulator.calc_restricted_global_corr_contribution(my_configdof,
  ///                                                       my_supercell.get_nlist(l_index).begin(),
  ///                                                       correlation_array.begin(),
  ///                                                       _corr_ind.begin(),
  ///                                                       _corr_ind.end());
  /// \endcode
  ///
  void calc_restricted_global_corr_contribution(
      ConfigDoF const &_input_configdof, long int const *_nlist_begin,
      long int const *_nlist_end, double *_corr_begin, double *_corr_end,
      size_type const *_corr_ind_begin, size_type const *_corr_ind_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_restricted_global_corr_contribution(_corr_begin, _corr_ind_begin,
                                              _corr_ind_end);
  }

  /// \brief Calculate point correlations about basis site 'neighbor_ind'
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point
  /// correlations int l_index = my_supercell.find(bijk); // Linear index of
  /// site in Configuration myclexulator.calc_point_corr(my_configdof,
  /// my_supercell.get_nlist(l_index).begin(), b, correlation_array.begin());
  /// \endcode
  ///
  void calc_point_corr(ConfigDoF const &_input_configdof,
                       long int const *_nlist_begin, long int const *_nlist_end,
                       int neighbor_ind, double *_corr_begin,
                       double *_corr_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_point_corr(neighbor_ind, _corr_begin);
  }

  /// \brief Calculate select point correlations about basis site 'neighbor_ind'
  ///
  /// \brief neighbor_ind Neighbor within specified neighborhood about which to
  /// calculate correlations \brief _corr_begin Pointer to beginning of data
  /// structure where correlations are written \param
  /// _corr_ind_begin,_corr_ind_end Pointers to range indicating which
  /// correlations should be calculated
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point
  /// correlations int l_index = my_supercell.find(bijk); // Linear index of
  /// site in Configuration std::vector<int> _corr_ind = {0, 2, 4, 6}; // Get
  /// contribution to correlations 0, 2, 4, and 6
  /// myclexulator.calc_restricted_point_corr(my_configdof,
  ///                                         my_supercell.get_nlist(l_index).begin(),
  ///                                         b,
  ///                                         correlation_array.begin(),
  ///                                         _corr_ind.begin(),
  ///                                         _corr_ind.end());
  /// \endcode
  ///
  void calc_restricted_point_corr(ConfigDoF const &_input_configdof,
                                  long int const *_nlist_begin,
                                  long int const *_nlist_end, int neighbor_ind,
                                  double *_corr_begin, double *_corr_end,
                                  size_type const *_corr_ind_begin,
                                  size_type const *_corr_ind_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_restricted_point_corr(neighbor_ind, _corr_begin, _corr_ind_begin,
                                _corr_ind_end);
  }

  /// \brief Calculate the change in point correlations due to changing an
  /// occupant
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief occ_i,occ_f Initial and final occupant variable
  /// \brief _corr_begin Pointer to beginning of data structure where difference
  /// in correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta
  /// point correlations int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration int occ_i=0, occ_f=1;  // Swap from occupant 0 to
  /// occupant 1 myclexulator.calc_delta_point_corr(my_configdof,
  /// my_supercell.get_nlist(l_index).begin(), b, occ_i, occ_f,
  /// correlation_array.begin()); \endcode
  ///
  void calc_delta_point_corr(ConfigDoF const &_input_configdof,
                             long int const *_nlist_begin,
                             long int const *_nlist_end, int neighbor_ind,
                             int occ_i, int occ_f, double *_corr_begin,
                             double *_corr_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_delta_point_corr(neighbor_ind, occ_i, occ_f, _corr_begin);
  }
  /// \brief Calculate the change in select point correlations due to changing
  /// an occupant
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief occ_i,occ_f Initial and final occupant variable
  /// \brief _corr_begin Pointer to beginning of data structure where difference
  /// in correlations are written \param _corr_ind_begin,_corr_ind_end Pointers
  /// to range indicating which correlations should be calculated
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta
  /// point correlations int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration int occ_i=0, occ_f=1;  // Swap from occupant 0 to
  /// occupant 1 std::vector<int> _corr_ind = {0, 2, 4, 6}; // Get contribution
  /// to correlations 0, 2, 4, and 6
  /// myclexulator.calc_restricted_delta_point_corr(my_configdof,
  ///                                               my_supercell.get_nlist(l_index).begin(),
  ///                                               b,
  ///                                               occ_i,
  ///                                               occ_f,
  ///                                               correlation_array.begin(),
  ///                                               _corr_ind.begin(),
  ///                                               _corr_ind.end());
  /// \endcode
  ///
  void calc_restricted_delta_point_corr(ConfigDoF const &_input_configdof,
                                        long int const *_nlist_begin,
                                        long int const *_nlist_end,
                                        int neighbor_ind, int occ_i, int occ_f,
                                        double *_corr_begin, double *_corr_end,
                                        size_type const *_corr_ind_begin,
                                        size_type const *_corr_ind_end) const {
    _set_configdof(_input_configdof);
    _set_nlist(_nlist_begin, _nlist_end);
    _calc_restricted_delta_point_corr(neighbor_ind, occ_i, occ_f, _corr_begin,
                                      _corr_ind_begin, _corr_ind_end);
  }

 private:
  /// \brief Clone the Clexulator
  virtual Base *_clone() const = 0;

  /// \brief Set internal pointers to correct DoFs
  void _set_configdof(ConfigDoF const &_input_configdof) const {
    if (m_config_ptr != &_input_configdof) {
      m_config_ptr = &_input_configdof;
      m_occ_ptr = _input_configdof.occupation().data();
      for (auto const &dof : m_local_dof_registry) {
        m_local_dof_ptrs[dof.second] = &(_configdof().local_dof(dof.first));
      }

      for (auto const &dof : m_global_dof_registry) {
        m_global_dof_ptrs[dof.second] = &(_configdof().global_dof(dof.first));
      }
    }
  }

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
    Index new_size = max(Index(_ind), Index(m_local_dof_ptrs.size())) + 1;
    m_local_dof_ptrs.resize(new_size, nullptr);
    m_global_dof_ptrs.resize(new_size, nullptr);
    m_local_dof_registry[_type_name] = _ind;
  }

  void _register_global_dof(std::string const &_type_name, Index _ind) {
    Index new_size = max(Index(_ind), Index(m_local_dof_ptrs.size())) + 1;
    m_local_dof_ptrs.resize(new_size, nullptr);
    m_global_dof_ptrs.resize(new_size, nullptr);
    m_global_dof_registry[_type_name] = _ind;
  }

  /// \brief access reference to internally pointed ConfigDoF
  ConfigDoF const &_configdof() const { return *m_config_ptr; }

  /// \brief access reference to internally pointed ConfigDoF
  Index const &_l(Index nlist_ind) const { return *(m_nlist_ptr + nlist_ind); }

  /// \brief access reference to internally pointed occupation list
  int const &_occ(Index nlist_ind) const {
    return *(m_occ_ptr + *(m_nlist_ptr + nlist_ind));
  }

  /// \brief Set pointer to neighbor list
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // UnitCellCoord of site in
  /// Configuration int l_index = my_supercell.find(bijk); // Linear index of
  /// site in Configuration \endcode
  ///
  void _set_nlist(const long int *_nlist_begin,
                  const long int *_nlist_end) const {
    m_nlist_ptr = _nlist_begin;
    // std::cout << "setting _nlist: ";
    // for(;_nlist_begin!=_nlist_end; ++_nlist_begin)
    // std::cout << *_nlist_begin << "  ";
    // std::cout << std::endl;
    return;
  }

  /// \brief The UnitCellCoord involved in calculating the basis functions,
  /// relative origin UnitCell
  std::set<UnitCell> m_neighborhood;

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  std::vector<std::set<UnitCell> > m_orbit_neighborhood;

  /// \brief The weight matrix used for ordering the neighbor list
  PrimNeighborList::Matrix3Type m_weight_matrix;

  mutable std::vector<LocalContinuousConfigDoFValues const *> m_local_dof_ptrs;

  mutable std::vector<GlobalContinuousConfigDoFValues const *>
      m_global_dof_ptrs;

 private:
  /// \brief Pointer to ConfigDoF for which evaluation is occuring
  mutable ConfigDoF const *m_config_ptr;

  /// \brief Pointer to neighbor list
  mutable long int const *m_nlist_ptr;

  /// \brief Pointer to occupation list
  mutable int const *m_occ_ptr;
};
}  // namespace Clexulator_impl

class RuntimeLibrary;

/// \brief Evaluates correlations
///
/// CASM generates code for very efficient calculation of basis functions via
/// the print_clexulator function. This source code may be compiled, linked,
/// and used at runtime via Clexulator.
///
/// \ingroup Clexulator
///
class Clexulator {
 public:
  typedef Clexulator_impl::Base::size_type size_type;

  Clexulator() {}

  /// \brief Construct a Clexulator
  Clexulator(std::string name, fs::path dirpath, PrimNeighborList &nlist,
             std::string compile_options, std::string so_options);

  /// \brief Copy constructor
  Clexulator(const Clexulator &B);

  /// \brief Move constructor
  Clexulator(Clexulator &&B);

  ~Clexulator();

  /// \brief Assignment operator
  Clexulator &operator=(Clexulator B);

  /// \brief Swap
  friend void swap(Clexulator &first, Clexulator &second) {
    using std::swap;

    swap(first.m_name, second.m_name);
    swap(first.m_clex, second.m_clex);
    swap(first.m_lib, second.m_lib);
  }

  /// \brief Is runtime library loaded?
  bool initialized() const { return m_lib.get() != nullptr; }

  /// \brief Name
  std::string name() const { return m_name; }

  /// \brief Neighbor list size
  size_type nlist_size() const { return m_clex->nlist_size(); }

  /// \brief Number of correlations
  size_type corr_size() const { return m_clex->corr_size(); }

  /// \brief Valid range for `neighbor_ind` argument to calc_point_corr
  ///
  /// - For periodic clex, this is the number of sublattices.
  /// - For local clex, this is the neighbor list size.
  size_type n_point_corr() const { return m_clex->n_point_corr(); }

  /// \brief Obtain const reference to abstract ClexParamPack object
  ClexParamPack const &param_pack() const { return m_clex->param_pack(); }

  /// \brief Obtain reference to abstract ClexParamPack object
  ClexParamPack &param_pack() { return m_clex->param_pack(); }

  /// \brief Obtain ClexParamKey for a particular parameter
  ClexParamKey const &param_key(std::string const &_param_name) const;

  /// \brief Alter evaluation of parameters specified by @param _param_key,
  /// using a custom double -> double function set
  void set_evaluation(
      ClexParamKey const _param_key,
      std::vector<std::function<double(ConfigDoF const &)> > const &_basis_set);

  /// \brief Alter evaluation of parameters specified by @param _param_key,
  /// using a custom int -> double function set
  void set_evaluation(
      ClexParamKey const _param_key,
      std::vector<std::function<double(std::vector<double> const &)> > const
          &_basis_set);

  /// \brief Alter evaluation of parameters specified by @param _param_key,
  /// using the string  @param _eval_type, which can be at least either "READ"
  /// (i.e., read from ClexParamPack) or "DEFAULT" (i.e., the Clexulator's
  /// default implementation)
  void set_evaluation(ClexParamKey const _param_key, std::string _eval_type);

  /// \brief Check evaluation mode of parameters specified by @param _param_key,
  /// which can be one of (at least) "READ" (i.e., read from ClexParamPack),
  /// "CUSTOM", or "DEFAULT" (i.e., the Clexulator's default implementation)
  std::string check_evaluation(ClexParamKey const _param_key) const;

  /// \brief The UnitCellCoord involved in calculating the basis functions,
  /// relative origin UnitCell
  const std::set<UnitCell> &neighborhood() const {
    return m_clex->neighborhood();
  }

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  const std::set<UnitCell> &neighborhood(size_type linear_orbit_index) const {
    return m_clex->neighborhood(linear_orbit_index);
  }

  /// \brief The weight matrix used for ordering the neighbor list
  const PrimNeighborList::Matrix3Type &weight_matrix() const {
    return m_clex->weight_matrix();
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get
  /// contribution from int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration
  /// myclexulator.calc_global_corr_contribution(my_configdof,
  /// my_supercell.get_nlist(l_index).begin(), correlation_array.begin());
  /// \endcode
  ///
  void calc_global_corr_contribution(ConfigDoF const &_input_configdof,
                                     long int const *_nlist_begin,
                                     long int const *_nlist_end,
                                     double *_corr_begin,
                                     double *_corr_end) const {
    m_clex->calc_global_corr_contribution(_input_configdof, _nlist_begin,
                                          _nlist_end, _corr_begin, _corr_end);
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get
  /// contribution from int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration
  /// myclexulator.calc_global_corr_contribution(my_configdof,
  ///                                            my_supercell.get_nlist(l_index).begin(),
  ///                                            my_supercell.get_nlist(l_index).end());
  /// \endcode
  ///
  void calc_global_corr_contribution(ConfigDoF const &_input_configdof,
                                     long int const *_nlist_begin,
                                     long int const *_nlist_end) const {
    m_clex->calc_global_corr_contribution(_input_configdof, _nlist_begin,
                                          _nlist_end);
  }

  /// \brief Calculate contribution to select global correlations from one unit
  /// cell
  ///
  /// \param _corr_begin Pointer to beginning of data structure where
  /// correlations are written \param _corr_ind_begin,_corr_ind_end Pointers to
  /// range indicating which correlations should be calculated
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get
  /// contribution from int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration std::vector<int> _corr_ind = {0, 2, 4, 6}; // Get
  /// contribution to correlations 0, 2, 4, and 6
  /// myclexulator.calc_restricted_global_corr_contribution(my_configdof,
  ///                                                       my_supercell.get_nlist(l_index).begin(),
  ///                                                       correlation_array.begin(),
  ///                                                       _corr_ind.begin(),
  ///                                                       _corr_ind.end());
  /// \endcode
  ///
  void calc_restricted_global_corr_contribution(
      ConfigDoF const &_input_configdof, long int const *_nlist_begin,
      long int const *_nlist_end, double *_corr_begin, double *_corr_end,
      size_type const *_corr_ind_begin, size_type const *_corr_ind_end) const {
    m_clex->calc_restricted_global_corr_contribution(
        _input_configdof, _nlist_begin, _nlist_end, _corr_begin, _corr_end,
        _corr_ind_begin, _corr_ind_end);
  }

  /// \brief Calculate point correlations about basis site 'neighbor_ind'
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief _corr_begin Pointer to beginning of data structure where
  /// correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point
  /// correlations int l_index = my_supercell.find(bijk); // Linear index of
  /// site in Configuration myclexulator.calc_point_corr(my_configdof,
  /// my_supercell.get_nlist(l_index).begin(), b, correlation_array.begin());
  /// \endcode
  ///
  void calc_point_corr(ConfigDoF const &_input_configdof,
                       long int const *_nlist_begin, long int const *_nlist_end,
                       int neighbor_ind, double *_corr_begin,
                       double *_corr_end) const {
    m_clex->calc_point_corr(_input_configdof, _nlist_begin, _nlist_end,
                            neighbor_ind, _corr_begin, _corr_end);
  }

  /// \brief Calculate select point correlations about basis site 'neighbor_ind'
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief _corr_begin Pointer to beginning of data structure where
  /// correlations are written \param _corr_ind_begin,_corr_ind_end Pointers to
  /// range indicating which correlations should be calculated
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point
  /// correlations int l_index = my_supercell.find(bijk); // Linear index of
  /// site in Configuration std::vector<int> _corr_ind = {0, 2, 4, 6}; // Get
  /// contribution to correlations 0, 2, 4, and 6
  /// myclexulator.calc_restricted_point_corr(my_configdof,
  /// my_supercell.get_nlist(l_index).begin(), b, correlation_array.begin(),
  /// _corr_ind.begin(), _corr_ind.end()); \endcode
  ///
  void calc_restricted_point_corr(ConfigDoF const &_input_configdof,
                                  long int const *_nlist_begin,
                                  long int const *_nlist_end, int neighbor_ind,
                                  double *_corr_begin, double *_corr_end,
                                  size_type const *_corr_ind_begin,
                                  size_type const *_corr_ind_end) const {
    m_clex->calc_restricted_point_corr(
        _input_configdof, _nlist_begin, _nlist_end, neighbor_ind, _corr_begin,
        _corr_end, _corr_ind_begin, _corr_ind_end);
  }

  /// \brief Calculate the change in point correlations due to changing an
  /// occupant
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief occ_i,occ_f Initial and final occupant variable
  /// \brief _corr_begin Pointer to beginning of data structure where difference
  /// in correlations are written
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta
  /// point correlations int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration int occ_i=0, occ_f=1;  // Swap from occupant 0 to
  /// occupant 1 myclexulator.calc_delta_point_corr(my_configdof,
  /// my_supercell.get_nlist(l_index).begin(), b, occ_i, occ_f,
  /// correlation_array.begin()); \endcode
  ///
  void calc_delta_point_corr(ConfigDoF const &_input_configdof,
                             long int const *_nlist_begin,
                             long int const *_nlist_end, int neighbor_ind,
                             int occ_i, int occ_f, double *_corr_begin,
                             double *_corr_end) const {
    m_clex->calc_delta_point_corr(_input_configdof, _nlist_begin, _nlist_end,
                                  neighbor_ind, occ_i, occ_f, _corr_begin,
                                  _corr_end);
  }

  /// \brief Calculate the change in select point correlations due to changing
  /// an occupant
  ///
  /// \brief neighbor_ind Basis site index about which to calculate correlations
  /// \brief occ_i,occ_f Initial and final occupant variable
  /// \brief _corr_begin Pointer to beginning of data structure where difference
  /// in correlations are written \param _corr_ind_begin,_corr_ind_end Pointers
  /// to range indicating which correlations should be calculated
  ///
  /// Call using:
  /// \code
  /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta
  /// point correlations int l_index = my_supercell.find(bijk); // Linear index
  /// of site in Configuration int occ_i=0, occ_f=1;  // Swap from occupant 0 to
  /// occupant 1 std::vector<int> _corr_ind = {0, 2, 4, 6}; // Get contribution
  /// to correlations 0, 2, 4, and 6
  /// myclexulator.calc_restricted_delta_point_corr(my_configdof,
  ///                                               my_supercell.get_nlist(l_index).begin(),
  ///                                               b,
  ///                                               occ_i,
  ///                                               occ_f,
  ///                                               correlation_array.begin(),
  ///                                               correlation_array.end(),
  ///                                               _corr_ind.begin(),
  ///                                               _corr_ind.end());
  /// \endcode
  ///
  void calc_restricted_delta_point_corr(ConfigDoF const &_input_configdof,
                                        long int const *_nlist_begin,
                                        long int const *_nlist_end,
                                        int neighbor_ind, int occ_i, int occ_f,
                                        double *_corr_begin, double *_corr_end,
                                        size_type const *_corr_ind_begin,
                                        size_type const *_corr_ind_end) const {
    m_clex->calc_restricted_delta_point_corr(
        _input_configdof, _nlist_begin, _nlist_end, neighbor_ind, occ_i, occ_f,
        _corr_begin, _corr_end, _corr_ind_begin, _corr_ind_end);
  }

 private:
  std::string m_name;
  std::unique_ptr<Clexulator_impl::Base> m_clex;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

/**  @} */
}  // namespace CASM

#endif
