#ifndef CLEXULATOR_HH
#define CLEXULATOR_HH
#include <boost/filesystem/path.hpp>
#include <cstddef>

#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clexulator/BaseClexulator.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
/** \defgroup Clexulator
 *  \ingroup ClexClex
 *  \brief Clexulator is a class for evaluating crystal basis functions and
 * cluster expansions
 *
 *  @{
 */

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
  typedef clexulator::BaseClexulator::size_type size_type;

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
  clexulator::ClexParamPack const &param_pack() const {
    return m_clex->param_pack();
  }

  /// \brief Obtain reference to abstract ClexParamPack object
  clexulator::ClexParamPack &param_pack() { return m_clex->param_pack(); }

  /// \brief Obtain ClexParamKey for a particular parameter
  clexulator::ClexParamKey const &param_key(
      std::string const &_param_name) const;

  /// \brief The UnitCellCoord involved in calculating the basis functions,
  /// relative origin UnitCell
  const std::set<xtal::UnitCell> &neighborhood() const {
    return m_clex->neighborhood();
  }

  /// \brief The UnitCellCoord involved in calculating the basis functions
  /// for a particular orbit, relative origin UnitCell
  const std::set<xtal::UnitCell> &neighborhood(
      size_type linear_orbit_index) const {
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_global_corr_contribution(_corr_begin);
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_global_corr_contribution();
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_restricted_global_corr_contribution(
        _corr_begin, _corr_ind_begin, _corr_ind_end);
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_point_corr(neighbor_ind, _corr_begin);
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_restricted_point_corr(neighbor_ind, _corr_begin,
                                       _corr_ind_begin, _corr_ind_end);
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_delta_point_corr(neighbor_ind, occ_i, occ_f, _corr_begin);
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
    m_clex->set_configdofvalues(_input_configdof.values());
    m_clex->set_nlist(_nlist_begin);
    m_clex->calc_restricted_delta_point_corr(neighbor_ind, occ_i, occ_f,
                                             _corr_begin, _corr_ind_begin,
                                             _corr_ind_end);
  }

 private:
  std::string m_name;
  std::unique_ptr<clexulator::BaseClexulator> m_clex;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

/**  @} */
}  // namespace CASM

#endif
