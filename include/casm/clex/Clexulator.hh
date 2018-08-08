#ifndef CLEXULATOR_HH
#define CLEXULATOR_HH
#include <cstddef>

#include "casm/external/boost.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  namespace Clexulator_impl {

    /// \brief Abstract base class for cluster expansion correlation calculations
    class Base {

    public:

      typedef unsigned int size_type;


      Base(size_type _nlist_size, size_type _corr_size) :
        m_nlist_size(_nlist_size),
        m_corr_size(_corr_size) {}

      virtual ~Base() {}

      /// \brief Neighbor list size
      size_type nlist_size() const {
        return m_nlist_size;
      }

      /// \brief Number of correlations
      size_type corr_size() const {
        return m_corr_size;
      }

      /// \brief Clone the Clexulator
      std::unique_ptr<Base> clone() const {
        return std::unique_ptr<Base>(_clone());
      }

      /// \brief The UnitCellCoord involved in calculating the basis functions,
      /// relative origin UnitCell
      const std::set<UnitCellCoord> &neighborhood() const {
        return m_neighborhood;
      }

      /// \brief The UnitCellCoord involved in calculating the basis functions
      /// for a particular orbit, relative origin UnitCell
      const std::set<UnitCellCoord> &neighborhood(size_type linear_orbit_index) const {
        return m_orbit_neighborhood[linear_orbit_index];
      }

      /// \brief The weight matrix used for ordering the neighbor list
      const PrimNeighborList::Matrix3Type &weight_matrix() const {
        return m_weight_matrix;
      }

      /// \brief Set pointer to data structure containing occupation variables
      ///
      /// \param _occ_ptr Pointer to beginning of data structure containing occupation variables
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// \endcode
      ///
      void set_config_occ(const int *_occ_ptr) {
        m_occ_ptr = _occ_ptr;
      }

      /// \brief Set pointer to neighbor list
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(b,i,j,k);           // UnitCellCoord of site in Configuration
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// \endcode
      ///
      void set_nlist(const long int *_nlist_ptr) {
        m_nlist_ptr = _nlist_ptr;
        return;
      };

      /// \brief Calculate contribution to global correlations from one unit cell
      ///
      /// \param corr_begin Pointer to beginning of data structure where correlations are written
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// myclexulator.calc_global_corr_contribution(correlation_array.begin());
      /// \endcode
      ///
      virtual void calc_global_corr_contribution(double *corr_begin) const = 0;

      /// \brief Calculate contribution to select global correlations from one unit cell
      ///
      /// \param corr_begin Pointer to beginning of data structure where correlations are written
      /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
      /// myclexulator.calc_restricted_global_corr_contribution(correlation_array.begin(), ind_list.begin(), ind_list.end());
      /// \endcode
      ///
      virtual void calc_restricted_global_corr_contribution(double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const = 0;

      /// \brief Calculate point correlations about basis site 'b_index'
      ///
      /// \brief b_index Basis site index about which to calculate correlations
      /// \brief corr_begin Pointer to beginning of data structure where correlations are written
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// myclexulator.calc_point_corr(b, correlation_array.begin());
      /// \endcode
      ///
      virtual void calc_point_corr(int b_index, double *corr_begin) const = 0;

      /// \brief Calculate select point correlations about basis site 'b_index'
      ///
      /// \brief b_index Basis site index about which to calculate correlations
      /// \brief corr_begin Pointer to beginning of data structure where correlations are written
      /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
      /// myclexulator.calc_restricted_point_corr(b, correlation_array.begin(), ind_list.begin(), ind_list.end());
      /// \endcode
      ///
      virtual void calc_restricted_point_corr(int b_index, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const = 0;

      /// \brief Calculate the change in point correlations due to changing an occupant
      ///
      /// \brief b_index Basis site index about which to calculate correlations
      /// \brief occ_i,occ_f Initial and final occupant variable
      /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
      /// myclexulator.calc_delta_point_corr(b, occ_i, occ_f, correlation_array.begin());
      /// \endcode
      ///
      virtual void calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const = 0;

      /// \brief Calculate the change in select point correlations due to changing an occupant
      ///
      /// \brief b_index Basis site index about which to calculate correlations
      /// \brief occ_i,occ_f Initial and final occupant variable
      /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
      /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
      ///
      /// Call using:
      /// \code
      /// myclexulator.set_config_occ(my_configdof.occupation().begin());
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
      /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
      /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
      /// myclexulator.calc_restricted_delta_point_corr(b, occ_i, occ_f, correlation_array.begin(), ind_list.begin(), ind_list.end());
      /// \endcode
      ///
      virtual void calc_restricted_delta_point_corr(int b_index,
                                                    int occ_i,
                                                    int occ_f,
                                                    double *corr_begin,
                                                    size_type const *ind_list_begin,
                                                    size_type const *ind_list_end) const = 0;


    private:

      /// \brief Clone the Clexulator
      virtual Base *_clone() const = 0;

      /// \brief The neighbor list size
      size_type m_nlist_size;

      /// \brief The number of correlations
      size_type m_corr_size;


    protected:

      /// \brief Pointer to beginning of data structure containing occupation variables
      const int *m_occ_ptr;

      /// \brief Pointer to neighbor list
      const long int *m_nlist_ptr;

      /// \brief The UnitCellCoord involved in calculating the basis functions,
      /// relative origin UnitCell
      std::set<UnitCellCoord> m_neighborhood;

      /// \brief The UnitCellCoord involved in calculating the basis functions
      /// for a particular orbit, relative origin UnitCell
      std::vector<std::set<UnitCellCoord> > m_orbit_neighborhood;

      /// \brief The weight matrix used for ordering the neighbor list
      PrimNeighborList::Matrix3Type m_weight_matrix;

    };
  }

  /// \brief Evaluates correlations
  ///
  /// CASM generates code for very efficient calculation of basis functions via
  /// the print_clexulator function. This source code may be compiled, linked,
  /// and used at runtime via Clexulator.
  ///
  /// \ingroup ClexClex
  ///
  class Clexulator {

  public:

    typedef Clexulator_impl::Base::size_type size_type;


    Clexulator() {}

    /// \brief Construct a Clexulator
    ///
    /// \param name Class name for the Clexulator, typically 'X_Clexulator', with X
    ///             referring to the system of interest (i.e. 'NiAl_Clexulator')
    /// \param dirpath Directory containing the source code and compiled object file.
    /// \param nlist, A PrimNeighborList to be updated to include the neighborhood
    ///        of this Clexulator
    /// \param logging Print messages to inform users that compilation is occuring
    /// \param compile_options Compilation options
    /// \param so_options Shared library compilation options
    ///
    /// If 'name' is 'X_Clexulator', and 'dirpath' is '/path/to':
    /// - Looks for '/path/to/X_Clexulator.so' and tries to load it.
    /// - If not found, looks for 'path/to/X_Clexulator.cc' and tries to compile and load it.
    /// - If unsuccesful, will throw std::runtime_error.
    ///
    /// The Clexulator has shared ownership of the loaded library,
    /// so it is preferrable to duplicate the Clexulator using it's copy constructor rather
    /// than construct another using this constructor which will re-load the library.
    ///
    Clexulator(std::string name,
               boost::filesystem::path dirpath,
               PrimNeighborList &nlist,
               const Logging &logging,
               std::string compile_options,
               std::string so_options) {

      namespace fs = boost::filesystem;

      // Construct the RuntimeLibrary that will store the loaded clexulator library
      try {
        m_lib = std::make_shared<RuntimeLibrary>(
                  (dirpath / name).string(),
                  compile_options,
                  so_options,
                  "compile time depends on how many basis functions are included");
      }
      catch(std::exception &e) {
        logging.log() << "Clexulator construction failed: could not construct runtime library." << std::endl;
        throw;
      }

      // Get the Clexulator factory function
      std::function<Clexulator_impl::Base* (void)> factory;
      factory = m_lib->get_function<Clexulator_impl::Base* (void)>("make_" + name);

      // Use the factory to construct the clexulator and store it in m_clex
      m_clex.reset(factory());

      // Check nlist has the right weight_matrix
      if(nlist.weight_matrix() != m_clex->weight_matrix()) {
        std::cerr << "Error in Clexulator constructor: weight matrix of neighbor "
                  "list does not match the weight matrix used to print the "
                  "clexulator." << std::endl;
        std::cerr << "nlist weight matrix: \n" << nlist.weight_matrix() << std::endl;
        std::cerr << "clexulator weight matrix: \n" << m_clex->weight_matrix() << std::endl;
        throw std::runtime_error(
          "Error in Clexulator constructor: weight matrix of neighbor list does "
          "not match the weight matrix used to print the clexulator. Try 'casm bset -uf'.");
      }

      // Expand the given neighbor list as necessary
      nlist.expand(neighborhood().begin(), neighborhood().end());

    }


    /// \brief Copy constructor
    Clexulator(const Clexulator &B) :
      m_name(B.name()),
      m_lib(B.m_lib) {

      if(B.m_clex.get() != nullptr) {
        m_clex.reset(B.m_clex->clone().release());
      }
    }

    /// \brief Move constructor
    Clexulator(Clexulator &&B) {
      swap(*this, B);
    }

    ~Clexulator() {
      // ensure Clexulator is deleted before library
      delete m_clex.release();
    }

    /// \brief Assignment operator
    Clexulator &operator=(Clexulator B) {

      swap(*this, B);

      return *this;
    }

    /// \brief Swap
    friend void swap(Clexulator &first, Clexulator &second) {

      using std::swap;

      swap(first.m_name, second.m_name);
      swap(first.m_clex, second.m_clex);
      swap(first.m_lib, second.m_lib);
    }

    /// \brief Is runtime library loaded?
    bool initialized() const {
      return m_lib.get() != nullptr;
    }

    /// \brief Name
    std::string name() const {
      return m_name;
    }

    /// \brief Neighbor list size
    size_type nlist_size() const {
      return m_clex->nlist_size();
    }

    /// \brief Number of correlations
    size_type corr_size() const {
      return m_clex->corr_size();
    }

    /// \brief The UnitCellCoord involved in calculating the basis functions,
    /// relative origin UnitCell
    const std::set<UnitCellCoord> &neighborhood() const {
      return m_clex->neighborhood();
    }

    /// \brief The UnitCellCoord involved in calculating the basis functions
    /// for a particular orbit, relative origin UnitCell
    const std::set<UnitCellCoord> &neighborhood(size_type linear_orbit_index) const {
      return m_clex->neighborhood(linear_orbit_index);
    }

    /// \brief The weight matrix used for ordering the neighbor list
    const PrimNeighborList::Matrix3Type &weight_matrix() const {
      return m_clex->weight_matrix();
    }

    /// \brief Set pointer to data structure containing occupation variables
    ///
    /// \param _occ_ptr Pointer to beginning of data structure containing occupation variables
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// \endcode
    ///
    void set_config_occ(const int *_occ_ptr) {
      return m_clex->set_config_occ(_occ_ptr);
    }

    /// \brief Set pointer to neighbor list
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(b,i,j,k);           // UnitCellCoord of site in Configuration
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// \endcode
    ///
    void set_nlist(const long int *_nlist_ptr) {
      return m_clex->set_nlist(_nlist_ptr);
    };

    /// \brief Calculate contribution to global correlations from one unit cell
    ///
    /// \param corr_begin Pointer to beginning of data structure where correlations are written
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// myclexulator.calc_global_corr_contribution(correlation_array.begin());
    /// \endcode
    ///
    void calc_global_corr_contribution(double *corr_begin) const {
      m_clex->calc_global_corr_contribution(corr_begin);
    }

    /// \brief Calculate contribution to select global correlations from one unit cell
    ///
    /// \param corr_begin Pointer to beginning of data structure where correlations are written
    /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
    /// myclexulator.calc_restricted_global_corr_contribution(correlation_array.begin(), ind_list.begin(), ind_list.end());
    /// \endcode
    ///
    void calc_restricted_global_corr_contribution(double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
      m_clex->calc_restricted_global_corr_contribution(corr_begin, ind_list_begin, ind_list_end);
    }

    /// \brief Calculate point correlations about basis site 'b_index'
    ///
    /// \brief b_index Basis site index about which to calculate correlations
    /// \brief corr_begin Pointer to beginning of data structure where correlations are written
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// myclexulator.calc_point_corr(b, correlation_array.begin());
    /// \endcode
    ///
    void calc_point_corr(int b_index, double *corr_begin) const {
      m_clex->calc_point_corr(b_index, corr_begin);
    }

    /// \brief Calculate select point correlations about basis site 'b_index'
    ///
    /// \brief b_index Basis site index about which to calculate correlations
    /// \brief corr_begin Pointer to beginning of data structure where correlations are written
    /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
    /// myclexulator.calc_restricted_point_corr(b, correlation_array.begin(), ind_list.begin(), ind_list.end());
    /// \endcode
    ///
    void calc_restricted_point_corr(int b_index, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
      m_clex->calc_restricted_point_corr(b_index, corr_begin, ind_list_begin, ind_list_end);
    }

    /// \brief Calculate the change in point correlations due to changing an occupant
    ///
    /// \brief b_index Basis site index about which to calculate correlations
    /// \brief occ_i,occ_f Initial and final occupant variable
    /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
    /// myclexulator.calc_delta_point_corr(b, occ_i, occ_f, correlation_array.begin());
    /// \endcode
    ///
    void calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const {
      m_clex->calc_delta_point_corr(b_index, occ_i, occ_f, corr_begin);
    }

    /// \brief Calculate the change in select point correlations due to changing an occupant
    ///
    /// \brief b_index Basis site index about which to calculate correlations
    /// \brief occ_i,occ_f Initial and final occupant variable
    /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
    /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
    ///
    /// Call using:
    /// \code
    /// myclexulator.set_config_occ(my_configdof.occupation().begin());
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.set_nlist(my_supercell.get_nlist(l_index).begin());
    /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
    /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
    /// myclexulator.calc_restricted_delta_point_corr(b, occ_i, occ_f, correlation_array.begin(), ind_list.begin(), ind_list.end());
    /// \endcode
    ///
    void calc_restricted_delta_point_corr(int b_index,
                                          int occ_i,
                                          int occ_f,
                                          double *corr_begin,
                                          size_type const *ind_list_begin,
                                          size_type const *ind_list_end) const {
      m_clex->calc_restricted_delta_point_corr(b_index, occ_i, occ_f, corr_begin, ind_list_begin, ind_list_end);
    }


  private:

    std::string m_name;
    std::unique_ptr<Clexulator_impl::Base> m_clex;
    std::shared_ptr<RuntimeLibrary> m_lib;

  };

}

#endif
