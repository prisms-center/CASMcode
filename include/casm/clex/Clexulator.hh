#ifndef CLEXULATOR_HH
#define CLEXULATOR_HH
#include <cstddef>

#include "casm/external/boost.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/ConfigDoF.hh"
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

      /// \brief Calculate contribution to global correlations from one unit cell
      ///
      /// \param corr_begin Pointer to beginning of data structure where correlations are written
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.calc_global_corr_contribution(my_configdof,
      ///                                            my_supercell.get_nlist(l_index).begin(),
      ///                                            correlation_array.begin());
      /// \endcode
      ///
      void calc_global_corr_contribution(ConfigDoF const &_config_dof,
                                         long int const *_n_list_begin,
                                         double *_corr_begin) const {
        m_config_ptr = &_config_dof;
        _set_nlist(_n_list_begin);
        _global_prepare();
        _calc_global_corr_contribution(corr_begin);
      }

      /// \brief Calculate contribution to select global correlations from one unit cell
      ///
      /// \param corr_begin Pointer to beginning of data structure where correlations are written
      /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
      /// myclexulator.calc_restricted_global_corr_contribution(my_configdof,
      ///                                                       my_supercell.get_nlist(l_index).begin(),
      ///                                                       correlation_array.begin(),
      ///                                                       ind_list.begin(),
      ///                                                       ind_list.end());
      /// \endcode
      ///
      void calc_restricted_global_corr_contribution(ConfigDoF const &_config_dof,
                                                    long int const *_n_list_begin,
                                                    double *corr_begin,
                                                    size_type const *ind_list_begin,
                                                    size_type const *ind_list_end) const {
        m_config_ptr = &_config_dof;
        _set_nlist(_n_list_begin);
        _global_prepare();
        _calc_restricted_global_corr_contribution(corr_begin, ind_list_begin, ind_list_end);
      }

      /// \brief Calculate point correlations about basis site 'sublat_ind'
      ///
      /// \brief sublat_ind Basis site index about which to calculate correlations
      /// \brief corr_begin Pointer to beginning of data structure where correlations are written
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// myclexulator.calc_point_corr(my_configdof, my_supercell.get_nlist(l_index).begin(), b, correlation_array.begin());
      /// \endcode
      ///
      void calc_point_corr(ConfigDoF const &_config_dof,
                           long int const *_n_list_begin,
                           int sublat_ind,
                           double *corr_begin) const {
        m_config_ptr = &_config_dof;
        _set_nlist(_n_list_begin);
        _point_prepare(sublat_ind);
        _calc_point_corr(sublat_ind, corr_begin);
      }

      /// \brief Calculate select point correlations about basis site 'sublat_ind'
      ///
      /// \brief sublat_ind Basis site index about which to calculate correlations
      /// \brief corr_begin Pointer to beginning of data structure where correlations are written
      /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
      /// myclexulator.calc_restricted_point_corr(my_configdof,
      ///                                         my_supercell.get_nlist(l_index).begin(),
      ///                                         b,
      ///                                         correlation_array.begin(),
      ///                                         ind_list.begin(),
      ///                                         ind_list.end());
      /// \endcode
      ///
      void calc_restricted_point_corr(ConfigDoF const &_config_dof,
                                      long int const *_n_list_begin,
                                      int sublat_ind,
                                      double *corr_begin,
                                      size_type const *ind_list_begin,
                                      size_type const *ind_list_end) const {

        m_config_ptr = &_config_dof;
        _set_nlist(_n_list_begin);
        _point_prepare(sublat_ind);
        _calc_restricted_point_corr(sublat_ind, corr_begin, ind_list_begin, ind_list_end);
      }

      /// \brief Calculate the change in point correlations due to changing an occupant
      ///
      /// \brief sublat_ind Basis site index about which to calculate correlations
      /// \brief occ_i,occ_f Initial and final occupant variable
      /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
      /// myclexulator.calc_delta_point_corr(my_configdof, my_supercell.get_nlist(l_index).begin(), b, occ_i, occ_f, correlation_array.begin());
      /// \endcode
      ///
      void calc_delta_point_corr(ConfigDoF const &_config_dof,
                                 long int const *_n_list_begin,
                                 int sublat_ind,
                                 int occ_i,
                                 int occ_f,
                                 double *corr_begin) const {
        m_config_ptr = &_config_dof;
        _set_nlist(_n_list_begin);
        _point_prepare(sublat_ind);
        _calc_delta_point_corr(sublat_ind, occ_i, occ_f, corr_begin);
      }
      /// \brief Calculate the change in select point correlations due to changing an occupant
      ///
      /// \brief sublat_ind Basis site index about which to calculate correlations
      /// \brief occ_i,occ_f Initial and final occupant variable
      /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
      /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
      /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
      /// myclexulator.calc_restricted_delta_point_corr(my_configdof,
      ///                                               my_supercell.get_nlist(l_index).begin(),
      ///                                               b,
      ///                                               occ_i,
      ///                                               occ_f,
      ///                                               correlation_array.begin(),
      ///                                               ind_list.begin(),
      ///                                               ind_list.end());
      /// \endcode
      ///
      void calc_restricted_delta_point_corr(ConfigDoF const &_config_dof,
                                            long int const *_n_list_begin,
                                            int sublat_ind,
                                            int occ_i,
                                            int occ_f,
                                            double *corr_begin,
                                            size_type const *ind_list_begin,
                                            size_type const *ind_list_end) const {
        m_config_ptr = &_config_dof;
        _set_nlist(_n_list_begin);
        _point_prepare(sublat_ind);
        _calc_delta_point_corr(sublat_ind, occ_i, occ_f, corr_begin, ind_list_begin, ind_list_end);
      }


    private:

      /// \brief Clone the Clexulator
      virtual Base *_clone() const = 0;

      /// \brief The neighbor list size
      size_type m_nlist_size;

      /// \brief The number of correlations
      size_type m_corr_size;


    protected:

      virtual void _calc_global_corr_contribution(double *corr_begin) const = 0;
      virtual void _calc_restricted_global_corr_contribution(double *corr_begin,
                                                             size_type const *ind_list_begin,
                                                             size_type const *ind_list_end) const = 0;

      virtual void _calc_point_corr(int sublat_ind,
                                    double *corr_begin) const = 0;

      virtual void _calc_restricted_point_corr(int sublat_ind,
                                               double *corr_begin,
                                               size_type const *ind_list_begin,
                                               size_type const *ind_list_end) const = 0;

      virtual void _calc_delta_point_corr(int sublat_ind,
                                          int occ_i,
                                          int occ_f,
                                          double *corr_begin) const = 0;

      virtual void _calc_delta_point_corr(int sublat_ind,
                                          int occ_i,
                                          int occ_f,
                                          double *corr_begin,
                                          size_type const *ind_list_begin,
                                          size_type const *ind_list_end) const = 0;

      virtual void _global_prepare() const = 0;

      virtual void _point_prepare(int sublat_ind) const = 0;


      /// \brief Set pointer to neighbor list
      ///
      /// Call using:
      /// \code
      /// UnitCellCoord bijk(b,i,j,k);           // UnitCellCoord of site in Configuration
      /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
      /// \endcode
      ///
      void _set_nlist(const long int *_nlist_ptr) {
        m_nlist_ptr = _nlist_ptr;
        return;
      }

      /// \brief Pointer to ConfigDoF for which evaluation is occuring
      ConfigDoF const *m_config_ptr;

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
  /// \ingroup Clex
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
    /// \param status_log Print a message to inform users that compilation is occuring
    /// \param compile_options Compilation options, by default "g++ -O3 -Wall -fPIC"
    /// \param so_options Shared library compilation options, by default "g++ -shared"
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
               Log &status_log = null_log(),
               std::string compile_options = RuntimeLibrary::default_compile_options(),
               std::string so_options = RuntimeLibrary::default_so_options()) {

      namespace fs = boost::filesystem;

      try {

        // Construct the RuntimeLibrary that will store the loaded clexulator library
        m_lib = std::make_shared<RuntimeLibrary>(compile_options, so_options);

        // If the shared library doesn't exist
        if(!fs::exists(dirpath / (name + ".so"))) {

          // But the library source code does
          if(fs::exists(dirpath / (name + ".cc"))) {

            // Compile it
            try {
              status_log.compiling<Log::standard>(name);
              status_log.begin_lap();
              status_log << "compile time depends on how many basis functions are included" << std::endl;
              m_lib->compile((dirpath / name).string());
              status_log << "compile time: " << status_log.lap_time() << " (s)\n" << std::endl;
            }
            catch(std::exception &e) {
              status_log << "Error compiling clexulator. To fix: \n";
              status_log << "  - Check compiler error messages.\n";
              status_log << "  - Check compiler options with 'casm settings -l'\n";
              status_log << "    - Update compiler options with 'casm settings --set-compile-options '...options...'\n";
              status_log << "    - Make sure the casm headers can be found by including '-I/path/to/casm'\n";
              status_log << "  - The default compiler is 'g++'. Override by setting the environment variable CXX\n" << std::endl;
              throw e;
            }
          }
          else {
            throw std::runtime_error(
              std::string("Error in Clexulator constructor\n") +
              "  Could not find '" + dirpath.string() + "/" + name + ".so' or '" + dirpath.string() + "/" + name + ".cc'");
          }
        }

        // If the shared library exists
        if(fs::exists(dirpath / (name + ".so"))) {

          // Load the library with the Clexulator
          m_lib->load((dirpath / name).string());

          // Get the Clexulator factory function
          std::function<Clexulator_impl::Base* (void)> factory;
          factory = m_lib->get_function<Clexulator_impl::Base* (void)>("make_" + name);

          // Use the factory to construct the clexulator and store it in m_clex
          m_clex.reset(factory());
        }
        else {
          throw std::runtime_error(
            std::string("Error in Clexulator constructor\n") +
            "  Did not find '" + dirpath.string() + "/" + name + ".so'");
        }

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
      catch(const std::exception &e) {
        std::cout << "Error in Clexulator constructor" << std::endl;
        std::cout << e.what() << std::endl;
        throw e;
      }

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

    /// \brief Calculate contribution to global correlations from one unit cell
    ///
    /// \param corr_begin Pointer to beginning of data structure where correlations are written
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.calc_global_corr_contribution(my_configdof, my_supercell.get_nlist(l_index).begin(), correlation_array.begin());
    /// \endcode
    ///
    void calc_global_corr_contribution(ConfigDoF const &_config_dof,
                                       long int const *_n_list_begin,
                                       double *corr_begin) const {
      m_clex->calc_global_corr_contribution(_config_dof,
                                            _n_list_begin,
                                            corr_begin);
    }

    /// \brief Calculate contribution to select global correlations from one unit cell
    ///
    /// \param corr_begin Pointer to beginning of data structure where correlations are written
    /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(0,i,j,k);           // i,j,k of unit cell to get contribution from
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
    /// myclexulator.calc_restricted_global_corr_contribution(my_configdof,
    ///                                                       my_supercell.get_nlist(l_index).begin(),
    ///                                                       correlation_array.begin(),
    ///                                                       ind_list.begin(),
    ///                                                       ind_list.end());
    /// \endcode
    ///
    void calc_restricted_global_corr_contribution(ConfigDoF const &_config_dof,
                                                  long int const *_n_list_begin,
                                                  double *corr_begin,
                                                  size_type const *ind_list_begin,
                                                  size_type const *ind_list_end) const {
      m_clex->calc_restricted_global_corr_contribution(_config_dof,
                                                       _n_list_begin,
                                                       corr_begin,
                                                       ind_list_begin,
                                                       ind_list_end);
    }

    /// \brief Calculate point correlations about basis site 'sublat_ind'
    ///
    /// \brief sublat_ind Basis site index about which to calculate correlations
    /// \brief corr_begin Pointer to beginning of data structure where correlations are written
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// myclexulator.calc_point_corr(my_configdof, my_supercell.get_nlist(l_index).begin(), b, correlation_array.begin());
    /// \endcode
    ///
    void calc_point_corr(ConfigDoF const &_config_dof,
                         long int const *_n_list_begin,
                         int sublat_ind,
                         double *corr_begin) const {
      m_clex->calc_point_corr(_config_dof,
                              _n_list_begin,
                              sublat_ind,
                              corr_begin);
    }

    /// \brief Calculate select point correlations about basis site 'sublat_ind'
    ///
    /// \brief sublat_ind Basis site index about which to calculate correlations
    /// \brief corr_begin Pointer to beginning of data structure where correlations are written
    /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
    /// myclexulator.calc_restricted_point_corr(my_configdof, my_supercell.get_nlist(l_index).begin(), b, correlation_array.begin(), ind_list.begin(), ind_list.end());
    /// \endcode
    ///
    void calc_restricted_point_corr(ConfigDoF const &_config_dof,
                                    long int const *_n_list_begin,
                                    int sublat_ind,
                                    double *corr_begin,
                                    size_type const *ind_list_begin,
                                    size_type const *ind_list_end) const {
      m_clex->calc_restricted_point_corr(_config_dof,
                                         _n_list_begin,
                                         sublat_ind,
                                         corr_begin,
                                         ind_list_begin,
                                         ind_list_end);
    }

    /// \brief Calculate the change in point correlations due to changing an occupant
    ///
    /// \brief sublat_ind Basis site index about which to calculate correlations
    /// \brief occ_i,occ_f Initial and final occupant variable
    /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
    /// myclexulator.calc_delta_point_corr(my_configdof, my_supercell.get_nlist(l_index).begin(), b, occ_i, occ_f, correlation_array.begin());
    /// \endcode
    ///
    void calc_delta_point_corr(ConfigDoF const &_config_dof,
                               long int const *_n_list_begin,
                               int sublat_ind,
                               int occ_i,
                               int occ_f,
                               double *corr_begin) const {
      m_clex->calc_delta_point_corr(_config_dof,
                                    _n_list_begin,
                                    sublat_ind,
                                    occ_i,
                                    occ_f,
                                    corr_begin);
    }

    /// \brief Calculate the change in select point correlations due to changing an occupant
    ///
    /// \brief sublat_ind Basis site index about which to calculate correlations
    /// \brief occ_i,occ_f Initial and final occupant variable
    /// \brief corr_begin Pointer to beginning of data structure where difference in correlations are written
    /// \param ind_list_begin,ind_list_end Pointers to range indicating which correlations should be calculated
    ///
    /// Call using:
    /// \code
    /// UnitCellCoord bijk(b,i,j,k);           // b,i,j,k of site to get delta point correlations
    /// int l_index = my_supercell.find(bijk); // Linear index of site in Configuration
    /// int occ_i=0, occ_f=1;  // Swap from occupant 0 to occupant 1
    /// std::vector<int> ind_list = {0, 2, 4, 6}; // Get contribution to correlations 0, 2, 4, and 6
    /// myclexulator.calc_restricted_delta_point_corr(my_configdof,
    ///                                               my_supercell.get_nlist(l_index).begin(),
    ///                                               b,
    ///                                               occ_i,
    ///                                               occ_f,
    ///                                               correlation_array.begin(),
    ///                                               ind_list.begin(),
    ///                                               ind_list.end());
    /// \endcode
    ///
    void calc_restricted_delta_point_corr(ConfigDoF const &_config_dof,
                                          long int const *_n_list_begin,
                                          int sublat_ind,
                                          int occ_i,
                                          int occ_f,
                                          double *corr_begin,
                                          size_type const *ind_list_begin,
                                          size_type const *ind_list_end) const {
      m_clex->calc_restricted_delta_point_corr(_config_dof,
                                               _n_list_begin,
                                               sublat_ind,
                                               occ_i,
                                               occ_f,
                                               corr_begin,
                                               ind_list_begin,
                                               ind_list_end);
    }


  private:
    std::string m_name;
    std::unique_ptr<Clexulator_impl::Base> m_clex;
    std::shared_ptr<RuntimeLibrary> m_lib;
  };

}

#endif

