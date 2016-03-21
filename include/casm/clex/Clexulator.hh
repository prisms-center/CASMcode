#ifndef CLEXULATOR_HH
#define CLEXULATOR_HH
#include <cstddef>

#define BOOST_NO_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include "casm/system/RuntimeLibrary.hh"

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
            m_lib->compile((dirpath / name).string());
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

