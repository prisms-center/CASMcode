#include "casm/clex/ECIContainer.hh"
namespace CASM {

  ECIContainer::ECIContainer(const fs::path &eci_fit_file) {
    if(!fs::exists(eci_fit_file))
      throw std::runtime_error("ECIContainer initialized with path " + eci_fit_file.string() + ", which does not exist.\n");


    ECIContainer_impl::populate_eci(eci_fit_file, m_eci_list, m_eci_index_list);
  }

  double operator*(const ECIContainer &_eci, const Correlation &_corr) {
    double result(0);
    auto ind_it(_eci.eci_index_list().cbegin()), ind_end(_eci.eci_index_list().cend());
    auto eci_it(_eci.eci_list().cbegin());
    for(; ind_it != ind_end; ++ind_it, ++eci_it)
      result += (*eci_it) * _corr[*ind_it];
    return result;
  }

  namespace ECIContainer_impl {
    /**
     * Read from specified path and fill up ECI values and indices.
     * Used primarily by MonteCarlo::populate_eci().
     *
     * I'm sure the format of eci.out is going to do nothing but change, so keep
     * an eye on this one. The current format is:
     *
     * Num clusters
     * Num structures
     * rms value
     * WCV structures
     * WCV score
     * some other rms value
     * #Header for ECI (we use these)   ECI/mult (ignore)   Cluster number (this is the index)
     * value                            value               value
     * value                            value               value
     * .
     * .
     * .
     * value                            value               value
     */

    void populate_eci(const fs::path &filepath, ECIContainer::ScalarECI &mc_eci, Array<ECIContainer::size_type> &mc_eci_index) {
      //Stop program if either of the arrays contains values already
      if(mc_eci.size() != 0 || mc_eci_index.size() != 0) {
        std::cerr << "ERROR in Monte::populate_eci" << std::endl;
        std::cerr << "Expected to be given empty arrays to fill but" << std::endl;
        std::cerr << "The eci array was size " << mc_eci.size() << std::endl;
        std::cerr << "and the eci index array was size " << mc_eci_index.size() << std::endl;
        exit(371);
      }

      //read in from eci,out file
      std::ifstream ecistream(filepath.string().c_str());

      //This is where we dump each line of the file as we read
      std::string stringdump;

      //skip all the lines until the actual eci values
      for(int i = 0; i < 7; i++) {
        std::getline(ecistream, stringdump);
      }

      //read the three columns of eci, eci/mult and index line by line
      while(std::getline(ecistream, stringdump)) {
        //make a stringstream to split the three values of the line we just read
        std::istringstream eciline_stream(stringdump);
        double eci_value, eci_over_mult;
        Index eci_index;

        eciline_stream >> eci_value;
        eciline_stream >> eci_over_mult;
        eciline_stream >> eci_index;

        //store the values we just read
        mc_eci.push_back(eci_value);
        mc_eci_index.push_back(eci_index);
      }
      return;
    }

  }

}
