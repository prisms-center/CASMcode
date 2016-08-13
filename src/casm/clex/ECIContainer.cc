#include "casm/clex/ECIContainer.hh"

namespace CASM {

  /// \brief Evaluate property given an ECIContainer and Correlation
  double operator*(const ECIContainer &_eci, const Correlation &_corr) {
    double result(0);
    auto ind_it(_eci.index().cbegin()), ind_end(_eci.index().cend());
    auto eci_it(_eci.value().cbegin());
    for(; ind_it != ind_end; ++ind_it, ++eci_it)
      result += (*eci_it) * _corr[*ind_it];
    return result;
  }

  /// \brief Evaluate property given an ECIContainer and pointer to beginning of range of correlation
  double operator*(const ECIContainer &_eci, double const *_corr_begin) {
    double result(0);
    auto ind_it(_eci.index().cbegin()), ind_end(_eci.index().cend());
    auto eci_it(_eci.value().cbegin());
    while(ind_it != ind_end) {
      result += (*eci_it) * (*(_corr_begin + *ind_it));
      ++ind_it;
      ++eci_it;
    }
    return result;
  }

  /// \brief Read eci.out file from specified path (deprecated)
  ///
  /// I'm sure the format of eci.out is going to do nothing but change, so keep
  /// an eye on this one. The current format is:
  /// \code
  /// Num clusters
  /// Num structures
  /// rms value
  /// WCV structures
  /// WCV score
  /// some other rms value
  /// #Header for ECI (we use these)   ECI/mult (ignore)   Cluster number (this is the index)
  /// value                            value               value
  /// value                            value               value
  /// .
  /// .
  /// .
  /// value                            value               value
  /// \endcode
  ///
  ECIContainer read_eci_out(const fs::path &filepath) {

    std::vector<double> value;
    std::vector<ECIContainer::size_type> index;

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
      value.push_back(eci_value);
      index.push_back(eci_index);
    }
    return ECIContainer(value.begin(), value.end(), index.begin());
  }


  /// \brief Read eci.json file from specified path
  ///
  /// Format:
  /// \code
  /// {
  ///   "site_functions":[
  ///     {
  ///       "asym_unit": X,
  ///       "sublat_indices: [2, 3],
  ///       "phi_b_0": {"Va":0.0, "O":1.0},
  ///       "phi_b_1": {"Va":0.0, "O":1.0},
  ///        ...
  ///     },
  ///     ...
  ///   ],
  ///   "cluster_functions":[
  ///     {
  ///       "eci": X.XXXXX,
  ///       "prototype_function": "\phi_b_i(s_j)...",
  ///       "orbit": [branch_index, orbit_index],
  ///       "linear_function_index": I,
  ///       "mult": X,
  ///       "prototype": [
  ///         [b, i, j, k],
  ///         ...
  ///       ]
  ///     },
  ///     ...
  ///   ]
  /// }
  /// \endcode
  ///
  ECIContainer read_eci(const fs::path &filepath) {

    std::vector<double> value;
    std::vector<ECIContainer::size_type> index;

    jsonParser json(filepath);
    for(auto it = json["cluster_functions"].begin(); it != json["cluster_functions"].end(); ++it) {
      auto eci = it->find("eci");
      if(eci != it->end()) {
        value.push_back(eci->get<double>());
        index.push_back(it->find("linear_function_index")->get<ECIContainer::size_type>());
      }
    }
    return ECIContainer(value.begin(), value.end(), index.begin());
  }

}
