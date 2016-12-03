#ifndef CASM_jsonIO_container
#define CASM_jsonIO_container

#include "casm/casm_io/jsonParser.hh"
#include "casm/CASM_global_definitions.hh"

//#include "casm/container/Array.hh"
//#include "casm/container/LinearAlgebra.hh"

namespace CASM {
  /*
    // --- Array<T> ------------------------------

    /// Converts to a JSON array
    template<typename T>
    jsonParser &to_json(const Array<T> &value, jsonParser &json) {
      return json.put_array(value.begin(), value.end());
    }

    /// This requires that 'T::T()' exists, if not, you must do this by hand
    template<typename T>
    void from_json(Array<T> &value, const jsonParser &json) {
      value.resize(json.size());
      for(int i = 0; i < json.size(); i++)
        from_json(value[i], json[i]);
    }

    // --- Matrix & Vector -----------------------

    /// Converts to a JSON array of arrays
    template <typename T, Index N>
    jsonParser &to_json(const Matrix<T, N> &value, jsonParser &json) {
      json.put_array();
      for(Index i = 0; i < value.num_rows(); i++) {
        jsonParser json_row;
        json_row.put_array();
        for(Index j = 0; j < value.num_cols(); j++) {
          json_row.push_back(value(i, j));
        }
        json.push_back(json_row);
      }
      return json;
    }

    template <typename T, Index N>
    void from_json(Matrix<T, N> &value, const jsonParser &json) {
      for(Index i = 0; i < value.num_rows(); i++) {
        for(Index j = 0; j < value.num_cols(); j++) {
          from_json(value.at(i, j), json[i][j]);
        }
      }
    }

    /// Converts to a JSON array
    template <typename T, Index N>
    jsonParser &to_json(const Vector<T, N> &value, jsonParser &json) {
      json.put_array();
      for(Index i = 0; i < value.size(); i++)
        json.push_back(value[i]);
      return json;
    }

    template <typename T, Index N>
    void from_json(Vector<T, N> &value, const jsonParser &json) {
      for(Index i = 0; i < value.size(); i++)
        from_json(value.at(i), json[i]);
    }


    // --- Permutation ---------------------------

    class Permutation;

    jsonParser &to_json(const Permutation &perm, jsonParser &json);

    template<>
    Permutation from_json(const jsonParser &json);

    void from_json(Permutation &perm, const jsonParser &json);
  */

  /**  \ingroup jsonParser
   *
   *  @{
  */

  // --- std::map<std::string, T> --------------

  /// Converts to a JSON object
  template<typename T>
  jsonParser &to_json(const std::map<std::string, T> &map, jsonParser &json) {
    return json.put_obj(map.begin(), map.end());
  }

  /// Read map from JSON
  ///
  /// Clears any previous contents
  template<typename T>
  void from_json(std::map<std::string, T> &map, const jsonParser &json) {
    map.clear();
    for(auto it = json.begin(); it != json.end(); ++it) {
      from_json(map[it.name()], *it);
    }
  }

  // --- std::vector<T> --------------

  /// Converts to a JSON array
  template<typename T>
  jsonParser &to_json(const std::vector<T> &vec, jsonParser &json) {
    return json.put_array(vec.begin(), vec.end());
  }

  /// Read map from JSON
  ///
  /// Clears any previous contents
  template<typename T>
  void from_json(std::vector<T> &vec, const jsonParser &json) {
    vec.resize(json.size());
    int i = 0;
    for(auto it = json.begin(); it != json.end(); ++it, ++i) {
      from_json(vec[i], *it);
    }
  }

  // --- std::array<T,std::size_t N> --------------

  /// Converts to a JSON array
  template<typename T, std::size_t N>
  jsonParser &to_json(const std::array<T, N> &arr, jsonParser &json) {
    return json.put_array(arr.begin(), arr.end());
  }

  /// Read map from JSON
  ///
  /// Clears any previous contents
  template<typename T, std::size_t N>
  void from_json(std::array<T, N> &arr, const jsonParser &json) {
    if(json.size() != N) {
      std::ostringstream ss;
      ss << "Attempting to to initialize std::array of size " << N << " from JSON array of size " << json.size() << ":\n    " << json << "\n\n";
      throw std::range_error(ss.str());
    }

    std::size_t i = 0;
    for(auto it = json.begin(); it != json.end(); ++it, ++i) {
      from_json(arr[i], *it);
    }
  }

  // --- std::set<T> --------------

  /// Converts to a JSON array
  template<typename T, typename Compare>
  jsonParser &to_json(const std::set<T, Compare> &set, jsonParser &json) {
    return json.put_array(set.begin(), set.end());
  }

  /// Read std::set from JSON array
  ///
  /// Clears any previous contents, constructs via jsonConstructor<T>
  template<typename T, typename Compare, typename...Args>
  void from_json(std::set<T, Compare> &set, const jsonParser &json, Args...args) {
    set.clear();
    for(auto it = json.begin(); it != json.end(); ++it) {
      set.insert(jsonConstructor<T>::from_json(*it, args...));
    }
  }
  //}
  //
  //namespace Eigen {

  /// \brief Write Eigen Matrix/Vector to JSON
  template <typename Derived>
  CASM::jsonParser &to_json(const Eigen::MatrixBase<Derived> &value, CASM::jsonParser &json) {
    json.put_array();
    for(int i = 0; i < value.rows(); i++) {
      CASM::jsonParser json_row;
      json_row.put_array();
      for(int j = 0; j < value.cols(); j++) {
        json_row.push_back(value(i, j));
      }
      json.push_back(json_row);
    }
    return json;
  }


  /// \brief Write Eigen Matrix with 1 row or 1 column to JSON array
  template <typename Derived>
  CASM::jsonParser &to_json_array(const Eigen::MatrixBase<Derived> &value, CASM::jsonParser &json) {
    json.put_array();
    if(value.rows() == 1) {
      for(int i = 0; i < value.cols(); ++i) {
        json.push_back(value(0, i));
      }
    }
    else if(value.cols() == 1) {
      for(int i = 0; i < value.rows(); ++i) {
        json.push_back(value(i, 0));
      }
    }
    else {
      throw std::runtime_error("Error in 'to_json_array': Not a vector");
    }
    return json;
  }

  /// \brief Read Eigen Matrix/Vector from JSON
  ///
  /// - Reads 1d JSON array as column vector
  template <typename Derived>
  void from_json(Eigen::MatrixBase<Derived>  &value, const CASM::jsonParser &json) {
    if(json.is_array() && !json[0].is_array()) {
      value.derived().resize(json.size(), 1);
      for(int i = 0; i < value.rows(); i++) {
        from_json(value(i, 0), json[i]);
      }
    }
    else {
      value.derived().resize(json.size(), json[0].size());
      for(int i = 0; i < value.rows(); i++) {
        for(int j = 0; j < value.cols(); j++) {
          from_json(value(i, j), json[i][j]);
        }
      }
    }
  }

  /** @} */
}

#endif
