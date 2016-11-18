#include "casm/casm_io/jsonParser.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

  /// Functions for converting basic types to/from json

  jsonParser &to_json(bool value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(value);
    return json;
  }

  jsonParser &to_json(int value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(value);
    return json;
  }

  jsonParser &to_json(unsigned int value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(boost::uint64_t(value));
    return json;
  }

  jsonParser &to_json(const long int value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(boost::int64_t(value));
    return json;
  }

  jsonParser &to_json(const unsigned long int value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(boost::uint64_t(value));
    return json;
  }

  jsonParser &to_json(double value, jsonParser &json) {
    if(value != value) {
      return to_json("nan", json);
    }
    else if(value == 1.0 / 0.0) {
      return to_json("inf", json);
    }
    else if(value == -1.0 / 0.0) {
      return to_json("-inf", json);
    }
    else {
      *((json_spirit::mValue *) &json) = json_spirit::mValue(value);
      return json;
    }

  }

  jsonParser &to_json(const std::string &value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(value);
    return json;
  }

  jsonParser &to_json(const char *value, jsonParser &json) {
    *((json_spirit::mValue *) &json) = json_spirit::mValue(value);
    return json;
  }

  jsonParser &to_json(const jsonParser &value, jsonParser &json) {
    return json = value;
  }


  template<>
  bool from_json<bool>(const jsonParser &json) {
    return json.get_bool();
  }

  template<>
  int from_json<int>(const jsonParser &json) {
    return json.get_int();
  }

  template<>
  unsigned int from_json<unsigned int>(const jsonParser &json) {
    return (unsigned int) json.get_int();
  }

  template<>
  long int from_json<long int>(const jsonParser &json) {
    return (long int) json.get_int64();
  }

  template<>
  unsigned long int from_json<unsigned long int>(const jsonParser &json) {
    return (unsigned long int) json.get_uint64();
  }

  template<>
  double from_json<double>(const jsonParser &json) {
    double d;
    from_json(d, json);
    return d;
  }

  template<>
  std::string from_json<std::string>(const jsonParser &json) {
    return json.get_str();
  }

  template<>
  jsonParser from_json<jsonParser>(const jsonParser &json) {
    return json;
  }


  void from_json(bool &value, const jsonParser &json) {
    value = json.get_bool();
  }

  void from_json(int &value, const jsonParser &json) {
    value = json.get_int();
  }

  void from_json(unsigned int &value, const jsonParser &json) {
    value = json.get_int();
  }

  void from_json(long int &value, const jsonParser &json) {
    value = json.get_int64();
  }

  void from_json(unsigned long int &value, const jsonParser &json) {
    value = json.get_uint64();
  }

  void from_json(double &value, const jsonParser &json) {
    if(json.is_string()) {
      std::string str = json.get_str();
      if(str == "nan") {
        value = sqrt(-1.0);
      }
      else if(str == "inf") {
        value = 1.0 / 0.0;
      }
      else if(str == "-inf") {
        value = -1.0 / 0.0;
      }
      else {
        throw std::runtime_error("Expected json real, received string other than 'nan', 'inf', or '-inf': '" + str + "'");
      }
    }
    else {
      value = json.get_real();
    }
  }

  void from_json(std::string &value, const jsonParser &json) {
    value = json.get_str();
  }

  void from_json(jsonParser &value, const jsonParser &json) {
    value = json;
  }

  void from_json(fs::path &value, const jsonParser &json) {
    value = fs::path(json.get_str());
  }



  // ---- Read/Print JSON  ----------------------------------

  bool jsonParser::read(std::istream &stream) {
    return json_spirit::read_stream(stream, (json_spirit::mValue &) * this);
  }

  bool jsonParser::read(const boost::filesystem::path &file_path) {
    boost::filesystem::ifstream stream(file_path);
    return read(stream);
  }

  std::istream &operator>>(std::istream &stream, jsonParser &json) {
    if(!json.read(stream)) {
      std::cerr << "ERROR: Unable to successfully parse JSON file.  File parsed as:\n"
                << json
                << "\nPlease correct input file and try again. Exiting...\n";
      exit(1);
    }
    return stream;
  }

  /// Writes json to stream
  void jsonParser::print(std::ostream &stream, unsigned int indent, unsigned int prec) const {
    json_spirit::write_stream((json_spirit::mValue &) *this, stream, indent, prec,
                              json_spirit::pretty_print | json_spirit::single_line_arrays);
  };

  /// Write json to file
  void jsonParser::write(const std::string &file_name, unsigned int indent, unsigned int prec) const {
    std::ofstream file(file_name.c_str());
    print(file, indent, prec);
    file.close();
    return;
  }

  /// Write json to file
  void jsonParser::write(const boost::filesystem::path &file_path, unsigned int indent, unsigned int prec) const {
    boost::filesystem::ofstream file(file_path);
    print(file, indent, prec);
    file.close();
    return;
  }

  std::ostream &operator<< (std::ostream &stream, const jsonParser &json) {
    json.print(stream);
    return stream;
  }

  bool jsonParser::almost_equal(const jsonParser &B, double tol) const {
    if(type() != B.type()) {
      return false;
    }

    if(is_array()) {
      auto f = [ = ](const jsonParser & _A, const jsonParser & _B) {
        return _A.almost_equal(_B, tol);
      };
      bool res = (size() == B.size() && std::equal(begin(), end(), B.begin(), f));
      return res;
    }
    else if(is_obj()) {
      if(size() != B.size()) {
        return false;
      }
      auto A_it = begin();
      auto A_end = end();
      auto B_it = B.begin();
      for(; A_it != A_end; ++A_it, ++B_it) {
        if(A_it.name() != B_it.name() || !A_it->almost_equal(*B_it, tol)) {
          return false;
        }
      }
      return true;
    }
    else if(is_number()) {
      bool res = CASM::almost_equal(this->get<double>(), B.get<double>(), tol);
      return res;
    }
    else {
      bool res = (*this == B);
      return res;
    }
  }


  // ------ Type Checking Methods ------------------------------------

  /// Check if null type
  bool jsonParser::is_null() const {
    return type() == json_spirit::null_type;
  }

  /// Check if bool type
  bool jsonParser::is_bool() const {
    return type() == json_spirit::bool_type;
  }

  /// Check if int type
  bool jsonParser::is_int() const {
    return type() == json_spirit::int_type;
  }

  /// Check if number type (not including int)
  bool jsonParser::is_number() const {
    return type() == json_spirit::real_type;
  }

  /// Check if string
  bool jsonParser::is_string() const {
    return type() == json_spirit::str_type;
  }

  /// Check if object type
  bool jsonParser::is_obj() const {
    return type() == json_spirit::obj_type;
  }

  /// Check if array type
  bool jsonParser::is_array() const {
    return type() == json_spirit::array_type;
  }


  // ---- Navigate the JSON data: ----------------------------

  /// Return a reference to the sub-jsonParser (JSON value) with 'name' if it exists
  ///   If it does not exist, create it with value == 'null' and return a reference
  jsonParser &jsonParser::operator[](const std::string &name) {

    json_spirit::mObject &obj = get_obj();
    json_spirit::mObject::iterator it = obj.find(name);

    // if 'name' not found, add it and with value 'null'
    if(it == obj.end()) {
      obj[name] = json_spirit::mValue(json_spirit::mObject());
    }
    return (jsonParser &) obj[name];
  }

  /// Return a reference to the sub-jsonParser (JSON value) with 'name' if it exists.
  ///   Will throw if the 'name' doesn't exist.
  const jsonParser &jsonParser::operator[](const std::string &name) const {

    const json_spirit::mObject &obj = get_obj();
    json_spirit::mObject::const_iterator it = obj.find(name);

    // if 'name' not found, add it and with value 'null'
    if(it == obj.end()) {
      throw std::runtime_error("Const operator[] access, but " + name + " does not exist");
    }
    return (const jsonParser &) it->second;
  }

  /// Return a reference to the sub-jsonParser (JSON value) with specified relative path
  ///   Will throw if the 'path' doesn't exist.
  ///
  /// - If 'path' is 'A/B/C', then json[path] is equivalent to json[A][B][C]
  /// - If any sub-jsonParser is an array, it will attempt to convert the filename to int
  jsonParser &jsonParser::at(const fs::path &path) {
    return const_cast<jsonParser &>(static_cast<const jsonParser *>(this)->at(path));
  }

  /// Return a reference to the sub-jsonParser (JSON value) with specified relative path
  ///   Will throw if the 'path' doesn't exist.
  ///
  /// - If 'path' is 'A/B/C', then json[path] is equivalent to json[A][B][C]
  /// - If any sub-jsonParser is an array, it will attempt to convert the filename to int
  const jsonParser &jsonParser::at(const fs::path &path) const {
    if(!path.is_relative()) {
      throw std::invalid_argument(
        "Error in jsonParser::operator[](const fs::path &path): path must be relative");
    }
    const jsonParser *curr = this;
    for(auto it = path.begin(); it != path.end(); ++it) {
      if(curr->is_array()) {
        int index = std::stoi(it->string());
        if(curr->size() > index) {
          curr = &((*curr)[index]);
        }
        else {
          std::string msg = "Error in jsonParser::at: attempted to access element outside of array range";
          std::cerr << "path: " << path << std::endl;
          std::cerr << "index: " << index << std::endl;
          std::cerr << "curr->size(): " << curr->size() << std::endl;
          throw std::invalid_argument(msg);
        }
      }
      else {
        auto res = curr->find(it->string());
        if(res != curr->end()) {
          curr = &((*curr)[it->string()]);
        }
        else {
          std::string msg = "Error in jsonParser::at: key not found";
          std::cerr << "path: " << path << std::endl;
          std::cerr << "key: " << it->string() << std::endl;
          throw std::invalid_argument(msg);
        }
      }
    }

    return *curr;
  }

  /// Return a reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
  jsonParser &jsonParser::operator[](const int &element) {

    return (jsonParser &) get_array()[element];
  }

  /// Return a const reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
  const jsonParser &jsonParser::operator[](const int &element) const {

    return (const jsonParser &) get_array()[element];
  }

  /// Return the location at which jsonParser 'A' != 'B' as a boost::filesystem::path
  boost::filesystem::path find_diff(const jsonParser &A, const jsonParser &B, boost::filesystem::path diff) {
    auto A_it = A.cbegin();
    auto B_it = B.cbegin();
    while(A_it != A.cend()) {
      if(*A_it != *B_it) {
        if(A.is_obj() && B.is_obj()) {
          return find_diff(*A_it, *B_it, diff / A_it.name());
        }
        else if(A.is_array() && B.is_array()) {
          std::stringstream ss;
          ss << std::distance(A.cbegin(), A_it);
          return find_diff(*A_it, *B_it, diff / ss.str());
        }
        return diff;
      }
      ++A_it;
      ++B_it;
    }
    return diff;
  }

  /// Return the location at which jsonParser !A.almost_equal(B, tol) as a boost::filesystem::path
  boost::filesystem::path find_diff(const jsonParser &A, const jsonParser &B, double tol, boost::filesystem::path diff) {
    auto A_it = A.cbegin();
    auto B_it = B.cbegin();
    while(A_it != A.cend()) {
      if(!A_it->almost_equal(*B_it, tol)) {
        if(A.is_obj() && B.is_obj()) {
          if(A.size() != B.size()) {
            return diff;
          }
          return find_diff(*A_it, *B_it, tol, diff / A_it.name());
        }
        else if(A.is_array() && B.is_array()) {
          if(A.size() != B.size()) {
            return diff;
          }
          std::stringstream ss;
          ss << std::distance(A.cbegin(), A_it);
          return find_diff(*A_it, *B_it, tol, diff / ss.str());
        }
        return diff;
      }
      ++A_it;
      ++B_it;
    }
    return diff;
  }

  /// Returns array size if *this is a JSON array, object size if *this is a JSON object, 1 otherwise
  jsonParser::size_type jsonParser::size() const {
    if(is_obj())
      return get_obj().size();
    else if(is_array())
      return get_array().size();
    else
      return 1;
  }

  /// Returns iterator to beginning of JSON object or JSON array
  jsonParser::iterator jsonParser::begin() {
    if(is_obj())
      return iterator(this, get_obj().begin());
    else if(is_array())
      return iterator(this, get_array().begin());
    else
      return iterator(this, 0);
  }

  /// Returns const_iterator to beginning of JSON object or JSON array
  jsonParser::const_iterator jsonParser::begin() const {
    return cbegin();
  }

  /// Returns const iterator to beginning of const JSON object or JSON array
  jsonParser::const_iterator jsonParser::cbegin() const {
    if(is_obj())
      return const_iterator(this, get_obj().cbegin());
    else if(is_array())
      return const_iterator(this, get_array().cbegin());
    else
      return const_iterator(this, 0);
  }

  /// Returns iterator to end of JSON object or JSON array
  jsonParser::iterator jsonParser::end() {
    if(is_obj())
      return iterator(this, get_obj().end());
    else if(is_array())
      return iterator(this, get_array().end());
    else
      return iterator(this, 0);
  }

  /// Returns iterator to end of JSON object or JSON array
  jsonParser::const_iterator jsonParser::end() const {
    return cend();
  }

  /// Returns const_iterator to end of JSON object or JSON array
  jsonParser::const_iterator jsonParser::cend() const {
    if(is_obj())
      return const_iterator(this, get_obj().cend());
    else if(is_array())
      return const_iterator(this, get_array().cend());
    else
      return const_iterator(this, 0);
  }

  /// Return iterator to JSON object value with 'name'
  jsonParser::iterator jsonParser::find(const std::string &name) {
    return iterator(this, get_obj().find(name));
  }

  /// Return const_iterator to JSON object value with 'name'
  jsonParser::const_iterator jsonParser::find(const std::string &name) const {
    return const_iterator(this, get_obj().find(name));
  }

  /// Return true if JSON object contains 'name'
  bool jsonParser::contains(const std::string &name) const {
    return find(name) != cend();
  }

  /// Erase key:value pair from an object
  ///   Returns the number of elements erased, which will be 0 or 1
  json_spirit::mObject::size_type jsonParser::erase(const std::string &name) {
    return get_obj().erase(name);
  }

}

