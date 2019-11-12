#include "casm/casm_io/json/jsonParser.hh"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
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

  /// Create a jsonParser by reading a file
  ///
  /// This function reads the contents of the file at 'file_path' as if it were JSON.
  /// Use 'to_json(file_path.string(), json)' if you only want the path as a string
  void to_json(fs::path file_path, jsonParser &json) {
    if(!json.read(file_path)) {
      throw std::runtime_error(
        std::string("ERROR: Could not read JSON file: '") + file_path.string() +
        "'.\n\nPlease check your formatting. For instance, try http://www.jsoneditoronline.org.");
    }
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

  bool jsonParser::read(const fs::path &file_path) {
    fs::ifstream stream(file_path);
    return read(stream);
  }

  std::istream &operator>>(std::istream &stream, jsonParser &json) {
    if(!json.read(stream)) {
      std::stringstream msg;
      msg << "Error: Unable to successfully parse JSON file.  File parsed as:\n"
          << json
          << "\nPlease correct input file and try again. Exiting...\n";
      throw std::invalid_argument(msg.str());
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
  void jsonParser::write(const fs::path &file_path, unsigned int indent, unsigned int prec) const {
    fs::ofstream file(file_path);
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
      throw std::runtime_error("JSON const operator[] access, but " + name + " does not exist");
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
          std::stringstream msg;
          msg << "Error in jsonParser::at: attempted to access element outside of array range. "
              << "path: '" << path << "' "
              << "index: " << index << " "
              << "curr->size(): " << curr->size();
          throw std::invalid_argument(msg.str());
        }
      }
      else {
        auto res = curr->find(it->string());
        if(res != curr->end()) {
          curr = &((*curr)[it->string()]);
        }
        else {
          std::stringstream msg;
          msg << "Error in jsonParser::at: key '" << it->string() << "' not found at '"
              << path << "'.";
          throw std::invalid_argument(msg.str());
        }
      }
    }

    return *curr;
  }

  /// Return a reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
  jsonParser &jsonParser::operator[](const size_type &element) {

    return (jsonParser &) get_array()[element];
  }

  /// Return a const reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
  const jsonParser &jsonParser::operator[](const size_type &element) const {

    return (const jsonParser &) get_array()[element];
  }

  /// Return a reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
  jsonParser &jsonParser::at(const size_type &element) {
    if(!is_array()) {
      throw std::invalid_argument("Error in jsonParser::at: attempting to access non-array with index");
    }
    if(!(element < size())) {
      throw std::out_of_range("Error in jsonParser::at: out of range");
    }
    return (jsonParser &) get_array()[element];
  }

  /// Return a const reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
  const jsonParser &jsonParser::at(const size_type &element) const {
    if(!is_array()) {
      throw std::invalid_argument("Error in jsonParser::at: attempting to access non-array with index");
    }
    if(!(element < size())) {
      throw std::out_of_range("Error in jsonParser::at: out of range");
    }
    return (const jsonParser &) get_array()[element];
  }

  namespace {
    /// Return the location at which jsonParser 'A' != 'B' as a fs::path
    fs::path find_diff(const jsonParser &A, const jsonParser &B, fs::path diff) {
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

    /// Return the location at which jsonParser !A.almost_equal(B, tol) as a fs::path
    fs::path find_diff(const jsonParser &A, const jsonParser &B, double tol, fs::path diff) {
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

  }

  /// Return the location at which jsonParser 'A' != 'B' as a fs::path
  fs::path find_diff(const jsonParser &A, const jsonParser &B) {
    return find_diff(A, B, fs::path());
  }

  /// Return the location at which jsonParser !A.almost_equal(B, tol) as a fs::path
  fs::path find_diff(const jsonParser &A, const jsonParser &B, double tol) {
    return find_diff(A, B, tol, fs::path());
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
      return iterator(this, 1);
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
      return const_iterator(this, 1);
  }

  /// Return iterator to JSON object value with 'name'
  jsonParser::iterator jsonParser::find(const std::string &name) {
    return iterator(this, get_obj().find(name));
  }

  /// Return const_iterator to JSON object value with 'name'
  jsonParser::const_iterator jsonParser::find(const std::string &name) const {
    return const_iterator(this, get_obj().find(name));
  }

  /// Return iterator to sub-object or element, or 'end' if not found
  ///
  /// - If path.empty(), return iterator that dereferences to this, and one increment results in end
  jsonParser::iterator jsonParser::find_at(const fs::path &path) {
    if(!path.is_relative()) {
      throw std::invalid_argument(
        "Error in jsonParser::operator[](const fs::path &path): path must be relative");
    }
    if(path.empty()) {
      return jsonParser::iterator(this, 0);
    }
    jsonParser *curr = this;
    jsonParser::iterator res = this->end();
    for(auto it = path.begin(); it != path.end(); ++it) {
      if(curr->is_array()) {
        int index = std::stoi(it->string());
        if(curr->size() > index) {
          res = curr->begin();
          for(int i = 0; i < index; ++i) {
            ++res;
          }
          curr = &(*res);
        }
        else {
          return this->end();
        }
      }
      else {
        res = curr->find(it->string());
        if(res != curr->end()) {
          curr = &((*curr)[it->string()]);
        }
        else {
          return this->end();
        }
      }
    }

    return res;
  }

  /// Return iterator to sub-object or element, or 'end' if not found
  ///
  /// - If path.empty(), return iterator that dereferences to this, and one increment results in end
  jsonParser::const_iterator jsonParser::find_at(const fs::path &path) const {
    return const_cast<jsonParser *>(this)->find_at(path);
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

  // ---- static Methods -------------------------------------

  /// Construct a jsonParser from a file containing JSON data
  jsonParser jsonParser::parse(const fs::path &path) {
    return jsonParser(path);
  }


  template<bool IsConst>
  jsonParserIterator<IsConst>::jsonParserIterator() {}

  template<bool IsConst>
  jsonParserIterator<IsConst>::jsonParserIterator(const jsonParserIterator &iter)
    : parser(iter.parser), type(iter.type), obj_iter(iter.obj_iter), array_iter(iter.array_iter), val_iter(iter.val_iter) {
  }

  template<bool IsConst>
  jsonParserIterator<IsConst> &jsonParserIterator<IsConst>::operator=(jsonParserIterator iter) {
    swap(*this, iter);
    return *this;
  }

  template<bool IsConst>
  jsonParserIterator<IsConst>::jsonParserIterator(
    typename jsonParserIterator<IsConst>::pointer j,
    const typename jsonParserIterator<IsConst>::object_iterator &iter)
    : parser(j), type(json_spirit::obj_type), obj_iter(iter) {
  }

  template<bool IsConst>
  jsonParserIterator<IsConst>::jsonParserIterator(
    typename jsonParserIterator<IsConst>::pointer j,
    const typename jsonParserIterator<IsConst>::array_iterator &iter)
    : parser(j), type(json_spirit::array_type), array_iter(iter) {
  }

  template<bool IsConst>
  jsonParserIterator<IsConst>::jsonParserIterator(
    typename jsonParserIterator<IsConst>::pointer j,
    const int &iter)
    : parser(j), type(json_spirit::null_type), val_iter(iter) {
  }

  template<bool IsConst>
  typename jsonParserIterator<IsConst>::reference jsonParserIterator<IsConst>::operator*() const {
    if(type == json_spirit::obj_type)
      return (reference) obj_iter->second;
    else if(type == json_spirit::array_type)
      return (reference) * array_iter;
    else
      return *parser;
  }

  template<bool IsConst>
  typename jsonParserIterator<IsConst>::pointer jsonParserIterator<IsConst>::operator->() const {
    if(type == json_spirit::obj_type)
      return (pointer) &obj_iter->second;
    else if(type == json_spirit::array_type)
      return (pointer) & (*array_iter);
    else
      return parser;
  }

  template<bool IsConst>
  bool jsonParserIterator<IsConst>::operator==(const jsonParserIterator &iter) const {
    if(parser != iter.parser) {
      return false;
    }

    bool this_is_end = this->is_end();
    bool that_is_end = iter.is_end();

    if(this_is_end && that_is_end) {
      return true;
    }
    else if(this_is_end != that_is_end) {
      return false;
    }
    else {
      if(type == json_spirit::obj_type) {
        return obj_iter == iter.obj_iter;
      }
      else if(type == json_spirit::array_type) {
        return array_iter == iter.array_iter;
      }
      else if(type == json_spirit::null_type) {
        return val_iter == iter.val_iter;
      }
    }

    return false;
  }

  template<bool IsConst>
  bool jsonParserIterator<IsConst>::is_end() const {
    if(type == json_spirit::obj_type && obj_iter == parser->get_obj().end()) {
      return true;
    }
    else if(type == json_spirit::array_type && array_iter == parser->get_array().end()) {
      return true;
    }
    else if(type == json_spirit::null_type && val_iter == 1) {
      return true;
    }
    return false;
  }

  template<bool IsConst>
  bool jsonParserIterator<IsConst>::operator!=(const jsonParserIterator &iter) const {
    return !(*this == iter);
  }

  template<bool IsConst>
  jsonParserIterator<IsConst> &jsonParserIterator<IsConst>::operator++() {
    if(type == json_spirit::obj_type) {
      ++obj_iter;
      return *this;
    }
    else if(type == json_spirit::array_type) {
      ++array_iter;
      return *this;
    }
    else {
      ++val_iter;
      return *this;
    }
  }

  template<bool IsConst>
  jsonParserIterator<IsConst> jsonParserIterator<IsConst>::operator++(int) {

    jsonParserIterator cp(*this);

    if(type == json_spirit::obj_type) {
      ++obj_iter;
      return cp;
    }
    else if(type == json_spirit::array_type) {
      ++array_iter;
      return cp;
    }
    else {
      ++val_iter;
      return cp;
    }
  }

  template<bool IsConst>
  jsonParserIterator<IsConst> &jsonParserIterator<IsConst>::operator--() {
    if(type == json_spirit::obj_type) {
      --obj_iter;
      return *this;
    }
    else if(type == json_spirit::array_type) {
      --array_iter;
      return *this;
    }
    else {
      --val_iter;
      return *this;
    }
  }

  template<bool IsConst>
  jsonParserIterator<IsConst> jsonParserIterator<IsConst>::operator--(int) {

    jsonParserIterator<IsConst> cp(*this);

    if(type == json_spirit::obj_type) {
      --obj_iter;
      return cp;
    }
    else if(type == json_spirit::array_type) {
      --array_iter;
      return cp;
    }
    else {
      --val_iter;
      return cp;
    }
  }

  template<bool IsConst>
  jsonParserIterator<IsConst>::operator jsonParser::const_iterator() const {
    if(type == json_spirit::obj_type)
      return jsonParser::const_iterator(parser, obj_iter);
    else if(type == json_spirit::array_type)
      return jsonParser::const_iterator(parser, array_iter);
    else
      return jsonParser::const_iterator(parser, val_iter);
  }

  /// When iterating over a JSON object, returns the 'name' of the 'name':value pair the iterator is pointing at
  template<bool IsConst>
  std::string jsonParserIterator<IsConst>::name() const {
    if(type == json_spirit::obj_type)
      return obj_iter->first;
    else
      throw std::runtime_error("Calling 'name' on non-object jsonParserIterator");
  }

  template<bool IsConst>
  void swap(jsonParserIterator<IsConst> &a, jsonParserIterator<IsConst> &b) {
    using std::swap;

    std::swap(a.parser, b.parser);
    swap(a.type, b.type);
    swap(a.obj_iter, b.obj_iter);
    swap(a.array_iter, b.array_iter);
    swap(a.val_iter, b.val_iter);
  }

  template class jsonParserIterator<true>;
  template class jsonParserIterator<false>;
  template void swap<true>(jsonParserIterator<true> &, jsonParserIterator<true> &);
  template void swap<false>(jsonParserIterator<false> &, jsonParserIterator<false> &);

}
