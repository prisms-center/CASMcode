#ifndef CASM_JSONPARSER_HH
#define CASM_JSONPARSER_HH

#include "casm/external/json_spirit/json_spirit_reader_template.h"
#include "casm/external/json_spirit/json_spirit_writer_template.h"

#include "casm/misc/CASM_TMP.hh"

#include <exception>
#include <complex>

#include "casm/external/boost.hh"

namespace CASM {

  template<bool IsConst>
  class jsonParserIterator;

  /**
   * \defgroup jsonParser
   *
   * \brief JSON input/output
   *
   * \ingroup casmIO
   *
   * @{
   */

  /// jsonParser allows for reading / writing JSON data
  ///
  /// JSON consists of values, which can be of type:
  ///     object, array, float, int, str, null, bool
  ///   object: map/dict of <string, JSON value>
  ///   array: vector/list of JSON values
  ///
  ///   value[string] -> value iff object
  ///   value[int] -> value iff array
  ///   value.get<type>() -> type must be correct type
  ///
  /// Assumes functions exist for each type T to read/write:
  ///   jsonParser& to_json( const T &value, jsonParser &json)
  ///   void from_json( T &value, const jsonParser &json);
  ///
  ///   These functions exist for basic types (bool, int, double, std::string),
  ///   and the containers:
  ///      Eigen::MatrixXd, CASM::Array, CASM::Matrix, CASM::Vector,
  ///   So if they also exist for the type being contained,
  ///   you can read/write the entire object as a JSON array or nested JSON arrays
  ///
  /// Simple Usage:
  ///
  ///   Reading data from a JSON file:
  ///
  ///     jsonParser json("myfile.json");
  ///
  ///   Three ways to get data of type T:
  ///
  ///     T t;
  ///     from_json(t, json["path"]["to"]["object"]["data"]);  <-- use [std::string] for JSON values of any type
  ///     t = json["path"]["to"]["object"]["data"].get<T>();
  ///     t = from_json<T>(json["path"]["to"]["array"][0]);   <-- use [int] for JSON arrays only
  ///
  ///   Check if 'name' in JSON object:
  ///
  ///     if( json.contains("other_data") )
  ///       from_json( array, json["other_data"]);
  ///
  ///   Writing data to a JSON file:
  ///
  ///     jsonParser json;
  ///     json["mydata"].put(mydata);        <-- if ["mydata"] does not exist, it is created, else overwrites
  ///     json["more_data"].put(more_data);
  ///
  ///     ofstream file("myfile.json");
  ///     json.write(file);
  ///     file.close();
  ///
  class jsonParser : public json_spirit::mValue {

  public:

    typedef json_spirit::mObject::size_type size_type;
    typedef jsonParserIterator<false> iterator;
    typedef jsonParserIterator<true> const_iterator;

    // ---- Constructors  ----------------------------------

    /// Create a new empty jsonParser
    jsonParser() :
      json_spirit::mValue(json_spirit::mObject()) {}

    /// Create a jsonParser from any other object for which 'to_json(t, json)' is defined
    ///
    /// To parse a std::string as JSON, rather than store it verbatim, use jsonParser::parse
    template<typename T>
    explicit jsonParser(T &t) {
      to_json(t, *this);
    }

    /// Create a jsonParser from any other object for which 'to_json(t, json)' is defined
    ///
    /// To parse a std::string as JSON, rather than store it verbatim, use jsonParser::parse
    template<typename T>
    explicit jsonParser(const T &t) {
      to_json(t, *this);
    }


    // ---- Read/Print JSON  ----------------------------------

    /// Reads json from the stream
    bool read(std::istream &stream);

    /// Reads json from a path
    bool read(const boost::filesystem::path &mypath);


    /// Print json to stream
    void print(std::ostream &stream, unsigned int indent = 2, unsigned int prec = 12) const;

    /// Write json to file
    void write(const std::string &file_name, unsigned int indent = 2, unsigned int prec = 12) const;

    /// Write json to file
    void write(const boost::filesystem::path &mypath, unsigned int indent = 2, unsigned int prec = 12) const;


    // ---- Value level printing options: ---------------------

    /// Force array printing as column
    using json_spirit::mValue::set_force_column;

    /// Force array printing as row
    using json_spirit::mValue::set_force_row;

    /// Force printing double using scientific notation
    using json_spirit::mValue::set_scientific;

    /// Remove trailing zeros for real (double) values
    using json_spirit::mValue::set_remove_trailing_zeros;

    /// Do not force array printing as column
    using json_spirit::mValue::unset_force_column;

    /// Do not force array printing as row
    using json_spirit::mValue::unset_force_row;

    /// Do not force printing double using scientific notation
    using json_spirit::mValue::unset_scientific;

    /// Do not remove trailing zeros for real (double) values
    using json_spirit::mValue::unset_remove_trailing_zeros;


    using json_spirit::mValue::operator==;

    bool operator!=(const jsonParser &json) const {
      return !(json_spirit::mValue::operator==(json));
    }

    bool almost_equal(const jsonParser &B, double tol) const;


    // ------ Type Checking Methods ------------------------------------

    bool is_null() const;
    bool is_bool() const;
    bool is_int() const;
    bool is_number() const;
    bool is_string() const;
    bool is_obj() const;
    bool is_array() const;


    // ---- Navigate the JSON data: ----------------------------

    /// Return a reference to the sub-jsonParser (JSON value) with 'name' if it exists
    ///   If it does not exist, create it with an empty JSON object and return a reference to it
    jsonParser &operator[](const std::string &name);

    /// Return a reference to the sub-jsonParser (JSON value) with 'name' if it exists.
    ///   Will throw if the 'name' doesn't exist.
    const jsonParser &operator[](const std::string &name) const;

    /// Return a reference to the sub-jsonParser (JSON value) with specified relative path
    ///   Will throw if the 'path' doesn't exist.
    ///
    /// - If 'path' is 'A/B/C', then json.at(path) is equivalent to json[A][B][C]
    /// - If any sub-jsonParser is an array, it will attempt to convert the filename to int
    jsonParser &at(const fs::path &path);

    /// Return a reference to the sub-jsonParser (JSON value) with specified relative path
    ///   Will throw if the 'path' doesn't exist.
    ///
    /// - If 'path' is 'A/B/C', then json.at(path) is equivalent to json[A][B][C]
    /// - If any sub-jsonParser is an array, it will attempt to convert the filename to int
    const jsonParser &at(const fs::path &path) const;

    /// Return a reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
    jsonParser &operator[](const int &element);

    /// Return a const reference to the sub-jsonParser (JSON value) from index 'element' iff jsonParser is a JSON array
    const jsonParser &operator[](const int &element) const;

    /// Returns array size if *this is a JSON array, object size if *this is a JSON object, 1 otherwise
    size_type size() const;

    /// Returns const_iterator to beginning of JSON object or JSON array
    iterator begin();

    /// Returns iterator to beginning of JSON object or JSON array
    const_iterator begin() const;

    /// Returns iterator to end of JSON object or JSON array
    iterator end();

    /// Returns const_iterator to end of JSON object or JSON array
    const_iterator end() const;

    /// Returns const_iterator to beginning of JSON object or JSON array
    const_iterator cbegin() const;

    /// Returns const_iterator to end of JSON object or JSON array
    const_iterator cend() const;


    /// Return iterator to JSON object value with 'name'
    iterator find(const std::string &name);

    /// Return const_iterator to JSON object value with 'name'
    const_iterator find(const std::string &name) const;

    /// Return true if JSON object contains 'name'
    bool contains(const std::string &name) const;

    /// Erase key:value pair from an object
    size_type erase(const std::string &name);

    // ---- Data-retrieval Methods -----------------------------------------

    /// Get data from json, using one of several alternatives
    template<typename T, typename...Args>
    T get(Args... args) const;

    /// Get data from json, for any type T for which 'void from_json( T &value, const jsonParser &json, Args... args)' is defined
    ///   Call using: T t; json.get(t);
    template<typename T, typename...Args>
    void get(T &t, Args... args) const;

    /// Get data from json, if 'this' contains 'key'
    ///   Returns true if 'key' found, else false
    template<typename T, typename...Args>
    bool get_if(T &t, const std::string &key, Args... args) const;

    /// Get data from json, if 'this' contains 'key', else set to 'default_value'
    ///   Returns true if 'key' found, else false
    template<typename T, typename...Args>
    bool get_else(T &t, const std::string &key, const T &default_value, Args... args) const;


    // ---- Data addition Methods (Overwrites any existing data with same 'name') ---

    /// Puts data of any type T for which 'jsonParser& to_json( const T &value, jsonParser &json)' is defined
    template<typename T>
    jsonParser &operator=(const T &value);

    /// Puts new valued element at end of array of any type T for which 'jsonParser& to_json( const T &value, jsonParser &json)' is defined
    template<typename T>
    jsonParser &push_back(const T &value);

    /// Puts data of any type T for which 'jsonParser& to_json( const T &value, jsonParser &json)' is defined (same as 'operator=')
    template<typename T>
    jsonParser &put(const T &value);

    /// Puts new empty JSON object
    jsonParser &put_obj() {
      return *this = object();
    }

    /// Puts new JSON object, from iterators over a range of values of type std::pair<std::string, T>
    template<typename Iterator>
    jsonParser &put_obj(Iterator begin, Iterator end);

    /// Puts new empty JSON array
    jsonParser &put_array() {
      return *this = array();
    }

    /// Puts new JSON array
    jsonParser &put_array(size_type N) {
      return *this = array(N);
    }

    /// Puts new JSON array, using the same value
    template<typename T>
    jsonParser &put_array(size_type N, const T &t);

    /// Puts new JSON array, from iterators
    template<typename Iterator>
    jsonParser &put_array(Iterator begin,
                          Iterator end,
                          typename CASM_TMP::enable_if_iterator<Iterator>::type * = nullptr);

    /// Puts 'null' JSON value
    jsonParser &put_null() {
      return *this = null();
    }

    // ---- static Methods -------------------------------------

    /// Construct a jsonParser from a string containing JSON data
    static jsonParser parse(const std::string &str) {
      std::stringstream ss;
      ss << str;
      return jsonParser(ss);
    }

    /// Construct a jsonParser from a file containing JSON data
    static jsonParser parse(const fs::path &path) {
      return jsonParser(path);
    }

    /// Construct a jsonParser from a stream containing JSON data
    static jsonParser parse(std::istream &stream) {
      return jsonParser(stream);
    }

    /// Returns an empty json object
    static jsonParser object() {
      jsonParser json;
      return json = json_spirit::mValue(json_spirit::mObject());
    }

    /// Puts new JSON object, from iterators over a range of values of type std::pair<std::string, T>
    template<typename Iterator>
    static jsonParser object(Iterator begin, Iterator end) {
      jsonParser json;
      return json.put_obj(begin, end);
    }

    /// Returns an empty json array
    static jsonParser array() {
      jsonParser json;
      return json = json_spirit::mValue(json_spirit::mArray());
    }

    /// Returns an empty json array
    static jsonParser array(size_type N) {
      jsonParser json;
      return json = json_spirit::mValue(json_spirit::mArray(N));
    }

    /// Puts new JSON array, using the same value
    template<typename T>
    static jsonParser array(size_type N, const T &t) {
      jsonParser json;
      return json.put_array(N, t);
    }

    /// Puts new JSON array, from iterators
    template<typename Iterator>
    static jsonParser array(Iterator begin,
                            Iterator end,
                            typename CASM_TMP::enable_if_iterator<Iterator>::type * = nullptr) {
      jsonParser json;
      return json.put_array(begin, end);
    }


    /// Returns a null JSON value
    static jsonParser null() {
      jsonParser json;
      return json = json_spirit::mValue();
    }

  private:

    jsonParser &operator=(const json_spirit::mValue &value) {
      this->json_spirit::mValue::operator=(value);
      return *this;
    }


  };

  std::ostream &operator<<(std::ostream &stream, const jsonParser &json);
  std::istream &operator>>(std::istream &stream, jsonParser &json);

  /// To JSON for basic types
  jsonParser &to_json(bool value, jsonParser &json);
  jsonParser &to_json(int value, jsonParser &json);
  jsonParser &to_json(unsigned int value, jsonParser &json);
  jsonParser &to_json(long int value, jsonParser &json);
  jsonParser &to_json(unsigned long int value, jsonParser &json);
  jsonParser &to_json(double value, jsonParser &json);
  jsonParser &to_json(const std::string &value, jsonParser &json);
  jsonParser &to_json(const char *value, jsonParser &json);
  jsonParser &to_json(const jsonParser &value, jsonParser &json);


  /// From JSON for basic types
  template<typename T>
  T from_json(const jsonParser &json);


  template<> bool from_json<bool>(const jsonParser &json);
  template<> int from_json<int>(const jsonParser &json);
  template<> unsigned int from_json<unsigned int>(const jsonParser &json);
  template<> long int from_json<long int>(const jsonParser &json);
  template<> unsigned long int from_json<unsigned long int>(const jsonParser &json);
  template<> double from_json<double>(const jsonParser &json);
  template<> std::string from_json<std::string>(const jsonParser &json);
  template<> jsonParser from_json<jsonParser>(const jsonParser &json);

  /// From JSON for basic types
  void from_json(bool &value, const jsonParser &json);
  void from_json(int &value, const jsonParser &json);
  void from_json(unsigned int &value, const jsonParser &json);
  void from_json(long int &value, const jsonParser &json);
  void from_json(unsigned long int &value, const jsonParser &json);
  void from_json(double &value, const jsonParser &json);
  void from_json(std::string &value, const jsonParser &json);
  void from_json(jsonParser &value, const jsonParser &json);
  void from_json(std::istream &stream, const jsonParser &json);
  void from_json(fs::path &value, const jsonParser &json);

  /// Create a jsonParser from a stream
  inline void to_json(std::istream &stream, jsonParser &json) {
    if(!json.read(stream)) {
      throw std::runtime_error(
        std::string("ERROR: Could not read JSON. Please check your formatting. "
                    "For instance, try http://www.jsoneditoronline.org."));
    }
  }

  /// Create a jsonParser by reading a file
  ///
  /// This function reads the contents of the file at 'file_path' as if it were JSON.
  /// Use 'to_json(file_path.string(), json)' if you only want the path as a string
  inline void to_json(fs::path file_path, jsonParser &json) {
    if(!json.read(file_path)) {
      throw std::runtime_error(
        std::string("ERROR: Could not read JSON file: '") + file_path.string() +
        "'.\n\nPlease check your formatting. For instance, try http://www.jsoneditoronline.org.");
    }
  }

  /// To JSON for complex
  template<typename T>
  jsonParser &to_json(const std::complex<T> &value, jsonParser &json) {
    json = jsonParser::object();
    json["real"] = value.real();
    json["imag"] = value.imag();
    return json;
  }

  /// From JSON for complex
  template<typename T>
  void from_json(std::complex<T> &value, const jsonParser &json) {
    try {
      value = std::complex<T>(json["real"].get<T>(), json["imag"].get<T>());
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


  /// To JSON for std::pair<std::string, T>
  template<typename Key, typename T>
  jsonParser &to_json(const std::pair<Key, T> &value, jsonParser &json);

  /// From JSON for std::pair<std::string, T>
  template<typename Key, typename T>
  void from_json(std::pair<Key, T> &value, const jsonParser &json);

  /// \brief Helper struct for constructing objects that need additional data
  ///
  /// \code jsonParser::get<T>(Args...args) \endcode is equivalent to:
  /// - \code jsonConstructor<T>::from_json(*this, args...) \endcode
  ///
  /// This struct can be specialized to create new jsonConstructor<T>::from_json
  /// as needed.
  template<typename ReturnType>
  struct jsonConstructor {

    /// \brief Default from_json is equivalent to \code CASM::from_json<ReturnType>(json) \endcode
    static ReturnType from_json(const jsonParser &json) {
      return CASM::from_json<ReturnType>(json);
    }
  };


  /// Default works if T::T(args...) and 'void from_json(T&, const jsonParser&)' exist
  template<typename T>
  T from_json(const jsonParser &json) {
    T value;
    from_json(value, json);
    return value;
  }


  /// Return the location at which jsonParser 'A' != 'B' as a boost::filesystem::path
  boost::filesystem::path find_diff(const jsonParser &A, const jsonParser &B, boost::filesystem::path diff = boost::filesystem::path());

  /// Return the location at which jsonParser !A.almost_equal(B, tol) as a boost::filesystem::path
  boost::filesystem::path find_diff(const jsonParser &A, const jsonParser &B, double tol, boost::filesystem::path diff = boost::filesystem::path());


  /// jsonParser bidirectional Iterator class
  ///   Can iterate over a JSON object or JSON array or JSON value (though this is only one value)
  ///   When iterating over a JSON object, can return the 'name' of the current 'name':value pair being pointed at
  template<bool IsConst>
  class jsonParserIterator {

  public:

    typedef std::forward_iterator_tag iterator_category;
    typedef typename std::conditional<IsConst, json_spirit::mObject::const_iterator, json_spirit::mObject::iterator>::type object_iterator;
    typedef typename std::conditional<IsConst, json_spirit::mArray::const_iterator, json_spirit::mArray::iterator>::type array_iterator;
    typedef int difference_type;
    typedef typename std::conditional<IsConst, const jsonParser, jsonParser>::type value_type;
    typedef value_type &reference;
    typedef value_type *pointer;


    jsonParserIterator() {}

    jsonParserIterator(const jsonParserIterator &iter)
      : parser(iter.parser), type(iter.type), obj_iter(iter.obj_iter), array_iter(iter.array_iter), val_iter(iter.val_iter) {
    }

    jsonParserIterator &operator=(jsonParserIterator iter) {
      swap(*this, iter);
      return *this;
    }

    jsonParserIterator(pointer j, const object_iterator &iter)
      : parser(j), type(json_spirit::obj_type), obj_iter(iter) {
    }

    jsonParserIterator(pointer j, const array_iterator &iter)
      : parser(j), type(json_spirit::array_type), array_iter(iter) {
    }

    jsonParserIterator(pointer j, const int &iter)
      : parser(j), type(json_spirit::null_type), val_iter(iter) {
    }

    reference operator*() {
      if(type == json_spirit::obj_type)
        return (reference) obj_iter->second;
      else if(type == json_spirit::array_type)
        return (reference) * array_iter;
      else
        return *parser;
    }

    pointer operator->() {
      if(type == json_spirit::obj_type)
        return (pointer) &obj_iter->second;
      else if(type == json_spirit::array_type)
        return (pointer) & (*array_iter);
      else
        return parser;
    }


    bool operator==(const jsonParserIterator &iter) {
      if(parser != iter.parser)
        return false;

      if(type == json_spirit::obj_type)
        return obj_iter == iter.obj_iter;
      else if(type == json_spirit::array_type)
        return array_iter == iter.array_iter;
      else
        return true;
    }

    bool operator!=(const jsonParserIterator &iter) {
      return !(*this == iter);
    }

    jsonParserIterator &operator++() {
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

    jsonParserIterator operator++(int) {

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

    jsonParserIterator &operator--() {
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

    jsonParserIterator operator--(int) {

      jsonParserIterator cp(*this);

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

    operator jsonParser::const_iterator() {
      if(type == json_spirit::obj_type)
        return jsonParser::const_iterator(parser, obj_iter);
      else if(type == json_spirit::array_type)
        return jsonParser::const_iterator(parser, array_iter);
      else
        return jsonParser::const_iterator(parser, val_iter);
    }

    /// When iterating over a JSON object, returns the 'name' of the 'name':value pair the iterator is pointing at
    std::string name() {
      if(type == json_spirit::obj_type)
        return obj_iter->first;
      else
        throw std::runtime_error("Calling 'name' on non-object jsonParserIterator");
    }

    friend void swap(jsonParserIterator &a, jsonParserIterator &b) {
      using std::swap;

      std::swap(a.parser, b.parser);
      swap(a.type, b.type);
      swap(a.obj_iter, b.obj_iter);
      swap(a.array_iter, b.array_iter);
      swap(a.val_iter, b.val_iter);
    }


  private:

    pointer parser;

    json_spirit::Value_type type;

    object_iterator obj_iter;

    array_iterator array_iter;

    int val_iter;

  };

  /// Puts new valued element at end of array of any type T for which 'jsonParser& to_json( const T &value, jsonParser &json)' is defined
  template<typename T>
  jsonParser &jsonParser::push_back(const T &value) {
    try {
      jsonParser json;
      get_array().push_back(to_json(value, json));
      return *this;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  /// Get data from json, using one of several alternatives
  ///
  /// Use for any type T for which the either of the following is specialized
  /// (they are called in the following order):
  /// - \code
  ///   template<typename T>
  ///   template<typename...Args>
  ///   T jsonConstructor<T>::from_json(const jsonParser& json, Args...args);
  ///   \endcode
  /// - \code
  ///   template<typename T>
  ///   T from_json(const jsonParser &json);
  ///   \endcode
  /// If neither is specialized, then this is equivalent to:
  /// - \code
  ///   T value;
  ///   from_json(value, *this);
  ///   return value;
  ///   \endcode
  ///
  template<typename T, typename...Args>
  T jsonParser::get(Args... args) const {
    return jsonConstructor<T>::from_json(*this, args...);
  }

  template<typename T, typename...Args>
  void jsonParser::get(T &t, Args... args) const {
    from_json(t, *this, args...);
  }


  template<typename T, typename...Args>
  bool jsonParser::get_if(T &t, const std::string &key, Args... args) const {
    if(find(key) != cend()) {
      from_json(t, (*this)[key], args...);
      return true;
    }
    return false;
  }

  template<typename T, typename...Args>
  bool jsonParser::get_else(T &t, const std::string &key, const T &default_value, Args... args) const {
    if(find(key) != cend()) {
      from_json(t, (*this)[key], args...);
      return true;
    }

    t = default_value;
    return false;
  }

  /// Puts data of any type T for which 'jsonParser& to_json( const T &value, jsonParser &json)' is defined (same as operator=)
  template <typename T>
  jsonParser &jsonParser::put(const T &value) {
    return to_json(value, *this);
  }

  /// Puts new JSON object, from iterators over a range of values of type std::pair<std::string, T>
  template<typename Iterator>
  jsonParser &jsonParser::put_obj(Iterator begin,
                                  Iterator end) {
    *this = object();
    for(auto it = begin; it != end; ++it) {
      to_json(it->second, (*this)[it->first]);
    }
    return *this;
  }

  /// Puts new JSON array, using the same value
  template<typename T>
  jsonParser &jsonParser::put_array(size_type N, const T &t) {
    *this = array();
    for(auto i = 0; i < N; ++i) {
      push_back(t);
    }
    return *this;
  }

  /// Puts new JSON array, from iterators
  template<typename Iterator>
  jsonParser &jsonParser::put_array(
    Iterator begin,
    Iterator end,
    typename CASM_TMP::enable_if_iterator<Iterator>::type *) {

    *this = array();
    for(auto it = begin; it != end; ++it) {
      push_back(*it);
    }
    return *this;
  }

  /// Puts data of any type T for which 'jsonParser& to_json( const T &value, jsonParser &json)' is defined
  template <typename T>
  jsonParser &jsonParser::operator=(const T &value) {
    return to_json(value, *this);
  }

  /// To JSON for std::pair<std::string, T> and other convertible types
  template<typename Key, typename T>
  jsonParser &to_json(const std::pair<Key, T> &value, jsonParser &json) {
    json = jsonParser::object();
    return json[value.first] = value.second;
  }

  /// From JSON for std::pair<std::string, T>
  template<typename Key, typename T>
  void from_json(std::pair<Key, T> &value, const jsonParser &json) {
    auto it = json.begin();
    value = std::make_pair<Key, T>(it.name(), *it);
  }

  /** @} */
}

#endif
