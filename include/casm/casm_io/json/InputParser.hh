#ifndef CASM_InputParser
#define CASM_InputParser

#include <boost/filesystem.hpp>
#include <set>
#include <map>
#include <memory>
#include <string>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/Validator.hh"

namespace CASM {

  class Log;

  /// Base data structure for parsing values from JSON and storing error and warning messages
  ///
  /// This contains two jsonParser references:
  /// - input: which is a reference to the top of the JSON document (technically, this can be any
  ///   object in the JSON document from which all subsequent paths with be relative to).
  /// - self: which is a reference to the JSON subobject `this->input.find_at(this->path)` that will
  ///   be parsed when member functions (require, optional, etc.) are called and which any error
  ///   or warning messages they create are referenced to.
  ///
  /// Basic example:
  /// \code
  /// jsonParser json = jsonParser::parse(R"({
  ///   "an_options_data_structure": {
  ///     "an_integer_parameter" : 20,
  ///     "another_integer_parameter": "a string",
  ///   }
  /// })";
  /// bool is_required = true;
  /// KwargsParser parser {json, fs::path {"an_options_data_structure"}, is_required};
  /// int an_integer_parameter;
  /// // this will assign 20 to `an_integer_parameter` without error
  /// parser.require(an_integer_parameter, "an_integer_parameter");
  ///
  /// int another_integer_parameter;
  /// // this will assign nothing `another_integer_parameter` and store an error message of type:
  /// // "Error: could not construct type 'int' from option 'another_integer_parameter'."
  /// parser.require(another_integer_parameter, "another_integer_parameter");
  ///
  /// // this will construct a KwargsParser with error message of type:
  /// // "Error: Required property 'another_options_data_structure' not found."
  /// KwargsParser other_parser {json, fs::path {"another_options_data_structure"}, is_required};
  ///
  /// // since other_parser was constructed with an error, other_parser.valid() will be false
  /// assert(!other_parser.valid()};
  /// \endcode
  ///
  /// Similar options allow returning parsed values instead of assigning, optionally parsing if
  /// values exists in the JSON, specifying default values if no value exists, etc.
  ///
  /// The `make_report`, `print_errors` and `print_warnings` methods allow formatted output of error
  /// and warning messages.
  struct KwargsParser : public Validator {

    /// Reference to the top of the JSON document being parsed
    jsonParser const &input;

    /// Path to the JSON component to be parsed from the opt of the JSON document.
    /// If it exists, `this->self = this->input.at(this->path))`
    fs::path path;

    /// Reference to the JSON component to be parsed, if it exists, otherwise set to this->input
    ///
    /// Note:
    /// - if path.empty(), or path not found in input, self = input;
    ///   - else: self = *input.find_at(path)
    /// - Use this->exists() to check if the JSON component exists in the document
    jsonParser const &self;

    /// If this->input.at(this->path) is required to exist
    bool required;

    /// Default empty, can be used to differentiate between parsers when multiple values are parsed
    /// from a single JSON object
    std::string type_name;


    /// Construct KwargsParser
    ///
    /// \param _input Reference to top of the JSON document
    /// \param _path Location of JSON component, relative to input, to be parsed
    /// \param _required If true, insert an error if JSON component refered to by _path does not exist
    ///
    /// If required to exist, but does not, the KwargsParser will be constructed and the following
    /// error will be inserted:
    /// - "Error: Required property '<_path>' not found."
    KwargsParser(jsonParser const &_input, fs::path _path, bool _required);

    virtual ~KwargsParser() {}


    // --- Accessing the JSON `input` ---

    /// Return a const reference to the parent JSON object of this->self
    ///
    /// If self==input, returns self.
    const jsonParser &parent() const;

    /// Return a fs::path from the top-level to the parent JSON object of this->self
    ///
    /// The result satisfies: `this->parent()==this->input.find_at(this->parent_path())`
    fs::path parent_path() const;

    /// Name of this->self, equivalent to `this->path.filename().string()`
    std::string name() const;

    /// Check if this->input.find_at(this->path) exists
    ///
    /// When KwargsParser is contructed with a path that does not exist, then this->self=this->input.
    bool exists() const;

    /// Return this->path / val, ensuring the result is a relative path
    fs::path relpath(const fs::path &val) const {
      return path.empty() ? val : path / val;
    }


    // --- Subparsers ---

    typedef std::multimap<fs::path, std::shared_ptr<KwargsParser>> map_type;

    /// Begin iterator over subparsers
    map_type::const_iterator begin() const;

    /// End iterator over subparsers
    map_type::const_iterator end() const;

    /// Return true if this and and all subparsers are valid
    bool valid() const;

    /// Return warning messages from this and all subparsers
    std::map<fs::path, std::set<std::string>> all_warnings() const;

    /// Return error messages from this and all subparsers
    std::map<fs::path, std::set<std::string>> all_errors() const;

    using Validator::insert;

    /// Insert a subparser
    ///
    /// Subparsers are stored in a map of path (from this->input) to the subparser. After being
    /// inserted, the subparser's errors and warnings are included in validation checks
    /// (`valid`) and in parser output (`make_report`, `print_warnings`, `print_errors`, etc.).
    void insert(fs::path path, const std::shared_ptr<KwargsParser> &subparser);

    /// Insert a subparser at location `option` with a single error `message`
    void insert_error(fs::path option, std::string message);

    /// Insert a subparser at location `option` with a single warning `message`
    void insert_warning(fs::path option, std::string message);

    /// Insert a warning if any unexpected JSON attributes are found in self
    bool warn_unnecessary(const std::set<std::string> &expected);


  private:

    typedef map_type::value_type PairType;

    /// Used to store sub-parsers.
    /// Allows code re-use w/ storage of all errors & warnings.
    /// Can use static_cast to get values from InputParser subparsers if they are included.
    map_type m_subparsers;

  };

  /// Constructs values from JSON and collects error and warning messages for easy printing, without throwing.
  ///
  /// To use InputParser for a type, T, you must write:
  ///    void parse(InputParser<T> &parser, ... any other required input ...);
  ///
  /// The `parse` function should call methods of InputParser (i.e. require, optional, subparse,
  /// etc.) to automate the handling of errors and capture of error and warning messages.
  ///
  /// The `parse` function is called when InputParser<T> is constructed from input json:
  ///    jsonParser json_input = ...;
  ///    InputParser<T> parser {json_input, ... any other required input ...};
  ///
  /// When parsing complex objects, subparsers may be called and they will store their errors and
  /// warnings in a map of (path to the JSON subobject being parsed) : (string with an
  /// error or warning message).
  ///
  /// To check if parsing was successful and conditionally print errors & warnings, with the path
  /// to the location in the json where the issue occurred, do:
  ///     Log& log = ...;
  ///     if(!parser.valid()) {
  ///        jsonParser report = make_report(parser);
  ///        print_errors(parser, log);
  ///        log << std::endl << report << std::endl << std::endl;
  ///        ... handle error or throw ...
  ///     }
  ///     if(parser.all_warnings().size()) {
  ///        jsonParser report = make_report(parser);
  ///        print_warnings(parser, log);
  ///        log << std::endl << report << std::endl << std::endl;
  ///     }
  ///
  /// To do the above and throw an exception if the parser has any errors, use:
  ///     MyExceptionType error {"... message ..."};
  ///     report_and_throw_if_invalid(parser, log, error);
  ///
  /// Get the set of errors and warnings, including all subparsers, use:
  ///     std::set<std::string> warnings = parser.all_warnings();
  ///     std::set<std::string> errors = parser.all_errors();
  ///
  /// If parsing could not occur, then the `parse` function may leave `parser.value` as an empty
  /// unique_ptr<T>, otherwise parser.value should hold the resulting value. So if parsing was
  /// successful you can get the constructed value with `*parser.value`.
  ///
  template<typename T>
  class InputParser : public KwargsParser {
  public:

    /// Store the object being read from JSON, use unique_ptr so default constructor not necessary
    std::unique_ptr<T> value;

    /// Construct parser and use `parse(*this)`
    template<typename...Args>
    InputParser(jsonParser const &_input, Args &&... args);

    /// Construct parser and use `parse(*this, std::forward<Args>(args)...)` if `_path` exists
    template<typename...Args>
    InputParser(jsonParser const &_input, fs::path _path, bool _required, Args &&... args);

    /// Construct parser and use custom parse function, `f_parse(*this)`
    template<typename CustomParse>
    InputParser(CustomParse f_parse, jsonParser const &_input);

    /// Construct parser and use custom parse function, `f_parse(*this)`, if `_path` exists
    template<typename CustomParse>
    InputParser(CustomParse f_parse, jsonParser const &_input, fs::path _path, bool _required);


    /// Require self.find_at(option) of type RequiredType, returning result in unique_ptr
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: missing required option '<option>'"
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require(fs::path option, Args &&...args);

    /// Require self.find_at(option) of type RequiredType, assigning result to value
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: missing required option '<option>'"
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    template<typename RequiredType, typename...Args>
    void require(RequiredType &value, fs::path option, Args &&...args);

    /// Check that if self.find_at(option) exists it can constructed as type RequiredType, returning result in unique_ptr
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find_at(option) does not exist, return empty unique_ptr
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional(fs::path option, Args &&...args);

    /// Check that if self.find_at(option) exists it can constructed as type RequiredType, assigning result to value
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find_at(option) does not exist, do not change `value` and do not insert an error.
    template<typename RequiredType, typename...Args>
    void optional(RequiredType &value, fs::path option, Args &&...args);

    /// Check for self.find_at(option), return value or default, error if cannot be constructed
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find_at(option) does not exist, return `_default` and do not insert an error.
    template<typename RequiredType, typename...Args>
    RequiredType optional_else(fs::path option, const RequiredType &_default, Args &&...args);

    /// Check for self.find_at(option), assign result or default, error if cannot be constructed
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find(option) does not exist, assign `_default` to `value` and do not insert an error.
    template<typename RequiredType, typename...Args>
    void optional_else(RequiredType &value, fs::path option, const RequiredType &_default, Args &&...args);


    /// Run an InputParser on the JSON subobject at this->path / option, collecting errors and warnings
    ///
    /// Will:
    /// - Subparser errors and warnings are stored using this->insert
    ///
    /// Equivalent to:
    /// \code
    /// auto subparser = std::make_shared<InputParser<RequiredType>>(
    ///   this->input, this->relpath(option), true, std::forward<Args>(args)...);
    /// this->insert(subparser->path, subparser);
    /// return subparser;
    /// \endcode
    template<typename RequiredType, typename...Args>
    std::shared_ptr<InputParser<RequiredType>> subparse(fs::path option, Args &&...args);

    /// Subparse, if `this->path / option` exists
    ///
    /// If the JSON subobject does not exist, the result->value will be empty
    template<typename RequiredType, typename...Args>
    std::shared_ptr<InputParser<RequiredType>> subparse_if(fs::path option, Args &&...args);

    /// Subparse, if `this->path / option` exists, else the result->value will be copy-constructed from `_default`
    template<typename RequiredType, typename...Args>
    std::shared_ptr<InputParser<RequiredType>> subparse_else(fs::path option, const RequiredType &_default, Args &&...args);

    /// Parse `this->self` as RequiredType
    ///
    /// \param args Arguments forwared to the `parse` method
    ///
    template<typename RequiredType, typename...Args>
    std::shared_ptr<InputParser<RequiredType>> parse_as(Args &&...args);

  };

  /// Use when the JSON document is not associated with a single resulting value but instead is
  /// parsed with multiple subparsers
  typedef InputParser<std::nullptr_t> ParentInputParser;

  /// Parse Log "verbosity" level from JSON
  ///
  /// Expected that parser.self.find("verbosity") returns iterator to:
  /// \code
  /// {
  ///    "verbosity": <int or string, default=10, range=[0,100], "none" = 0, "quiet"=5,
  ///                  "standard"=10, "verbose"=20, "debug"=100>
  /// }
  /// \endcode
  int parse_verbosity(KwargsParser &parser, int default_verbosity = 10);

  /// Temporary -- enables compilation of legacy code
  void parse(InputParser<std::nullptr_t> &parser);

  template<typename T>
  void parse(InputParser<T> &parser);


  // --- Methods for formatting error and warning messages ---

  /// Formatted print warning messages, including all subparsers
  void print_warnings(KwargsParser const &parser, Log &log, std::string header = "Warnings");

  /// Formatted print error messages, including all subparsers
  void print_errors(KwargsParser const &parser, Log &log, std::string header = "Errors");

  /// Return parser.input with error and warning messages added in place, including all subparsers
  jsonParser make_report(KwargsParser const &parser);

  /// Print errors and warnings, throwing as specified if any errors exist in parser (and subparsers)
  template<typename ErrorType>
  void report_and_throw_if_invalid(KwargsParser const &parser, Log &log, ErrorType error);

}

#endif
