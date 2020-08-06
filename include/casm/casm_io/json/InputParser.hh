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
  /// The `report`, `print_errors` and `print_warnings` members allow formatted output of error and
  /// warning messages.
  struct KwargsParser : public Validator {

    /// Reference to the top of the JSON document being parsed
    jsonParser &input;

    /// Path to the JSON component to be parsed from the opt of the JSON document.
    /// If it exists, `this->self = this->input.at(this->path))`
    fs::path path;

    /// Reference to the JSON component to be parsed, if it exists, otherwise set to this->input
    ///
    /// Note:
    /// - if path.empty(), or path not found in input, self = input;
    ///   - else: self = *input.find_at(path)
    /// - Use this->exists() to check if the JSON component exists in the document
    jsonParser &self;

    /// If this->input.at(this->path) is required to exist
    bool required;

    /// Construct KwargsParser
    ///
    /// \param _input Reference to top of the JSON document
    /// \param _path Location of JSON component, relative to input, to be parsed
    /// \param _required If true, insert an error if JSON component refered to by _path does not exist
    ///
    /// If required to exist, but does not, the KwargsParser will be constructed and the following
    /// error will be inserted:
    /// - "Error: Required property '<_path>' not found."
    KwargsParser(jsonParser &_input, fs::path _path, bool _required);

    virtual ~KwargsParser() {}

    /// Formatted print warning messages
    virtual void print_warnings(Log &log, std::string header = "Warnings") const;

    /// Formatted print error messages
    virtual void print_errors(Log &log, std::string header = "Errors") const;


    /// Require self.find(option) of type RequiredType, returning result in unique_ptr
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: missing required option '<option>'"
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require(std::string option, Args &&...args);

    /// Require self.find_at(option) of type RequiredType, returning result in unique_ptr
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: missing required option '<option>'"
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require_at(fs::path option, Args &&...args);


    /// Require self.find(option) of type RequiredType, assigning result to value
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: missing required option '<option>'"
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    template<typename RequiredType, typename...Args>
    void require(RequiredType &value, std::string option, Args &&...args);

    /// Require self.find_at(option) of type RequiredType, assigning result to value
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: missing required option '<option>'"
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    template<typename RequiredType, typename...Args>
    void require_at(RequiredType &value, fs::path option, Args &&...args);


    /// Check that if self.find(option) exists it can constructed as type RequiredType, returning result in unique_ptr
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find(option) does not exist, return empty unique_ptr
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional(std::string option, Args &&...args);

    /// Check that if self.find_at(option) exists it can constructed as type RequiredType, returning result in unique_ptr
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find_at(option) does not exist, return empty unique_ptr
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional_at(fs::path option, Args &&...args);


    /// Check that if self.find(option) exists it can constructed as type RequiredType, assigning result to value
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find(option) does not exist, do not change `value` and do not insert an error.
    template<typename RequiredType, typename...Args>
    void optional(RequiredType &value, std::string option, Args &&...args);

    /// Check that if self.find_at(option) exists it can constructed as type RequiredType, assigning result to value
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find_at(option) does not exist, do not change `value` and do not insert an error.
    template<typename RequiredType, typename...Args>
    void optional_at(RequiredType &value, fs::path option, Args &&...args);


    /// Check for self.find(option), return value or default, error if cannot be constructed
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find(option) does not exist, return `_default` and do not insert an error.
    template<typename RequiredType, typename...Args>
    RequiredType optional_else(std::string option, const RequiredType &_default, Args &&...args);

    /// Check for self.find_at(option), return value or default, error if cannot be constructed
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find_at(option) does not exist, return `_default` and do not insert an error.
    template<typename RequiredType, typename...Args>
    RequiredType optional_at_else(fs::path option, const RequiredType &_default, Args &&...args);

    /// Check for self.find(option), assign result or default, error if cannot be constructed
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find(option) does not exist, assign `_default` to `value` and do not insert an error.
    template<typename RequiredType, typename...Args>
    void optional_else(RequiredType &value, std::string option, const RequiredType &_default, Args &&...args);

    /// Check for self.find_at(option), assign result or default, error if cannot be constructed
    ///
    /// Failing to parse the component will result in errors of type:
    /// - "Error: could not construct type '<type_name<RequiredType>()>' from option '<option>'. <
    ///    singleline_help<RequiredType>()>"
    ///
    /// If self.find(option) does not exist, assign `_default` to `value` and do not insert an error.
    template<typename RequiredType, typename...Args>
    void optional_at_else(RequiredType &value, fs::path option, const RequiredType &_default, Args &&...args);


    /// Insert a warning if any unexpected JSON attributes are found in self
    bool warn_unnecessary(const std::set<std::string> &expected);


    /// Return true if this has no errors (may have warnings)
    virtual bool valid() const;

    /// Modifies input JSON document to include error and warning messages
    ///
    /// If this has any errors, they are inserted as a JSON array at this->parants() with name
    /// `this->name() + ".ERROR"`.
    /// If this has any warnings, they are inserted as a JSON array at this->parent() with name
    /// `this->name() + ".WARNING"`
    virtual jsonParser &report();

    /// Return a reference to the parent JSON object of this->self
    ///
    /// If self==input, returns self.
    jsonParser &parent();

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

  };

  /// Constructs values from JSON and collects error and warning messages for easy printing, without throwing.
  ///
  /// To use InputParser for a type, T, you must write:
  ///    void parse(InputParser<T> &parser, ... any other required input ...);
  ///
  /// The `parse` function should call methods of parser (i.e. require, optional, make, subparse,
  /// etc.) to automate the handling of errors and capture of error and warning messages.
  ///
  /// That function is called when InputParser<T> is constructed from input json:
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
  ///        parser.print_errors(log);
  ///        log << std::endl << parser.report() << std::endl << std::endl;
  ///        ... handle error or throw ...
  ///     }
  ///     if(parser.all_warnings().size()) {
  ///         parser.print_warnings(log);
  ///         log << std::endl << parser.report() << std::endl << std::endl;
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

    template<typename...Args>
    InputParser(jsonParser &_input, Args &&... args);

    template<typename...Args>
    InputParser(jsonParser &_input, fs::path _path, bool _required, Args &&... args);

    /// Return true if this and and all subparsers are valid
    bool valid() const override;

    /// Modifies input JSON document to include error and warning messages from this and all subparsers
    jsonParser &report() override;

    /// Formatted print warning messages from this and all subparsers
    void print_warnings(Log &log, std::string header = "Warnings") const override;

    /// Formatted print error messages from this and all subparsers
    void print_errors(Log &log, std::string header = "Errors") const override;

    /// Return warning messages from this and all subparsers
    std::set<std::string> all_warnings() const;

    /// Return error messages from this and all subparsers
    std::set<std::string> all_errors() const;

    /// If exists(), make this->value from JSON; else add error
    ///
    /// Note:
    /// - If `this->exists()==false`, will insert an error of type:
    ///   - "Error: Required input at '<this->path>' does not exist."
    /// - Otherwise, make from JSON: `this->value = self.make<T>(std::forward<Args>(args)...)`
    ///   - If make from JSON fails, will insert an error of type:
    ///     - "Error: could not construct type '<type_name<T>()>' from option '<option>'.
    ///       <singleline_help<T>()>"
    template<typename...Args>
    void make(Args &&...args);

    /// If exists(), make this->value from JSON; else no error
    ///
    /// Note:
    /// - If `this->exists()==false`: do nothing
    /// - Otherwise, equivalent to `this->make(std::forward<Args>(args)...)`
    template<typename...Args>
    void make_if(Args &&...args);

    /// If exists(), make this->value from JSON; else set this->value to _default
    ///
    /// Note:
    /// - If `this->exists()==false`, set this->value to _default: `this->value = std::move(_default)`
    /// - Otherwise, equivalent to `this->make(std::forward<Args>(args)...)`
    template<typename...Args>
    void make_else(std::unique_ptr<T> _default, Args &&...args);

    /// If exists(), make this->value from current json; else make with default constructor
    ///
    /// Note:
    /// - If `this->exists()==false`, construct default value: `this->value = notstd::make_unique<T>()`
    /// - Otherwise, equivalent to `this->make(std::forward<Args>(args)...)`
    template<typename...Args>
    void make_else_construct_default(Args &&...args);

    /// Run an InputParser on the JSON subobject at this->path / option, collecting errors and warnings
    ///
    /// Will:
    /// - If the subparser constructs a value, it will be move-assigned to _value
    /// - Subparser errors and warnings are stored using this->insert
    ///
    /// Equivalent to:
    /// \code
    /// auto subparser = std::make_shared<InputParser<RequiredType>>(
    ///   this->input, this->relpath(option), true, std::forward<Args>(args)...);
    /// this->insert(subparser->path, subparser);
    /// if(subparser->value) {
    ///   _value = std::move(*(subparser->value));
    /// }
    /// \endcode
    template<typename RequiredType, typename...Args>
    void subparse(RequiredType &_value, std::string option, Args &&...args);

    /// Subparse, if `this->path / option` exists
    ///
    /// If the JSON subobject does not exist, `_value` is unchanged, and no errors or warnings are
    /// inserted.
    template<typename RequiredType, typename...Args>
    void subparse_if(RequiredType &_value, std::string option, Args &&...args);

    /// Subparse, if `this->path / option` exists, else assign `_default`
    template<typename RequiredType, typename...Args>
    void subparse_else(RequiredType &_value, std::string option, const RequiredType &_default, Args &&...args);

    /// Parse `this->self` as a type derived from type `T`
    ///
    /// This can be used to construct a type derived from type `T` and assign it to `this->value`
    template<typename RequiredType, typename...Args>
    void parse_as(std::unique_ptr<RequiredType> &_value, Args &&...args);

    /// Parse a JSON subobject as a type derived from type `T`
    ///
    /// This can be used to construct a type derived from type `T` from a JSON subobject and assign
    /// it to `this->value`.
    template<typename ParseAsType, typename...Args>
    void subparse_as(std::string option, Args &&...args);

    using Validator::insert;

    typedef std::map<fs::path, std::shared_ptr<KwargsParser>> map_type;

    /// Insert a subparser
    ///
    /// Subparsers are stored in a map of path (from this->input) to the subparser. After being
    /// inserted, the subparser's errors and warnings are included in validation checks
    /// (`valid`) and in parser output (`report`, `print_warnings`, `print_errors`, etc.).
    void insert(fs::path path, const std::shared_ptr<KwargsParser> &subparser);

    /// Begin iterator over subparsers
    map_type::const_iterator begin() const;

    /// End iterator over subparsers
    map_type::const_iterator end() const;

    /// Find a subparser by path (from this->input)
    map_type::const_iterator find(fs::path path) const;

  private:

    typedef map_type::value_type PairType;

    /// Used to store sub-parsers.
    /// Allows code re-use w/ storage of all errors & warnings.
    /// Can use static_cast to get values from InputParser subparsers if they are included.
    map_type kwargs;

  };

  /// Print errors and warnings, throwing as specified if any errors exist in parser (and subparsers)
  template<typename T, typename ErrorType>
  void report_and_throw_if_invalid(InputParser<T> &parser, Log &log, ErrorType error);

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

}

#endif
