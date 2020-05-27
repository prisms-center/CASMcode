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

  /// Base struct for creating individual input parsers
  struct KwargsParser : public Validator {

    // top-level
    jsonParser &input;

    // path to self from input
    fs::path path;

    // component to be parsed
    // if path.empty(), or path not found in input, self = input;
    // else: self = *input.find_at(path)
    jsonParser &self;

    // if required to exist
    bool required;

    // std::set<std::string> warning;
    // std::set<std::string> error;

    KwargsParser(jsonParser &_input, fs::path _path, bool _required);

    virtual ~KwargsParser() {}

    virtual void print_warnings(Log &log, std::string header = "Warnings") const;

    virtual void print_errors(Log &log, std::string header = "Errors") const;

    /// equivalent to require_at fs::path(it.name()) / option
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require(jsonParser::const_iterator it, std::string option, Args &&...args);

    /// require option self.find(option) of type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require(std::string option, Args &&...args);

    /// require option self.find_at(option) of type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require_at(fs::path option, Args &&...args);


    /// equivalent to require_at fs::path(it.name()) / option
    template<typename RequiredType, typename...Args>
    void require(RequiredType &value, jsonParser::const_iterator it, std::string option, Args &&...args);

    /// require option self.find(option) of type RequiredType
    template<typename RequiredType, typename...Args>
    void require(RequiredType &value, std::string option, Args &&...args);

    /// require option self.find_at(option) of type RequiredType
    template<typename RequiredType, typename...Args>
    void require_at(RequiredType &value, fs::path option, Args &&...args);


    /// equivalent to optional_at fs::path(it.name()) / option
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional(jsonParser::const_iterator it, std::string option, Args &&...args);

    /// check that if option self.find(option) exists it can constructed as type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional(std::string option, Args &&...args);

    /// check that if option self.find_at(option) exists it can constructed as type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional_at(fs::path option, Args &&...args);


    /// equivalent to optional_at fs::path(it.name()) / option
    template<typename RequiredType, typename...Args>
    void optional(RequiredType &value, jsonParser::const_iterator it, std::string option, Args &&...args);

    /// check that if option self.find(option) exists it can constructed as type RequiredType
    template<typename RequiredType, typename...Args>
    void optional(RequiredType &value, std::string option, Args &&...args);

    /// check that if option self.find_at(option) exists it can constructed as type RequiredType
    template<typename RequiredType, typename...Args>
    void optional_at(RequiredType &value, fs::path option, Args &&...args);


    /// check for option, return value or default, error if cannot be constructed
    template<typename RequiredType, typename...Args>
    RequiredType optional_else(std::string option, const RequiredType &_default, Args &&...args);

    /// check for option, return value or default, error if cannot be constructed
    template<typename RequiredType, typename...Args>
    RequiredType optional_at_else(fs::path option, const RequiredType &_default, Args &&...args);

    /// check for option, return value or default, error if cannot be constructed
    template<typename RequiredType, typename...Args>
    void optional_else(RequiredType &value, std::string option, const RequiredType &_default, Args &&...args);

    /// check for option, return value or default, error if cannot be constructed
    template<typename RequiredType, typename...Args>
    void optional_at_else(RequiredType &value, fs::path option, const RequiredType &_default, Args &&...args);


    /// add warning if unrecognized settings are found in self
    bool warn_unnecessary(const std::set<std::string> &expected);

    /// add warning if unrecognized settings are found in obj located at path
    bool warn_unnecessary(const jsonParser &obj, fs::path path, const std::set<std::string> &expected);

    virtual bool valid() const;

    virtual jsonParser &report();

    jsonParser &parent();

    const jsonParser &parent() const;

    fs::path parent_path() const;

    std::string name() const;

    bool exists() const;

    /// Return this->path / val, ensuring the result is a relative path
    fs::path relpath(const fs::path &val) const {
      return path.empty() ? val : path / val;
    }

  };

  /// Base struct for creating input parsers that include multiple KwargsParser
  ///
  /// The input parsers are supposed to construct C++ object from JSON and collect errors and
  /// warning messages for easy printing, without throwing. It calls
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
  /// warnings in a map of fs::path to the subobject -> std::string error or warning message.
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
  /// To do the above and just throw error, use:
  ///     report_and_throw_if_invalid(parser, log, error);
  ///
  /// Get the set of errors and warnings with:
  ///     std::set<std::string> warnings = parser.all_warnings();
  ///     std::set<std::string> errors = parser.all_errors();
  ///
  /// If parsing could not occur, then the `parse` function may leave parser.value as an empty
  /// unique_ptr<T>, otherwise parser.value should hold the resulting value. So if parsing was
  /// successful you can get the constructed value with:
  ///     return *parser.value;
  ///
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

    /// \brief Return true if all parsers in kwargs are valid
    bool valid() const override;

    /// \brief Modifies this->input to include error and warning messages from all parsers in kwargs
    jsonParser &report() override;

    void print_warnings(Log &log, std::string header = "Warnings") const override;

    void print_errors(Log &log, std::string header = "Errors") const override;

    std::set<std::string> all_warnings() const;

    std::set<std::string> all_errors() const;

    /// \brief If exists(), make value from current json; else add error
    template<typename...Args>
    void make(Args &&...args);

    /// \brief If exists(), make value from current json; else no error
    template<typename...Args>
    void make_if(Args &&...args);

    /// \brief If exists(), make value from current json; else set value to default
    template<typename...Args>
    void make_else(std::unique_ptr<T> _default, Args &&...args);

    /// \brief If exists(), make value from current json; else make with default constructor
    template<typename...Args>
    void make_else_construct_default(Args &&...args);

    /// \brief Run a subparser for a required subobject
    template<typename RequiredType, typename...Args>
    void subparse(RequiredType &_value, std::string option, Args &&...args);

    /// \brief Run a subparser for an optional subobject, else leave as is
    template<typename RequiredType, typename...Args>
    void subparse_if(RequiredType &_value, std::string option, Args &&...args);

    /// \brief Run a subparser for an optional subobject, else use default value
    template<typename RequiredType, typename...Args>
    void subparse_else(RequiredType &_value, std::string option, const RequiredType &_default, Args &&...args);

    /// \brief Parse self as a different type
    template<typename RequiredType, typename...Args>
    void parse_as(std::unique_ptr<RequiredType> &_value, Args &&...args);

    /// \brief Run a subparser for a derived type, ParseAsType, and store in base unique_ptr
    template<typename ParseAsType, typename...Args>
    void subparse_as(std::string option, Args &&...args);

    using Validator::insert;

    typedef std::map<fs::path, std::shared_ptr<KwargsParser>> map_type;

    /// \brief Insert subparser
    void insert(fs::path path, const std::shared_ptr<KwargsParser> &subparser);

    /// \brief Iterator over subparser
    map_type::const_iterator begin() const;

    /// \brief Iterator over subparser
    map_type::const_iterator end() const;

    /// \brief Find subparser
    map_type::const_iterator find(fs::path path) const;

  private:

    typedef map_type::value_type PairType;

    /// Can be used to store sub-parsers.
    /// Allows code re-use w/ storage of all errors & warnings.
    /// Can use static_cast to get values from sub-parsers.
    map_type kwargs;

  };

  template<typename T, typename ErrorType>
  void report_and_throw_if_invalid(InputParser<T> &parser, Log &log, ErrorType error);

  /// Parse Log "verbosity" level from JSON
  int parse_verbosity(KwargsParser &parser, int default_verbosity = 10);

  /// Temporary
  void parse(InputParser<std::nullptr_t> &parser);

}

#endif
