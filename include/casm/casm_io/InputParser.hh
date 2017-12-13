#ifndef CASM_InputParser
#define CASM_InputParser

#include <boost/filesystem.hpp>
#include <set>
#include <map>
#include <memory>
#include <string>

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/Log.hh"
#include "casm/misc/TypeInfo.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/app/AppIO_impl.hh"

namespace CASM {

  /// Base struct for creating individual input parsers
  struct KwargsParser {

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

    std::set<std::string> warning;
    std::set<std::string> error;

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


    /// equivalent to optional_at fs::path(it.name()) / option
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional(jsonParser::const_iterator it, std::string option, Args &&...args);

    /// check that if option self.find(option) exists it can constructed as type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional(std::string option, Args &&...args);

    /// check that if option self.find_at(option) exists it can constructed as type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> optional_at(fs::path option, Args &&...args);


    /// check for option, return value or default, error if cannot be constructed
    template<typename RequiredType, typename...Args>
    RequiredType optional_else(std::string option, const RequiredType &_default, Args &&...args);

    /// check for option, return value or default, error if cannot be constructed
    template<typename RequiredType, typename...Args>
    RequiredType optional_at_else(fs::path option, const RequiredType &_default, Args &&...args);


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


    // --- Optional, check CLI and JSON options ---

    template<typename OptHandlerType>
    int parse_verbosity(const OptHandlerType &opt);

    template<typename OptHandlerType>
    bool parse_dry_run(const OptHandlerType &opt);

    template<typename OptHandlerType>
    COORD_TYPE parse_coord_type(const OptHandlerType &opt);

    template<typename OptHandlerType>
    ORBIT_PRINT_MODE parse_orbit_print_mode(const OptHandlerType &opt);

    template<typename OptHandlerType>
    std::vector<std::string> parse_filter_expr(const OptHandlerType &opt);

  };

  /// Base struct for creating input parsers that include multiple KwargsParser
  struct InputParser : public KwargsParser {

    static std::string dry_run_help();
    static std::string coordinate_mode_help();
    static std::string orbit_print_mode_help();
    static std::string verbosity_help();


    typedef std::map<fs::path, std::shared_ptr<KwargsParser>> map_type;
    map_type kwargs;
    typedef map_type::value_type PairType;

    InputParser(jsonParser &_input, fs::path _path, bool _required);

    /// \brief Return true if all parsers in kwargs are valid
    bool valid() const override;

    /// \brief Modifies this->input to include error and warning messages from all parsers in kwargs
    jsonParser &report() override;

    void print_warnings(Log &log, std::string header = "Warnings") const override;

    void print_errors(Log &log, std::string header = "Errors") const override;

    std::set<std::string> all_warnings() const;

    std::set<std::string> all_errors() const;
  };

  /// \brief Return path `base / val`, ensuring result is a relative path
  struct Relpath {

    fs::path base;

    Relpath(const fs::path &_base) :
      base(_base) {}

    fs::path operator()(const fs::path &val) const {
      return base.empty() ? val : base / val;
    };
  };
}

#endif
