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

namespace CASM {

  /// Base struct for creating individual input parsers
  struct KwargsParser {

    jsonParser &input;
    fs::path path;

    /// result of 'parent().find(name())'; compare with parent().end()
    jsonParser::const_iterator self_it;
    bool required;

    std::set<std::string> warning;
    std::set<std::string> error;

    KwargsParser(jsonParser &_input, fs::path _path, bool _required);

    void print_warnings(Log &log, std::string header = "Warnings") const;

    void print_errors(Log &log, std::string header = "Errors") const;

    /// equivalent to require_at fs::path(it.name()) / option
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require(jsonParser::const_iterator it, std::string option, Args &&...args);

    /// require option self_it->find(option) of type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require(std::string option, Args &&...args);

    /// require option self_it->find(option) of type RequiredType
    template<typename RequiredType, typename...Args>
    std::unique_ptr<RequiredType> require_at(fs::path option, Args &&...args);

    /// add warning if setting in JSON object is unnecessary or unrecognized
    bool warn_unnecessary(const jsonParser &obj, fs::path path, const std::set<std::string> &expected);

    bool valid() const;

    void report();

    jsonParser &parent();

    const jsonParser &parent() const;

    fs::path parent_path() const;

    std::string name() const;

  };

  /// Base struct for creating input parsers that include multiple KwargsParser
  struct InputParser {

    jsonParser input;

    typedef std::map<fs::path, std::unique_ptr<KwargsParser>> map_type;
    map_type kwargs;
    typedef map_type::value_type PairType;

    InputParser(const jsonParser &_input);

    /// \brief Return true if all parsers in kwargs are valid
    bool valid() const;

    /// \brief Modifies this->input to include error and warning messages from all parsers in kwargs
    jsonParser &report();

    void print_warnings(Log &log, std::string header = "Warnings") const;

    void print_errors(Log &log, std::string header = "Errors") const;
  };
}

#endif
