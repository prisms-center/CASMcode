#include "casm/casm_io/InputParser_impl.hh"
#include "casm/casm_io/Help.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  namespace {
    jsonParser &_self(jsonParser &_input, fs::path _path) {
      if(_path.empty()) {
        return _input;
      }
      auto it = _input.find_at(_path);
      if(it == _input.end()) {
        return _input;
      }
      return *it;
    }
  }

  KwargsParser::KwargsParser(jsonParser &_input, fs::path _path, bool _required) :
    input(_input),
    path(_path),
    self(_self(_input, _path)),
    required(_required) {}

  void KwargsParser::print_warnings(Log &log, std::string header) const {
    bool top = false;
    if(!header.empty()) {
      log.custom(header);
      header = "";
      top = true;
    }
    for(const auto &msg : warning) {
      log << msg << std::endl;
    }
    if(top) log << std::endl;
  }

  void KwargsParser::print_errors(Log &log, std::string header) const {
    bool top = false;
    if(!header.empty()) {
      log.custom(header);
      header = "";
      top = true;
    }
    for(const auto &msg : error) {
      log << msg << std::endl;
    }
    if(top) log << std::endl;
  }

  /// add warning if unrecognized settings are found in self
  bool KwargsParser::warn_unnecessary(const std::set<std::string> &expected) {
    return warn_unnecessary(self, fs::path(), expected);
  }

  /// add warning if unrecognized settings are found in obj located at path
  bool KwargsParser::warn_unnecessary(const jsonParser &obj, fs::path path, const std::set<std::string> &expected) {
    jsonParser json;
    bool all_necessary = true;
    for(auto opt_it = obj.begin(); opt_it != obj.end(); ++opt_it) {
      if(expected.find(opt_it.name()) == expected.end()) {
        std::string _path = path.empty() ? opt_it.name() : (path / opt_it.name()).string();
        warning.insert(
          std::string("Warning: ") + "Ignoring setting '" + _path + "' (it is unrecognized or unncessary).");
        all_necessary = false;
      }
    }
    return all_necessary;
  }

  bool KwargsParser::valid() const {
    return !error.size();
  }

  jsonParser &KwargsParser::report() {
    if(warning.size()) {
      parent()[name() + ".WARNING"] = warning;
      parent()[name() + ".WARNING"].set_force_column();
    }
    if(error.size()) {
      parent()[name() + ".ERROR"] = error;
      parent()[name() + ".ERROR"].set_force_column();
    }
    return input;
  }

  jsonParser &KwargsParser::parent() {
    if(parent_path().empty()) {
      return input;
    }
    return input[parent_path().string()];
  }

  const jsonParser &KwargsParser::parent() const {
    if(parent_path().empty()) {
      return input;
    }
    return input[parent_path().string()];
  }

  fs::path KwargsParser::parent_path() const {
    return path.parent_path();
  }

  std::string KwargsParser::name() const {
    return path.filename().string();
  }

  bool KwargsParser::exists() const {
    return input.find_at(path) != input.end();
  }


  std::string InputParser::dry_run_help() {
    return
      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n";
  }

  std::string InputParser::indent_space_help() {
    return
      "  indent_space: string (optional, default=6)\n"
      "    Number of spaces to indent for pretty-printing.\n\n";
  }

  std::string InputParser::coordinate_mode_help() {
    return
      "  coordinate_mode: string (optional, default=FRAC)\n"
      "    Coordinate mode (FRAC, CART, INTEGRAL) for printing orbits.\n\n";
  }

  std::string InputParser::orbit_print_mode_help() {
    return
      "  orbit_print_mode: string (optional, default=\"PROTO\")\n"
      "    Mode (FULL, PROTO) to select printing full orbits or just orbit prototypes.\n\n";
  }

  std::string InputParser::prec_help(std::string what, int default_prec) {
    std::stringstream ss;
    ss <<  "  prec: int (optional, default=" << default_prec << ")\n"
       <<  "    Precision for printing " << what << "\n\n";
    return ss.str();
  }

  std::string InputParser::verbosity_help() {
    return
      "  verbosity: string or int (optional, default=\"standard\")\n"
      "    Verbosity of output. Options are 'none', 'quiet', 'standard', 'verbose', 'debug',\n"
      "    or an integer 0-100 (0: 'none', 100: 'debug').\n\n";
  }


  InputParser::InputParser(jsonParser &_input, fs::path _path, bool _required):
    KwargsParser(_input, _path, _required) {}

  /// \brief Return true if all parsers in kwargs are valid
  bool InputParser::valid() const {
    auto lambda = [](const PairType & pair) {
      return pair.second->valid();
    };
    return KwargsParser::valid() && std::all_of(kwargs.begin(), kwargs.end(), lambda);
  }

  /// \brief Modifies this->input to include error and warning messages from all parsers in kwargs
  jsonParser &InputParser::report() {
    KwargsParser::report();
    std::for_each(kwargs.begin(), kwargs.end(), [](const PairType & pair) {
      pair.second->report();
    });
    return input;
  }

  void InputParser::print_warnings(Log &log, std::string header) const {
    auto lambda = [&](const PairType & pair) {
      for(const auto &msg : pair.second->warning) {
        log << msg << std::endl;
      }
    };

    bool top = false;
    if(!header.empty()) {
      log.custom(header);
      header = "";
      top = true;
    }
    KwargsParser::print_warnings(log, header);
    std::for_each(kwargs.begin(), kwargs.end(), lambda);
    if(top) log << std::endl;
  }

  void InputParser::print_errors(Log &log, std::string header) const {
    auto lambda = [&](const PairType & pair) {
      for(const auto &msg : pair.second->error) {
        log << msg << std::endl;
      }
    };

    bool top = false;
    if(!header.empty()) {
      log.custom(header);
      header = "";
      top = true;
    }
    KwargsParser::print_errors(log, header);
    std::for_each(kwargs.begin(), kwargs.end(), lambda);
    if(top) log << std::endl;
  }

  std::set<std::string> InputParser::all_warnings() const {
    std::set<std::string> res = this->warning;
    for(const auto &val : kwargs) {
      const auto &parser = *val.second;
      res.insert(parser.warning.begin(), parser.warning.end());
    }
    return res;
  }

  std::set<std::string> InputParser::all_errors() const {
    std::set<std::string> res = this->error;
    for(const auto &val : kwargs) {
      const auto &parser = *val.second;
      res.insert(parser.error.begin(), parser.error.end());
    }
    return res;
  }
}
