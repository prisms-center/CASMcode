#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/Help.hh"
#include "casm/casm_io/container/json_io.hh"

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
    required(_required) {
    if(required && !exists()) {
      error.insert(std::string("Error: ") + "Required property '" + _path.string() + "' not found.");
    }
  }

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

  bool KwargsParser::warn_unnecessary(const std::set<std::string> &expected) {
    bool all_necessary = true;
    for(auto opt_it = self.begin(); opt_it != self.end(); ++opt_it) {
      if(expected.find(opt_it.name()) == expected.end()) {
        warning.insert(
          std::string("Warning: ") + "Ignoring setting '" + opt_it.name() + "' (it is unrecognized or unncessary).");
        all_necessary = false;
      }
    }
    return all_necessary;
  }

  bool KwargsParser::valid() const {
    return Validator::valid();
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
    return path.empty() || input.find_at(path) != input.end();
  }

  int parse_verbosity(KwargsParser &parser, int default_verbosity) {
    auto it = parser.self.find("verbosity");
    if(it != parser.self.end()) {
      std::string verbosity_string;
      if(it->is_string()) {
        verbosity_string = it->get<std::string>();
      }
      else if(it->is_int()) {
        verbosity_string = std::to_string(it->get<int>());
      }
      else {
        parser.error.insert(Log::invalid_verbosity_msg(verbosity_string));
        return default_verbosity;
      }

      auto res = Log::verbosity_level(verbosity_string);
      if(res.first) {
        return res.second;
      }
      else {
        parser.error.insert(Log::invalid_verbosity_msg(verbosity_string));
        return default_verbosity;
      }
    }
    else {
      return default_verbosity;
    }
  }

  void parse(InputParser<std::nullptr_t> &parser) {
    return;
  }

}
