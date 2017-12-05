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
    log.custom(header);
    for(const auto &msg : warning) {
      log << msg << std::endl;
    }
    log << std::endl;
  }

  void KwargsParser::print_errors(Log &log, std::string header) const {
    log.custom(header);
    for(const auto &msg : error) {
      log << msg << std::endl;
    }
    log << std::endl;
  }

  /// add warning if setting in JSON object is unnecessary or unrecognized
  bool KwargsParser::warn_unnecessary(const jsonParser &obj, fs::path path, const std::set<std::string> &expected) {
    bool all_necessary = true;
    for(auto opt_it = obj.begin(); opt_it != obj.end(); ++opt_it) {
      if(expected.find(opt_it.name()) == expected.end()) {
        warning.insert(
          std::string("Warning: ") + "Ignoring setting '" + (path / opt_it.name()).string() + "' (it is unrecognized or unncessary).");
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


  InputParser::InputParser(jsonParser &_input, fs::path _path, bool _required):
    KwargsParser(_input, _path, _required) {}

  /// \brief Return true if all parsers in kwargs are valid
  bool InputParser::valid() const {
    auto lambda = [](const PairType & pair) {
      return pair.second->valid();
    };
    return std::all_of(kwargs.begin(), kwargs.end(), lambda);
  }

  /// \brief Modifies this->input to include error and warning messages from all parsers in kwargs
  jsonParser &InputParser::report() {
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

    log.custom(header);
    std::for_each(kwargs.begin(), kwargs.end(), lambda);
    log << std::endl;
  }

  void InputParser::print_errors(Log &log, std::string header) const {
    auto lambda = [&](const PairType & pair) {
      for(const auto &msg : pair.second->error) {
        log << msg << std::endl;
      }
    };

    log.custom(header);
    std::for_each(kwargs.begin(), kwargs.end(), lambda);
    log << std::endl;
  }

  std::set<std::string> InputParser::all_warning() const {
    std::set<std::string> res = this->warning;
    for(const auto &val : kwargs) {
      const auto &parser = *val.second;
      res.insert(parser.warning.begin(), parser.warning.end());
    }
    return res;
  }

  std::set<std::string> InputParser::all_error() const {
    std::set<std::string> res = this->error;
    for(const auto &val : kwargs) {
      const auto &parser = *val.second;
      res.insert(parser.error.begin(), parser.error.end());
    }
    return res;
  }
}
