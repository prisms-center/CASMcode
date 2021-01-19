#include "casm/casm_io/Help.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

namespace {
jsonParser const &_self(jsonParser const &_input, fs::path _path) {
  if (_path.empty()) {
    return _input;
  }
  auto it = _input.find_at(_path);
  if (it == _input.end()) {
    return _input;
  }
  return *it;
}
}  // namespace

KwargsParser::KwargsParser(jsonParser const &_input, fs::path _path,
                           bool _required)
    : input(_input),
      path(_path),
      self(_self(_input, _path)),
      required(_required),
      type_name("") {
  if (required && !exists()) {
    error.insert(std::string("Error: ") + "Required property '" +
                 _path.string() + "' not found.");
  }
}

const jsonParser &KwargsParser::parent() const {
  if (parent_path().empty()) {
    return input;
  }
  return input[parent_path().string()];
}

fs::path KwargsParser::parent_path() const { return path.parent_path(); }

std::string KwargsParser::name() const { return path.filename().string(); }

bool KwargsParser::exists() const {
  return path.empty() || input.find_at(path) != input.end();
}

KwargsParser::map_type::const_iterator KwargsParser::begin() const {
  return m_subparsers.begin();
}

KwargsParser::map_type::const_iterator KwargsParser::end() const {
  return m_subparsers.end();
}

bool KwargsParser::valid() const {
  auto lambda = [](const PairType &pair) { return pair.second->valid(); };
  return Validator::valid() && std::all_of(begin(), end(), lambda);
}

std::map<fs::path, std::set<std::string>> KwargsParser::all_warnings() const {
  std::map<fs::path, std::set<std::string>> result;
  if (this->warning.size()) {
    result[this->path] = this->warning;
  }
  for (const auto &val : *this) {
    const auto &parser = *val.second;
    auto all_subparser_warnings = parser.all_warnings();
    if (all_subparser_warnings.size()) {
      result.insert(all_subparser_warnings.begin(),
                    all_subparser_warnings.end());
    }
  }
  return result;
}

std::map<fs::path, std::set<std::string>> KwargsParser::all_errors() const {
  std::map<fs::path, std::set<std::string>> result;
  if (this->error.size()) {
    result[this->path] = this->error;
  }
  for (const auto &val : *this) {
    const auto &parser = *val.second;
    auto all_subparser_errors = parser.all_errors();
    if (all_subparser_errors.size()) {
      result.insert(all_subparser_errors.begin(), all_subparser_errors.end());
    }
  }
  return result;
}

void KwargsParser::insert(fs::path path,
                          const std::shared_ptr<KwargsParser> &subparser) {
  m_subparsers.emplace(path, subparser);
}

void KwargsParser::insert_error(fs::path option, std::string message) {
  auto subparser =
      std::make_shared<KwargsParser>(this->input, this->relpath(option), false);
  subparser->error.insert(message);
  insert(subparser->path, subparser);
}

void KwargsParser::insert_warning(fs::path option, std::string message) {
  auto subparser =
      std::make_shared<KwargsParser>(this->input, this->relpath(option), false);
  subparser->warning.insert(message);
  insert(subparser->path, subparser);
}

bool KwargsParser::warn_unnecessary(const std::set<std::string> &expected) {
  bool all_necessary = true;
  for (auto opt_it = self.begin(); opt_it != self.end(); ++opt_it) {
    if (expected.find(opt_it.name()) == expected.end()) {
      warning.insert(std::string("Warning: ") + "Ignoring setting '" +
                     opt_it.name() + "' (it is unrecognized or unncessary).");
      all_necessary = false;
    }
  }
  return all_necessary;
}

int parse_verbosity(KwargsParser &parser, int default_verbosity) {
  auto it = parser.self.find("verbosity");
  if (it != parser.self.end()) {
    std::string verbosity_string;
    if (it->is_string()) {
      verbosity_string = it->get<std::string>();
    } else if (it->is_int()) {
      verbosity_string = std::to_string(it->get<int>());
    } else {
      parser.error.insert(Log::invalid_verbosity_msg(verbosity_string));
      return default_verbosity;
    }

    auto res = Log::verbosity_level(verbosity_string);
    if (res.first) {
      return res.second;
    } else {
      parser.error.insert(Log::invalid_verbosity_msg(verbosity_string));
      return default_verbosity;
    }
  } else {
    return default_verbosity;
  }
}

void parse(InputParser<std::nullptr_t> &parser) { return; }

void print_warnings(KwargsParser const &parser, Log &log, std::string header) {
  bool top = false;
  if (!header.empty()) {
    log.custom(header);
    header = "";
    top = true;
  }
  for (const auto &pair : parser.all_warnings()) {
    log << std::endl;
    std::string location = "/" + pair.first.string();
    log.custom(std::string("Warnings at location: ") + location);
    for (const auto &msg : pair.second) {
      log.indent() << msg << std::endl;
    }
  }
  if (top) log << std::endl;
}

void print_errors(KwargsParser const &parser, Log &log, std::string header) {
  bool top = false;
  if (!header.empty()) {
    log.custom(header);
    header = "";
    top = true;
  }
  for (const auto &pair : parser.all_errors()) {
    log << std::endl;
    std::string location = "/" + pair.first.string();
    log.custom(std::string("Errors at location: ") + location);
    for (const auto &msg : pair.second) {
      log.indent() << msg << std::endl;
    }
  }
  if (top) log << std::endl;
}

namespace make_report_impl {

jsonParser &get_parent(jsonParser &report, KwargsParser const &parser) {
  if (parser.parent_path().empty()) {
    return report;
  } else {
    return report[parser.parent_path().string()];
  }
}

void add_to_report(jsonParser &report, KwargsParser const &parser) {
  std::string sep = ".";
  if (parser.name().empty()) {
    sep = "";
  }
  jsonParser &parent = get_parent(report, parser);
  if (parser.warning.size()) {
    parent[parser.name() + sep + "WARNING"] = parser.warning;
    parent[parser.name() + sep + "WARNING"].set_force_column();
  }
  if (parser.error.size()) {
    parent[parser.name() + sep + "ERROR"] = parser.error;
    parent[parser.name() + sep + "ERROR"].set_force_column();
  }
}
}  // namespace make_report_impl

/// Return parser.input with error and warning messages added in place
jsonParser make_report(KwargsParser const &parser) {
  jsonParser report = parser.input;
  make_report_impl::add_to_report(report, parser);
  for (auto const &subparser_pair : parser) {
    make_report_impl::add_to_report(report, *subparser_pair.second);
  }
  return report;
}

}  // namespace CASM
