#ifndef CASM_InputParser_impl
#define CASM_InputParser_impl

#include "casm/casm_io/json/InputParser.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/casm_io/Help.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::require(std::string option, Args &&...args) {
    return require_at<RequiredType>(fs::path(option), std::forward<Args>(args)...);
  }

  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::require_at(fs::path option, Args &&...args) {
    auto it = self.find_at(option);
    std::unique_ptr<RequiredType> res;
    if(it == self.end()) {
      error.insert(std::string("Error: missing required option '") + option.string() + "'.");
      return res;
    }

    try {
      return it->make<RequiredType>(std::forward<Args>(args)...);
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return res;
    }
  }


  template<typename RequiredType, typename...Args>
  void KwargsParser::require(RequiredType &value, std::string option, Args &&...args) {
    require_at<RequiredType>(value, fs::path(option), std::forward<Args>(args)...);
  }

  template<typename RequiredType, typename...Args>
  void KwargsParser::require_at(RequiredType &value, fs::path option, Args &&...args) {
    auto it = self.find_at(option);
    std::unique_ptr<RequiredType> res;
    if(it == self.end()) {
      error.insert(std::string("Error: missing required option '") + option.string() + "'.");
      return;
    }

    try {
      it->get<RequiredType>(value, std::forward<Args>(args)...);
      return;
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return;
    }
  }


  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::optional(std::string option, Args &&...args) {
    return optional_at<RequiredType>(fs::path(option), std::forward<Args>(args)...);
  }

  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::optional_at(fs::path option, Args &&...args) {
    jsonParser *ptr;
    if(option.empty()) {
      ptr = &self;
    }
    else {
      auto it = self.find_at(option);
      if(it == self.end()) {
        return std::unique_ptr<RequiredType>();
      }
      else {
        ptr = &(*it);
      }
    }

    try {
      return ptr->make<RequiredType>(std::forward<Args>(args)...);
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return std::unique_ptr<RequiredType>();
    }
  }


  template<typename RequiredType, typename...Args>
  void KwargsParser::optional(RequiredType &value, std::string option, Args &&...args) {
    optional_at<RequiredType>(value, fs::path(option), std::forward<Args>(args)...);
  }

  template<typename RequiredType, typename...Args>
  void KwargsParser::optional_at(RequiredType &value, fs::path option, Args &&...args) {
    jsonParser *ptr;
    if(option.empty()) {
      ptr = &self;
    }
    else {
      auto it = self.find_at(option);
      if(it == self.end()) {
        return;
      }
      else {
        ptr = &(*it);
      }
    }

    try {
      ptr->get<RequiredType>(value, std::forward<Args>(args)...);
      return;
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return;
    }
  }


  template<typename RequiredType, typename...Args>
  RequiredType KwargsParser::optional_else(std::string option, const RequiredType &_default, Args &&...args) {
    return optional_at_else<RequiredType>(fs::path(option), _default, std::forward<Args>(args)...);
  }

  template<typename RequiredType, typename...Args>
  RequiredType KwargsParser::optional_at_else(fs::path option, const RequiredType &_default, Args &&...args) {
    jsonParser *ptr;
    if(option.empty()) {
      ptr = &self;
    }
    else {
      auto it = self.find_at(option);
      if(it == self.end()) {
        return _default;
      }
      else {
        ptr = &(*it);
      }
    }

    try {
      return ptr->get<RequiredType>(std::forward<Args>(args)...);
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return _default;
    }
  }


  template<typename RequiredType, typename...Args>
  void KwargsParser::optional_else(RequiredType &value, std::string option, const RequiredType &_default, Args &&...args) {
    optional_at_else<RequiredType>(value, fs::path(option), _default, std::forward<Args>(args)...);
  }

  template<typename RequiredType, typename...Args>
  void KwargsParser::optional_at_else(RequiredType &value, fs::path option, const RequiredType &_default, Args &&...args) {
    jsonParser *ptr;
    if(option.empty()) {
      ptr = &self;
    }
    else {
      auto it = self.find_at(option);
      if(it == self.end()) {
        value = _default;
        return;
      }
      else {
        ptr = &(*it);
      }
    }

    try {
      ptr->get<RequiredType>(value, std::forward<Args>(args)...);
      return;
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return;
    }
  }


  template<typename T>
  template<typename... Args>
  InputParser<T>::InputParser(jsonParser &_input, Args &&... args):
    KwargsParser(_input, "", true) {
    parse(*this, std::forward<Args>(args)...);
  }

  template<typename T>
  template<typename... Args>
  InputParser<T>::InputParser(jsonParser &_input, fs::path _path, bool _required, Args &&... args):
    KwargsParser(_input, _path, _required) {
    parse(*this, std::forward<Args>(args)...);
  }

  template<typename T>
  bool InputParser<T>::valid() const {
    auto lambda = [](const PairType & pair) {
      return pair.second->valid();
    };
    return KwargsParser::valid() && std::all_of(kwargs.begin(), kwargs.end(), lambda);
  }

  template<typename T>
  jsonParser &InputParser<T>::report() {
    KwargsParser::report();
    std::for_each(kwargs.begin(), kwargs.end(), [](const PairType & pair) {
      pair.second->report();
    });
    return input;
  }

  template<typename T>
  void InputParser<T>::print_warnings(Log &log, std::string header) const {
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

  template<typename T>
  void InputParser<T>::print_errors(Log &log, std::string header) const {
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

  template<typename T>
  std::set<std::string> InputParser<T>::all_warnings() const {
    std::set<std::string> res = this->warning;
    for(const auto &val : kwargs) {
      const auto &parser = *val.second;
      res.insert(parser.warning.begin(), parser.warning.end());
    }
    return res;
  }

  template<typename T>
  std::set<std::string> InputParser<T>::all_errors() const {
    std::set<std::string> res = this->error;
    for(const auto &val : kwargs) {
      const auto &parser = *val.second;
      res.insert(parser.error.begin(), parser.error.end());
    }
    return res;
  }

  template<typename T>
  template<typename...Args>
  void InputParser<T>::make(Args &&...args) {
    if(!exists()) {
      error.insert(std::string("Error: ") + "Required input at '" + path.string() + "' does not exist.");
    }
    else {
      try {
        this->value = self.make<T>(std::forward<Args>(args)...);
      }
      catch(std::exception &e) {
        error.insert(std::string("Error: could not construct type '")
                     + type_name<T>() + "' from input at '" + path.string() + "'. "
                     + singleline_help<T>());
      }
    }
  }

  template<typename T>
  template<typename...Args>
  void InputParser<T>::make_if(Args &&...args) {
    if(exists()) {
      this->make(std::forward<Args>(args)...);
    }
  }

  template<typename T>
  template<typename...Args>
  void InputParser<T>::make_else(std::unique_ptr<T> _default, Args &&...args) {
    if(!exists()) {
      this->value = std::move(_default);
    }
    else {
      this->make(std::forward<Args>(args)...);
    }
  }

  template<typename T>
  template<typename...Args>
  void InputParser<T>::make_else_construct_default(Args &&...args) {
    if(!exists()) {
      this->value = notstd::make_unique<T>();
    }
    else {
      this->make(std::forward<Args>(args)...);
    }
  }

  template<typename T>
  template<typename RequiredType, typename...Args>
  void InputParser<T>::subparse(RequiredType &_value, std::string option, Args &&...args) {

    auto subparser = std::make_shared<InputParser<RequiredType>>(
                       this->input, this->relpath(option), true, std::forward<Args>(args)...);
    kwargs[subparser->path] = subparser;
    if(subparser->value) {
      _value = std::move(*(subparser->value));
    }
  }

  template<typename T>
  template<typename RequiredType, typename...Args>
  void InputParser<T>::subparse_if(RequiredType &_value, std::string option, Args &&...args) {
    if(self.contains(option)) {
      subparse(_value, option, std::forward<Args>(args)...);
    }
  }

  template<typename T>
  template<typename RequiredType, typename...Args>
  void InputParser<T>::subparse_else(RequiredType &_value, std::string option, const RequiredType &_default, Args &&...args)  {
    if(self.contains(option)) {
      subparse(_value, option, std::forward<Args>(args)...);
    }
    else {
      _value = _default;
    }
  }

  template<typename T>
  template<typename RequiredType, typename...Args>
  void InputParser<T>::parse_as(std::unique_ptr<RequiredType> &_value, Args &&...args) {

    auto subparser = std::make_shared<InputParser<RequiredType>>(
                       this->input, this->path, true, std::forward<Args>(args)...);
    kwargs[subparser->path] = subparser;
    _value = std::move(subparser->value);
  }

  template<typename T>
  template<typename ParseAsType, typename...Args>
  void InputParser<T>::subparse_as(std::string option, Args &&...args) {

    auto subparser = std::make_shared<InputParser<ParseAsType>>(
                       this->input, this->relpath(option), true, std::forward<Args>(args)...);
    kwargs[subparser->path] = subparser;
    if(subparser->value) {
      this->value = std::move(subparser->value);
    }
  }

  template<typename T>
  void InputParser<T>::insert(fs::path path, const std::shared_ptr<KwargsParser> &subparser) {
    kwargs[path] = subparser;
  }

  template<typename T>
  InputParser<T>::map_type::const_iterator InputParser<T>::begin() const {
    return kwargs.begin();
  }

  template<typename T>
  InputParser<T>::map_type::const_iterator InputParser<T>::end() const {
    return kwargs.end();
  }

  template<typename T>
  InputParser<T>::map_type::const_iterator InputParser<T>::find(fs::path path) const {
    return kwargs.find(path);
  }

  template<typename T, typename ErrorType>
  void report_and_throw_if_invalid(InputParser<T> &parser, Log &log, ErrorType error) {
    if(!parser.valid()) {
      parser.print_errors(log);
      log << std::endl;
      log.indent() << parser.report() << std::endl << std::endl;
      throw error;
    }
    if(parser.all_warnings().size()) {
      parser.print_warnings(log);
      log << std::endl;
      log.indent() << parser.report() << std::endl << std::endl;
    }
  }

}

#endif
