#ifndef CASM_InputParser_impl
#define CASM_InputParser_impl

#include "casm/casm_io/InputParser.hh"
#include "casm/casm_io/Help.hh"

namespace CASM {

  /// equivalent to require_at fs::path(it.name()) / option
  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::require(jsonParser::const_iterator it, std::string option, Args &&...args) {
    return require_at<RequiredType>(fs::path(it.name()) / option, std::forward<Args>(args)...);
  }

  /// require option self_it->find(option) of type RequiredType
  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::require(std::string option, Args &&...args) {
    return require_at<RequiredType>(fs::path(option), std::forward<Args>(args)...);
  }

  /// require option self_it->find(option) of type RequiredType
  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::require_at(fs::path option, Args &&...args) {
    auto it = self.find_at(option);
    std::unique_ptr<RequiredType> res;
    if(it == self.end()) {
      error.insert(std::string("Error: missing required option '") + option.string() + "'.");
      return res;
    }

    try {
      res = notstd::make_unique<RequiredType>(it->get<RequiredType>(std::forward<Args>(args)...));
      return res;
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return res;
    }
  }

  /// equivalent to optional_at fs::path(it.name()) / option
  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::optional(jsonParser::const_iterator it, std::string option, Args &&...args) {
    return optional_at<RequiredType>(fs::path(it.name()) / option, std::forward<Args>(args)...);
  }

  /// check that if option self.find(option) exists it can constructed as type RequiredType
  template<typename RequiredType, typename...Args>
  std::unique_ptr<RequiredType> KwargsParser::optional(std::string option, Args &&...args) {
    return optional_at<RequiredType>(fs::path(option), std::forward<Args>(args)...);
  }

  /// check that if self.find_at(option) exists, it can constructed as type RequiredType
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
      return notstd::make_unique<RequiredType>(ptr->get<RequiredType>(std::forward<Args>(args)...));
    }
    catch(std::exception &e) {
      error.insert(std::string("Error: could not construct type '")
                   + type_name<RequiredType>() + "' from option '" + option.string() + "'. "
                   + singleline_help<RequiredType>());
      return std::unique_ptr<RequiredType>();
    }
  }

  /// check that if option self.find(option) exists it can constructed as type RequiredType
  template<typename RequiredType, typename...Args>
  RequiredType KwargsParser::optional_else(std::string option, const RequiredType &_default, Args &&...args) {
    return optional_at_else<RequiredType>(fs::path(option), _default, std::forward<Args>(args)...);
  }

  /// check that if self.find_at(option) exists, it can constructed as type RequiredType
  template<typename RequiredType, typename...Args>
  RequiredType KwargsParser::optional_at_else(fs::path option, const RequiredType &_default, Args &&...args) {
    auto res = optional_at(option, std::forward<Args>(args)...);
    return res ? *res : _default;
  }

  template<typename OptHandlerType>
  int KwargsParser::parse_verbosity(const OptHandlerType &opt) {
    std::string verbosity_str = opt.vm().count("verbosity") ?
                                opt.verbosity_str() :
                                optional_else<std::string>("verbosity", "standard");
    auto val = Log::verbosity_level(verbosity_str);
    if(val.first) {
      return val.second;
    }
    else {
      error.insert(Log::invalid_verbosity_msg(verbosity_str));
      return Log::standard;
    }
    return 0;
  }

  template<typename OptHandlerType>
  bool KwargsParser::parse_dry_run(const OptHandlerType &opt) {
    return opt.vm().count("dry-run") ?
           true :
           optional_else<bool>("dry_run", false);
  }

  template<typename OptHandlerType>
  COORD_TYPE KwargsParser::parse_coord_type(const OptHandlerType &opt) {
    return opt.vm().count("coord") ?
           opt.coordtype_enum() :
           optional_else<COORD_TYPE>(traits<COORD_TYPE>::name, COORD_TYPE::FRAC);
  }

  template<typename OptHandlerType>
  ORBIT_PRINT_MODE KwargsParser::parse_orbit_print_mode(const OptHandlerType &opt) {
    return optional_else<ORBIT_PRINT_MODE>(traits<ORBIT_PRINT_MODE>::name, ORBIT_PRINT_MODE::PROTO);
  }

}

#endif
