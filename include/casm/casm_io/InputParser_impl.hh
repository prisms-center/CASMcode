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

}

#endif
