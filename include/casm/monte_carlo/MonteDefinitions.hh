#ifndef CASM_MonteDefinitions_HH
#define CASM_MonteDefinitions_HH

#include <string>
#include <stdexcept>
#include <iostream>
#include "casm/external/boost.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  namespace Monte {

    template <typename T>
    struct traits {};

    template<typename ENUM>
    std::string help() {
      std::stringstream ss;
      ss << "Options are:\n";
      for(auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end(); ++it) {
        ss << "  ";
        for(auto sit = it->second.begin(); sit != it->second.end(); sit++) {
          if(sit != it->second.begin()) {
            ss << " or " << *sit;
          }
          else {
            ss << *sit;
          }
        }
        ss << "\n";
      }
      return ss.str();
    }

    template<typename ENUM>
    void invalid_enum_string(std::string val, std::ostream &serr) {
      std::stringstream s;
      s << "Invalid " << traits<ENUM>::name << ": " << val;
      serr << s.str() << "\n";
      serr << help<ENUM>();
      throw std::invalid_argument(std::string("ERROR: ") + s.str());
    }

    template<typename ENUM>
    std::string to_string(ENUM val) {
      return traits<ENUM>::strval.find(val)->second[0];
    }

    template<typename ENUM>
    ENUM from_string(const std::string &val) {
      for(auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end(); ++it) {
        for(auto sit = it->second.begin(); sit != it->second.end(); sit++) {
          if(*sit == val) {
            return it->first;
          }
        }
      }

      invalid_enum_string<ENUM>(val, std::cerr);
      return traits<ENUM>::strval.begin()->first;
    }

#define ENUM_DEF(ENUM) \
    template<> \
    struct traits<ENUM> { \
    \
      static const std::string name; \
    \
      static const std::multimap<ENUM, std::vector<std::string> > strval; \
    \
    }; \
    \
    inline std::ostream &operator<<(std::ostream &sout, const ENUM& val) { \
      sout << to_string<ENUM>(val); \
      return sout; \
    } \
    \
    inline std::istream &operator>>(std::istream &sin, ENUM& val) { \
      std::string s; \
      sin >> s; \
      val = from_string<ENUM>(s); \
      return sin; \
    } \
    \
    inline jsonParser &to_json(const ENUM &val, jsonParser &json) { \
      return to_json(to_string<ENUM>(val), json); \
    } \
    \
    inline void from_json(ENUM& val, const jsonParser& json) { \
      val = from_string<ENUM>(json.get<std::string>()); \
    } \
 

    /// \brief Monte Carlo ensemble type
    enum class ENSEMBLE {
      GrandCanonical
    };

    ENUM_DEF(ENSEMBLE)


    /// \brief Monte Carlo method type
    enum class METHOD {
      Metropolis, LTE1
    };

    ENUM_DEF(METHOD)


    ///How often to sample runs
    enum class SAMPLE_MODE {
      STEP, PASS
    };

    ENUM_DEF(SAMPLE_MODE)


    ///How to change conditions
    enum class DRIVE_MODE {
      INCREMENTAL, CUSTOM
    };

    ENUM_DEF(DRIVE_MODE)


  }

}
#endif


